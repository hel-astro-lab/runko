.. _entity-component-system:

Tiles as entities
#################

There is no ``Tile`` class that has a Yee lattice or particles as members.
Instead, a tile is only a number, and its data is stored elsewhere, grouped by type.
This is the *entity component system* (ECS) pattern,
and runko uses the `EnTT <https://github.com/skypjack/entt>`_ library
(``external/entt``) to implement it.

This page explains what that pattern is and how runko maps tiles onto it.


Problem ECS solves
==================

A runko simulation is split into *tiles*.
Each MPI rank owns some tiles, and every tile covers a small, regular part of the grid.
But not all tiles are alike:

* a *local* tile is owned by this rank and holds fields, and in PIC runs also particles,
* a *boundary* tile is a local tile with at least one neighbor on another rank,
  and it needs to know where to send its data,
* a *virtual* tile is a stand-in for a neighbor that lives on another rank.
  It holds only receive buffers and knows where its data comes from.

On top of that, what a tile holds depends on the run.
An FDTD-only run has no particles.
A run with the atomic current depositer needs a scratch grid for the current,
and a run with reflecting walls needs a correction current on some tiles.

The classic object-oriented answer is a hierarchy:
a ``Tile`` base class, an ``emf::Tile`` that adds a Yee lattice,
a ``pic::Tile`` that inherits from ``emf::Tile`` and adds particles,
and so on.
Earlier versions of runko were built like this, on top of the corgi grid library.
This works as long as the kinds of tiles form a tree,
but it breaks down when they do not:

* **Features do not form a tree.**
  "Is virtual", "is on a boundary", "has particles" and "has a current cache"
  are independent of each other.
  A hierarchy has to either multiply classes for every combination,
  or put everything in a base class and leave most of it unused.
* **Every new feature changes the base.**
  Adding a new kind of per-tile data means changing a class
  that every other module depends on.
* **Changing a tile's kind is awkward.**
  An object cannot change its class.
  Turning a tile into a boundary tile, or giving it particles later,
  means constructing a new object and copying its data over.
* **Loops are about data, not about classes.**
  The field pusher only cares about tiles with a Yee lattice.
  With a hierarchy it has to loop over all tiles and ask each one,
  with ``dynamic_cast`` or a type flag, whether it is the right kind.

ECS turns this around.
Instead of asking "what class is this tile?",
code asks "which tiles have these pieces of data?".


What is ECS?
============

ECS has three concepts.

**Entity.**
An entity is just an identifier, an integer.
It has no data and no behavior of its own.
In runko every tile is an entity,
and ``simulation_context::tile_id_type`` is EnTT's ``entt::registry::entity_type``.

**Component.**
A component is a plain C++ value attached to an entity.
Any movable type can be a component:
``emf::YeeLattice``, ``runko::cartesian_index<3>``,
a ``std::map`` of particle containers, or an empty struct.
An entity has at most one component of each type.
The *type* is the key, so asking an entity for its ``emf::YeeLattice``
is unambiguous.

**System.**
A system is ordinary code that loops over all entities
that have a given set of components, and works on them.
In runko, systems are the functions behind the actions:
``push_e``, ``deposit_current``, ``comm_local`` and so on.

The container that holds all of this is the *registry*.
``simulation_context`` has exactly one, and it is called ``tiles``:

.. code:: c++

   // src/runko/simulation_context.h
   struct simulation_context {
     // ...
     const toolbox::ConfigParser config;
     entt::registry tiles;
     using tile_id_type = entt::registry::entity_type;
   };

Internally the registry keeps one densely packed storage, a *pool*, per component type.
All ``emf::YeeLattice`` values are in one pool,
all ``runko::cartesian_index<3>`` values in another.
An entity "has" a component if the pool for that type contains an entry for its id.
There is no per-tile object that ties the pieces together.

One way to picture it is as a sparse table
in which rows are entities and columns are component types:

.. code:: text

             cartesian   local  boundary  virtual   YeeLattice  particle     comm
             _index<3>   _tag   _tag      _tag                  _containers  _buffs
   tile 0    (0,0,0)     x                          ...         ...
   tile 1    (1,0,0)     x      ...                 ...         ...          ...
   tile 2    (2,0,0)                      ...                                ...
   tile 3    (-1,0,0)                     ...                                ...

Tile 0 is an interior local tile, tile 1 is a local tile on a rank boundary,
and tiles 2 and 3 are virtual tiles standing in for neighbors on other ranks.
The "kind" of a tile is nothing more than which columns are filled in.

This way of representing program state is desirable for Runko,
but we don't want to maintain our own implementation of it.
Many simpler alternatives are just badly implemented ECS in disguise.


How runko uses the registry
===========================

Creating tiles
--------------

A tile is created by making a new entity and attaching components to it.
This is the whole of ``add_tile`` in ``src/runko/bindings/pyactions.c++``:

.. code:: c++

   const auto id = sim.tiles.create();
   sim.tiles.emplace<runko::cartesian_index<3>>(id, idx);
   sim.tiles.emplace<runko::local_tile_tag>(id);

A new tile only knows where it is and that it is local.
Everything else is added later, by the actions that need it.

Tags
----

``runko::local_tile_tag`` is an empty struct.
Its only purpose is to be present or absent.
Such components are called *tags*.
EnTT stores no values for empty types, only which entities have them,
so a tag costs almost nothing.

Not all "tags" in runko are empty.
``boundary_tile_tag`` holds the list of ranks and MPI tags to send to,
and ``virtual_tile_tag`` holds where the data comes from
(``src/runko/comm/cartesian_grid.h``).
They mark a kind of tile *and* carry the data that only that kind needs.

Views
-----

A system finds its tiles with a *view*.
``simulation_context::view_tiles<Ts...>()`` is a thin wrapper around
``tiles.view<Ts...>().each()``.
It iterates over every entity that has *all* of ``Ts``,
and yields a tuple of the entity id followed by references to its components:

.. code:: c++

   for(auto&& [id, yee, particles, idx]: sim.view_tiles<
                                        emf::YeeLattice,
                                        pic::particle_containers,
                                        const runko::cartesian_index<3>,
                                        runko::local_tile_tag>()) {
     // yee:       emf::YeeLattice&
     // particles: pic::particle_containers&
     // idx:       const runko::cartesian_index<3>&
   }

A few details are easy to miss when reading such loops:

* The entity id always comes first. Most loops ignore it with ``_``.
* **Empty types do not appear in the tuple.**
  ``local_tile_tag`` filters the tiles but gives no binding.
  That is why ``view_tiles<const index_type, local_tile_tag>()`` binds ``[_, index]``,
  while ``view_tiles<const index_type, virtual_tile_tag>()`` binds ``[_, index, virt]``,
  because ``virtual_tile_tag`` has a member.
* ``const T`` in the list gives a ``const T&`` and documents that the loop only reads it.
  The ``const`` overload of ``view_tiles`` adds ``const`` to every type.

To select tiles that *lack* a component,
use the registry directly with ``entt::exclude``:

.. code:: c++

   for(auto&& id: sim.tiles.view<runko::local_tile_tag>(entt::exclude<emf::YeeLattice>)) {
     // local tiles that do not have a Yee lattice yet
   }

Single-tile access
------------------

A tile id, for example a neighbor id from ``runko::cartesian_neighbors<3>``,
can be used to access the tile directly:

* ``tiles.try_get<T>(id)`` returns a pointer, or ``nullptr`` if the tile has no ``T``.
  This is the ECS way to ask "what kind of tile is this?".
  ``comm_local`` in ``src/runko/comm/emf.c++`` copies from the neighbor's
  ``YeeLattice`` if it has one, and from its ``comm_buffs`` otherwise,
  which is exactly the local-versus-virtual distinction.
* ``tiles.emplace<T>(id, args...)`` constructs a component.
  It is an error if the tile already has one.
* ``tiles.emplace_or_replace<T>(id, args...)`` constructs or overwrites it.
* ``tiles.get<T>(id)`` returns a reference and assumes the component exists.
  Prefer ``try_get`` and a clear error message.

Tile ids leave C++ as plain integers.
``simulation_context.get_local_tile_ids()`` returns them to Python,
and ``ProxyTile`` in ``runko/tiles.py`` holds one together with its context.
Methods such as ``get_EBJ(id)`` then look up the component with ``try_get``.
An id is only meaningful for the registry it came from, and it is local to one rank.
The global identity of a tile is its ``cartesian_index<3>``.

Adding data lazily: the ``ensure_constructed`` pattern
------------------------------------------------------

Because components can be added at any time,
runko does not build all of a tile's data up front.
Instead, the actions that need a component make sure it exists.
``emf::ensure_constructed_yee_lattices`` is a typical example:

.. code:: c++

   auto where_to_construct = std::vector<runko::simulation_context::tile_id_type> {};

   for(auto&& id: sim.get().tiles.view<runko::local_tile_tag>(entt::exclude<emf::YeeLattice>)) {
     where_to_construct.push_back(id);
   }
   if(where_to_construct.empty()) { return; }

   // ... read the configuration once ...
   std::ranges::for_each(where_to_construct, [&](const auto x) {
     sim.get().tiles.template emplace<emf::YeeLattice>(x, args);
   });

It is *idempotent*: calling it twice does nothing the second time.
That makes it safe to put in front of any program that needs Yee lattices.
The same pattern is used for particle containers and for communication buffers.
Note that it first collects the ids and only then emplaces.
The view excludes ``YeeLattice``,
so adding a ``YeeLattice`` while iterating would change the set the view is iterating over.

A component can also be a private cache.
``pic::deposit_current`` defines a struct ``J_cache`` inside the function,
and attaches one to each tile the first time the atomic depositer runs.
No other code knows the type exists, so no other code can depend on it.

Context variables
-----------------

Some data belongs to the simulation, not to any tile:
the chosen field propagator, the stencil coefficients,
registered antennas, edge boundary conditions and reflector walls.
EnTT stores such singletons as *context variables* in ``tiles.ctx()``,
which is again keyed by type.
``simulation_context`` implements a few wrappers for the context variables.

For example,
``register_reflector_wall`` appends to a ``pic::reflectors`` context variable,
and ``reflect_particles`` later reads it, returning early if nothing has been registered.


Type is the name
----------------

The registry identifies components and context variables by their C++ type.
Defining a new type is the way to name things.
For a new kind of data, give it its own ``struct``,
even if it only wraps a single member.

.. note::

   C++ type aliases are indistinguishable from the underlying type:

   .. code:: c++

      using vec = std::vector<int>;

      // The following are equivalent:
      tiles.get<vec>(id);
      tiles.get<std::vector<int>>(id);


Filter with components, not with flags
--------------------------------------

If a system should only touch some tiles, express that as a component in the view,
not as an ``if`` inside the loop.
Almost every system over fields or particles includes ``local_tile_tag``,
because virtual tiles are not evolved.


Structural changes and iteration
--------------------------------

Creating or destroying entities, and adding or removing components,
are *structural* changes to the registry.
Modifying the *values* of the components a view yields is always fine.
Structural changes to pools that the view iterates over, including excluded types, are not.
Either collect ids first, as ``ensure_constructed_yee_lattices`` does,
or make sure the type added is unrelated to the view.
``set_cartesian_neighbors`` adds ``boundary_tile_tag`` while iterating over
``cartesian_index<3>`` and ``local_tile_tag``, which is fine because the tag is not part of that view.

Also do not keep references or pointers to components across code
that removes components of the same type.
EnTT keeps each pool packed, so removing one element may move another.

Iteration order is unspecified
------------------------------

A view visits entities in the order of its pools, which depends on the history of insertions
and removals.
It is not sorted by ``cartesian_index`` and is not the same on every rank.
Do not rely on it for anything that has to match across ranks, such as the order of MPI messages.
Runko's communication matches messages by MPI tag, computed from the cartesian index,
for this reason.

The registry is not thread-safe
-------------------------------

EnTT does no locking.
If concurrent work on pika's thread pool needs to access the same data,
there has to be some kind of synchronization. Otherwise there will be a data race.
Views and reads of *different* pools from different threads are fine as long as nobody
changes the structure at the same time.

Components are touched when the sender runs
-------------------------------------------

As explained in :ref:`senders-and-receivers`, code in an action runs at two different times.
Reading context variables that come from the configuration, such as the field propagator,
can happen when the sender is built.
Iterating over tiles and emplacing components must happen inside a sender,
so that it sees the tiles as earlier steps of the program left them.

Adding new per-tile data
------------------------

Adding a new kind of per-tile data usually looks like this:

#. define a ``struct`` for it next to the code that uses it,
#. write an idempotent ``ensure_constructed_*`` function that emplaces it
   on the tiles that need it, and expose it as an action
   (see :ref:`actions-language`),
#. write systems that ``view_tiles`` over it together with the tags that select the right tiles,
#. if Python needs to read it, add a method to ``SimulationContext``
   in ``src/runko/bindings/pyactions.c++`` that takes a tile id and uses ``try_get``.


Further reading
===============

* In Runko: ``src/runko/simulation_context.h`` (the registry and the context-variable
  wrappers), ``src/runko/comm/cartesian_grid.h`` (tile kinds and how neighbors and virtual tiles
  are created), ``src/runko/actions/emf.c++`` and ``src/runko/actions/pic.c++``
  (systems and ``ensure_constructed`` functions), and ``src/runko/bindings/pyactions.c++``
  (creating tiles and accessing them from Python).
* :ref:`actions-language`, for how the systems are called from Python.
* `EnTT wiki: Entity Component System <https://github.com/skypjack/entt/wiki/Entity-Component-System>`_,
  the reference for the registry, views, ``exclude`` and context variables.
* `EnTT Doxygen documentation <https://skypjack.github.io/entt/>`_ for API reference.
* `ECS back and forth <https://skypjack.github.io/2019-02-14-ecs-baf-part-1/>`_,
  a series by the author of EnTT on how ECS libraries are designed.
* `The ECS FAQ <https://github.com/SanderMertens/ecs-faq>`_,
  a general introduction to the pattern that is not tied to any library.
