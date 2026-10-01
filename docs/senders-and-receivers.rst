.. _senders-and-receivers:

Senders and receivers
#####################

C++ code in ``src/runko/actions/`` or ``src/runko/comm/``
contains functions that return a ``tyvi::actions::sexpr_sender``
and bodies that look like this:

.. code:: c++

   return te::just(std::ref(x)) | te::then([prop, cfl](runko::simulation_context& sim) {
            for(auto&& [_, yee]: sim.view_tiles<emf::YeeLattice, runko::local_tile_tag>()) {
              yee.push_e_fdtd2(static_cast<vt>(cfl));
            }
            return ta::null;
          });

This page explains what this style of code is and why runko uses it.


Problem senders solve
=====================

A simulation time step in runko is a series of operations:
push the fields, push the particles, deposit current, and exchange data with neighboring tiles,
both between tiles in the same process and between MPI ranks.
Most of these operations are ordinary loops over tiles.
Some of them, like MPI communication, are *asynchronous*.

The simplest way to write this is a list of function calls that each finish before returning.
That works, but it has two drawbacks:

* **No overlap.** A blocking ``MPI_Recv`` makes the CPU core sit idle while it waits.
  If we have dozens of neighbors to talk to,
  we want all of those messages to be in flight at the same time,
  and we want to continue as soon as all of them have arrived.
* **No common shape.** Some operations return a value, some return nothing, some might fail,
  and some are asynchronous.
  If every operation has its own calling convention
  (callbacks, futures, ``MPI_Request`` handles, plain return values),
  it becomes hard to build generic tools that combine them.

*Senders and receivers* are a C++ model of asynchronous work that addresses both issues.
It was standardized as ``std::execution`` in C++26
(`proposal P2300: std::execution <https://wg21.link/p2300>`_).
Runko uses the implementation provided by the `pika <https://pikacpp.org/>`_ library,
exposed through the tyvi portability library:

.. code:: c++

   // external/tyvi/src/tyvi/execution.h
   namespace tyvi {
   namespace exec        = pika::execution::experimental;
   namespace this_thread = pika::this_thread::experimental;
   }

Inside runko the convention is ``namespace te = tyvi::exec;``,
so ``te::just``, ``te::then`` and so on come from pika.
The names match ``std::execution``,
so most material written about the standard version also applies here.

.. note::
   Pika is exposed via tyvi because the standard library support
   for senders and receivers is poor at the moment.
   Ideally, ``tyvi::actions`` would only depend on the standard library
   and users of tyvi would be able to choose any scheduler (e.g. pika).


Senders are descriptions of lazy work
=====================================

A **sender** is an object that *describes* some work.
Creating a sender does not do the work.
Every piece of work finishes in exactly one of three ways,
called *completion channels*:

* **value**: the work succeeded and produced zero or more values,
* **error**: the work failed, usually with a C++ exception,
* **stopped**: the work was cancelled before it finished.

A **receiver** is the other half.
It is the object that is told how the work ended.
A receiver has three entry points that match the three channels:
``set_value(...)``, ``set_error(e)`` and ``set_stopped()``.
The work calls exactly one of them, exactly once.

To actually run the work, a sender is *connected* to a receiver.
This produces an *operation state*, an object that holds everything the work needs while it runs.
The operation state is then *started*.

.. code:: c++

   // This operation state holds everything the work needs while it runs.
   auto op = std::exec::connect(some_sender, some_receiver);
   std::exec::start(op);
   // After the async operation completes, one of the Receiver's
   // set_value / set_error / set_stopped will be called.


Developers of runko will **almost never write a receiver or call connect/start explicitly**.
The standard library is responsible for that.
The developer's job is the **composition of senders**.
Receivers are still worth knowing about,
because they explain why senders behave the way they do.
A sender cannot "return" its result like a function.
It can only pass it on to whatever receiver it eventually gets connected to.
This is also why the model works for asynchronous work:
the result can arrive on another thread, at any later time, and still reach the right place.


Building senders out of smaller senders
=======================================

Senders become useful when combined.
The library provides *sender factories*, which create senders from nothing,
and *sender adaptors*, which take a sender and return a new, bigger sender.
Adaptors are usually chained with ``|``, which reads left to right like a Unix pipe.
``a | f(...)`` means "when ``a`` completes with a value, continue with ``f``".

These are the building blocks that appear in runko's source.

``te::just(args...)``
   A factory: a sender that immediately completes with ``args...`` as its values.
   It is the usual starting point of a chain.
   ``te::just()`` with no arguments completes with no values.
   Useful when there is nothing to pass on yet
   and the following steps require a starting point.

``sender | te::then(f)``
   When ``sender`` completes with values ``v...``, call ``f(v...)``.
   Whatever ``f`` returns becomes the value of the new sender.
   If ``f`` throws, the new sender completes on the error channel instead.
   Most of runko's physics actions are a single ``then``
   that loops over the local tiles, like ``push_e`` above.

``sender | te::let_value(f)``
   Like ``then``, but ``f`` returns *another sender*,
   and the result is whatever that inner sender produces.
   Use it when the next step is itself asynchronous work,
   or when the value is required to decide which work to do next.
   ``runko::comm_local`` is an example:
   it picks one of several communication routines based on the ``comm_mode``,
   and each routine returns a sender.

   .. code:: c++

      auto f = [sim = std::ref(sim), mode] {
        switch(mode) {
          case runko::comm_mode::emf_E: /* ... */ return emf::comm_local(sim, mode);
          case runko::comm_mode::pic_particle:    return pic::comm_local_particles(sim);
          // ...
        }
      };
      return te::just() | te::let_value(f);


   ``let_value`` also keeps the values it received alive
   for as long as the inner sender runs, so the inner work can safely refer to them.

   .. tip::

      If the lambda returns a plain value, use ``then``. If it returns a sender, use ``let_value``.

``te::when_all(s1, s2, ...)`` and ``te::when_all_vector(vec)``
   Run several senders and complete when *all* of them have completed.
   The values of all inputs are combined into one completion.
   If any input fails, the combined sender fails.
   This is how runko expresses "do these things concurrently, then continue".
   ``emf::comm_external`` builds one sender per MPI message and joins them:

   .. code:: c++

      return te::just() |
             te::then(std::bind_front(&ensure_constructed_emf_comm_buffs, std::ref(x))) |
             te::then(std::bind_front(&update_emf_send_buff, std::ref(x), mode)) |
             te::let_value([sim = std::ref(x), sends, recvs] {
               return te::when_all(sends(sim), recvs(sim));
             }) |
             te::then([] { return tyvi::actions::null; });

   Read it top to bottom: make sure the buffers exist, fill the send buffers,
   then start all sends and receives at once and wait for them all, then return ``null``.
   ``when_all_vector`` is the pika extension used when the number of senders
   is only known at run time, for example one per neighbor tile.

``sender | te::continues_on(scheduler)``
   Continue the rest of the chain on the execution resource given by ``scheduler``.
   In runko this is ``te::thread_pool_scheduler{}``, pika's pool of worker threads.
   Without it, all work would run on whichever thread happened to complete the previous step.
   With it, independent branches of a ``when_all`` can run in parallel.
   The ``mt_showcase`` action in ``src/runko/actions/env.c++`` exists only to show this.

``sender | pmpi::transform_mpi(MPI_Isend)``
   A pika adaptor (``pika::mpi::experimental``) that turns a non-blocking MPI call into a sender.
   It calls the MPI function with the incoming values plus an ``MPI_Request``,
   and completes the sender once that request has finished.
   Explicit ``MPI_Wait`` or ``MPI_Test`` is not required.
   Pika polls outstanding requests in the background.
   The ``RuntimeActivator`` in ``src/runko/runtime.h`` turns that polling on.

``tyvi::this_thread::sync_wait(sender)``
   The bridge back to ordinary, blocking code.
   It connects the sender to a receiver of its own, starts it, blocks the calling thread
   until the receiver is called, and then returns the value or rethrows the error.

Put together, a chain like

.. code:: c++

   te::just(ptr, count, MPI_FLOAT, rank, tag, comm)
     | te::continues_on(te::thread_pool_scheduler {})
     | pmpi::transform_mpi(MPI_Irecv)

reads as: "take these arguments, move to a worker thread, post a non-blocking receive,
and complete when the message has arrived".
Nothing happens when this expression is evaluated.
It only builds the description.


Where senders fit in runko
==========================

Runko's C++ code is driven from Python through a small interpreter called *actions*
(``tyvi/actions*.h``).
The Python ``Simulation`` class builds nested tuples such as
``(actions.push_e, actions.current_context)``.
These are converted to *s-expressions* (``tyvi::actions::sexpr``),
a simple tree of atoms and lists, and evaluated.
:ref:`actions-language` explains this language in detail.

The important design choice is that **evaluating an s-expression produces a sender**,
not a result.
Every procedure the interpreter can call has the same signature:

.. code:: c++

   // external/tyvi/src/tyvi/actions_ast.h
   using sexpr_sender = exec::unique_any_sender<sexpr>;
   using procedure    = std::function<exec::unique_any_sender<sexpr>(sexpr)>;

A procedure takes its (already evaluated) arguments as an ``sexpr``
and returns a sender that will eventually complete with an ``sexpr``.
The symbol table that maps names to procedures is built in
``src/runko/actions/env.c++`` and ``src/runko/actions/sim_env.c++``.
A typical entry parses its arguments into a sender
(``runko::parse_atom_args`` returns a ``te::just`` of the parsed values)
and pipes it into the real implementation:

.. code:: c++

   ta::cons(
     runko::symbol::ensure_constructed_yee_lattices,
     ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
       return parse_atom_args<std::reference_wrapper<runko::simulation_context>>(args) |
              te::then(&emf::ensure_constructed_yee_lattices) |
              te::then([] { return ta::null; });
     } }),

The interpreter itself (``tyvi/actions_eval.h``) is written with the same tools.
To evaluate a call ``(f a b)`` it evaluates ``f``, ``a`` and ``b`` with
``when_all``, then uses ``let_value`` to invoke the procedure on the results.
Finally, ``src/runko/bindings/pyactions.c++`` calls ``sync_wait`` on the whole tree.
This is the single point where Python waits for C++ to finish.

So there is one uniform contract all the way through:
*every runko action is a function that returns a sender*.
That uniformity is what lets the interpreter compose physics kernels, MPI communication and I/O
without knowing anything about them.


Type erasure: why ``unique_any_sender``
---------------------------------------

Each adaptor produces a new, distinct and often very long type.
``te::just(x) | te::then(f)`` has a different type from ``te::just(x) | te::then(g)``.
This is good for performance, because the compiler sees the whole chain,
but it makes it impossible to store different senders in one ``std::vector``
or to return them from a ``std::function``.

``unique_any_sender<T...>`` is a *type-erased* wrapper.
It can hold any sender that completes with ``T...``, in the same way that ``std::function``
can hold any callable with a given signature.
Runko uses it at the boundaries: the ``procedure`` signature, and the vectors passed to
``when_all_vector``. ``unique_`` means it is move-only. You can pass it on, but not copy it,
which is why you see ``std::move(senders)`` around these vectors.

Two phases: building and running
--------------------------------

Code in an action function runs at two different times.
Code *outside* the lambdas runs when the sender is **built**.
Code *inside* ``then``/``let_value`` lambdas runs when the sender is **executed**.
Look at ``emf::push_half_b``:

.. code:: c++

   const auto cfl  = x.config.template get_or_throw<double>("cfl");                 // built
   const auto prop = x.get_n_set_config<emf::FieldPropagator>(&parse_field_propagator);

   return te::just(std::ref(x)) | te::then([prop, cfl](runko::simulation_context& sim) {
            for(auto&& [_, yee]: sim.view_tiles<...>()) { /* ... */ }                  // executed
            return ta::null;
          });

Reading the configuration early is fine, because it does not change during a step.
Touching the tiles must happen inside the lambda.
Otherwise it would run before earlier steps in the chain have finished.
When you are unsure, put the work inside a ``then``.
That way it runs in the order the chain describes.

Lifetimes and captures
----------------------

A sender may run after the function that created it has returned,
so lambdas must not capture local variables by reference.
Runko code captures small values by copy (``[prop, cfl]``)
and the large, long-lived ``simulation_context`` through ``std::ref``/``std::reference_wrapper``.
That is safe because the context outlives any evaluation.
If data is needed that lives during the operation,
pass it as a value through the chain (for example ``te::just(buffer) | te::let_value(...)``),
so that the operation state owns it.

Concurrency is opt-in and explicit
----------------------------------

Steps joined with ``|`` always run one after another.
Only branches of a ``when_all`` or ``when_all_vector`` may run at the same time,
and only in parallel if they have been moved to the thread pool with ``continues_on``.
If two branches touch the same tile data, it is a data race unless other means
of synchronization are used.

A related detail: the interpreter evaluates the arguments of a call with ``when_all``,
so there is no guarantee that they are evaluated in any particular order.
The ``sequence`` procedure in ``src/runko/actions/env.c++`` exists for this reason.
It evaluates its arguments one at a time, calling ``sync_wait`` on each one in turn.

Errors travel as exceptions
---------------------------

Throwing inside a ``then`` lambda does not crash the process.
The exception is caught and sent down the error channel,
``when_all`` passes it on, and ``sync_wait`` rethrows it.
``pyactions.c++`` catches it there and wraps it in a ``std::runtime_error``,
which pybind11 turns into a Python exception.
So in an action, plain ``throw std::runtime_error{...}`` is the correct way to report a problem.

Not every "work" or ``when_all`` is a sender
--------------------------------------------

The GPU code in ``src/runko/emf/`` and ``src/runko/pic/`` uses ``tyvi::mdgrid_work``
and a function also called ``tyvi::when_all(w1, w2, ...)``.
Despite the name, these are not senders.
They are tyvi's own mechanism for ordering and waiting on device kernels.
In practice, sender chains decide *which* operation runs *when*,
while the tile methods they call (for example ``YeeLattice::push_e_fdtd2``)
launch and wait for GPU work internally with ``mdgrid_work``.

.. note::
   The long-term plan is to move away from ``tyvi::mdgrid_work`` and replace it with senders.


Why not something simpler?
==========================

Some alternatives, and why they fit runko less well:

* **Plain blocking calls.** These are the easiest to read,
  but they give up overlapping MPI messages with each other and with computation.
  They also leave no generic way for the Python-driven interpreter
  to treat computation and communication the same way.
* **Futures** (``std::future``, ``std::async``) represent one value that is already being computed.
  They start work eagerly, usually allocate, and chaining them
  (``.then`` continuations) was never standardized.
  Senders are lazy, so a whole pipeline can be built, inspected and optimized before anything runs.
* **Hand-written ``MPI_Request`` bookkeeping.** This can be efficient,
  but every communication routine then has to manage arrays of requests and wait calls,
  and none of it composes with the rest of the code.

Senders add verbosity and produce longer compilation errors
while providing one composable model for synchronous work, threaded work and MPI.
The same model is the one the C++ standard is adopting.


Further reading
===============

* In Runko: ``src/runko/actions/env.c++`` (small examples),
  ``src/runko/comm/emf.c++`` (MPI with senders) and ``external/tyvi/src/tyvi/actions_eval.h``
  (the interpreter).
* `pika documentation <https://pikacpp.org/>`_, especially the sections on
  ``pika::execution::experimental`` and ``pika::mpi::experimental``.
* `cppreference: Execution control library <https://en.cppreference.com/w/cpp/execution>`_
* `Documentation <https://nvidia.github.io/stdexec/>`_ of the `reference implementation <https://github.com/NVIDIA/stdexec>`_
  of ``std::exec`` by NVIDIA.
* `Documentation <https://intel.github.io/cpp-baremetal-senders-and-receivers/>`_ of the
  `C++ Baremetal Senders & Receivers <https://github.com/intel/cpp-baremetal-senders-and-receivers>`_ by Intel.
* `P2300: std::execution <https://wg21.link/p2300>`_, the proposal that introduced senders and receivers.
  The introduction and "motivation" sections are readable without deep C++ knowledge.
