.. _actions-language:

Actions language
################

Often Python code in ``runko/simulation.py`` or ``runko/tiles.py``
does not call C++ functions directly.
Instead, nested tuples are built and handed to ``eval``:

.. code:: python

   import runko_cpp_bindings.actions as actions

   set_neighs     = (actions.set_cartesian_neighbors, actions.current_context)
   set_comm_infos = (actions.set_cartesian_comm_infos, actions.current_context)

   self._simulation_context.eval((actions.sequence,
                                  (actions.quote, set_neighs),
                                  (actions.quote, set_comm_infos)))

These tuples are programs in a very small, Lisp-like language.
The language is implemented in C++ by ``tyvi::actions``
(``external/tyvi/src/tyvi/actions*.h``).
This page explains what that language is and how runko uses it between Python and C++.

This page builds on :ref:`senders-and-receivers`.
The main points are that a *sender* is a lazy description of work,
and ``sync_wait`` starts the sender and blocks until it finishes.


Problem the language solves
===========================

Runko is split in two halves.
Python is where a simulation is *configured and driven*.
It reads the configuration, creates tiles, sets initial conditions from Python functions,
and decides what happens in each lap.
C++ is where the *work* happens:
field and particle pushers on the GPU, and communication between tiles and MPI ranks.
Much of that C++ work is asynchronous and is written with senders,
so that many MPI messages can be in flight at once and independent work can run on
pika's thread pool.

This raises the question: how should Python ask C++ to do something?

The obvious answer is to bind every C++ function with
`pybind11 <https://pybind11.readthedocs.io/en/stable/>`_ and call it.
That works for synchronous functions, but it breaks down for asynchronous work:

* **Senders cannot cross the boundary.**
  A sender's type depends on the whole chain that built it,
  for example ``te::just(x) | te::then(f) | te::let_value(g)``.
  Such types cannot be meaningfully exposed to Python one by one.
  If each binding has to return a finished result,
  every call has to ``sync_wait`` before returning,
  and all the concurrency the senders describe is lost at the boundary.
* **Python threads do not help.**
  The obvious workaround is to start the blocking calls from several Python threads,
  or to wrap them in ``asyncio``.
  Because of the global interpreter lock (GIL), Python threads interleave rather than run in parallel.
  They also add a second threading system that has to cooperate with pika's worker threads, MPI,
  and GPU streams, and that makes reasoning about ordering and deadlocks much harder.
  It also puts the scheduling decisions in Python,
  where the C++ side cannot see or optimize them.
* **Many small calls are slow and scattered.**
  If a lap is a long list of Python-to-C++ calls, every call pays the boundary cost,
  and the structure of the lap exists only as Python control flow.
  There is no single object that describes "what this lap does".

Runko's answer is to have Python send *data*, not calls.
Python builds a description of the work as a binary tree.
C++ turns that tree into a single sender, and then waits for it once.
The tree is the program, and ``tyvi::actions`` is the Lisp-inspired language it is written in.


S-expressions
=============

The language is built on *s-expressions*, the data structure that Lisp uses for both code and data.
In ``tyvi/actions_ast.h`` an s-expression is one of three things:

.. code:: c++

   using sexpr = std::variant<null_type, cons, atom>;

``null``
   The empty list, written ``()``.

``atom``
   A single value of *any* C++ type.
   An atom is type-erased: it stores the value in a way that hides the actual type.
   Most types are implicitly converted to ``tyvi::actions::atom`` when used in a place
   that expects an atom.
   Thus, for example, a ``long``, a ``double``, a ``std::string``, an enum value,
   a ``std::reference_wrapper<simulation_context>``, a ``pybind11::function``
   and a ``tyvi::actions::procedure`` all work as atoms.
   The contained value is obtained with ``atom_cast<T>``,
   which returns a ``std::optional<T>`` that holds the value only if ``T`` is *exactly* the stored type.

``cons``
   A pair of two s-expressions, traditionally called ``car`` (the first element)
   and ``cdr`` (the rest). Denoted using ``(a . b)``, where ``a`` is the ``car``
   and ``b`` is the ``cdr``.


S-expressions form binary trees in which cons cells
are the nodes and atoms are the children. Null denotes missing children.
For example, the s-expression ``(a . ((b . (() . c)) . d))`` represents the following binary tree:

.. code::

         .
        / \
       a   .
          / \
         .   d
        / \
       b   .
          / \
        ()   c


Many data structures can be represented using only s-expressions.

Lists
-----

A list is a chain of ``cons`` cells that ends in ``null``.
For example, ``(a . (b . (c . ())))`` is a list that holds elements ``a``, ``b`` and ``c``.

.. code::

         .
        / \
       a   .
          / \
         b   .
            / \
           c   ()


By convention, lists are usually written without the cons dots, nested parentheses
and terminating null. So, ``(a . (b . (c . ())))`` is equivalent to ``(a b c)``.
As syntactic sugar, ``ta::list(a, b, c)`` builds such a chain in C++.
In Python there are no explicit syntax constructs.
Python tuples are interpreted as lists of their elements.
``parse_element`` in ``src/runko/bindings/pyactions.c++``
converts Python objects into s-expressions:

* a Python ``tuple`` becomes a list, converted element by element,
* ``int`` becomes an atom holding ``long``,
  ``float`` an atom holding ``double``, and ``str`` an atom holding ``std::string``,
* members of the bound enums (``actions.symbol``, ``actions.intrinsic``,
  ``comm_mode``, ``antenna_mode``, ``edge_bc``) become atoms holding that enum value,
* Python callables and ``list`` objects become atoms that hold the Python object itself,
  and a few bound C++ classes (for example ``reflector_wall``) are stored by value,
* anything else raises ``"Trying to parse unsupported type."``.

So ``(actions.push_e, actions.current_context)`` in Python is the two-element list
``(push_e current_context)`` in C++,
whose elements are atoms that hold ``runko::symbol`` enum values
which are exposed to Python from the runko C++ bindings module:

.. code:: python

   import runko_cpp_bindings.actions as actions
   # ...
   sim_ctx.eval((actions.emf_snapshot, actions.current_context, self.lap))
   actions.empty_context_eval((actions.println, (actions.format, "running", actions.version)))
   prog = (actions.sequence,
           (actions.quote, (actions.comm_external, actions.current_context, mode_a)),
           (actions.quote, (actions.comm_external, actions.current_context, mode_b)))
   sim_ctx.eval(prog)


Improper lists
--------------

A list-like s-expression that does not end in a null is called an *improper list*.
For example, ``(a . (b . c))`` is an improper list:

.. code::

     .
    / \
   a   .
      / \
     b   c


Facilities in Runko try to make it hard to create an improper list
and thus they are rarely encountered.


Association lists
-----------------

An *association list*, or *alist* for short, is a structure that holds key-value pairs.
It is a special case of a list in which each element of the list is a cons cell
whose ``car`` is the key and ``cdr`` is the value.
For example, ``((A . 1) (B . 2) (C . 3))`` is an alist that maps ``A`` to ``1``,
``B`` to ``2`` and ``C`` to ``3``.

.. code::

         .
        / \
       .   \
      / \   \
     A   1   .
            / \
           .   \
          / \   \
         B   2   .
                / \
               .   \
              / \   \
             C   3   ()


How a program is evaluated
==========================

``tyvi::actions::eval<Symbols...>(body, env)`` in
``tyvi/actions_eval.h`` is the entire interpreter.
It has only a handful of rules.

**Symbols are looked up.**
An atom whose type is in ``Symbols...`` is treated as a symbol.
When the interpreter encounters a symbol it looks up the corresponding s-expression
from the *environment*, which is an alist that maps symbols to s-expressions.
In Runko, ``runko::symbol`` is used as the symbol type.

**Other atoms evaluate to themselves.**
``42``, ``1.5``, a string or a ``comm_mode`` value is simply returned.

**A list is a call.**
To evaluate ``(f a b)``, the interpreter recursively evaluates ``f``, ``a`` and ``b``.
The value of ``f`` must be an atom holding a procedure.
The procedure is then called with the list of evaluated arguments ``(a' b')``.

**Quoting stops evaluation.**
``(quote x)`` evaluates to ``x`` itself, unevaluated.
This is the one *special form* built into the interpreter.
``quote`` is recognised only in the head position of a list.
It is represented by an atom containing ``tyvi::actions::intrinsic::quote``.
This is the trick that allows Lisp to treat code as data.

**The empty list cannot be evaluated.**
Evaluating ``()`` is an error.


The language has no user-defined variables, no ``lambda``, no ``if`` and no loops.
This is by design. It is not meant for writing algorithms.
It is meant for *wiring together* operations that are written in C++.
Any logic that needs loops or conditions stays either in Python,
which decides which program to build, or in C++, inside a procedure.

Environments
------------

The environment is an *association list*:
a list of ``(symbol . value)`` pairs that is searched from front to back.
Runko builds two of them in ``src/runko/actions/``:

* ``build_std_env()`` in ``env.c++`` contains procedures that do not need a simulation,
  such as ``print``, ``println``, ``format``, ``sequence`` and ``mt_showcase``,
  and the plain value ``version``.
* ``build_sim_env(sim)`` in ``sim_env.c++`` contains the physics, communication and I/O
  procedures, and appends the standard environment at the end.

Not every entry is a procedure.
``current_context`` is bound to ``std::ref(sim)``.
Evaluating the symbol ``current_context`` therefore produces an atom holding a reference to
the ``simulation_context`` that ``eval`` was called on.
This is how ``(push_e current_context)`` passes the simulation to ``push_e``
without Python ever having to hold a C++ reference.

Python reaches the two environments through two entry points:
``actions.empty_context_eval(prog)`` evaluates in the standard environment,
and ``SimulationContext.eval(prog)`` evaluates in the simulation environment of that context.

Procedures
----------

A procedure is an ordinary C++ callable with one fixed signature:

.. code:: c++

   using procedure = std::function<sexpr_sender(sexpr)>;

It receives its already evaluated arguments as one list,
and returns a *sender* that will eventually complete with an s-expression.
Most procedures in runko are thin adaptors.
They unpack the argument list into typed C++ values with ``runko::parse_atom_args<T...>``,
which checks the number and exact types of the arguments and returns ``te::just(values...)``.
They then pipe the values into the real implementation:

.. code:: c++

   ta::cons(
     runko::symbol::comm_external,
     ta::procedure { [](const ta::sexpr& args) -> ta::sexpr_sender {
       return parse_atom_args<std::reference_wrapper<simulation_context>,
                              runko::comm_mode>(args) |
              te::let_value(&runko::comm_external);
     } }),

A second kind, ``procedure_with_eval``, also receives the evaluator itself:

.. code:: c++

   using procedure_with_eval = std::function<sexpr_sender(sexpr, procedure)>;

This lets a procedure decide *when* and *whether* to evaluate parts of its input,
so that new special forms can be written in C++ without changing the interpreter
(e.g. ``sequence`` and ``parallel``).


Evaluation produces a sender
----------------------------

The most important property of ``eval`` is its return type.
It does not return a result. It returns a ``sexpr_sender``.
Every rule above is written with sender adaptors:

* a symbol lookup or a self-evaluating atom becomes ``te::just(value)``,
* a call ``(f a b)`` evaluates the head and every argument into senders,
  joins them with ``te::when_all`` (the arguments through ``tyvi::actions::map``, which is
  itself a ``when_all`` over the list),
  and then uses ``te::let_value`` to call the procedure on the results,
* the procedure returns a sender, and that sender becomes the result of the call.

So a whole program such as

.. code:: text

   (sequence (quote (comm_external current_context emf_E))
             (quote (comm_local current_context emf_E)))

is turned into one sender tree that mirrors the program's structure.
Nothing runs while this tree is built.
``pyactions.c++`` then calls ``tyvi::this_thread::sync_wait`` on it exactly once.
While it waits, the C++ side is free to run the MPI requests,
the thread pool tasks and the GPU kernels described by the tree, in whatever order the
senders allow.
Python sees one blocking call that returns when all of it is done,
or raises an exception if any part failed.

This is the reason the interpreter is written in terms of senders.
It lets the *structure* of the work come from Python
while the *execution* of the work stays entirely inside C++.


Arguments are evaluated concurrently
------------------------------------

The arguments of a call are joined with ``when_all``,
so there is no guaranteed order between them.
They may even run at the same time.
``(f (g ctx) (h ctx))`` does *not* mean "``g`` then ``h``".
If the order matters, or if both touch the same tile data,
do not put them next to each other as arguments.


One ``eval`` call is one synchronisation point
----------------------------------------------

Everything inside one ``eval`` can overlap.
Nothing overlaps *across* two ``eval`` calls, because Python waits for the first one to
finish before it can build the second one.
If two pieces of work need to overlap, they have to be in the same program.
This is also a reason to prefer one larger program over many small ``eval`` calls.


Python does not run while C++ evaluates
---------------------------------------

The thread that called ``eval`` blocks inside ``sync_wait``,
and runko does not release the GIL there.
This is exactly what we want: there are no Python threads and no Python code
running concurrently with the simulation.

It has one consequence.
Python callables passed as atoms, such as the field functions given to ``set_EBJ``,
are called from C++ while the calling thread is still waiting.
This works because those procedures run their loop in a plain ``te::then``,
on the thread that is already waiting in ``sync_wait``.
Code that calls back into Python must not be moved to pika's thread pool
with ``continues_on``, because a worker thread cannot acquire the GIL that the
waiting thread holds.


Types must match exactly
------------------------

``parse_atom_args<T...>`` uses ``atom_cast``, which compares types exactly.
A Python ``int`` arrives as ``long`` and a Python ``float`` as ``double``.
A procedure that asks for ``int`` or ``float`` will reject them with an
"atom at argument n does not hold object of type ..." error.


Errors become Python exceptions
-------------------------------

Parsing errors, unknown symbols, calling something that is not a procedure,
and any exception thrown inside a procedure's sender travel down the sender's error channel.
``sync_wait`` rethrows them, and ``pyactions.c++`` wraps them in a ``std::runtime_error``
that pybind11 raises in Python.
A failing action therefore shows up as an ordinary Python traceback at the ``eval`` call.


Adding a new action
===================

A new action touches the same few places every time:

#. add a value to ``enum class symbol`` in ``src/runko/actions/env.h``,
#. give it a Python name in the ``py::enum_<runko::symbol>`` block of
   ``src/runko/bindings/pyactions.c++``,
#. bind the symbol to a ``procedure`` in ``build_std_env`` or ``build_sim_env``.
   The procedure parses its arguments and returns a sender, usually one that finishes with
   ``ta::null`` when there is nothing meaningful to return,
#. call it from Python as ``(actions.your_action, actions.current_context, ...)``.

The interpreter itself does not need to change.

Further reading
===============

* In Runko: ``external/tyvi/src/tyvi/actions_eval.h`` (the interpreter),
  ``external/tyvi/src/tyvi/actions_ast.h`` (``atom``, ``cons`` and ``sexpr``),
  ``src/runko/actions/env.c++`` (small procedures, including ``sequence``),
  ``src/runko/actions/sim_env.c++`` (the simulation environment) and
  ``src/runko/bindings/pyactions.c++`` (Python conversion and ``sync_wait``).
* :ref:`senders-and-receivers`, for the execution model the interpreter is built on.
* `A Scheme Primer <https://files.spritely.institute/papers/scheme-primer.html>`_,
  a Scheme (Lisp dialect) tutorial whose chapter 12 was the inspiration for ``tyvi::actions::eval``.


Extra: But why choose Lisp as the DSL between Python and C++?
=============================================================

Something has to define the semantics between Python and C++.
One option is to invent a new domain-specific language (DSL).
However, this was deemed unnecessary when one could just use Lisp
instead of designing a syntax, a data model and an evaluator,
because those questions were settled decades ago.

As `Greenspun's tenth rule <https://en.wikipedia.org/wiki/Greenspun%27s_tenth_rule>`_
states:

   Any sufficiently complicated C or Fortran program contains an ad hoc,
   informally-specified, bug-ridden, slow implementation of half of Common Lisp.

By intentionally adding Lisp to Runko, we make sure that it stays small,
robust, and extensible. The design is inspired by
`GNU Emacs <https://www.gnu.org/software/emacs/>`_ which uses
`elisp <https://en.wikipedia.org/wiki/Emacs_Lisp>`_
and `GNU Guix <https://guix.gnu.org/>`_ which uses `Guile <https://www.gnu.org/software/guile/>`_.


There is no parser, because code is already data
------------------------------------------------

An s-expression is both the source code and the syntax tree.
Python tuples are already nested lists, so ``parse_element``
is a short structural conversion, not a grammar.
A newly invented DSL would need one of two things:

- A text syntax. This requires a lexer, a parser, error recovery and error locations,
  and Python would have to build strings, with the quoting and injection problems that brings.
- A custom abstract syntax tree (AST) for representing programs.
  This requires node classes on both the Python and the C++ side, and a binding for each node type.

With s-expressions, the Python side's "AST builder" is the tuple literal.


Very little is built in
-----------------------

The language has three data forms (null, atom, cons) and a handful of eval rules.
Adding an action means binding a symbol in the environment, and the interpreter never changes.
Invented DSLs tend to grow a new syntax rule for each feature ("comm blocks", "io statements", …),
and each new rule touches the parser, the AST and the evaluator.

"It's a tiny Lisp" is a sentence that carries a lot of meaning,
and there are decades of material about it. A bespoke DSL needs its own documentation
from scratch. Runko's version is small enough that the whole language fits on one page.

Semantics that are already known
--------------------------------

Many hard questions already have well-understood Lisp answers:

- What happens to arguments before a call? They are evaluated first.
- How to delay evaluation? With quote.
- How are names resolved? By looking them up in an association list (the environment).
- How to add control flow without changing the core?
  Through special forms, which here are ``procedure_with_eval`` functions like ``sequence``.

An invented language usually rediscovers these, often inconsistently.


It maps naturally onto senders
------------------------------

Lisp evaluation is a recursive tree walk: evaluate the head and the arguments, then apply.
This matches sender composition almost exactly:

- evaluating the children becomes ``when_all``,
- applying the procedure becomes ``let_value``,
- a constant becomes ``just``.

That is why the interpreter is so short.
It also means the concurrency meaning of a program is clear from its tree shape:
siblings may overlap, and nesting implies a dependency.


Other programs can build and change it
--------------------------------------

Because programs are plain data, Python can build, store, reuse and change them with ordinary code.
The same applies to C++: a ``procedure_with_eval`` receives unevaluated code
as a list it can inspect, rewrite, repeat or run in another order.
This is Lisp's macro idea, and it comes for free.
