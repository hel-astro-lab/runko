Intro
#####

Runko module is a python module provided by the Runko framework,
which is completely new API in order to facilitate GPU based computing.
It is separate from previous `pyrunko` module and should not be mixed with it.
Later references to runko are about `runko` module if not stated otherwise.


Motivation
==========

Previous Runko API exposed some of the implementation details to it's public API.
For example, emf tile initial conditions are set by directly modifying the underlying buffers.
In order to do GPU computations on these buffers there is two options.

A) Use unified GPU memory and let data transfers be implicit.
B) Add new functions for explicit memory transfers.

We did not want to do A) because it would limit us to unified memory
and we did not want to do B) either, because it would complicate the already complex API.
So we ended up redesinging the API from ground up such that it would be agnostic of
the underlying implementation. This way we could also simplify many aspects of the API at the same time.


Status
======

Runko module is at MVP stage.
This means that a simple emf and pic simulation is possible, but many of the features
present in the previous API have not been implemented.
