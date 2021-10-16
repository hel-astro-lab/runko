SIMD/GPU vectorization
======================


Many of the Runko algorithms are vectorized to run faster with modern CPU architectures with SIMD support and with GPUs that require multi-core execution.


Vectorization control flags
---------------------------

In some platforms these vectorizations can become a performance bottleneck; typical reasons are either launching of too many micro-kernels (that would be more efficiently calculated by grouping them) or memory-write-blocking because of atomic write operations (that would be performed faster by just using serial operations).

In these cases, we can revert back to non-vectorized operations by un-defining special compile time flags. These are:

- `VEC_FLD2D` use vectorization for 2D electromagnetic tile mesh copy operations.
- `VEC_FLD3D` use vectorization for 3D electromagnetic tile mesh copy operations.
- `VEC_CUR2D` use vectorization for 2D electromagnetic tile mesh addition operations.
- `VEC_CUR3D` use vectorization for 3D electromagnetic tile mesh addition operations.

Typical cause for bad performance are the last two flags that require atomic additions.



.. note::

    TODO: expand this section.





