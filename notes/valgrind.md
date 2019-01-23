# Using valgrind with mpi and python


## in serial mode

```
valgrind --tool=memcheck --suppressions=../../tools/valgrind-python.supp python3 pic.py
```
Important part is the suppression file located in tools.



## with MPI

```
mpirun -n 2 valgrind --leak-check=full --show-reachable=yes --log-file=pic.vg.%p --suppressions=../../tools/valgrind-python.supp python3 pic.py
```
Its good idea to use log-file instead of screen to make sense of the output.



## Refs:
- https://stackoverflow.com/questions/3982036/how-can-i-use-valgrind-with-python-c-extensions
- https://stackoverflow.com/questions/34851643/using-valgrind-to-spot-error-in-mpi-code






