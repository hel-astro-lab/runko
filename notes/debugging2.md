# Debugging mpi+python processes with GDB



run 

```
mpirun -n <NP> xterm -hold -e gdb -ex run --args ./program [arg1] [arg2] [...]
```



Yes, it's just `--` instead of `--args`. From the help:

```
lldb -v [[--] <PROGRAM-ARG-1> [<PROGRAM_ARG-2> ...]]
```

Thus:

```
$ lldb -- exe --lots --of --flags -a -b -c d e
```



```
(lldb) br s -n main
(lldb) r
(lldb) pr h -p true -n true -s true SIGSEGV
NAME        PASS   STOP   NOTIFY
==========  =====  =====  ======
SIGSEGV     true   true   true 
(lldb) c
```



set a breakpoint into malloc error:

```
br set --name malloc_error_break
```



## Case of segfault

`lldb -- python2 ~/projects/plasma/projects/rad-reconn/pic.pyf`

- run

- bt

- frame select <X>

- frame var [XX]

  ​