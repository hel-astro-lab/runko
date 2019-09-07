# Debugging with `LLDB`

After some tweaking, it is possible to debug also a python library directly.
This can be done with the following tricks:

- make sure the library is compiled with `-g`

- launch `python2`
- see what is the PID of the process: `ps aux | grep python`
- launch `lldb`
- inside it use `attach --pid XXX`
- follow this by `continue` to release the process.
- Voila! Import the library and start debugging normally



# Debugging MPI programs

install tmux
```
brew install tmux
```

Install tmpi (https://github.com/Azrael3000/tmp://github.com/Azrael3000/tmpi) and put it somewhere along the common search path so it can be executed from bash.
Note that in MacOS you have to modify the script where it says mktemp and change --dry-run to -u

Need to be inside tmux to run these, run:
```
tmux
```

run `tmpi` as
```
tmpi 4 lldb -- python3 pic.py 
```
as an example.

This should launch 4 windows and the keyboard should be multiplexed to each simultaneously.






