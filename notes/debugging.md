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



