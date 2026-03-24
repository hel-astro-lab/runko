import contextlib
import unittest
import io
import contextlib
import runko

# Source - https://stackoverflow.com/a/22434262
# Posted by jfs, modified by community. See post 'Timeline' for change history
# Retrieved 2026-03-23, License - CC BY-SA 4.0
try:
    import ctypes
    from ctypes.util import find_library
except ImportError:
    libc = None
else:
    try:
        libc = ctypes.cdll.msvcrt # Windows
    except OSError:
        libc = ctypes.cdll.LoadLibrary(find_library('c'))

def flush(stream):
    try:
        libc.fflush(None)
        stream.flush()
    except (AttributeError, ValueError, IOError):
        pass # unsupported


# Source - https://stackoverflow.com/a/22434262
# Posted by jfs, modified by community. See post 'Timeline' for change history
# Retrieved 2026-03-22, License - CC BY-SA 4.0
import os
import sys
from contextlib import contextmanager

def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            flush(stdout_fd)
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

ra = runko.actions

class actions_eval(unittest.TestCase):

    def test_print(self):
        fd = os.memfd_create("test.txt")

        with stdout_redirected(fd) as mem:
            strings = (ra.quote, ("foo", "bar"))
            body = (ra.print, (ra.car, strings))
            ra.print_context_eval(body)

        correct = "foo"
        s = os.pread(fd, len(correct) + 1, 0).decode("utf-8")
        os.close(fd)
        self.assertEqual(s, correct)


    def test_println(self):
        fd = os.memfd_create("test.txt")

        with stdout_redirected(fd) as mem:
            strings = (ra.quote, ("foo", "bar"))
            body = (ra.println, (ra.car, (ra.cdr, strings)))
            ra.print_context_eval(body)

        correct = "bar\n"
        s = os.pread(fd, len(correct) + 1, 0).decode("utf-8")
        os.close(fd)
        self.assertEqual(s, correct)


    def test_print_version(self):
        fd = os.memfd_create("test.txt")

        with stdout_redirected(fd) as mem:
            body = (ra.print, ra.version)
            ra.print_context_eval(body)

        correct = "runko v6.x"
        s = os.pread(fd, len(correct) + 1, 0).decode("utf-8")
        os.close(fd)
        self.assertEqual(s, correct)


if __name__ == "__main__":
    unittest.main()
