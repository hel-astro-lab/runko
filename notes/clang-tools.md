# clang-tools

## creating compile_commands.json

https://sarcasm.github.io/notes/dev/compilation-database.html

https://eli.thegreenplace.net/2014/05/21/compilation-databases-for-clang-based-tools

First we need compile_commands.json. 

- from `cmake`

- or using `makefile` and bear

  - https://github.com/rizsotto/Bear
  - does not work in MacOS exepct in special folder: /usr/local, run it there 
  - does not work!

- try compiledb instead: `pip2 install compiledb` and then run make with `compiledb make`

- or by hand:

  - Clang’s [-MJ option](https://clang.llvm.org/docs/ClangCommandLineReference.html#cmdoption-clang-mj) generates a compilation database entry per input (requires `Clang >= 5.0`).

    Usage:

    ```
    clang++ -MJ a.o.json -Wall -std=c++11 -o a.o -c a.cpp
    clang++ -MJ b.o.json -Wall -std=c++11 -o b.o -c b.cpp
    ```

    To merge the compilation database entries into a valid compilation database, it is possible to use `sed`:

    ```
    sed -e '1s/^/[\n/' -e '$s/,$/\n]/' *.o.json > compile_commands.json
    ```

## running clang tools

make sure llvm is on your path by appending your `.bash_profile`: `export PATH="/usr/local/opt/llvm/bin/:$PATH"`


https://clang.llvm.org/extra/clang-rename.html

http://clang.llvm.org/extra/clang-tidy/checks/list.html

`python2 tools/run-clang-tidy.py -header-filter='.*' -checks='-*,modernize-*'`

`python2 tools/run-clang-tidy.py -header-filter='.*' -checks='-*,performance-*`'

`python2 tools/run-clang-tidy.py -header-filter='.*' -checks='-*,readability-*'`



/usr/local/opt/llvm/share/clang/run-clang-tidy.py -header-filter=.* -checks=performance-*


## bug detection

can be used to detect bugs during compile-time.
`scan-build make`


## optimization with Polly

http://polly.llvm.org/docs/





