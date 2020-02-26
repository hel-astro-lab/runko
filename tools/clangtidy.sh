
# all headers
# python3 tools/run-clang-tidy.py -p=build/ \
#	-header-filter=.* \
#	-checks=-*,modernize-* 

# pick only runko headers
#python3 tools/run-clang-tidy.py -p=build/ \
#	-header-filter=/home/natj/runko \
#	-fix \
#	-export-fixes=cfixes.yaml \
#	-target-filter=/home/natj/runko \
#	-checks=-*,\
#modernize-use-using,\

python3 tools/run-clang-tidy.py -p=build/ \
	-export-fixes=fixits \
	-target-filter='pyrunko|pyprtcls|pygol|pycorgitest' \
	-exclude-files='Eigen|optional-lite|pybind|cppitertools' \
	-header-filter=/home/natj/runko \
	-checks=-*,\
performance-* \
| tee ctidy.hist


#clang analyzer:
#hilbert 496 localh

#modernize-use-default-member-init
#modernize-use-equals-default
#modernize-use-trailing-return-type,\

#modernize-use-using,\
# -checks=-*,modernize-*,-modernize-use-trailing-return-type,-modernize-avoid-c-arrays

#checks not to be run
# modernize-loop-convert
# modernize-pass-by-value
# modernize-use-equals-default
# modernize-use-default-member-init
# modernize-use-using

# checks already run
#modernize-avoid-bind,\
#modernize-concat-nested-namespaces,\
#modernize-deprecated-headers,\
#modernize-deprecated-ios-base-aliases,\
#modernize-make-shared,\
#modernize-make-unique,\
#modernize-raw-string-literal,\
#modernize-redundant-void-arg,\
#modernize-replace-auto-ptr,\
#modernize-replace-random-shuffle,\
#modernize-return-braced-init-list,\
#modernize-shrink-to-fit,\
#modernize-unary-static-assert,\
#modernize-use-auto,\
#modernize-use-default-member-init,\
#modernize-use-bool-literals,\
#modernize-use-emplace,\
#modernize-use-equals-default,\
#modernize-use-equals-delete,\
#modernize-use-nodiscard,\
#modernize-use-noexcept,\
#modernize-use-nullptr,\
#modernize-use-override,\
#modernize-use-transparent-functors,\


