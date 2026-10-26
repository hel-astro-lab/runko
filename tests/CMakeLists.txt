configure_file("test_wrapper.py.in" "test_wrapper.py" @ONLY)
set(py_test_wrapper "${CMAKE_CURRENT_BINARY_DIR}/test_wrapper.py")

add_subdirectory(py)
add_subdirectory(py-multirank)
add_subdirectory(cpp)
add_subdirectory(cpp-multirank)
