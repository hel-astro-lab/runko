# FindRunkoHDF5.cmake
# Finds HDF5 with fallback to pkg-config, and creates hdf5::hdf5 imported target.

find_package(HDF5 QUIET)

if(NOT HDF5_FOUND OR NOT HDF5_LIBRARIES)
    find_package(PkgConfig QUIET)
    pkg_check_modules(HDF5_PKG QUIET hdf5)
    if(HDF5_PKG_FOUND)
        find_library(HDF5_LIBRARY NAMES hdf5 PATHS ${HDF5_PKG_LIBRARY_DIRS} NO_DEFAULT_PATH)
        find_library(HDF5_HL_LIBRARY NAMES hdf5_hl PATHS ${HDF5_PKG_LIBRARY_DIRS} NO_DEFAULT_PATH)
        if(HDF5_LIBRARY)
            set(HDF5_LIBRARIES ${HDF5_LIBRARY})
            list(APPEND HDF5_LIBRARIES ${HDF5_HL_LIBRARY})
            set(HDF5_INCLUDE_DIRS ${HDF5_PKG_INCLUDE_DIRS})
            set(HDF5_FOUND TRUE)
        endif()
    endif()
endif()

if(NOT HDF5_FOUND)
    message(FATAL_ERROR "Could NOT find HDF5. Please check that HDF5_ROOT is set correctly or that the cray-hdf5 module is loaded.")
endif()

if(HDF5_FOUND AND NOT TARGET hdf5::hdf5)
    list(GET HDF5_LIBRARIES 0 HDF5_MAIN_LIB)
    add_library(hdf5::hdf5 UNKNOWN IMPORTED)
    set_target_properties(hdf5::hdf5 PROPERTIES
        IMPORTED_LOCATION "${HDF5_MAIN_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
    )
    set_target_properties(hdf5::hdf5 PROPERTIES INTERFACE_LINK_LIBRARIES "${HDF5_HL_LIBRARY}")
endif()
