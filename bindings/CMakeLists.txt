set (TOOLS_FILES
     pytools.c++
     ../tools/config_parser.c++
     ../tools/system.c++
     )

set (EMF_FILES
     pyemf.c++
     ../core/emf/tile.c++
     ../core/emf/virtual_tile.c++
     ../core/emf/yee_lattice.c++
     ../core/emf/yee_lattice_FDTD2.c++
     ../core/emf/yee_lattice_stencil.c++
     ../core/emf/yee_lattice_current_filter_binomial2.c++
     )

set (PIC_FILES
     pypic.c++
     ../core/pic/tile.c++
     ../core/pic/tile_communication.c++
     ../core/pic/virtual_tile.c++
     ../core/pic/particle.c++
     ../core/pic/particle_current_zigzag_1st.c++
     ../core/pic/reflector_wall.c++
     )

set (IO_FILES
     ../io/snapshots/hdf5_fields.c++
     ../io/snapshots/mpiio_fields.c++
     ../io/snapshots/mpiio_particles.c++
     ../io/snapshots/mpiio_spectra.c++
    )


pybind11_add_module(runko_cpp_bindings
    runko_cpp_bindings.c++
    ${TOOLS_FILES}
    ${EMF_FILES}
    ${PIC_FILES}
    ${IO_FILES}
)

target_link_libraries(runko_cpp_bindings
    PRIVATE
    runko_includes
    runko_warnings
    corgi
    tyvi::tyvi
    hdf5::hdf5
)

# Cray MPICH GPU Transport Layer: link explicitly so LD_PRELOAD is no
# longer needed at runtime. CRAY_MPICH_ROOTDIR is set by the cray-mpich
# module; module must be loaded before configure on LUMI.
if(DEFINED ENV{CRAY_MPICH_ROOTDIR})
    find_library(MPI_GTL_HSA_LIB
        NAMES mpi_gtl_hsa
        PATHS "$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib"
        NO_DEFAULT_PATH
        REQUIRED
    )
    target_link_libraries(runko_cpp_bindings PRIVATE "${MPI_GTL_HSA_LIB}")
    set_property(TARGET runko_cpp_bindings APPEND PROPERTY
        INSTALL_RPATH "$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib")
endif()

if(RUNKO_CODE_COVERAGE)
    target_link_libraries(runko_cpp_bindings PUBLIC runko_coverage_config)
endif()

set_target_properties(runko_cpp_bindings PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}"
)
