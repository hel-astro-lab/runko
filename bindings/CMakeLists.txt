project (pyrunko)


#find_package(PYBIND11 REQUIRED) #loaded externally from corgi
#INCLUDE_DIRECTORIES ( "$ENV{pybind11_INCLUDE_DIR}" )

set (TOOLS_FILES 
     pytools.c++
     ../tools/hilbert.c++
     )

set (FIELDS_FILES 
     pyfields.c++
     ../em-fields/tile.c++ 
     ../em-fields/damping_tile.c++
     ../em-fields/propagator/fdtd2.c++ 
     ../em-fields/propagator/fdtd4.c++ 
     ../em-fields/filters/digital.c++
     )

set (VLV_FILES 
     pyvlv.c++
     ../vlasov/tile.c++
     ../vlasov/boundaries/piston.c++
     ../vlasov/boundaries/outflow.c++
     ../vlasov/momentum-solvers/amr_momentum_solver.c++
     ../vlasov/momentum-solvers/bwd_lagrangian.c++
     )


set (PIC_FILES 
     pypic.c++
     ../pic/tile.c++
     ../pic/particle.c++
     ../pic/boundaries/wall.c++
     ../pic/boundaries/piston.c++
     ../pic/pushers/boris.c++
     ../pic/pushers/boris_drag.c++
     ../pic/pushers/vay.c++
     ../pic/interpolators/linear.c++
     ../pic/depositers/zigzag.c++
     )

set (RAD_FILES 
     pyrad.c++
     )

set (FFE_FILES 
     pyffe.c++
     ../ffe/tile.c++
     ../ffe/currents/driftcurrent.c++
     )

set (IO_FILES 
     ../io/writer.c++
     ../io/reader.c++
     ../io/quick_writer.c++
     ../io/test_prtcl_writer.c++
    )



# python module creation
#--------------------------------------------------
pybind11_add_module(pyrunko 
                    runko.c++ 
                    ${TOOLS_FILES}
                    ${FIELDS_FILES}
                    ${VLV_FILES}
                    ${PIC_FILES}
                    ${RAD_FILES}
                    ${FFE_FILES}
                    ${IO_FILES}
                    )

#target_link_libraries(pyrunko PRIVATE -lhdf5)
target_link_libraries(pyrunko PRIVATE ${HDF5_C_LIBRARIES})
target_include_directories(pyrunko PRIVATE ${HDF5_INCLUDE_DIRS})

#target_link_libraries(pyrunko PRIVATE -lfftw3)
#target_link_libraries(pyrunko PRIVATE -lfftw3f)
target_link_libraries(pyrunko PUBLIC coverage_config)

target_compile_options(pyrunko PRIVATE ${WARNING_FLAGS})




# installation
install (TARGETS pyrunko DESTINATION lib)
