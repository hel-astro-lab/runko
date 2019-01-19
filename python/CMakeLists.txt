project (pyplasmabox)


#find_package(PYBIND11 REQUIRED) #loaded externally from corgi
#INCLUDE_DIRECTORIES ( "$ENV{pybind11_INCLUDE_DIR}" )

set (TOOLS_FILES 
     pytools.c++
     )

set (FIELDS_FILES 
     pyfields.c++
     ../em-fields/tile.c++ 
     ../em-fields/damping_tile.c++
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
     ../pic/solvers/boris.c++
     ../pic/interpolators/linear.c++
     ../pic/depositers/zigzag.c++
     )

set (RAD_FILES 
     pyrad.c++
     )

set (IO_FILES 
     ../io/writer.c++
     ../io/reader.c++
    )


# python module creation
#--------------------------------------------------
pybind11_add_module(pyplasmabox 
                    pyplasmabox.c++ 
                    ${TOOLS_FILES}
                    ${FIELDS_FILES}
                    ${VLV_FILES}
                    ${PIC_FILES}
                    ${RAD_FILES}
                    #${IO_FILES}
                    )


target_link_libraries(pyplasmabox PRIVATE -lhdf5)
#target_link_libraries(pyplasmabox PRIVATE -lfftw3)
#target_link_libraries(pyplasmabox PRIVATE -lfftw3f)

target_compile_options(pyplasmabox PRIVATE ${WARNING_FLAGS})




# installation
install (TARGETS pyplasmabox DESTINATION lib)
