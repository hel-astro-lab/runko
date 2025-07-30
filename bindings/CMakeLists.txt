set (TOOLS_FILES 
     pytools.c++
     ../tools/hilbert.c++
     ../tools/config_parser.c++
     ../tools/system.c++
     )

set (FIELDS_FILES 
     pyemf.c++
     ../core/emf/tile.c++ 
     ../core/emf/propagators/fdtd2.c++ 
     ../core/emf/propagators/fdtd2_pml.c++ 
     ../core/emf/propagators/fdtd4.c++ 
     ../core/emf/propagators/fdtd_general.c++ 
     ../core/emf/filters/binomial2.c++
     ../core/emf/filters/compensator.c++
     ../core/emf/filters/general_binomial.c++
     ../core/emf/filters/strided_binomial.c++
    #../core/emf/filters/sweeping_binomial.c++
     ../core/emf/boundaries/damping_tile.c++
     ../core/emf/boundaries/conductor.c++
     )

set (FIELDS2_FILES
     pyemf2.c++
     ../core/emf2/tile.c++
     ../core/emf2/yee_lattice.c++
     ../core/emf2/yee_lattice_FDTD2.c++
     ../core/emf2/yee_lattice_interpolate_linear_1st.c++
     )

set (VLV_FILES
      pyvlv.c++
      ../core/vlv/tile.c++
      ../core/vlv/boundaries/piston.c++
      ../core/vlv/boundaries/outflow.c++
      ../core/vlv/momentum-solvers/amr_momentum_solver.c++
      ../core/vlv/momentum-solvers/bwd_lagrangian.c++
     )


set (PIC_FILES 
     pypic.c++
     ../core/pic/tile.c++
     ../core/pic/particle.c++
     ../core/pic/boundaries/wall.c++
     ../core/pic/boundaries/piston.c++
     ../core/pic/boundaries/piston_z.c++
     ../core/pic/boundaries/star_surface_injector.c++
     ../core/pic/pushers/boris.c++
     ../core/pic/pushers/boris_drag.c++
     ../core/pic/pushers/boris_rad.c++
     ../core/pic/pushers/boris_grav.c++
     ../core/pic/pushers/vay.c++
     ../core/pic/pushers/higuera_cary.c++
     ../core/pic/pushers/rgca.c++
     ../core/pic/pushers/photon.c++
     ../core/pic/pushers/pulsar.c++
     ../core/pic/interpolators/linear_1st.c++
     ../core/pic/interpolators/quadratic_2nd.c++
     ../core/pic/interpolators/cubic_3rd.c++
     ../core/pic/interpolators/quartic_4th.c++
     ../core/pic/depositers/zigzag.c++
     ../core/pic/depositers/zigzag_2nd.c++
     ../core/pic/depositers/zigzag_3rd.c++
     ../core/pic/depositers/zigzag_4th.c++
     ../core/pic/depositers/esikerpov_2nd.c++
     ../core/pic/depositers/esikerpov_4th.c++
    #../core/pic/depositers/esikerpov_4th_vec.c++
     )

set (PIC2_FILES
     pypic2.c++
     ../core/pic2/tile.c++
     ../core/pic2/tile_communication.c++
     ../core/pic2/particle.c++
     ../core/pic2/particle_boris.c++
     ../core/pic2/particle_current_zigzag_1st.c++
     )

set (QED_FILES 
     pyqed.c++
     ../core/qed/interactions/pair_ann.c++
     ../core/qed/interactions/phot_ann.c++
     ../core/qed/interactions/compton.c++
     ../core/qed/interactions/synchrotron.c++
     ../core/qed/interactions/multi_phot_ann.c++
     ../core/qed/pairing.c++
     ../core/qed/pairing_all2all.c++
     )

set (FFE_FILES 
     pyffe.c++
     ../core/ffe/tile.c++
     ../core/ffe/slim_grids.c++
     ../core/ffe/currents/rffe2.c++
     ../core/ffe/currents/rffe4.c++
     ../core/ffe/currents/ffe2.c++
     ../core/ffe/currents/ffe4.c++
     )



set (IO_FILES 
     ../io/writers/writer.c++
     ../io/readers/reader.c++
     ../io/snapshots/fields.c++
     ../io/snapshots/fields2.c++
     ../io/snapshots/test_prtcls.c++
     ../io/snapshots/pic_moments.c++
     ../io/snapshots/field_slices.c++
     ../io/snapshots/master_only_fields.c++
     ../io/snapshots/master_only_moments.c++
    )

set (ACC_FILES 
    )


set (MISC_FILES 
    )

set (RUNKO_NEXT_FILES
     runko_next.c++
    )

pybind11_add_module(pyrunko
    runko.c++
    ${TOOLS_FILES}
    ${FIELDS_FILES}
    ${FIELDS2_FILES}
    ${VLV_FILES}
    ${PIC_FILES}
    ${PIC2_FILES}
    ${QED_FILES}
    ${FFE_FILES}
    ${IO_FILES}
    ${ACC_FILES}
    ${MISC_FILES}
    ${RUNKO_NEXT_FILES}
)

target_link_libraries(pyrunko
    PRIVATE
    runko_includes
    runko_warnings
    corgi
    tyvi::tyvi
    hdf5::hdf5
    OpenMP::OpenMP_CXX
)

if(RUNKO_CODE_COVERAGE)
    target_link_libraries(pyrunko PUBLIC runko_coverage_config)
endif()

set_target_properties(pyrunko PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
)
