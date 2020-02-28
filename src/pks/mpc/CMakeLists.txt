# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories

add_subdirectory(constitutive_relations)
#add_subdirectory(morphology_pk)

include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${ATS_SOURCE_DIR}/pks)
include_directories(${ATS_SOURCE_DIR}/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/pks/flow)
include_directories(${ATS_SOURCE_DIR}/pks/energy)
include_directories(${ATS_SOURCE_DIR}/pks/transport)
include_directories(${ATS_SOURCE_DIR}/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/operators/advection)
include_directories(${ATS_SOURCE_DIR}/pks/mpc/reactive_transport)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/constitutive_relations)

set(ats_mpc_src_files
  weak_mpc.cc
  operator_split_mpc.cc
  weak_mpc_semi_coupled.cc
  weak_mpc_semi_coupled_deform.cc
  mpc_coupled_cells.cc
  pk_mpcsubcycled_ats.cc
  mpc_surface_subsurface_helpers.cc
  weak_mpc_semi_coupled_helper.cc
  mpc_weak_subgrid.cc
  mpc_delegate_ewc.cc
  mpc_reactivetransport_pk.cc
  #mpc_flowreactivetransport_pk.cc
  mpc_delegate_ewc_subsurface.cc
  mpc_delegate_water.cc
  mpc_coupled_water.cc
  mpc_coupled_water_split_flux.cc
  mpc_coupled_transport.cc
  mpc_permafrost_split_flux.cc
  mpc_permafrost_split_flux_columns.cc
  mpc_permafrost_split_flux_columns_subcycled.cc
  mpc_subsurface.cc
  mpc_surface.cc
  mpc_permafrost.cc
  mpc_morphology_pk.cc
  biomass_evaluator.cc
  )

set(ats_mpc_inc_files
  mpc.hh
  weak_mpc.hh
  strong_mpc.hh
  operator_split_mpc.hh
  weak_mpc_semi_coupled.hh
  weak_mpc_semi_coupled_deform.hh
  mpc_coupled_cells.hh
  pk_mpcsubcycled_ats.hh
  mpc_surface_subsurface_helpers.hh
  weak_mpc_semi_coupled_helper.hh
  mpc_weak_subgrid.hh
  mpc_delegate_ewc.hh
  mpc_reactivetransport_pk.hh
  #mpc_flowreactivetransport_pk.hh
  mpc_delegate_ewc_subsurface.hh
  mpc_delegate_water.hh
  mpc_coupled_water.hh
  mpc_coupled_transport.hh
  mpc_coupled_water_split_flux.hh
  mpc_permafrost_split_flux.hh
  mpc_permafrost_split_flux_columns.hh
  mpc_permafrost_split_flux_columns_subcycled.hh
  mpc_subsurface.hh
  mpc_surface.hh
  mpc_permafrost.hh
  biomass_evaluator.hh
  )

set(ats_mpc_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  operators
  solvers
  time_integration
  state
  pks
  chemistry_pk
  ats_operators
  ats_eos
  ats_generic_evals
  ats_pks
  ats_transport
  ats_flow
  ats_mpc_relations
  )

add_amanzi_library(ats_mpc
                   SOURCE ${ats_mpc_src_files}
                   HEADERS ${ats_mpc_inc_files}
		   LINK_LIBS ${ats_mpc_link_libs})

# register factories
register_evaluator_with_factory(
  HEADERFILE weak_mpc_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE strong_mpc_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_predictor_corrector_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE operator_split_mpc_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE pk_mpcsubcycled_ats_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE reactive_transport/mpc_flowreactivetransport_pk_reg.hh
  LISTNAME   MPC_MPC_FAC
  )

register_evaluator_with_factory(
  HEADERFILE reactive_transport/mpc_reactivetransport_pk_reg.hh
  LISTNAME   MPC_MPC_FAC
  )

register_evaluator_with_factory(
  HEADERFILE reactive_transport/mpc_coupled_reactive_transport_pk_reg.hh
  LISTNAME   MPC_MPC_FAC
  )


register_evaluator_with_factory(
   HEADERFILE mpc_morphology_reg.hh
   LISTNAME   ATS_MPC_REG
   )
    
register_evaluator_with_factory(
  HEADERFILE mpc_subsurface_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_surface_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_cells_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_split_flux_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_split_flux_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_split_flux_columns_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_split_flux_columns_subcycled_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
    HEADERFILE mpc_coupled_transport_reg.hh
    LISTNAME   ATS_MPC_REG
   )


register_evaluator_with_factory(
  HEADERFILE weak_mpc_semi_coupled_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_weak_subgrid_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE weak_mpc_semi_coupled_deform_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE biomass_evaluator_reg.hh
  LISTNAME ATS_MPC_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_mpc_registration.hh
  LISTNAME   ATS_MPC_REG
  INSTALL    True
)


include_evaluators_directories(LISTNAME ATS_POROSITY_EVAL_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_EOS_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_GENERIC_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_CONSTITUTIVE_RELATIONS_POROSITY_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME WRM_MODELS_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME WMR_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME RICHARDS_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME PERMAFROST_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_INTERNAL_ENERGY_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_SOURCE_TERMS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_TWO_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_THREE_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME MPC_PK_FAC_INCLUDES)


#add_subdirectory(coupled_transport)
add_subdirectory(reactive_transport)

#if ( BUILD_TESTS )
if (0)
  # Add UnitTest includes
  set(amanzi_libs atk geometry mesh mesh_audit simple_mesh state whetstone mesh_functions)

  include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

  include_directories(${ATS_SOURCE_DIR}/src/pks/energy/advection_diffusion)

  # Test: miscaleneous transport routines
  add_executable(test_mpc_subsurface test/Main.cc test/test_mpc_subsurface.cc)
#  add_executable(test_mpc_subsurface test/Main.cc )
  target_link_libraries(test_mpc_subsurface 
    mpc pk_energy_two_phase 
    pk_energy_three_phase 
    pk_energy_base
    pk_flow_permafrost pk_flow_richards
    pk_transport
    pk_coupled_transport
    pk_bases
    energy_relations_thermal_conductivity energy_relations_internal_energy
    energy_relations_enthalpy energy_relations_source_terms
    flow_relations
    generic_evaluators relations_eos relations_ewc bc_factory divgrad
    amanzi_solvers amanzi_time_integration amanzi_operators amanzi_pks
    amanzi_state amanzi_whetstone amanzi_mesh amanzi_geometry 
    amanzi_mesh_functions amanzi_mesh_factory
    amanzi_output amanzi_data_structures amanzi_atk
    ${Amanzi_TPL_UnitTest_LIBRARIES})

endif()
