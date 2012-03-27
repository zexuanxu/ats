/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "air_water_rock.hh"

namespace Amanzi {
namespace Energy {

void AirWaterRock::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // get needed variables
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> mol_frac_gas = S->GetFieldData("mol_frac_gas");

  // and the secondary variables to be calculated
  Teuchos::RCP<CompositeVector> int_energy_gas =
    S->GetFieldData("internal_energy_gas", "energy");
  Teuchos::RCP<CompositeVector> int_energy_liquid =
    S->GetFieldData("internal_energy_liquid", "energy");
  Teuchos::RCP<CompositeVector> int_energy_rock =
    S->GetFieldData("internal_energy_rock", "energy");

  // update secondary variables
  InternalEnergyGas_(*temp, *mol_frac_gas, int_energy_gas);
  InternalEnergyLiquid_(*temp, int_energy_liquid);
  InternalEnergyRock_(*temp, int_energy_rock);

};

void AirWaterRock::UpdateSpecificEnthalpyLiquid_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");

  Teuchos::RCP<const CompositeVector> dens_liq;
  if (internal_energy_liquid_model_->IsMolarBasis()) {
    dens_liq = S->GetFieldData("molar_density_liquid");
  } else {
    dens_liq = S->GetFieldData("density_liquid");
  }

  Teuchos::RCP<const CompositeVector> int_energy_liquid =
    S->GetFieldData("internal_energy_liquid");

  Teuchos::RCP<CompositeVector> spec_enthalpy_liq =
    S->GetFieldData("specific_enthalpy_liquid", "energy");

  // update enthalpy of liquid
  SpecificEnthalpyLiquid_(*int_energy_liquid, *pres, *dens_liq, spec_enthalpy_liq);
};

void AirWaterRock::UpdateThermalConductivity_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> poro =
    S->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> sat_liq =
    S->GetFieldData("saturation_liquid");
  Teuchos::RCP<CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity", "energy");

  ThermalConductivity_(*poro, *sat_liq, thermal_conductivity);
};

void AirWaterRock::AddAccumulation_(Teuchos::RCP<CompositeVector> f) {
  Teuchos::RCP<const CompositeVector> poro0 =
    S_inter_->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> poro1 =
    S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> density_gas0;
  Teuchos::RCP<const CompositeVector> density_gas1;
  if (internal_energy_gas_model_->IsMolarBasis()) {
    density_gas0 = S_inter_->GetFieldData("molar_density_gas");
    density_gas1 = S_next_->GetFieldData("molar_density_gas");
  } else {
    density_gas0 = S_inter_->GetFieldData("density_gas");
    density_gas1 = S_next_->GetFieldData("density_gas");
  }

  Teuchos::RCP<const CompositeVector> density_liq0;
  Teuchos::RCP<const CompositeVector> density_liq1;
  if (internal_energy_liquid_model_->IsMolarBasis()) {
    density_liq0 = S_inter_->GetFieldData("molar_density_liquid");
    density_liq1 = S_next_->GetFieldData("molar_density_liquid");
  } else {
    density_liq0 = S_inter_->GetFieldData("density_liquid");
    density_liq1 = S_next_->GetFieldData("density_liquid");
  }

  Teuchos::RCP<const CompositeVector> sat_liq0 =
    S_inter_->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> sat_liq1 =
    S_next_->GetFieldData("saturation_liquid");

  Teuchos::RCP<const CompositeVector> sat_gas0 =
    S_inter_->GetFieldData("saturation_gas");
  Teuchos::RCP<const CompositeVector> sat_gas1 =
    S_next_->GetFieldData("saturation_gas");

  Teuchos::RCP<const CompositeVector> int_energy_gas0 =
    S_inter_->GetFieldData("internal_energy_gas");
  Teuchos::RCP<const CompositeVector> int_energy_gas1 =
    S_next_->GetFieldData("internal_energy_gas");

  Teuchos::RCP<const CompositeVector> int_energy_liq0 =
    S_inter_->GetFieldData("internal_energy_liquid");
  Teuchos::RCP<const CompositeVector> int_energy_liq1 =
    S_next_->GetFieldData("internal_energy_liquid");

  Teuchos::RCP<const CompositeVector> int_energy_rock0 =
    S_inter_->GetFieldData("internal_energy_rock");
  Teuchos::RCP<const CompositeVector> int_energy_rock1 =
    S_next_->GetFieldData("internal_energy_rock");

  Teuchos::RCP<const CompositeVector> cell_volume0 =
    S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cell_volume1 =
    S_next_->GetFieldData("cell_volume");

  Teuchos::RCP<const double> density_rock =
    S_next_->GetScalarData("density_rock");

  double dt = S_next_->time() - S_inter_->time();

  // NOTE: gas and liquid are done in a ?? basis, but rock is done in a mass basis

  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c != c_owned; ++c) {
    // calculte the energy density at the old and new times
    double edens_liq1 = (*density_liq1)(c) * (*sat_liq1)(c) *
      (*int_energy_liq1)(c);
    double edens_gas1 = (*density_gas1)(c) * (*sat_gas1)(c) *
      (*int_energy_gas1)(c);
    double edens_rock1 = (*density_rock) * (*int_energy_rock1)(c);
    double energy1 = ((*poro1)(c) * (edens_gas1 + edens_liq1) +
                      (1-(*poro1)(c)) * (edens_rock1)) * (*cell_volume1)(c);

    double edens_liq0 = (*density_liq0)(c) * (*sat_liq0)(c) *
      (*int_energy_liq0)(c);
    double edens_gas0 = (*density_gas0)(c) * (*sat_gas0)(c) *
      (*int_energy_gas0)(c);
    double edens_rock0 = (*density_rock) * (*int_energy_rock0)(c);
    double energy0 = ((*poro0)(c) * (edens_gas0 + edens_liq0) +
                      (1-(*poro0)(c)) * (edens_rock0)) * (*cell_volume0)(c);

    // add the time derivative of energy density to the residual
    (*f)("cell",0,c) += (energy1 - energy0)/dt;
  }
};

void AirWaterRock::AddAdvection_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> f, bool negate) {
  advection_->set_flux(S->GetFieldData("darcy_flux"));
  Teuchos::RCP<CompositeVector> field = advection_->field();

  // stuff density_liquid * enthalpy_liquid into the field cells
  Teuchos::RCP<const CompositeVector> density_liq =
    S->GetFieldData("density_liquid");

  UpdateSpecificEnthalpyLiquid_(S);
  Teuchos::RCP<const CompositeVector> enthalpy_liq =
    S->GetFieldData("specific_enthalpy_liquid");

  field->ViewComponent("cell")->PutScalar(0);
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    (*field)("cell",0,c) = (*density_liq)("cell",0,c)*(*enthalpy_liq)("cell",0,c);
  }

  // apply the advection operator and add to residual
  advection_->Apply();
  if (negate) {
    for (int c=0; c!=c_owned; ++c) {
      (*f)("cell",c) -= (*field)("cell",c);
    }
  } else {
    for (int c=0; c!=c_owned; ++c) {
      (*f)("cell",c) = (*field)("cell",c);
    }
  }
};

void AirWaterRock::ApplyConduction_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> f) {
  // compute the stiffness matrix at the new time
  Teuchos::RCP<const CompositeVector> temp =
    S->GetFieldData("temperature");

  // get conductivity, and push it into whetstone tensor
  UpdateThermalConductivity_(S);
  Teuchos::RCP<CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity", "energy");

  for (int c=0; c != Ke_.size(); ++c) {
    Ke_[c](0,0) = (*thermal_conductivity)("cell", c);
  }

  // calculate the div-grad operator, apply it to temperature, and add to residual
  matrix_->CreateMFDstiffnessMatrices(Ke_, *thermal_conductivity);
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();
  matrix_->ComputeNegativeResidual(*temp, f);
};

void AirWaterRock::InternalEnergyGas_(const CompositeVector& temp,
        const CompositeVector& mol_frac_gas,
        const Teuchos::RCP<CompositeVector>& int_energy_gas) {
  // just a single model for now -- ignore blocks
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c != c_owned; ++c) {
    (*int_energy_gas)("cell",0,c) = internal_energy_gas_model_->
      InternalEnergy(temp("cell",0,c), mol_frac_gas("cell",0,c));
  }
};

void AirWaterRock::InternalEnergyLiquid_(const CompositeVector& temp,
        const Teuchos::RCP<CompositeVector>& int_energy_liquid) {
  // just a single model for now -- ignore blocks
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c != c_owned; ++c) {
    (*int_energy_liquid)("cell",0,c) = internal_energy_liquid_model_->
      InternalEnergy(temp("cell",0,c));
  }
};

void AirWaterRock::InternalEnergyRock_(const CompositeVector& temp,
        const Teuchos::RCP<CompositeVector>& int_energy_rock) {
  // just a single model for now -- ignore blocks
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c != c_owned; ++c) {
    (*int_energy_rock)("cell",0,c) = internal_energy_rock_model_->
      InternalEnergy(temp("cell",0,c));
  }
};

void AirWaterRock::SpecificEnthalpyLiquid_(const CompositeVector& int_energy_liquid,
        const CompositeVector& pres, const CompositeVector& dens_liq,
        const Teuchos::RCP<CompositeVector>& spec_enthalpy_liq) {

  // just a single model for now -- ignore blocks
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c != c_owned; ++c) {
    (*spec_enthalpy_liq)("cell",0,c) = dens_liq("cell",0,c)*int_energy_liquid("cell",0,c)
                                              + pres("cell",0,c);
  }
};

void AirWaterRock::ThermalConductivity_(const CompositeVector& porosity,
        const CompositeVector& sat_liq,
        const Teuchos::RCP<CompositeVector>& thermal_conductivity) {

  // just a single model for now -- ignore blocks
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c != c_owned; ++c) {
    (*thermal_conductivity)("cell",0,c) = thermal_conductivity_model_->
      CalculateConductivity(porosity("cell",0,c), sat_liq("cell",0,c));
  }
};

} //namespace Energy
} //namespace Amanzi
