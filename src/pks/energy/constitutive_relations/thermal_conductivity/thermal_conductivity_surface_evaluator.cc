/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "dbc.hh"
#include "thermal_conductivity_surface_factory.hh"
#include "thermal_conductivity_surface_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

ThermalConductivitySurfaceEvaluator::ThermalConductivitySurfaceEvaluator(
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("thermal conductivity key",
          "surface_thermal_conductivity");
  setLinePrefix(my_key_+std::string(" evaluator"));

  uf_key_ = plist_.get<std::string>("unfrozen fraction key", "unfrozen_fraction");
  dependencies_.insert(uf_key_);

  height_key_ = plist_.get<std::string>("height key", "ponded_depth");
  dependencies_.insert(height_key_);

  ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  K_liq_ = sublist.get<double>("thermal conductivity of water");
  K_ice_ = sublist.get<double>("thermal conductivity of ice");
}


ThermalConductivitySurfaceEvaluator::ThermalConductivitySurfaceEvaluator(
      const ThermalConductivitySurfaceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    uf_key_(other.uf_key_),
    height_key_(other.height_key_),
    K_liq_(other.K_liq_),
    K_ice_(other.K_ice_) {}
    

Teuchos::RCP<FieldEvaluator>
ThermalConductivitySurfaceEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivitySurfaceEvaluator(*this));
}


void ThermalConductivitySurfaceEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> uf = S->GetFieldData(uf_key_);
  Teuchos::RCP<const CompositeVector> height = S->GetFieldData(height_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& uf_v = *uf->ViewComponent(*comp,false);
    const Epetra_MultiVector& height_v = *height->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = height[0][i] * (K_liq_ * eta[0][i] + K_ice_ * (1. - eta[0][i]));
    }
  }
}


void ThermalConductivitySurfaceEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0); // not implemented, not yet needed
}

} //namespace
} //namespace
} //namespace
