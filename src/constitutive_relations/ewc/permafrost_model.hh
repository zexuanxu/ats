/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_PERMAFROST_MODEL_
#define AMANZI_PERMAFROST_MODEL_

#include "Teuchos_ParameterList.hpp"

#include "Tensor.hh"
#include "Point.hh"

#include "wrm_partition.hh"
#include "compressible_porosity_model_partition.hh"
#include "compressible_porosity_leijnse_model_partition.hh"
#include "ewc_model_base.hh"

namespace Amanzi {

namespace Flow {
namespace Flow {
class PCIceWater;
class PCLiqAtm;
}
}

namespace Energy { namespace Energy { class IEM; class IEMWaterVapor;} }
namespace Relations { class EOS; class VaporPressureRelation; }


class PermafrostModel : public EWCModelBase {

 public:
  PermafrostModel() {}

  virtual void InitializeModel(const Teuchos::Ptr<State>& S,
                               Teuchos::ParameterList& plist);
  virtual void UpdateModel(const Teuchos::Ptr<State>& S, int c);
  virtual bool Freezing(double T, double p);
  virtual int EvaluateSaturations(double T, double p,
                                  double& s_gas, double& s_liq, double& s_ice);
  
 protected:
  bool IsSetUp_();

  int EvaluateEnergyAndWaterContent_(double T, double p,
          AmanziGeometry::Point& result);

 protected:
  Teuchos::RCP<Flow::Flow::WRMPermafrostModelPartition> wrms_;
  Teuchos::RCP<Flow::Flow::WRMPermafrostModel> wrm_;
  Teuchos::RCP<Relations::EOS> liquid_eos_;
  Teuchos::RCP<Relations::EOS> gas_eos_;
  Teuchos::RCP<Relations::EOS> ice_eos_;
  Teuchos::RCP<Flow::Flow::PCIceWater> pc_i_;
  Teuchos::RCP<Flow::Flow::PCLiqAtm> pc_l_;
  Teuchos::RCP<Relations::VaporPressureRelation> vpr_;
  Teuchos::RCP<Energy::Energy::IEM> liquid_iem_;
  Teuchos::RCP<Energy::Energy::IEMWaterVapor> gas_iem_;
  Teuchos::RCP<Energy::Energy::IEM> ice_iem_;
  Teuchos::RCP<Energy::Energy::IEM> rock_iem_;
  Teuchos::RCP<Flow::Flow::CompressiblePorosityModelPartition> poro_models_;
  Teuchos::RCP<Flow::Flow::CompressiblePorosityModel> poro_model_;

  Teuchos::RCP<Flow::Flow::CompressiblePorosityLeijnseModelPartition> poro_leij_models_;
  Teuchos::RCP<Flow::Flow::CompressiblePorosityLeijnseModel> poro_leij_model_;

  double p_atm_;
  double poro_;
  double rho_rock_;
  bool poro_leij_;
  Key domain;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};



} //namespace


#endif
