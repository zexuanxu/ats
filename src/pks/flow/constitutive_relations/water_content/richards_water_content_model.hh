/*
  The richards water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RICHARDS_WATER_CONTENT_MODEL_HH_
#define AMANZI_FLOW_RICHARDS_WATER_CONTENT_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class RichardsWaterContentModel {

 public:
  explicit
  RichardsWaterContentModel(Teuchos::ParameterList& plist);

  double WaterContent(double phi, double sl, double nl, double cv) const;

  double DWaterContentDPorosity(double phi, double sl, double nl, double cv) const;
  double DWaterContentDSaturationLiquid(double phi, double sl, double nl, double cv) const;
  double DWaterContentDMolarDensityLiquid(double phi, double sl, double nl, double cv) const;
  double DWaterContentDCellVolume(double phi, double sl, double nl, double cv) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:



};

} //namespace
} //namespace
} //namespace

#endif