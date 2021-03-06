/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*

  EOS for salt water (does not implement viscosity at this point!)

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_SW_HH_
#define AMANZI_RELATIONS_EOS_SW_HH_

#include "Teuchos_ParameterList.hpp"
#include "eos.hh"
#include "dbc.hh"

namespace Amanzi {
namespace Relations {

class EOS_SW : public EOS {

 public:
  explicit EOS_SW(Teuchos::ParameterList& eos_plist);
  virtual ~EOS_SW() {};

  // Virtual methods that form the EOS
  virtual double MassDensity(std::vector<double>& params);
  virtual double DMassDensityDC(std::vector<double>& params);  
  virtual double DMassDensityDT(std::vector<double>& params);
  virtual double DMassDensityDp(std::vector<double>& params);


  virtual double MolarDensity(std::vector<double>& params);
  virtual double DMolarDensityDC(std::vector<double>& params);  
  virtual double DMolarDensityDT(std::vector<double>& params);
  virtual double DMolarDensityDp(std::vector<double>& params);


  // If molar mass is constant, we can take some shortcuts if we need both
  // molar and mass densities.  MolarMass() is undefined if
  // !IsConstantMolarMass()
  virtual bool IsConstantMolarMass(){return false;}
  virtual double MolarMass() {AMANZI_ASSERT(0); return 0.0;}

protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double E_;
  double rho_f_;
  double M_salt_;
  double M_water_;

 private:
  static Utils::RegisteredFactory<EOS, EOS_SW> factory_;

  
};

} // namespace
} // namespace

#endif
