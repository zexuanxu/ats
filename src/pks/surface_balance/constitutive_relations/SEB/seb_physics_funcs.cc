/*
  Functions for calculating the snow / surface energy balance.
*/

#include <cmath>
#include <algorithm>
#include <boost/math/tools/roots.hpp>

#include "dbc.hh"

#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {

void UpdateIncomingRadiation(const SEB& seb, EnergyBalance& eb, bool debug) {
  // Calculate incoming short-wave radiation
  eb.fQswIn = (1 - seb.in.surf.albedo) * seb.in.met.QswIn;

  // Calculate incoming long-wave radiation
  const ThermoProperties& vp_air = seb.in.met.vp_air;
  double e_air = std::pow(10*vp_air.actual_vaporpressure, vp_air.temp / 2016.);
  e_air = 1.08 * (1 - std::exp(-e_air));
  eb.fQlwIn = e_air * seb.params.stephB * std::pow(vp_air.temp,4);

  // Calculate D_h, D_e, sqig
  eb.Dhe = std::pow(seb.params.VKc,2) * seb.in.met.Us
                       / std::pow(std::log(seb.params.Zr / seb.in.surf.Zo), 2);

  if (debug) {
    std::cout << "Incoming Radiation Energy Terms:" << "\n"
              << "  fQswIn   = " << eb.fQswIn << "\n"
              << "  fQlwIn   = " << eb.fQlwIn << std::endl;
  }

}


void UpdateEnergyBalance(const SEB& seb, const ThermoProperties& vp_surf, EnergyBalance& eb, bool debug) {
  const ThermoProperties& vp_air = seb.in.met.vp_air;

  // Calculate outgoing long-wave radiation
  eb.fQlwOut = -seb.in.surf.emissivity*seb.params.stephB*std::pow(vp_surf.temp,4);

  double Sqig;          // special constant for use in e and h, precalculated for efficiency
  double air_temp = seb.in.met.vp_air.temp;
  double Ri  = seb.params.gravity * seb.params.Zr * (air_temp - vp_surf.temp)
      / (air_temp * std::pow(seb.in.met.Us,2));
  if (Ri >= 0.) {
    // stable condition or snow?
    Sqig = 1 / (1 + 10*Ri);
  } else {
    // Unstable condition
    Sqig = (1-10*Ri);
  }

  // Calculate sensible heat flux
  eb.fQh = seb.params.density_air * seb.params.Cp * eb.Dhe * Sqig * (vp_air.temp - vp_surf.temp);

  // Calculate latent heat flux
  eb.fQe = vp_surf.porosity * seb.params.density_air * seb.params.Ls * eb.Dhe * Sqig * 0.622
      * (vp_air.actual_vaporpressure - vp_surf.actual_vaporpressure) / seb.params.Apa;

  // Calculate heat conducted to ground, if snow
  if (seb.out.snow_new.ht > 0.) {
    double Ks = 2.9e-6 * std::pow(seb.out.snow_new.density,2);
    eb.fQc = Ks * (vp_surf.temp - seb.in.vp_ground.temp) / seb.out.snow_new.ht;
  }


  if (debug) {
    std::cout << "Energy Balance Terms (ht_snow = " << seb.out.snow_new.ht << "):" << "\n"
              << "  fQlwOut  = " << eb.fQlwOut << "\n"
              << "  fQh      = " << eb.fQh << "\n"
              << "  fQe      = " << eb.fQe << "\n"
              << "  fQc      = " << eb.fQc << std::endl;
  }

}


void UpdateMassBalance(const SEB& seb, MassBalance& mb, EnergyBalance& eb, SnowProperties& snow_new, bool debug) {
  // Melt rate given by available energy rate divided by heat of fusion.
  mb.Mm = eb.fQm / (seb.in.vp_ground.density_w * seb.params.Hf);

  // Condensation rate
  mb.Me = eb.fQe / (seb.in.vp_ground.density_w * seb.params.Ls);

  // Snow balance
  if (seb.in.snow_old.ht > 0.) {
    double swe_old = seb.in.snow_old.ht * seb.in.snow_old.density / seb.in.vp_ground.density_w;
    double swe_new = swe_old + (seb.in.met.Ps - mb.Mm + mb.Me)*seb.in.dt;

    // First do a pass to ensure we are not melting or sublimation ALL
    // of the available snow.  If so, adjust...
    if (swe_new < 0.) {
      // Must adjust... we are melting or sublimating all of the available snow.
      if (mb.Mm > 0.) {
	// -- stop melting, and push the extra heat off into the ground
	double swe_new_nomelting = swe_old + (seb.in.met.Ps + mb.Me)*seb.in.dt;

	double swe_melt = 0.;
	if (swe_new_nomelting >= 0) {
	  swe_melt = swe_new_nomelting;  // melt all we have
	}
	mb.Mm = swe_melt / seb.in.dt;
	swe_new = swe_old + (seb.in.met.Ps - mb.Mm + mb.Me)*seb.in.dt;

	// fix energy balance
	double fQm_new = mb.Mm * (seb.in.vp_ground.density_w * seb.params.Hf);
	eb.fQc += eb.fQm - fQm_new;
	eb.fQm = fQm_new;
      }

      if (swe_new < 0.) {
	// -- Either there was no melting (sublimating it all) or we have
	// turned off all melting and are still in the negative
	// (sublimating too much of it).  Turn down sublimation.
	// Note we do not fix the energy balance, as I'm not sure how to do so...
	double swe_subl = -mb.Me * seb.in.dt;
	ASSERT(swe_subl > 0.);
	swe_subl += swe_new;
	ASSERT(swe_subl >= 0.);
	mb.Me = -swe_subl / seb.in.dt;

	swe_new = swe_old + (seb.in.met.Ps - mb.Mm + mb.Me)*seb.in.dt;
      }
      ASSERT(swe_new > -1.e-20);
    }

    // age the old snow
    double age_settled = seb.in.snow_old.age + seb.in.dt / 86400.;
    double dens_settled = seb.params.density_freshsnow
        * std::max(std::pow(age_settled, 0.3), 1.);

    // Match frost age with assigned density -- Calculate which day frost
    // density matched snow defermation function from (Martinec, 1977)
    double age_frost = std::pow((seb.params.density_frost / seb.params.density_freshsnow),(1/0.3)) - 1 + seb.in.dt / 86400.;

    // precip
    double age_precip = seb.in.dt / 86400.;

    // determine the new height
    // -- sources
    double swe_settled = swe_old;
    double swe_frost = mb.Me > 0. ? mb.Me*seb.in.dt : 0.;
    double swe_precip = seb.in.met.Ps*seb.in.dt;

    // -- sinks
    double swe_subl =  mb.Me < 0. ? -mb.Me*seb.in.dt : 0.;
    double swe_melt = mb.Mm * seb.in.dt;

    // -- sublimate precip first
    ASSERT(swe_subl >= 0.);
    if (swe_subl > 0.) {
      if (swe_subl > swe_precip) {
	swe_precip = 0.;
	swe_subl -= swe_precip;
      } else {
	swe_precip -= swe_subl;
	swe_subl = 0.;
      }
    }

    // -- next sublimate settled snow
    if (swe_subl > 0.) {
      if (swe_subl > swe_settled) {
	ASSERT(false);
      } else {
	swe_settled -= swe_subl;
	swe_subl = 0.;
      }
    }

    // -- melt settled snow first
    ASSERT(swe_melt >= 0.);
    if (swe_melt > 0.) {
      if (swe_melt > swe_settled) {
        swe_settled = 0.;
        swe_melt -= swe_settled;
      } else {
        swe_settled -= swe_melt;
        swe_melt = 0.;
      }
    }

    // -- now melt frost, precip by even amounts
    if (swe_melt > 0.) {
      double swe_melt_from_frost = swe_melt * (swe_frost / (swe_frost + swe_precip));
      double swe_melt_from_precip = swe_melt - swe_melt_from_frost;

      swe_frost -= swe_melt_from_frost;
      swe_precip -= swe_melt_from_precip;
    }

    // -- check we didn't screw up
    ASSERT(swe_settled >= 0.);
    ASSERT(swe_frost >= 0.);
    ASSERT(swe_precip >= 0.);

    // -- convert these to heights
    double ht_settled = swe_settled * seb.in.vp_ground.density_w / dens_settled;
    double ht_frost = swe_frost * seb.in.vp_ground.density_w / seb.params.density_frost;
    double ht_precip = swe_precip * seb.in.vp_ground.density_w / seb.params.density_freshsnow;

    // set the snow properties
    snow_new.ht = ht_settled + ht_frost + ht_precip;
    snow_new.age = (swe_settled*age_settled + swe_frost*age_frost + swe_precip*age_precip)
        / (swe_settled + swe_frost + swe_precip);
    snow_new.density = swe_new * seb.in.vp_ground.density_w / snow_new.ht;

    // set the water properties
    // -- water source to ground is (corrected) melt and rainfall
    mb.MWg = mb.Mm + seb.in.met.Pr;
    mb.MWg_subsurf = 0.;
    mb.MWg_temp = mb.MWg > 0. ? (mb.Mm * 273.15 + seb.in.met.Pr * seb.in.met.vp_air.temp) / mb.MWg : seb.in.met.vp_air.temp;


  } else {
    // set the snow properties
    snow_new.ht = seb.in.met.Ps * seb.in.dt
        * seb.in.vp_ground.density_w / seb.params.density_freshsnow;
    snow_new.age = seb.in.dt / 86400.;
    snow_new.density = seb.params.density_freshsnow;

    // set the water properties
    // -- water source to ground is rainfall + condensation
    // -- evaporation is taken from ground if ponded water, from cell source if not (with transition)
    mb.MWg_temp = seb.in.met.vp_air.temp;
    mb.MWg = seb.in.met.Pr;
    mb.MWg_subsurf = 0.;
    if (mb.Me > 0.) {
      mb.MWg += mb.Me;

    } else {
      double surf_p = seb.in.vp_ground.pressure;
      double trans_factor = surf_p > seb.params.Apa*1000. ? 0. :
          surf_p < seb.params.Apa*1000. - seb.params.evap_transition_width ? 1. :
          (seb.params.Apa*1000. - surf_p) / seb.params.evap_transition_width;

      std::cout << "ground pres = " << surf_p << std::endl;
      std::cout << "trans factor = " << trans_factor << std::endl;

      mb.MWg += (1-trans_factor) * mb.Me;
      mb.MWg_subsurf += trans_factor * mb.Me;
    }
  }


  if (debug) {
    std::cout << "Mass Balance:\n"
              << "  Mm   = " << mb.Mm << "\n"
              << "  Me   = " << mb.Me << "\n"
              << "  Snow Melt:\n"
              << "    new ht   = " << snow_new.ht << "\n"
              << "    new age  = " << snow_new.age << "\n"
              << "    new dens = " << snow_new.density << "\n"
              << "  Water Balance:\n"
              << "    surf src = " << mb.MWg << "\n"
              << "    sub src  = " << mb.MWg_subsurf << std::endl;
  }
}


// Snow temperature calculation.
double DetermineSnowTemperature(const SEB& seb, ThermoProperties& vp_snow,
        EnergyBalance& eb, std::string method) {
  SnowTemperatureFunctor_ func(&seb, &vp_snow, &eb);
  Tol_ tol(1.e-6);
  uintmax_t max_it(50);
  double left, right;
  double res_left, res_right;

  double res_init = func(vp_snow.temp);
  if (res_init < 0.) {
    right = vp_snow.temp;
    res_right = res_init;

    left = vp_snow.temp - 1.;
    res_left = func(left);
    while (res_left < 0.) {
      left = left - 1.;
      res_left = func(left);
    }
  } else {
    left = vp_snow.temp;
    res_left = res_init;

    right = vp_snow.temp + 1.;
    res_right = func(right);
    while (res_right > 0.) {
      right = right + 1.;
      res_right = func(right);
    }
  }

  std::pair<double,double> result;
  if (method == "bisection") {
    result = boost::math::tools::bisect(func, left, right, tol, max_it);
  } else if (method == "toms") {
    result = boost::math::tools::toms748_solve(func, left, right, res_left, res_right, tol, max_it);
  }

  return (result.first + result.second)/2.;
}

// master driver
void CalculateSurfaceBalance(SEB& seb, bool debug) {
  // initialize the data
  seb.in.met.vp_air.UpdateVaporPressure();

  // Energy balance
  UpdateIncomingRadiation(seb, seb.out.eb, debug);
  if (seb.in.snow_old.ht > 0.) {
    // snow on the ground, solve for snow temperature
    double T_snow = DetermineSnowTemperature(seb, seb.in.vp_snow, seb.out.eb);
    if (T_snow > 273.15) {
      // limit snow temp to 0, then melt with the remaining energy
      seb.in.vp_snow.temp = 273.15;
      seb.in.vp_snow.UpdateVaporPressure();
      UpdateEnergyBalance(seb, seb.in.vp_snow, seb.out.eb, debug);
      seb.out.eb.BalanceViaMelt();
    } else {
      // snow not melting
      seb.in.vp_snow.temp = T_snow;
      seb.in.vp_snow.UpdateVaporPressure();
      UpdateEnergyBalance(seb, seb.in.vp_snow, seb.out.eb, debug);
      seb.out.eb.fQm = 0.;
    }

  } else {
    // no snow on the ground, balance given by conduction
    seb.in.vp_ground.UpdateVaporPressure();
    UpdateEnergyBalance(seb, seb.in.vp_ground, seb.out.eb, debug);
    seb.out.eb.BalanceViaConduction();
  }

  // Mass balance
  UpdateMassBalance(seb, seb.out.mb, seb.out.eb, seb.out.snow_new, debug);
}

double CalcAlbedoSnow(double density_snow) {
  double AlSnow;
  if (density_snow <= 432.23309912785146) {
    AlSnow = 1.0 - 0.247 * std::pow(0.16 + 110*std::pow(density_snow/1000, 4), 0.5);
  } else {
    AlSnow = 0.6 - density_snow / 4600;
  }
  return AlSnow;
}

double CalcRoughnessFactor(double air_temp) {
  double Zsmooth = 0.005;
  double Zrough = 0.04;

  double Zfraction = air_temp < 270. ? 1. : air_temp > 280. ? 0. : -0.1*air_temp + 28;
  return Zsmooth * Zfraction + Zrough * (1. - Zfraction);
}

} // namespace
} // namespace
} // namespace