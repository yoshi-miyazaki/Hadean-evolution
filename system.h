/*
  system.h
  - created by Yoshi Miyazaki (November 2019)
  
  `surface_mantle_system`
  ... contains information of two components:
  .   1. surface (ocean + atmosphere)
  .      - Ts:    surface temperature
  .      - natms: volatile amount & composition
  .
  .   2. mantle
  .      - Tm:    mantle potential temperature
  .      - nmant: volatile composition

 */

#include "../GibbsFE_minimization/CGgibbsmin.h"
#include "../tensor.h"

class surface_mantle_system{
 public:
    surface_mantle_system(double, double, double, double, bool, double);
    void    next(double, double, double);
    void    set_Tm(double Tm0){ Tm = Tm0; }
    
    /* atmosphere evolution */
    double  surface_T(double);
    void    volatile_loss(double, double);
    
    /* mantle cooling */
    tuple<double, double, double> Nu_heat_scaling(double, double, double);
    tuple<double, double>         lithosphere_cooling(double, double);
    double                        radioactive_heating(double);
    double                        melt_fraction(double, bool);
    
    /* core cooling */
    tuple<double, double, double> core_heat_flux(double, double);
    
 private:
    /* atmosphere + ocean component */
    double           Ts;          // surface temperatBure
    tensor1d<double> natms;       // volatile composition @surface
    
    /* mantle component */
    double Tm;                    // mantle potential temperature
    double D;                     // the mantle depth
    bool   homo;                  // is the mantle homogeneous?
    tensor1d<double> nmant;       // volatile composition @interior
    
    /* core comonent */
    double Tc;
    
    /* constants */
    const double alpha = 3.3e-5;  // thermal expansivity
    const double kappa = 1e-6;    // thermal diffusivity
    const double kT    = 3.8;     // thermal conductivity
    const double gv    = 9.807;   // gravity
    const double rhom  = 3300;    // mantle density
    
    /* reference parameters */
    const double T_now   = 1350 + (Tref - 25);
    const double eta_now = 1e19;  // mantle viscosity at T = T_now for a wet mantle
    const double E_diff  = 300e3; // activation energy for olivine

};
