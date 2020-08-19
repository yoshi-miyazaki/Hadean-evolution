/*
  atmos1d.h
  ... construct 1-D atmosphere model object.

  Inputs:
  - Psurf0: surface pressure of non-greenhouse gases
  - xCO20 : molar ratio of CO2 of total atmosphere (0: no CO2, 1: pure CO2)
  
  Jun 5, 2020 by Yoshi Miyazaki
*/

#include <iomanip>
#include "../GibbsFE_minimization/tensor.h"
#include "../GibbsFE_minimization/const.h"

class atmos1d{
 public:
    atmos1d(double, double);

    void   record_result(string);
    
    void             calc_Tprofile(double);
    tensor1d<double> calc_tot_vapor();
    
    double get_Tsurf(){ return T[T.size()-1]; };
    double get_Ttp()  { return Ttp;           };
    double get_Ptp()  { return Ptp;           };
    double get_topxv(){ return q[0];          };
    double get_tautp(){ return tautp;         };
    
    tensor1d<double> get_T(){   return T;   };
    tensor1d<double> get_q(){   return q;   };
    tensor1d<double> get_tau(){ return tau; };
    tensor1d<double> get_Fup(){ return Fup; };
    
 private:
    double           Psurf;
    double           xCO2;
    double           Ttp;
    double           Ptp;
    double           tautp;
    
    tensor1d<double> T;     /* temperature grid   */
    tensor1d<double> P;     /* pressure    grid   */
    tensor1d<double> q;     /* water mixing ratio */
    tensor1d<double> tau;   /* optical thickness  */
    tensor1d<double> Fup;   /* upward radiative flux */
    
    /* functions used inside */
    double                 water_saturation(double);
    double                 moist_adiabat(double, double);
    double                 radflux_up   (double, double, double, double);
    
    tuple<double, double>  tropopause(double, double, double);
    double                 delta_tau  (double, double&, double&, double&, double&);
    double                 d_delta_tau(double, double&, double&, double&, double&);

    void                   calc_Tprofile_3layers(double, double, double, double, double, double);
    double                 radflux_up_3layers(double, double, double, double, double, double, double);
    tuple<double, double, double>  minxv_for_tropopause(double, double, double, double);
    
    /* constants */
    double muv   = 18e-3;
    double mun   = 28e-3;
    double muCO2 = 44e-3;
    
    double Cpv   = 4*8.31;
    double Cpn   = 3.5*8.31;
    
    double kpav   = 1.0e-2;
    double kpaCO2 = 1.3e-4;
    
    double L      = 43655; // 2442e3*18e-3;
    double gv     = 9.8;

};
