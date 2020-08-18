/*
  mantle.cpp
  - created by Yoshi Miyazaki (November 2019)
  
  ... solves for the thermal evolution of the mantle after the solidification of magma ocean
  
 */

#include "./system.h"

tuple<double, double, double> surface_mantle_system::Nu_heat_scaling(double T0, double T1, double Deta_stiff){
    /* 
       the mantle heat loss is calculated based on the scaling of Korenaga (2010a)
       :
       surface heat loss: Q  = (k*dT/D) * Nu
       Nusselt number   : Nu = (Rai/Rac)^(1/3) / (DetaL)^(1/3)
       
       Rai: internal Rayleigh number

       [ INPUT ]
       T0:       surface temperature
       T1:       mantle potential temperature
       Deta_dry: dehydration stiffening
     */
    
    /* calc the Rayleigh number
       :
       Ra = alpha*g*rho0*(T1-T0)*D^3 / (kappa*eta)
       eta:    viscosity
       E_diff: activation energy for viscosity
       --- reference viscosity is 10^19 at 1350 C.
     */
    double al_phi = 26;
    double phi = 0.0; //melt_fraction(T1, homo);
    double eta = eta_now * exp(E_diff/R*(1/T1 - 1/T_now)) * exp(-al_phi * phi);
    
    /* the effect of water
       : T > 1600C corresponds to early Earth, when Earth was completely dry.
       . The mantle should have been more viscous.

       But, for the purpose of this paper, we assume that the water was as hydrated as present.
    */
    eta *= 1.; // exp(log(Deta_stiff) * 1.);
    
    double Rai = alpha*gv*rhom*T1*(D*D*D)/(kappa*eta);
    // cout << "phi: "<< T1 << "\t" << exp(-al_phi * phi) << "\t" << eta_now  << endl;
    
    /* calc the "effective" viscosity contrast w/in the lithosphere (DetaL)
       : Eq.(37-38) of Korenaga (2010a) 
       
       DetaL = exp( 0.327 * gamma^0.647 * FK )
       gamma: normalized friction coefficient
       FK   : Frank-Kamenetskii parameter     (Ediff*DT)/(R*T^2)
    */
    double mu_brit   = 0.03;
    double gam_brit  = mu_brit/(alpha*(T1-T0));
    double FK        = E_diff*(T1-T0)/(R*T1*T1);
    double DetaL_eff = exp(0.327*pow(gam_brit, 0.647) * FK);
    
    /* heat loss scaling
       : Eq. (35) of Korenaga (2010a) */
    double Rac = 1000;
    double Nu_nonstiff = 2. * pow( (Rai/Rac)/DetaL_eff, 1.0/3);
    
    /* the effect of dehydration stiffening 
       : Eq. (44) of Korenaga (2010a)
       
       DetaL = DetaL_eff * exp( ln(Deta_stiff) * minv( 1.0, H_deyd/(chi*H_therm) ) )
       chi   = 6.0
     */
    double H_therm = D/Nu_nonstiff;
    double H_dehyd, t_cooling;  tie(H_dehyd, t_cooling) = lithosphere_cooling(Ts, Tm);
    double r_stiff   = minv(1.0, H_dehyd/(6*H_therm));
    
    double DetaL      = DetaL_eff * exp( log(Deta_stiff) * r_stiff);
    if (!homo){ DetaL = DetaL_eff; }
    
    /* heat loss scaling (again)
       : Eq. (35) of Korenaga (2010a) */
    double Nu = 2. * pow( (Rai/Rac)/DetaL, 1.0/3 );
    double W  = 1. * Nu * kT * (T1-T0)/D;
    double Q  = W * (4*pi*square(rE));

    /* plate velocity scaling
       : 
       1. vs = 16. * kappa * aspect * Nu^2 / D
       . To match the present-day mantle, the factor 16 is ignored here.
       
       2. compare with the modern
       . vs = 5 * ((Q/DT)/(Q0/DT0))^2
     */
    double aspect = 4.0;
    // double vs = aspect * Nu*Nu*kappa/D;
    double Q0 = 38.3e12, DT0 = 1300;
    double vs = 0.05 * square((Q/(T1-T0))/(Q0/DT0)) / yr2sec;
    // cout << Q << "\t" << T1-T0 << " -- \t " << Q << "\t" << DT0 << endl;
    // exit(2);
    
    /*
      for heterogeneous mantle, we scale the plate velocity using modern.
      considering that limited dehyration takes place, it should be faster, but we will ignore such effect.
     */
    if (0){
        double Ra_ref = alpha*gv*rhom*1300*(D*D*D)/(kappa*1e19);
        vs = 0.05/yr2sec * pow(Rai/Ra_ref, 2./3);
    }

    /* output */
    // cout << Q/1e12 << " TW \t" << vs << "\t" << square((Q/(T1-T0))/(Q0/DT0)) << "\t" << Nu << "\t" << DetaL << "\t" << DetaL_eff << "\t" << exp( log(Deta_stiff) * r_stiff) << endl;
    
    /* the effect of chemical buoyancy
       ... calc critical cooling time for lithosphere to be negatively buoyant.
       .   t_now (in Myr) shows the critical time for modern (20 Myr)           */

    //cout << "thermal: " << H_therm << "\t" << rH << " \t eff_dry: " << pow(eff_dry, 1.0/3) << "\t" << endl;

    // cout << "Nu: " <<  Nu << "\t DetaL_eff: " << DetaL_eff << "\t T1 [K]: " << T1 << endl;
    // cout << "\t DetaL_eff: " << DetaL_eff << endl;
    
    return forward_as_tuple(W, Q, vs);
}
tuple<double, double> surface_mantle_system::lithosphere_cooling(double T0, double T1){
    
    /* lithosphere thickness from Korenaga (2006) */
    double dphidP = 0.13;
    double Pmelt  = (T1-(1150+273.15))/100;                       // pressure where melting starts (in GPa)
    double Pcrus  = Pmelt + 1/dphidP*(1- sqrt(1+2*Pmelt*dphidP)); // pressure at the crust bottom  (in GPa)
    double phi    = 0.5*Pmelt*dphidP;                             // average melt fraction
    double rhol   = rhom - 15*8*phi;                              // lithosphreic density
    double dl    = (Pmelt-Pcrus)*(rhom*gv);                       // thickness of lithosphere
    
    /* crustal thickness from Korenaga (2006) */
    double Pa = (Pmelt+Pcrus)/2.0;
    double Pa_t  = 1.0;
    double phi_t = 0.05;
    double WL    = 1./4 * (1 - tanh(0.6*(Pa-Pa_t))) * (1 - tanh(8.4*(phi-phi_t)));
    double WH    = 1./4 * (1 + tanh(0.6*(Pa-Pa_t))) * (1 + tanh(8.4*(phi-phi_t)));
    double Vp    = 7.52 + WL*(-1.73 - 0.55*Pa + 7.71*phi - 0.11*Pa*Pa + 8.87*Pa*phi - 146.11*phi*phi);
    Vp += WH*(-0.35 + 0.034*Pa + 0.51*phi + 0.0016*Pa*Pa - 0.04*Pa*phi + 0.046*phi*phi);
    
    double rhoc  = (0.77 + 0.3*Vp)*1000;
    double dc    = dl*phi;

    /* calculate thermal boundary layer (TBL) from thermal buoanycy
       ... 
       .   drho_c (chemical buoanycy) = drho_t (thermal negative buoanycy)
       .
       .   drho_c = (rhom-rhomc)*dc + (rhom-rhol)+dl  
       .   drho_t = alpha*rho* average(dT) * h
     */
    double thermal_b  = alpha*rhom*(T1-T0)*0.52;
    double TBL        = ((rhom-rhoc)*dc + (rhom-rhol)*dl) / thermal_b;

    /* calculate cooling timescale assuming half-space cooling [in Myr!!]
       ...
       TBL = 2. * sqrt(kappa * tau) 
    */
    double t_cooling = square(TBL/2.0) / kappa / yr2sec / 1e6;
    
    return forward_as_tuple(dl, t_cooling);
    
}
double surface_mantle_system::radioactive_heating(double time){

    /* Radiogenic heat production from Korenaga (2006)
       hn: hear contribution of each isotope
       hn[0]: 238U
       hn[1]: 235U
       hn[2]: 232Th
       hn[3]: 40K
    */
    tensor1d<double>  hn(0.0, 4);
    hn[0] = 0.372;  hn[1] = 0.0164;  hn[2] = 0.430;  hn[3] = 0.181;

    /* time decay from Korenaga (2006)
       lb: exponential decay factor
    */
    tensor1d<double>  lb(0.0, 4);
    lb[0] = 0.155;  lb[1] = 0.985;   lb[2] = 0.0495; lb[3] = 0.555;

    /* heat production */
    double H = 0.0, H_now = 38e12*0.31;
    for (int i=0; i<(int)hn.size(); i++){
        H += hn[i] * exp(lb[i] * (4.5 - time/yr2sec/1e9));
    }
    H *= H_now;
    
    return H;
}
tuple<double, double, double> surface_mantle_system::core_heat_flux(double T0, double T1){
    /* constants */
    const double CpFe = 800;                     /* heat capacity of Fe metal */
    const double Mcore  = (1.835+0.09675)*1e24;  /* core mass                 */
    const double rcore  = 3485e3;                /* core radius               */

    /* core heat capacity */
    const double Cpcore = Mcore*CpFe;
    
    /* heat loss scaling */
    double Rai = alpha*gv*rhom*T0*(D*D*D)/(kappa*eta_now*1e2);
    double Rac = 1000;
    
    double Wc = kT * (T1-T0)/D * pow(Rai/Rac,1.0/3);
    double Qc = Wc * 4*M_PI*(rcore*rcore);

    return forward_as_tuple(Wc, Qc, Cpcore);
}
double surface_mantle_system::melt_fraction(double T1, bool homo0){
    /* 
       Because solidus & liquidus change with Mg#,
       two cases are calculated separately.
       
       For homogeneous mantle, we assume 
     */
    double phi;
    if (homo0){
        phi = maxv(0., (T1 - 1800)/1000);
    }else{
        phi = maxv(0., (T1 - 2000)/1000);
    }

    return minv(phi, 1.);
}
