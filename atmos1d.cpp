/*
  atmos1d.cpp 
  ... construct 1-D atmosphere model object.
  
  (see atmos1d.h for information)

  Jun 5, 2020 by Yoshi
*/

#include "./system.h"

atmos1d::atmos1d(double Psurf0, double xCO20){
    /* 1D atmosphere model
       ... prepare grid and atmos composition */
    Psurf = Psurf0;
    xCO2  = xCO20;
    
    // grid size
    int N = 500;
    T.resize(0.0, N);
    P.resize(0.0, N);
    q.resize(0.0, N);
    tau.resize(0.0, N);
    Fup.resize(0.0, N);
}
tensor1d<double> atmos1d::calc_tot_vapor(){
    /* calc total mass of each atmospheric component 
       Matm[0]: water vapor
       Matm[1]: carbon dioxide */
    tensor1d<double> Matm(0.0, 2);

    int N = (int)T.size();
    for (int i=1; i<N; i++){
        /* grid size */
        double mu   = q[i]*muv + (1-q[i])*(muCO2*xCO2 + mun*(1-xCO2));
        double rho  = P[i]*mu/(R*T[i]);
        double dz   = (P[i]-P[i-1])/(rho*gv);
        
        /* water vapor */
        double Pv   = q[i]*P[i];
        double rhov = Pv*muv/(R*T[i]);
        Matm[0] += rhov * 4*pi*square(rE) * dz;
        
        /* carbon dioxide */
        double Pc   = (1.-q[i])*P[i]*xCO2;
        double rhoc = Pc*muCO2/(R*T[i]);
        Matm[1] += rhoc * 4*pi*square(rE) * dz;
    }
    // cout << Psurf*4*pi*rE*rE/gv * (44./29) * xCO2 << "\t" << Matm[1] << endl;
    
    return Matm;
}
void atmos1d::calc_Tprofile(double Fout){
    
    /* search for correct xv (water vapor mixing ratio at the top of the atmosphere */
    double xv0 = 1e-10,  xv1 = 0.99;
    
    double rh = 1.0;
    double df0 = radflux_up(Fout, Psurf, xv0, rh);
    double df1 = radflux_up(Fout, Psurf, xv1, rh);
    
    /* determine search range:
       before performimng bisection search, check the signs of `df0` and `df1`
       ...
       The standard case is:
       -> df0 < 0 (low Tsurf; too little flux) & df1 > 0 (high Tsurf; too much flux)
    */
    // cout << setprecision(6) << "beg: " << xv0 << "\t" << xv1 << "\t" << df0 << "\t" << df1 << endl;
    
    /* 1. df0 == 0: tropopause cannot be calclulated. 
       This means too little vapor (or too much CO2), so greenhouse effect by CO2 shoudl be reduced 
       by increasing x_vapor. */
    while (df0 == 0){
        xv0 *= 1.05;
        df0 = radflux_up(Fout, Psurf, xv0, rh);
        // cout << "f0-1: " << xv0 << "\t" << df0 << "\t" << P[99] <<  "\t" << tau[99] << " " << df1 << endl;
        
        if (xv0>xv1){
            if (df1==0){ cout << "df0-Komabayashi-Ingersoll limit" << endl; }
            else       { cout << "A. ??? limit " << xv0 << "\t" << xv1 << "\t" << df0 << "\t" << df1 << endl; }
            exit(1);
        }
    }

    /* 1-2. If df0 > 0: it is possible that the initial xv0 is too large.
       But, what is more likely is that the standard case does not satisfy any radiatively stable atmosphere. 
       In such a case, the atmosphere should be consisted of three layers, 
       - stratosphere
       - radiative troposphere
       - convective troposphere
    */
    while (df0 > 0 and df1 > 0){
        // cout << "search for thermal profile of 3 layers: " << flush;
        double xvlim;
        tie(xvlim, Ttp, tautp) = minxv_for_tropopause(Fout, rh, xv0*0.1, xv1); 
        // cout << xvlim << "\t" << Ttp << "\t" << tautp << "\t df: " << df0 << "\t" << df1 << endl;       
        calc_Tprofile_3layers(Fout, Psurf, rh, xvlim, Ttp, tautp);
        return;
        /*
          xv0 *= 1.01;
          df0 = radflux_up(Fout, Psurf, xv0, rh);
          cout << "f0-2: " << xv0 << "\t" << df0 << "\t" << P[99] <<  "\t" << tau[99] << " " << df1 << endl;
          
          if (xv0>xv1){
          if (df1==0){ cout << "B. Komabayashi-Ingersoll limit" << endl; }
          else       { cout << "B. ??? limit " << xv0 << "\t" << xv1 << "\t" << df0 << "\t" << df1 << endl; }
          exit(2);
          }
        */
        
    }
    
    while (df0 < 0 and df1 <= 0){
        double df1_old = df1;
        
        /*  */
        xv1 *= 0.95;
        df1 = radflux_up(Fout, Psurf, xv1, rh);
        
        // cout << "f1: " << xv1 << "\t" << df1 << "\t" << P[99] <<  "\t" << tau[99] << " " << endl;
        if (xv1 < xv0){
            if (df0 == 0 and df1==0.){
                cout << "df1-Komabayashi-Ingersoll " << xv0 << "\t" << xv1 << "\t" << df0 << "\t" << df1 << endl;
            }
            else{ cout << "C. Troposphere limit" << xv0 << "\t" << xv1 << "\t" << df0 << "\t" << df1 << endl; }
            exit(2); 
        }
    }
    // cout << setprecision(6) << xv0 << "\t" << xv1 << "\t" << df0 << "\t" << df1 << endl;
    
    /* bisection search */
    double eps = 1e-3;
    while (log(xv1/xv0) > eps){
        /* new point -- half way */
        double xvA = (xv0+xv1)*0.5;
        double dfA = radflux_up(Fout, Psurf, xvA, rh);
        
        if (df0*dfA < 0){
            xv1 = xvA;
            df1 = dfA;
        }else{
            xv0 = xvA;
            df0 = dfA;
        }
        // cout << "* " << flush;
    }
    
    // cout << endl << setprecision(6) << xv0 << "\t" << xv1 << "\t" << df0 << "\t" << df1 << endl;

}
double atmos1d::water_saturation(double T){
    double P0 = 6.103e-3*1.013e5;
    double Ps  = P0 * exp( - L/R*(1/T - 1/273.15));
    Ps  = 1.4e11 * exp( - L/R*(1/T));
    
    return Ps;
}
tuple<double, double> atmos1d::tropopause(double Fout, double xv, double rh){
    /*
      calculate tau of tropopause for given Fout and xv (vapor mixing ratio)
    */
    double mu = muv*xv + (1-xv)*(muCO2*xCO2 + mun*(1-xCO2));
    
    /* bisection search for the tropopause temperature
       ...
       A. when the tropopause temperature is too low:
       -    (optical depth from radiation) < (that from absorbant amount)
       -    i.e., f0 < 0
       B. when the tropopause temperature is too high:
       -    (optical depth from radiation) > (that from absorbant amount)
       -    i.e., f1 > 0
       
       we want to find a search range so that A is met at T=T0, and B is met at T=T1
    */
    double T0 = sqrt(sqrt(0.5*Fout/sb)),T1 = 320;
    double f0 = delta_tau(T0, mu, xv, rh, Fout);
    double f1 = delta_tau(T1, mu, xv, rh, Fout);
    // cout << "trop init: " << setprecision(6) << f0 << "\t" << f1 << " @" << Fout << "\t" << xv << endl;
    while (f0>0){
        T0 -= 1.;  f0 = delta_tau(T0, mu, xv, rh, Fout);
    }
    while (0){ //f1<0){
        //cout << T1 << "\t" << f1 << "\t" << xv << "\t" << Fout << "\t";
        T1 -= 1.;  f1 = delta_tau(T1, mu, xv, rh, Fout);
        
        if (T1<T0){
            cout << "trop fnd : " << setprecision(6) << T1 << "\t" << f1 << " @" << Fout << "\t" << xv << endl;
            return forward_as_tuple(0., 0.);
        }
    }

    if (f1<0){
        /* to find the maximum f1, conduct a bisection search */
        double t0  = T0, t1 = T1;
        double ff0 = d_delta_tau(t0, mu, xv, rh, Fout);
        double ff1 = d_delta_tau(t1, mu, xv, rh, Fout);
        
        /* bisection search */
        double eps = 1e-2;
        double dt  = t1-t0;
        while (dt > eps){
            double tA  = (t0+t1)*0.5;
            double ffA = d_delta_tau(tA, mu, xv, rh, Fout);
            
            if (ff0*ffA<0){
                t1 = tA;   ff1 = ffA;
            }else{
                t0 = tA;   ff0 = ffA;
            }
            
            dt = t1-t0;
        }
        // cout << "d-trop end: " << t0 << "\t" << t1 << "\t" << ff0 << "\t" << ff1 << endl;
        T1 = t0;
        f1 = delta_tau(T1, mu, xv, rh, Fout);
        
        /* not enough absorbant in stratosphere */
        if (f1<0){
            // cout << "no tropo." << xv << "\t df = " << f0 << "\t" << f1 << "\t" << T0 << "\t" << T1 << endl;
            return forward_as_tuple(0., 0.);
        }
    }
    
    //cout << "trop end : " << setprecision(6) << f0 << "\t" << f1 << " @" << Fout << "\t" << xv << endl;
    
    double eps = 1e-2;
    double dT  = T1-T0;
    while (dT > eps){
        double TA = (T0+T1)*0.5;
        double fA = delta_tau(TA, mu, xv, rh, Fout);
        
        if (f0*fA<0){
            T1 = TA;   f1 = fA;
        }else{
            T0 = TA;   f0 = fA;
        }
        
        dT = T1-T0;
    }
    //cout << "tropo: " << setprecision(6) << T0 << "\t" << T1 << "\t" << f0 << "\t" << f1 << " @ " << Fout << " W" << endl;
    
    double T_tp    = (T0+T1)*0.5;
    double P_vp    = rh*water_saturation(T_tp);
    double tau_vp  = kpav  * muv  * P_vp;
    double tau_CO2 = kpaCO2* muCO2* P_vp*(1./xv-1) * xCO2;
    double tau_tp  = (tau_vp + tau_CO2) / (mu*gv);
    
    return forward_as_tuple(T_tp, tau_tp);
}
double atmos1d::delta_tau(double T, double& mu, double& xv, double& rh, double& Fout){
    /* equation 1:
       pi*B(tau) = 0.5*Fout*(1+1.5tau) */
    double tau_rad = (2*sb*(T*T*T*T)/Fout - 1)*2./3;
    
    /* equation 2:
       tau = (kv*muv*Pv + kCO2*muCO2*PCO2) / (mu*g) 
       
       Pv is humidity * saturation vapor. This eq. requires that there is enough absorbant. */
    double Pvp = rh*water_saturation(T);
    double Pnc = Pvp*(1./xv - 1.);
    double tau_sat = (kpav*muv*Pvp + kpaCO2*muCO2*(xCO2*Pnc)) / (mu * gv);
    //cout << tau_rad << "\t" << tau_sat << endl;
    
    return tau_rad - tau_sat;
}
double atmos1d::d_delta_tau(double T, double& mu, double& xv, double& rh, double& Fout){
    /* original equation:
       delta_tau = 2/3*(2*sb/Fout * T^4 - 1) - (kpav*muv + kpaC*muC*xC*(1/xv-1))*rH*Psat/(mu*g)
       
       derivative:
       d(delta_tau)/dT = 16/3*sb/Fout * T^3 - (kpav*muv + kpaC*muC*xC*(1/xv-1))/(mu*g)* rH*dPsat/dT
    */

    double dPsatdT = L/(R*T*T)*water_saturation(T);
    double ddtaudT = 16.*sb/(3*Fout) * cube(T) - (kpav*muv + kpaCO2*muCO2*xCO2*(1./xv-1))/(mu*gv) * (rh*dPsatdT);
    
    return ddtaudT;
}
double atmos1d::moist_adiabat(double Tnow, double Pnow){
    double Pvp  = water_saturation(Tnow);
    double Pnc  = (Pnow - Pvp)*(1.-xCO2);
    double PCO2 = (Pnow - Pvp)*xCO2;
    
    // average per 1 mole
    double Cp = (Cpn*Pnc + Cpn  *PCO2 + Cpv*Pvp)/Pnow;
    double mu = (mun*Pnc + muCO2*PCO2 + muv*Pvp)/Pnow;
    
    // mixing ratio
    double q = Pvp/(Pnc+PCO2);
    
    // equation of state
    double rho = Pnow/(R*Tnow)*mu;
    
    // moist lapse rate
    double Gamma = mu/(rho*Cp) * (1. + q*L/(R*Tnow))/(1. + q*L*L/(Cp*R*Tnow*Tnow));
    return Gamma;
}
double atmos1d::radflux_up(double Fout, double Psurf, double xv, double rh){
    /* calculate thermal structure for given
       -- Fout (net solar flux absorbed)
       -- xv0  (water mixing ratio at the top)
       -- rh   (relative humidity at the tropopause)
       
       and determine the net outgoing flux
    */
    
    /* determine tropopause first */
    // double Ttp, tautp;
    tie(Ttp, tautp) = tropopause(Fout, xv, rh);
    if (Ttp ==0.0){
        /* Komabayashi-Ingersoll limit */
        return 0.0;
    }
    
    /* radiation regime */
    int    N   = T.size();
    int    sep = int(N/2);
    double mu  = muv*xv + (1-xv)*(mun*(1-xCO2) + muCO2*xCO2);
    
    for (int i=0; i<sep; i++){
        /* set up a tau-grid */
        double ds = 6;
        tau[i] = pow(10., -ds + (ds+log10(tautp))*i/(sep-1));
        
        /* calc pressure */
        P[i]  = (mu*gv) / (kpav*xv*muv + kpaCO2*(1-xv)*xCO2*muCO2) * tau[i];
        q[i]  = xv;
        
        /* radiative equilibrium */
        T[i] = sqrt(sqrt((0.5*Fout*(1.5*tau[i] + 1) / sb)));
        
        //cout << i << "\t" << T[i] << "\t" << tau[i] << "\t" << tautp << "\t" << P[i]*xv << "\t" << water_saturation(T[i]) << endl;
    }
    
    /* compare estimated tropopause and the amonut of non-condensible atmosphere */
    Ptp = (mu*gv) / (kpav*xv*muv + kpaCO2*(1-xv)*xCO2*muCO2) * tautp;
    double Pnc = Ptp - rh*water_saturation(T[sep-1]);
    if (Pnc > Psurf){
        // cout << Ptp << "\t" << P[99] << "\t" << Pnc << "\t" << tautp << "\t" << tau[99] << endl;
        return -Fout;
    }
    
    /* convective regime */
    double Th   =   T[sep-1];
    double Ph   =   P[sep-1];
    double tauh = tau[sep-1];
    
    for (int i=sep; i<N; i++){
        Pnc = Ph - rh*water_saturation(Th);
        if (xv<1.e-10){ cout << i << "\t -- " << Ph << "\t" << Pnc << endl; }
        if (Pnc<0){ return Fout; }  // return value here (op.2)
        
        double dP  = Pnc * (pow(Psurf/Pnc, 1.0/(N-i))-1.) * (Ph/Pnc);
        
        double div = 100, ddP = dP/div;
        for (int j=0; j<div; j++){
            /* moist lapse rate
               using 2nd Runge-Kutta */
            double Gamm0 = moist_adiabat(Th, Ph);
            double Thalf = Th + Gamm0*(ddP/2);
            
            double Gamma = moist_adiabat(Thalf, Ph+ddP/2);
            Th   += Gamma * ddP;
            Ph   += ddP;
            
            /* optical depth */
            double Psath = water_saturation(Th);
            double xvh   = rh*Psath/Ph;
            double mu    = muv*xv + (1-xv)*(mun*(1-xCO2)+muCO2*xCO2);
            tauh += (kpav*muv*xvh + kpaCO2*muCO2*(1-xvh)*xCO2) * ddP / (mu*gv);
            
            double Pnc = Ph-rh*Psath;
        }
        
        T[i] = Th;
        P[i] = Ph;
        tau[i] = tauh;
    }
    
    for (int i=0; i<N; i++){
        /* calc radiation from the troposphere */
        // Fup[i] = sb*square(square(T[N-1])) * exp(-1.5 * (tau[N-1]-tau[i]));
        Fup[i] = sb*cuatro(T[i]);
        for (int j=(i+1); j<N; j++){
            double trap0 = exp(-1.5*(tau[j  ]-tau[i]));
            double trap1 = exp(-1.5*(tau[j-1]-tau[i]));
            Fup[i] -= 0.5*(trap0+trap1) * sb * (cuatro(T[j-1]) - cuatro(T[j]));
            
            //double trap0 = exp(-1.5*(tau[j  ]-tau[i])) * sb * cuatro(T[j]);
            //double trap1 = exp(-1.5*(tau[j-1]-tau[i])) * sb * cuatro(T[j-1]);
            //Fup[i] += 0.5*(trap0+trap1) * (tau[j] - tau[j-1]);
        }
        
        /* water amount */
        double Ps = water_saturation(T[i]);
        
        if (i<sep){ q[i] = xv; }
        else{       q[i] = rh*Ps/P[i]; }
        
    }
    
    return Fup[0]-Fout;
}



/*
  for CO2-rich atmospheres, where atmoshpere can be divided into
  --
  1. stratosphere
  2. radiative  troposphere 
  3. convective troposphere

 */
tuple<double, double, double> atmos1d::minxv_for_tropopause(double Fout, double rh, double xv0, double xv1){
    /*
      find a minimum water mixing ratio (xv) that has tropopause
      ...
      bisection search for the presence of tropopause
      
      A. when the mixing ratio is too small:
      -    Ttp = 0
      B. when the mixing ratio is too large:
      -    Ttp > 0
    */

    double x0 = xv0, x1 = xv1;
    double Ttp0, Ttp1, TtpA, tau0, tau1, tauA;
    tie(Ttp0, tau0) = tropopause(Fout, x0, rh);
    tie(Ttp1, tau1) = tropopause(Fout, x1, rh);
    
    if (Ttp0!=0.0 or Ttp1==0.0){
        cout << "minxv_for: not properly bound: " << Ttp0 << "\t" << Ttp1 << endl;
        exit(2);
    }
    
    double eps = 1e-4;
    while (log(x1/x0)>eps){
        /* new point - bisection */
        double xA = (x0+x1)*0.5;
        tie(TtpA, tauA) = tropopause(Fout, xA, rh);
        //cout << x0 << "\t" << x1 << "\t" << Ttp0 << "\t" << Ttp1 << " -> " << TtpA << endl;

        if (TtpA==0.0){ x0 = xA;  Ttp0 = TtpA;  tau0 = tauA; }
        else          { x1 = xA;  Ttp1 = TtpA;  tau1 = tauA; }
        
    }

    /* check results */
    // cout << "check minxv" << endl;
    
    return forward_as_tuple((x0+x1)*0.5, Ttp1, tau1);
}
void atmos1d::calc_Tprofile_3layers(double Fout, double Psurf, double rh,
                                    double xv, double Ttp, double tautp){
    double mu  = muv*xv + (1-xv)*(mun*(1-xCO2) + muCO2*xCO2);
    
    /* search for correct boundary between radiative and convetive tropopause layers */
    double tau0 = tautp*1.001;
    double tau1 = (kpav*muv*xv + kpaCO2*muCO2*(1-xv)*xCO2) * (Psurf/(1-xv)) / (mu*gv*10);
    double df0 = radflux_up_3layers(tau0, Fout, Psurf, rh, xv, Ttp, tautp);
    double df1 = radflux_up_3layers(tau1, Fout, Psurf, rh, xv, Ttp, tautp);
    // cout << endl << "3l: " << setprecision(6) << tau0 << "\t" << tau1 << " / \t" << df0 << "\t" << df1 << flush;
    
    /* bisection search */
    double eps = 1e-3;
    while (log(tau1/tau0) > eps){
        /* new point -- half way */
        double tauA = (tau0+tau1)*0.5;
        double dfA = radflux_up_3layers(tauA, Fout, Psurf, rh, xv, Ttp, tautp);
        
        if (df0*dfA < 0 or dfA == Fout){
            tau1 = tauA;
            df1  = dfA;
        }else{
            tau0 = tauA;
            df0  = dfA;
        }
        // cout << tauA << "\t" << dfA << " \t <- " << df0 << "\t" << df1 << endl;
        // cout << " *" << flush;
    }
    
    // cout << setprecision(6) << tau0 << "\t" << tau1 << "\t" << df0 << "\t" << df1 << endl;

}
double atmos1d::radflux_up_3layers(double taub, double Fout, double Psurf, double rh,
                                   double xv,   double Ttp,  double tautp){
    /* calculate thermal structure for given
       -- taub (optical depth where the thermal profile of troposphere determines by moist lapse rate)
       
       -- Fout  (net solar flux absorbed)
       -- Psurf (non-condensible atmospheric pressure)
       -- rh    (relative humidity at the tropopause)
       -- xv    (water mixing ratio at the top)
       .   -- Ttp   &
       .   -- tautp : calculated in minxv_for_tropopause
       
       and determine the net outgoing flux
    */
    
    int    N   = T.size();
    int    sep = int(N/4), tro = int(N/2);
    double mu  = muv*xv + (1-xv)*(mun*(1-xCO2) + muCO2*xCO2);    
    
    /* 1. -----
       the structure of stratosphere:

       thermal structure is in radiative equilibrium, and the water mixing ratio is under-saturated.
       `q` is a constant, governved by the mixing ratio at the top of the troposphere.
    */
    for (int i=0; i<sep; i++){
        /* set up a tau-grid */
        double ds = 6;
        tau[i] = pow(10., -ds + (ds+log10(tautp))*i/(sep-1));
        
        /* calc pressure */
        P[i]  = (mu*gv) / (kpav*xv*muv + kpaCO2*(1-xv)*xCO2*muCO2) * tau[i];
        q[i]  = xv;
        
        /* radiative equilibrium */
        T[i] = sqrt(sqrt((0.5*Fout*(1.5*tau[i] + 1) / sb)));

        //cout << i << "\t" << T[i] << "\t" << tau[i] << "\t" << tautp << "\t" << P[i]*xv << "\t" << water_saturation(T[i]) << endl;
    }
    
    /* compare estimated tropopause and the amonut of non-condensible atmosphere */
    Ptp = (mu*gv) / (kpav*xv*muv + kpaCO2*(1-xv)*xCO2*muCO2) * tautp;
    double Pnc = Ptp - rh*water_saturation(T[sep-1]);
    if (Pnc > Psurf){
        // cout << Ptp << "\t" << P[99] << "\t" << Pnc << "\t" << tautp << "\t" << tau[99] << endl;
        return -Fout;
    }
    
    /* 2. -----
       the structure of radiative troposphere:

       thermal structure is in radiative equilibrium, but the water mixing ratio is always saturated.
    */
    double Ta   =   T[sep-1];
    double Pa   =   P[sep-1];
    double taus = tau[sep-1];
    
    for (int i=sep; i<tro; i++){
        /* raditaive equilibrium */
        tau[i] = taus * pow(taub/taus, (double)(i+1-sep)/(tro-sep));
        T[i]   = sqrt(sqrt(0.5*Fout*(1.5*tau[i] + 1)/sb));
        
        /* molar mass */
        double xvh = rh*water_saturation(T[i-1])/P[i-1];
        mu  = muv*xvh + (1-xvh)*(mun*(1-xCO2) + muCO2*xCO2);

        /* calc pressure */
        double dtau  = tau[i] - tau[i-1];
        double dPv   = rh*( water_saturation(T[i])-water_saturation(T[i-1]) );
        double dPCO2 = (dtau*(mu*gv) - kpav*muv*dPv)/(kpaCO2*muCO2);
        
        P[i] = P[i-1] + (dPCO2/xCO2) + dPv;

        /* non-condensible amount */
        // Pnc = P[i] - rh*water_saturation(T[i]);
        // cout << i << "\t" << tau[i] << "\t" << P[i] << "\t" << Pnc << "\t" << endl;
    }
    
    /* 3. -----
       the structure of convective troposphere:
    
       theermal structure is in convective lapse rate */
    double Th   =   T[tro-1];
    double Ph   =   P[tro-1];
    double tauh = tau[tro-1];
    
    /* compare tropopause and the amonut of non-condensible atmosphere */
    double Pb = (mu*gv) / (kpav*xv*muv + kpaCO2*(1-xv)*xCO2*muCO2) * tauh;
    Pnc = Pb - rh*water_saturation(Th);
    if (Pnc > Psurf){
        // cout << Ptp << "\t" << P[99] << "\t" << Pnc << "\t" << tautp << "\t" << tau[99] << endl;
        return -Fout;
    }

    for (int i=tro; i<N; i++){
        Pnc = Ph - rh*water_saturation(Th);
        if (Pnc<0){ return Fout; }  // return value here (op.2)
        
        double dP  = Pnc * (pow(Psurf/Pnc, 1.0/(N-i))-1.) * (Ph/Pnc);
        
        double div = 100, ddP = dP/div;
        for (int j=0; j<div; j++){
            /* moist lapse rate
               using 2nd Runge-Kutta */
            double Gamm0 = moist_adiabat(Th, Ph);
            double Thalf = Th + Gamm0*(ddP/2);
            
            double Gamma = moist_adiabat(Thalf, Ph+ddP/2);
            Th   += Gamma * ddP;
            Ph   += ddP;
            
            /* optical depth */
            double Psath = water_saturation(Th);
            double xvh   = rh*Psath/Ph;
            double mu    = muv*xv + (1-xv)*(mun*(1-xCO2)+muCO2*xCO2);
            tauh += (kpav*muv*xvh + kpaCO2*muCO2*(1-xvh)*xCO2) * ddP / (mu*gv);
            
            double Pnc = Ph-rh*Psath;
        }
        
        T[i] = Th;
        P[i] = Ph;
        tau[i] = tauh;
    }

    
    /* calc radiation from the troposphere */
    for (int i=0; i<N; i++){
        // Fup[i] = sb*square(square(T[N-1])) * exp(-1.5 * (tau[N-1]-tau[i]));
        Fup[i] = sb*cuatro(T[i]);
        for (int j=(i+1); j<N; j++){
            double trap0 = exp(-1.5*(tau[j  ]-tau[i]));
            double trap1 = exp(-1.5*(tau[j-1]-tau[i]));
            Fup[i] -= 0.5*(trap0+trap1) * sb * (cuatro(T[j-1]) - cuatro(T[j]));
            
            //double trap0 = exp(-1.5*(tau[j  ]-tau[i])) * sb * cuatro(T[j]);
            //double trap1 = exp(-1.5*(tau[j-1]-tau[i])) * sb * cuatro(T[j-1]);
            //Fup[i] += 0.5*(trap0+trap1) * (tau[j] - tau[j-1]);
        }
        
        /* water amount */
        double Ps = water_saturation(T[i]);
        
        if (i<sep){ q[i] = xv; }
        else{       q[i] = rh*Ps/P[i]; }

        // cout << i << "\t" << T[i] << "\t" << tau[i] << "\t" << P[i] << "\t" << Ps << "\t" << xv*P[i] << "\t" << q[i] << "\t" << Ps/P[i] << "\t" << Fup[i]-Fout << endl;
        
    }
    
    return Fup[0]-Fout;
}


/*
  to record output
 */
void atmos1d::record_result(string filename){
    ofstream fout(filename, ios::out);

    for (int i=0; i<(int)T.size(); i++){
        fout << P[i] << "\t" << tau[i] << "\t" << T[i] << "\t" << q[i] << endl;
    }

    fout.close();
}
