/*
  system.cpp
  
  ... see system.h for class `surface_mantle_system`

 */

#include "./system.h"

surface_mantle_system::surface_mantle_system(double atm, double Ts0, double Tm0, double D0, bool homo0,
                                             double Tc0){
    
    /* 1. Initialize atmosphere conditions
       
       ... natms/nmant/nbase vector has 4 components (0: H2O, 1: H2, 2: CO2, 3: CO)
       .   the initial atmosphere is assumed to be oxidized, only consisted of H2O and CO2.
       . 
       .   nbase corresponds to an atomsphere of 1 atm (1atm <--> 2.92e20 mole of H2O, 1.198e20 mole of CO2).
       .   the unit is in 10^18 mole. */
    
    int n_comp = 5;
    tensor1d<double> nbase(0.0, n_comp);
    nbase[0] = 292;   nbase[2] = 119.8;    // <--- Change here for different CO2/H2O ratio
    
    // initial condition
    natms = nbase*atm;
    nmant = nbase*0.0;  // the mantle initially contains no volatile
    
    /* 2. Initialize mantle conditions */
    D    = D0;       /* the depth of the mantle     */
    homo = homo0;    /* is the mantle homogeneous?  */
    
    // initial condition
    Ts = Ts0;        /* initial surface temperature */
    Tm = Tm0;        /* initial mantle  temperature */
    Tc = Tc0;        /* initial core    temperature */
}
void surface_mantle_system::next(double time, double dt, double Deta_stiff){
    /*
      Solve for the thermal evolutions of core, mantle, and atmosphere.
      The evolution is solved from the bottom to the top
      :
      1. Core evolution.
      Core temperature is cooled and does not consider radiogenic heating
    */
    double Wc, Qc, Cc;
    tie(Wc, Qc, Cc) = core_heat_flux(Tm, Tc);
    
    /* 
       2. evolution of mantle
       :
       Cm dTm/dt = H - Q + Qc
       H:  internal heat production (radioactive heat)
       Q:  surface heat loss
       Qc: core heat flux, calculated above

       vs: plate velocity
    */
    // surface heat loss
    double W, Q, vs;
    tie(W, Q, vs) = Nu_heat_scaling(Ts, Tm, Deta_stiff);
    
    // internal heating
    double H = radioactive_heating(time);
    
    // thermal evolution of the mantle  
    double Cm = 7.0e27;      // the heat capacity of solid Earth
    Tm += (H-Q+Qc)*dt/Cm;
    
    cout << time/yr2sec/1e6 << "\t" << Ts << "\t" << Tm << "\t";
    //cout << "vs [m/yr]: " << vs * yr2sec << "\t Q [TW]: " << Q/1e12 << "\t" << H/1e12 << "\t" << Qc/1e12 << endl;
    //cout << Tm-273.15 << "\t" << Q/1e12 << "\t" << vs*100*yr2sec << " \t";
    
    /* 3. 
       evolution of volatile, and solve for updated compostions
    */
    volatile_loss(dt, vs);
        
    /* Gibbs free energy minimization */
    MoleculeData_G   gibbs_list("./spec_atmos.txt", Ts, 1.0);
    gibbsminCG       min(gibbs_list, natms);
    tensor1d<double> nequil = min.getnbest();
    natms = nequil;

    // updated atmos composition
    for (int j=0; j<(int)natms.size(); j++){ cout << natms[j] << "\t"; }   cout << vs*yr2sec << endl;
    
}
double surface_mantle_system::surface_T(double W){
    /* constants */
    double S0 = 960;                 // solar radiation, assumed to be 30% weaker than present
    double kappa0 = 0.01;
    
    /* calculate opacity */
    double kappaR  = sqrt(kappa0*gv/3/Patm);
    double P_H2O   = natms[0]/170 * Patm * 0;
    double P_CO2   = natms[2]/170 * Patm;
    double tau     = 1.5*kappaR* (P_H2O+P_CO2) /gv;
    
    /* From Abe and Matsui (1985) */
    double F0   = 140;
    double Tnow = pow( (F0*(1+tau/2) + S0/4)/sb, 0.25);
    
    return Tnow;
}
void surface_mantle_system::volatile_loss(double dt, double vs){
    /*
      Calculate the amount of volatile removed from the atmosphere based on plate velocity.
     */
    
    /* reference state:
       use 1cm/yr (vs_mod) & 100Myr (renew_mod) for modern plate speed & plate age
       
       r_react: the ratio of seafloor reacted within `dt`. 
       (it is assumed that the entire seafloor is altered in one cycle)
    */
    double vs_mod = 0.01/yr2sec, renew_mod = 100e6*yr2sec;
    double r_react = dt/(renew_mod * vs_mod/vs);
    
    double S_ocean = 4*pi*rE*rE * 0.7;
    double d_react = 500;               // the depth of reaction
    double V_ocean = S_ocean*d_react;   // reaction volume

    double V_react = V_ocean * r_react;
    double M_react = V_react * rhom;
    
    /* calculate stoichiometry
       
       r_olivine: mass ratio of olivine included in seafloor
       .   (for homogeneous mantle, the seafloor is basaltic and thus contains smaller amount of olivine)
       M_react_olivine: the amount of olivine reacted within `dt`
     */
    double r_olivine = 0.60;
    if (homo){ r_olivine = 0.04; }
    double M_react_olivine   = M_react * r_olivine;
    // cout << "react ratio: " << r_react << "\t " << renew_mod * vs_mod/vs / yr2sec /1e6 << " [Myr]" << endl;
    
    /* The mamga ocean model predicts  Mg# = 95 for the early Earth
       reactions
       1. 3 Fe-olivine + 2H2O          <-->  2 magnesite + 2 H2 + 3 SiO2  --> produce_H2
       2. 3 Mg-olivine + SiO2 + 4 H2O  <-->  2 serpentine                 --> consume_Mgol_byH2
       
       3. 2 Mg-olivine + 2 H2O + CO2   <-->  serpentine + magnesite       --> consume_CO2
     */
    double Mgn = 0.95;
    if (homo){ Mgn = 0.88; }
    double mol_react_olivine = M_react_olivine / (140.693e-3*Mgn + 203.77e-3*(1-Mgn));
    double produce_H2        = mol_react_olivine * (1-Mgn) * 2.0/3;
    double consume_Mgol_byH2   = produce_H2 * 4.5;
    double consume_Mgol_bycarb = mol_react_olivine - consume_Mgol_byH2;
    double consume_CO2         = consume_Mgol_bycarb * 0.5;
    // double consume_H2O         = produce_H2 * 1.0;
    
    // cout << "react ol, CO2, H2O: " << mol_react_olivine/dt*yr2sec << "\t" << consume_CO2/dt*yr2sec << "\t" << consume_H2O/dt*yr2sec << endl;
    //cout << "comp: " << natms[0] << "\t" << natms[1] << "\t" << natms[2] << "\t" << natms[3] << endl;
    
    /* updated composition */
    double unit = 1e18;
    natms[0] -= consume_H2O/unit;
    natms[1] += produce_H2/unit;
    natms[2] -= consume_CO2/unit;

}
