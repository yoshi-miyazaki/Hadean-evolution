/*
  main.cpp
  ... 
  Tracks the coupled evolution of mantle and atmosphere,
  starting with a thick atmosphere of CO2.
  Plate velocity 
  
  by Yoshi Miyazaki (November 2019)
 */

#include "./system.h"

int main(){
    /* 
       consturct mantle-ocean-atmosphere coupled model 
    */

    /* ---------- [ input parameters ] ----------
       
       1. the surface (ocean+atmosphere) component 
       
       Ts  : surface temperature [K]  (Tref = 298.15K)
       atm : initial atmospheric pressure [atm]
    */
    double Ts  = 200 + (Tref-25);
    double atm = 177;
    
    /* 2. mantle component 
       
       Tm : initial mantle potential temperature
       D  : mantle depth       
     */
    double Tm         = 1600 + (Tref-25);   // initial potential temperature 
    double D          = 2900e3;             // mantle  depth               
    
    /* 2-2. mantle compontent (dehydration stiffening)

       homogeneous : whether mantle has a homogeneous or heterogeneous structure.
       -- if the mantle if homogeneous:
       .     the mantle produces depleted lithospheric mantle, which limits plate motion

       Deta_stiff : viscosicty contrast by dehydration
       -- which only matters when homogeneous is `true`
       .  Korenaga (2011) suggests 1000
    */
    bool   homogeneous = true;              // is mantle homongeneous?
    double Deta_stiff  = 1000.;
    
    /* 3. core component */
    double Tc = 2450;


    /* ---------- [ model ] ----------
       Do not change below unless necessary 
    */
    surface_mantle_system sys(atm, Ts, Tm, D, homogeneous, Tc);
    
    double dt = 5e6*yr2sec; 
    double time = 0.;
    for (int i=0; i<20000; i++){
        // sys.set_Tm(T0);
        sys.next(time, dt, Deta_stiff);
        time += dt;

        if (time/yr2sec > 5e9){
            dt = 2e6*yr2sec;
        }
    }        
    
}
