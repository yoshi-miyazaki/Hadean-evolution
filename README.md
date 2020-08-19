# Hadean-evolution
Models the coupled evolution of the surface and interior during the Hadean.
This is the code to reproduce Figures 3b,c and 4 of 
"A new mode of geodynamics in the Hadean facilitates the emergence of early life"
by Miyazaki and Korenaga (2020)

-- Yoshi Miyazaki (yoshinori.miyazaki@yale.edu)

Makefile  
  - This file can be used to compile all source files.
    - CXX needs to be changed to your C++ compiler
    - Additional files (tensor.h, tensor_s.h, element.h/cpp, gibbse.h/cpp, massb.h/cpp, solution.h/cpp, and CGgibbsmin.h/cpp) are also necessary to compile this program (stored in my `GibbsE-minimization` repository), and they need to be stored in `../GibbsFE_minimization/`
  - To reproduce Figure 4,
    1. `make`
    2. `./early > filename`
  - To reproduce Figure 3, run the function `Nu_heat_scaling`

main.cpp
  - The main script to run the model. Initial conditions need to be set in this file.
    - mantle potential temperature (`Tm`)
    - atmospheric mass/composition (`natms`)
    - chemical heterogeneity of the mantle (`homo`)
  
system.cpp / system.h
  - This constructs `surface_mantle_system` class
    - contains information of the interior and surface
  - Functions to solve for the thermochemical evolution of the surface
  
mantle.cpp
  - This stores functions to solve for the thermal evolution of the mantle
  - Function to calculate plate velocity `Nu_heat_scaling`
    - based on heat flow scaling of Korenega (2010, JGR)
  
atmos1d.cpp / atmos1d.h
  - This constructs `atmos1d` object, which models the 1-D atmospheric strcture
    - The model assumes that the existence of ocean and the troposphere is saturated with water vapor
  - Input
    - The surface pressure of non-condensible species (`Psurf`)
    - The molar ratio of CO2 (`xCO2`)
    - The effect of other greenhouse gases can be included by calculating a equivalent `xCO2` value.
  - Output
    - `get_Tsurf()` returns the surface temperature

spec_atmos.txt
  - This file lists all the species considered in Gibbs energy minimization
    - Species can be added or removed, but their thermodynamic information needs to be included in `../GibbsFE_minimization/Gibbs_data/` to successfully run the model.
