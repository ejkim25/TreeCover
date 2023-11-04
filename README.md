
This code is for the simulation of the models as part of the master's thesis with title
"Multimodal Tropical Tree Cover States in Monostable Dynamical Systems - Reconsidering theory and observations"
at Technical University of Munich

This code yields the results presented in Chapter 3:
3.1 Simulation of Hirota's and Good's models
3.2 Results from Variant models in a monostable system
    3.2.1 Conceptual monostable model
    3.2.2 Extended Good's model
3.3 Analysis of observed background variables

1) Hirota 2011:
Simulation of the paper Hirota 2011.
"Global Resilience of Tropical Forest and Savanna to Critical Transitions" (Hirota et al., 2011)

2) Good 2016:
Simulation of the paper Good 2016.
"Are strong fire–vegetation feedbacks needed to explain the spatial distribution of tropical tree cover?" (Good et al., 2016)

3) Good 2016 model extension:  
Simulation for the model extension of the Good 2016.
It involves the bimodality check using different method instead of the "large scale" of productivity and mortality of Good 2016:
Precipitation in function of Tree cover -> Productivity in function of Trecipitation -> Tree cover in function of Productivity
The code contains two different sections. Thus, it's needed to consider them separately:
  1. Fixed P0 and k (manually varying P0 and k)
  2. Multiple P0 and k  

4) Good 2016 model extension with the mortality depending on the precipitation, temperature, and wood density
