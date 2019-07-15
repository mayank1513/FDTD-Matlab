# Finite Difference Time Domain Analysis of Metal and Plasma-Metallic Photonic Crystals

Photonic Crystals also referred to as Photonic Band Gap (PBG) structures are engineered materials designed to exhibit desired optical properties. Their properties depend on their design and structure rather than the base materials. The metal PBG structures have higher power handling capabilities which can be dynamically tuned by introducing plasma as the background medium. These kind of tuneable photonic bandgap structures may find very important applications. The analysis of plasma metallic photonic crystal has been reported in two research papers communicated to Journal of Physics of Plasmas and Journal of Progress in Electrical and Electronics Research. The later reports the non-conventional grid structure and formalism for simulating plasma metallic photonic crystal under application of eternal magnetic field.

## About the codes
At first I apologize that due to lack of time I was unable to well comment the methods used herein.

This is a FDTD based simulator that simulates the evolution of electric and magnetic fields with time. In the codes above the visualization of fields is disabled to improve the runtime as the simulation is run for 10^12 to 10^14 time steps for each configuration of the photonic crystal reported (by varying radious and plasma frequency).

Depending on your system the simulatoin may run overnight or for two or three days.

### What's really happenning
Ok! so in brief, the computation domain is setup by assigning material properties accross computation domain (Physical space represented by matrices). In this implementation we have used point sources with gausian profile (to get the broadest spectrum) however in your case you may chose different sources or different profiles. 

Next, the engine of the simulation defined in time_march_PMPC_TE.m and time_march_PMPC_TM.m is activated to let the fields evolve over time.

We have set up probes in stratagic points to probe the field profile with respect to time which is later analysed to obtain the frequency response. It's good time to mention here that the coputaion domain is terminated with Bloch boundaries to simulate an infinite lattice. The simulation needs to run for each k-vector and in the post process we grt the resonant modes. We combine all these to generate an E-k diagram for our photonic crystal.

*This code automatically generates and saves E-K diagrams for different values or r/a ratio and plasma frequencies. This is a very valuable feature as it practically automates all your work. To speed up the process, simulations are deviced in three files Band_Diagram_main.m, Band_Diagram_main_1.m and Band_Diagram_main_2.m which you can run on separate machines and then combine the results.*

For details you may refer to my research papers quoted below. Thank you.

## What if I ant to visualize the fields as they are evolving.
- Go ahead and edit the time_march_PMPC_TE.m and time_march_PMPC_TM. Hopefully you just need to uncoment a few lines.

## How do I use these codes
Please download all the files in one folder and run Band_Diagram_main.m, Band_Diagram_main1.m, or Band_Diagram_main2.m to generate photonic bnd diagrams for different parameters. All figure and .mat files are generated and stored in systematic directories.

You can modifiy Band_Diagram_main.m file to obtained desired results.

Please Quote 

M. K. Chaudhari, “Non-conventional discretization grid based FDTD for EM wave propagation in magnetized plasma metallic photonic crystal,” Prog. Electromagn. Res. M, vol. 49, 2016.

M. K. Chaudhari and S. Chaudhari, “Tuning photonic bands in plasma metallic photonic crystals,” Phys. Plasmas, vol. 23, no. 11, p. 112118, 2016.

For any asssistance reach out to me at mayank.srmu@gmail.com
