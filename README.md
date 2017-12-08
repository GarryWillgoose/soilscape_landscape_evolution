# Chemical Weathering

This directory contains the codes that I used for the chemical weathering simulations in Chapters 8 and 11. There are two codes. Note that versioning was for my benefit to keep track of what codes was used while I was writing the book and extending the codes as I realised more figures/simulations were required. The main difference between V3 (used for the two minerals simulations) and V4 (all other simulations) is that V4 has a more comprehensive set of processes, and the V4 has been cleaned up a bit more that V3.

1. weather_sens_v4.py: this is the single mineral code used for the majority of simulations. See details below of how to run this code.

2. weather_sens_v3_2minerals.py: This is the two minerals code used for the simulations in Figure 8.26. 

# Running weather_sens_v4

1. The simplest way to run this code is to run it at the command line to replicate the figures in the book (e.g. to replicate figure 8.4)

python weather_sens_v4.py fig8.4

If you look in the code what this command does is construct the command line arguments and run the simulation that were used to create the figures in the book.

2. The more flexible way to run this code is to set inside the code the command that is used by the simulation. If you look to the commented lines that set a value for NEWDIR you will see a range of command lines that can be easily commented out. Also just above this code is an explanation for some of the NEWDIR arguments. There are also global variables to allow you to change the depth of the soil, the values of parameters used in the sensitivity studies. Beyond that you'll need to fiddle with the code directly.

# Running weather_sens_v3_2minerals

This code has not been as thoroughly tested as the single mineral version. These code was developed in response to a question from Alex McBratney at the Pedometrics 2015 conference in Cordoba. It was used to generate Figure 8.26 but I can't guarantee that it will work perfectly for other applications. Maybe it will, maybe it won't. To generate Figure 8.26 run the command

python weather_sens_v3_2minerals.py

Note that unlike the single mineral version of the code there are no command line options.
