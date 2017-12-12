# Tectonics
These are the codes used to generate figures and analysis in Chapter 12 of Willgoose G R (2018) "Principles of Soilscape and Landscape Evolution", Cambridge University Press

This github repository include the sample landscape evolution codes that were used to create figures in Willgoose G R (2018) “ Principles of Soilscape and Landscape Evolution”, Cambridge University Press. You may download these codes and play with them, modify them and distribute them to others under only two conditions. These two conditions are that (i.e. these codes are distributed under the GPL licence)

(1) if you publish any results using them you provide a reference to my Cambridge University Press book and this github site, and 

(2) if you modify these codes you and distribute them to others that you (a) do not remove the existing header block at the top of the Python files that identify that Garry Willgoose was the original developer of the code (which includes a reference to the Cambridge University Press book), and (b) you add a subsequent header block that indicates that you have modified the code from the original (and summarises how you modified it). The first (a) condition ensures that I continue to get credit for the original work of developing the code, while the (b) condition means that I don’t get flooded with bug fix requests for bits of code that others have added to my original code.

As a courtesy it would be good that after you have enhanced the code you send it to me so that it can be added to this github site, so that all researchers may benefit from your enhancements.

Note that the solver in the code for plate-ice.py is identical to plate.py. The difference between plate-ice.py and plate.py is the application of the ice sheet loading on one half of the domain (i.e. plate-ice.py) versus centred around the middle of the domain (i.e. plate.py)

# plate.py

This code was used to generate Figures 12.6 and 12.7 in the book.There are no command line arguments. If you wish to change any of the parameter values you will need to change them in the global parameters at the top of the file. 

# plate-ice.py

This code was used to generate Figures 12.8 and 12.9 in the book. There are no command line arguments. If you wish to change any of the parameter values you will need to change them in the global parameters at the top of the file. 
