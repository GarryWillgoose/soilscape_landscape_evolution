# Cliff Retreat

These are the codes used to generate Figures 13.10 to 13.13 in Chapter 13 " High Slope Gravity Processes".

This github repository include the sample landscape evolution codes that were used to create figures in Willgoose G R (2018) “ Principles of Soilscape and Landscape Evolution”, Cambridge University Press. You may download these codes and play with them, modify them and distribute them to others under only two conditions. These two conditions are that (i.e. these codes are distributed under the GPL licence)

(1) if you publish any results using them you provide a reference to my Cambridge University Press book and this github site, and

(2) if you modify these codes you and distribute them to others that you (a) do not remove the existing header block at the top of the Python files that identify that Garry Willgoose was the original developer of the code (which includes a reference to the Cambridge University Press book), and (b) you add a subsequent header block that indicates that you have modified the code from the original (and summarises how you modified it). The first (a) condition ensures that I continue to get credit for the original work of developing the code, while the (b) condition means that I don’t get flooded with bug fix requests for bits of code that others have added to my original code.

As a courtesy it would be good that after you have enhanced the code you send it to me so that it can be added to this github site, so that all researchers may benefit from your enhancements.

# cliff-trans.py

This file does cliff retreat by transport limited erosion (Figure 13.11). There are no command line parameters so to generate the various figures you need to change the global constants at the top of the Python files. Note that the code is entirely unoptimised. THis is intentional so that the physics in the code and its solver is completely transparent. Note that about midway down the code you change the way slope is calculated to explore some of the effects of numerical approximations. To run the code

python cliff-trans.py

# cliff-detach2.py

This file does cliff retreat by detachment limited erosion (Figure 13.10, 13.12, 13.13). There are no command line parameters so to generate the various figures you need to change the global constants at the top of the Python files. Note that the code is entirely unoptimised. THis is intentional so that the physics in the code and its solver is completely transparent. Note that about midway down the code you change the way slope is calculated to explore some of the effects of numerical approximations. To run the code

python cliff-detach2.py
