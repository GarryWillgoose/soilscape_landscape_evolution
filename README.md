# soilscape_landscape_evolution
Codes used to generate figures and analysis in Willgoose G R (2018) "Principles of Soilscape and Landscape Evolution", Cambridge University Press

This github repository include the sample landscape evolution codes that were used to create figures in Willgoose G R (2018) “ Principles of Soilscape and Landscape Evolution”, Cambridge University Press. You may download these codes and play with them, modify them and distribute them to others under only two conditions. These two conditions are that 

(1) if you publish any results using them you provide a reference to my Cambridge University Press book and this github site, and 

(2) if you modify these codes you and distribute them to others that you (a) do not remove the existing header block at the top of the Python files that identify that Garry Willgoose was the original developer of the code (which includes a reference to the Cambridge University Press book), and (b) you add a subsequent header block that indicates that you have modified the code from the original (and summarises how you modified it). The first (a) condition ensures that I continue to get credit for the original work of developing the code, while the (b) condition means that I don’t get flooded with bug fix requests for bits of code that others have added to my original code.

As a courtesy it would be good that after you have enhanced the code you send it to me so that it can be added to this github site, so that all researchers may benefit from your enhancements.

The organisation of this site roughly corresponds to the chapter organisation of the book. Each directory has information about how to use the codes in that chapter. 

One note about the codes themselves. The codes are not particularly well optimised (either in appearance or performance). My focus was to make the codes as simple to understand as practical, while minimising the amount of time I spent generating the figures in the book. 

PYTHON

Almost all the codes provided are written in python (I’ve tested them against V2.7 and they won’t work on Python V3). Linux and OSX machines come with a python implementation pre-installed but it is missing many useful science library extensions so I suggest you install your own python implementation. There are a number of “science-heavy” python distributions out there. I have used the one from www.enthought.com called Canopy. The Canopy python is available for free for people associated with a university.

If you wish to use another python distribution or if you wish to construct your own science python using the public domain python (www. Python.org) or the one already installed on your computer you will need to add the following (free) extension libraries

numpy
matplotlib

SIBERIA

While it is only used lightly in the book this site also has the latest source code for SIBERIA. The source code for SIBERIA is also available at a number of other places (e.g. www.csdms.org) but those versions are typically not the latest version. The intent here is to have the latest version of SIBERIA available here at all times. The terms of use of SIBERIA are different for the other codes on this site. You may download and use SIBERIA. You may not modify it, or redistribute it in any form to other users. If you wish to do some modifications you would like to share with others then send it to me and after quality assuring it I will add it to this site (either as the next version of SIBERIA or a standalone code separate from the main SIBERIA code).

The different usage conditions for SIBERIA are important because SIBERIA is heavily used in industrial applications and I need to be able to assure users that it correctly operates for all applications.

The code for SIBERIA is written in fortran 90 and while I’ve not tested it against all fortran 90 compilers experience is that it should work on all compilers. There are some OPENMP commands in the codes used for monte-carlo simulation but they haven’t been used for a long time so may or may not work these days. These days for monte-carlo simulation we find it easier and more flexible to use a Python script that calls SIBERIA to do monte-carlo simulation.

The code for SIBERIA has been highly optimised for performance (at the expense in many cases of readability). The code also includes many alternative solvers that can be selected by parameters. Read the manual carefully about the use of these parameters, but all of the solvers that we’ve found useful at some time are still in the codes as standalone subroutines, so when you are looking at the codes you will see many subroutines that have the same calling convention and very similar names. This is the reason why. It makes the code longer and the logic more difficult to follow but it means that we can always go back to and rerun past simulations without having to go back and find old source code.
