# SIBERIA

While it is only used lightly in the book this site also has the latest source code for SIBERIA. The source code for SIBERIA is also available at a number of other places (e.g. www.csdms.org) but those versions are typically not the latest version. The intent here is to have the latest version of SIBERIA available here at all times. The terms of use of SIBERIA are different for the other codes on this site. You may download and use SIBERIA. You may not modify it, or redistribute it in any form to other users. If you wish to do some modifications you would like to share with others then send it to me and after quality assuring it I will add it to this site (either as the next version of SIBERIA or a standalone code separate from the main SIBERIA code).

The different usage conditions for SIBERIA are important because SIBERIA is heavily used in industrial applications and I need to be able to assure users that it correctly operates for all applications.

The code for SIBERIA is written in fortran 90 and while I’ve not tested it against all fortran 90 compilers experience is that it should work on all compilers. There are some OPENMP commands in the codes used for monte-carlo simulation but they haven’t been used for a long time so may or may not work these days. These days for monte-carlo simulation we find it easier and more flexible to use a Python script that calls SIBERIA to do monte-carlo simulation.

The code for SIBERIA has been highly optimised for performance (at the expense in many cases of readability). The code also includes many alternative solvers that can be selected by parameters. Read the manual carefully about the use of these parameters, but all of the solvers that we’ve found useful at some time are still in the codes as standalone subroutines, so when you are looking at the codes you will see many subroutines that have the same calling convention and very similar names. This is the reason why. It makes the code longer and the logic more difficult to follow but it means that we can always go back to and rerun past simulations without having to go back and find old source code.

COMPILATION

The user should look carefully at the makefile. The order of compilation of the various .f90 is important because the fortran modules in many files are importaed into other files, and the module headers need to be created before depdent modules can be compiled.

If you can't figure out the makefiule then a cruide fallback is to repeatedly run the command
fortran -c *.f90

until all the files finally compile. After then run the command
fortran -o SIBERIA *.o

OTHER FILES

siberia-revision-history.txt: has a brief history of revisions to the code since V8.00
siberia.setup: This file should be in the same directory as the executable of siberia as it provides run time options for how Siberia will execute.
