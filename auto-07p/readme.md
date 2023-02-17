# Bifurcation diagram of the self-coupled Jansen system with `auto-07p`

The file `jrnmm.f90` contains the Jansen-Rit equations and the Jacobian.
Also, we have made used of "parameter overspecificaction" to compute the extrema of
periodic solutions using $y_1-y_2$ (you do not need to understand this).

The `c.jrnmm` contains the auto constants to run the program.

Scritpts `bifur_eps0.auto.py` and `bifur_eps50.auto.py` contain the necessary calls
to obtain the 1-parameter bifurcation diagrams displayied in Fig. 1(a) and 1(b) of the paper
for $\epsilon=0$ (standard Jansen-Rit model) and $\epsilon=50$.

If you have `auto-07p` installed, just go to the interactive Python terminal and run the scripts.

In order to plot the results I suggest using the external tool `gnuplot`, which is more updated
and produce nicer results. 
Commands to plot the diagrams are provided in the commends of the `.auto.py` scripts.

Finally, `bifur_eps50.auto.py` also contains a script to export the initial
conditions to run the Master Stability Function (`floquet.c`).



