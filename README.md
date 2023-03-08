# Transverse-instabilities

This software provides the main tools to reproduce the results of the paper: 
[Clusella et al. , *Complex spatiotemporal oscillations emerge from transverse instabilities in large-scale brain networks* (2022)](https://doi.org/10.1101/2022.12.02.518809).
Cite this work if you use this, or parts of this, repository.


## Summary and dependencies

Mainly, this repository provides 3 tools that could be used independently:

- Integration of a network of coupled Jansen-Rit NMMs and computation of their largest Lyapunov exponents. This is provided by the `src/netdynamics.c` and `src/system.c` C codes.
- Scripts to compute 1-parameter bifurcation diagrams as those in Fig. 1(a-d) of the paper in `auto-07p`.
- Computation of the MSF of the coupled Jansen model, provided by the `src/floquet.c` code.

The C codes in this repository require the [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/index.html) (tested using version 2.7).
Auto-07p can be obtained from the [project github](https://github.com/auto-07p/auto-07p/releases).

## Integration of the Jansen-Rit network

The C-codes to simulate a network of Jansen-Rit elements can be compiled using `gcc`
with the provided `makefile`. From the `src` directory just run:

```
make
```

If the compilation is successful, the `netdyn` executable will be created in the main workspace.
This can be executed with:
`
./netdyn <output identifier> <file to network topology> <p> <epsilon> <number of Lyap. exp. to compute>
`

where:
- `<output identifier>` is a string to identify the output files (no file extension or path should be included).
- `<network topology>` is the path to the file containing the network topology. Notice that this needs to be in the appropiate format (see **Network data format** below).
- `<p>` is the parameter $p$ of the Jansen model (see paper).
- `<epsilon>` is the parameter $\epsilon$ of the Jansen model (see paper).
- `<nLE>` is the number of Lyapunov exponents to be computed, can be from 0 up to the number of system dimensions. Notice that computing Lyap. exp. increases the computation time.

An example call to the software could be:

```
./netdyn testing data/network.dat 230 50 2
```

The results can be plot, for instance, using [gnuplot](http://www.gnuplot.info/):

```
set xlabel 'Node index'
set ylabel 'time [s]'
set cblabel 'Activity [mV]'
plot 'outputs/pattern_testing.bin'  binary format='%double' array=(90,1000) with image

```

### Network data format

The network files must adhere to the following format:

- First line of the file should be the number of nodes (`N`).
- The following lines must contain three columns separated by spaces: `i j w`, where `A[i,j]=w` is the structural connectivity matrix of the system.
Indices should start from 0: `i,j=0,...,N-1`.


The file should then contain `L+1` lines, where `L` are the number of non-zero entries in `A`.
Notice that the simulation algorithm performs the row-normalization used in the paper by default.

## Bifurcation diagrams with `auto-07p`

The directory `auto-07p` contains the necessary files to compute 1-dimensional bifurcation
diagrams with [auto-07p](https://github.com/auto-07p/auto-07p).
The files `bifur_epsXX.auto.py` contain step-by-step explanations to obtain the diagrams using
the Python interface. The also contain the instruction to obtain initial conditions
to be used for computation of the Floquet exponents `floquet.c` (see **Master Stability Function** below).

## Master Stability Function

The standalone C code  `floquet.c` in `src` computes the largest Floquet exponent
of the self-coupled Jansen model.
The code can be compiled with:

```
gcc src/floquet.c -o floquet -lm -lgsl -lgslcblas -lblas -O3
```

The simulations are called with:
```
./floquet <p> <eps> <lambda> <ic file>
```

where 

- `p` is the $p$ parameter of the Jansen model.
- `eps` is the coupling strength.
- `lambda` is the structural connectivity eigenvalue that you want to use. To obtain the entire dispersion relation of the system one needs to loop over all the eigenvalues.
- `ic file` is the path to the initial conditions of limit cycle, as provided by the auto scripts.

Upon success the algorithm provides, in standard output 4 numbers: `<p> <eps> <lambda> <largest Floquet exponent>`.
For instance, a call to

```
./floquet 120 50 0.7 auto-07p/msf_initial_conditions.dat 
```
gives

```
230.000000 50.000000 0.700000 0.6125043913964573

```
## Structural Connectivity Data

The connectome used in the article is provided in `data/sc2017.dat`. 
The data is given in the format specified in the **Network data format** section above.
This structural connectivity corresponds to the average connectome
of 16 healthy subjects obtained from Diffusion Tensor Imaging and parcelled with the AAL90 atlas.  
For details about data collection and postprocessing users should refer to the
original publication by [Deco et al. (2017)](https://doi.org/10.1016/j.neuroimage.2017.03.023).
Any works using this data must cite that paper as well as [Clusella et. al (2022)](https://doi.org/10.1101/2022.12.02.518809).



