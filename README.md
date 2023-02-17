# Transverse-instabilities

The C codes in this repository require the [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/index.html) (tested using version 2.7).

## Summary 

The `src` folder contains the source code:

- `netdynamics`: contains the main functions for monitoring the integration of the large-scale brain model. Can work with either `system` or `system_lyap_exp`.
- `system`: contains the main functions for integration of the large-scale brain model, including the Runge-Kutta algorithm, the velocity fields, and the algorithms to compute Lyapunov Exponents
- `floquet`: independent code to compute the Floquet exponents using the Master Stability Function formalism. 

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
