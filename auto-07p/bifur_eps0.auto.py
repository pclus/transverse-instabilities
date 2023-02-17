# Initializing and exploring one-parameter line bifurcations from a fixed point:
init=run('jrnmm',IPS=-2,NMX=100000,PAR={'p':0.0, 'eps' : 0.0})                  # Use Euler integrator up to steady state to initialize the diagram
ic=init(201)                                                
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'p' : 500.0, 'p' : -200})     # Continue along p an uncover the bifurcations (2 saddle-node and 3 Hopf)
save(fp,'fp')                                                                   # Save results to a file

# Limit cycles:
# Continue the limit cycle solution (IPS=2) emerging from the first Hopf.
# This branch dies on a SNIC bifurcation, characterized by a infinite period.
# Therefore, we stop the bifurcation at PERIOD = 100.0 (it can be continued further if needed)
hb=fp('HB1')                                                                                   
lc1=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={'PERIOD' : 100.0}) 
save(lc1,'lc1')

# The second Hopf produces a limit cycle that vanishes at another Hopf (which auto labels as 'BP' when continuing limit-cycles)
hb=fp('HB2')
lc2=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, STOP=['BP1'])
save(lc2,'lc2')

# The results can be plot in gnuplot with:
# plot 'b.fp' skip 20 u 5:($8-$9):($2<0) lc var, 'b.lc1' skip 20 u 5:11:(2*($2<0)) lc var, '' skip 20 u 5:12:(2*($2<0)) lc var, 'b.lc2' skip 20 u 5:11:(2*($2<0)) lc var, '' skip 20 u 5:12:(2*($2<0)) lc var

