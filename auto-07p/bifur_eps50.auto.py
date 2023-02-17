# Initializing and exploring one-parameter line bifurcations from a fixed point:
init=run('jrnmm',IPS=-2,NMX=100000,PAR={'p':0.0, 'eps' : 50.0})                  # Use Euler integrator up to steady state to initialize the diagram
ic=init(201)                                                
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'p' : 500.0, 'p' : -300})     # Continue along p an uncover the bifurcations (2 saddle-node and 1 Hopf)
save(fp,'fp')                                                                   # Save results to a file

# Limit cycles:
# Continue the limit cycle solution (IPS=2) emerging from the Hopf.
# This branch dies on a SNIC bifurcation, characterized by a infinite period.
# Therefore, we stop the bifurcation at PERIOD = 100.0 (it can be continued further if needed)
hb=fp('HB1')                                                                                   
lc1=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=30000,ISW=1, UZSTOP={'PERIOD' : 100.0}) 
save(lc1,'lc1')

# The results can be plot in gnuplot with:
# plot 'b.fp' skip 20 u 5:($8-$9):($2<0) lc var, 'b.lc1' skip 20 u 5:11:(2*($2<0)) lc var, '' skip 20 u 5:12:(2*($2<0)) lc var

# In order to get initial conditions for the Master Stability Function we can iterate over the
# limit cycle and print the solution at the desired points:
filename = "msf_initial_conditions.dat"
fout=open(filename,"w");
icloop=hb;
for fl_p in range(329,85,-1):
    loop=run(icloop,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=10000,ISW=1,UZSTOP={'p':fl_p});
    fout.write("%lf %lf %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n"%(round( loop('UZ1')['p'] ), round( loop('UZ1')['eps'] ) , loop('UZ1')['PAR(11)'] , loop('UZ1')['y0'][0] ,loop('UZ1')['y1'][0],loop('UZ1')['y2'][0],loop('UZ1')['y3'][0],loop('UZ1')['y4'][0],loop('UZ1')['y5'][0]));
    icloop=loop('UZ1')


fout.close();




