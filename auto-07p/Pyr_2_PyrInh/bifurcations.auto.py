# Exploring one-parameter bifurcation diagram:

# In order to make things easy we choose eps=20 (larger eps causes a new lc to appear through a hopf and die to homoclinic)
# To get the lc, we first take eps=0 and large p, and then increase eps to find a Hopf. We stop at our desired eps=20

init=run('jrnmm',IPS=-2,NMX=100000,PAR={'p': 400.0, 'eps' : 0.0})   # Euler integration to steady state
ic=init(201)                                                        # Take final integration point   
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[2,1], UZSTOP={'eps' : 20.0})   # 1-param bifurcations for the steady state

# Limit cycle from a Hopf (along eps, fixed p=400)
hb=fp('HB1')                                                        # Hopf bifurcation
lc0=run(hb,IPS=2,ISP=2,ICP=[2,11,1,3,4],NMX=10000,ISW=1)             # Continuation of limit-cycle increasing eps until eps=20

# Limit cycle from user defined point (along p, fixed eps=20.0)
lc_ic=lc0('UZ1')
lc1=run(lc_ic,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=10000,ISW=1,DS="-",STOP=['LP1'])    # Continuation of limit-cycle through decreasing p
lc2=run(lc_ic,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=10000,ISW=1,UZSTOP={'p' : 600})     # Continuation of limit-cycle through increasing p
lc=merge(lc1+lc2)
save(lc,'lc')

# Steady states along p with fixed eps=20.0
init=run('jrnmm',IPS=-2,NMX=100000,PAR={'p': -100.0, 'eps' : 20.0})   
ic=init(201)                                                        
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2])   
save(fp,'fp')

# Plot the results in gnuplot:
# plot 'b.fp' skip 20 u 5:($8-$9):($2<0) lc var, 'b.lc' skip 20 u 5:11 lc 2, 'b.lc' skip 20 u 5:12 lc 2

# Notice that there is a bit of overlap between the fp and the lc close to the SN,
# this indicates that the onset of oscillations is through a simple homoclinic rather than a SNIC

