# Initializing and exploring one-parameter line bifurcations

init=run('jrnmm',IPS=-2,NMX=100000)
ic=init(201)
fp=run(ic,IPS=1,NMX=10000,ISW=1)
save(fp,'fp')

# Limit cycle:

hb=init('HB1')
lc=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=1000,ISW=1)
save(lc,'lc')

