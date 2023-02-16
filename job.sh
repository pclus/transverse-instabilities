#!/bin/bash
#
#for i in {136..600..5};
#do 
#	./floquet $i 50.0 0.7393487564226234 >> transverse_Pyr2PyrInh_eps50.dat; 
#done

for p in {150,200,250};
do
	for i in $(seq -1.0 0.1 1.0);
	do
	        ./floquet $p 50.0 $i >> results/Pyr2PyrInh/disprel_p${p}_eps50.dat; 
	done
	
	for i in {1..90};
	do
		EIG=$(awk -v k=$i 'NR==k{print $1}' eigenvalues.dat )
	        ./floquet $p 50.0 $EIG >> results/Pyr2PyrInh/disprel_eigs_p${p}_eps50.dat; 
	done
done







