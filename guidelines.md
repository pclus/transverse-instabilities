for i in {120..600..10};
do 
        ./floquet $i 10.0 0.7393487564226234 >> transverse_Pyr2PyrInh_eps10.dat; 
done


for i in {131..600..5};
do 
	./floquet $i 30.0 0.7393487564226234 >> transverse_Pyr2PyrInh_eps30.dat; 
done

for i in {126..600..5};
do 
	./floquet $i 20.0 0.7393487564226234 >> transverse_Pyr2PyrInh_eps20.dat; 
done


