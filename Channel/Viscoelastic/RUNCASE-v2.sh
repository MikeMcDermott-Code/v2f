# Define parameters

DR=(0 Re180-19 Re180-31 Re180-39 Re180-54
    Re300-33 Re300-47
    Re395-18 Re395-30 Re395-38 Re395-37 Re395-48
    Re590-39
    Re1000-30)
    
Wi=(0 25 50 100 100
    36 60
    25 50 50 100 100 100
    50
    50)
    
L=(0 900 900 900 3600 
   3600 3600
   900 900 3600 900 3600
   3600
   900)

lambda=(0 0.1388888 0.27777776 0.55555555 0.55555555 
        0.12 0.2
        0.06329 0.12658 0.12658 0.25316 0.25316
        0.08475
        0.05)

nus=(0 0.005 0.005 0.005 0.005
     0.003 0.003
     0.002278481 0.002278481 0.002278481 0.002278481 0.002278481
     0.001525424
     0.0009)

UNewt=(0 18.3 18.3 18.3 18.3
       19.5 19.5
       20.13 20.13 20.13 20.13 20.13
       21.3
       22)
# For loop running each case
#for ((i=1; i<=13; i++));
for i in 7;
do
	rm -r 0;
	mv 0-DR%/0-${DR["$i"]}% 0;
	sed -i "37s/.*/L2 ${L["$i"]};/" constant/RASProperties;
	sed -i "38s/.*/lambda ${lambda["$i"]};/" constant/RASProperties;
	sed -i "39s/.*/UNewt ${UNewt["$i"]};/" constant/RASProperties;
	sed -i "18s/.*/nu   nu [0 2 -1 0 0 0 0]   ${nus["$i"]};/" constant/transportProperties;
	TurbFoam >log;
	paraFoam -touch;
	paraview --script=DATA.py;
	var=$(ls -v | grep -E '^[0-9]+$' | tail -n 1);
	var2=$(cat $var/DR | sed -n '20p' | cut -c25-29);
	sed -i ""$i"s/.*/${DR["$i"]}%    $var2%    ${Wi["$i"]}      ${L["$i"]}      ${kappa["$i"]}    /" DR-Values;
	mv $var 0-DR%/0-${DR["$i"]}%; 
	unset var;
	mv DATA.csv Data/${DR["$i"]}.csv;
done
