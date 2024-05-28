	#
	#	Script for calculation of far-field pattern in the liquid case (sphere)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Physical parameters
	a = 0.1 ; # Major semiaxis spheroid
	b = a/10 ; # Minor semiaxis spheroid
	theta_inc = 0 ; # Incidence angle in radians
	theta_incdeg= round(Int,180*theta_inc/pi); 
	# Software parameters
	M = 8 ;
	method = 2 ;
	delta_eta = 2.0 ; # tuve que hacer lo de poner un paso como 1/2^8, de lo contrario da NaN

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	c = 1500 # speed of sound in water
	eta_inc = cos( theta_inc ) ; 
	Pattern_all = [;;]
	# Coefficient evaluation - replicate results from Andronov-paper
	kb_range = range(start=60, stop=80, step=0.5)
 	for kb = kb_range
 	    k_0 = kb / a  # (major semiaxis is called b in the Andronov paper)

        println("k_0", k_0)
	    AmnD=Coeff_ProDN( k_0, a, b, M, method, eta_inc, "D")
        println("After coeff")
	    # Far-Field pattern calculation
 	    Pattern = Pattern_LiquidPro( k_0, a, b, M, method, eta_inc, AmnD, delta_eta ) ;
 	    kb_vec = fill(kb, size(Pattern)[1])
 	    f_vec = fill(k_0 * c /(2 * pi), size(Pattern)[1])
 	    Pattern2 = [kb_vec f_vec Pattern]
 	    Pattern_pi = Pattern2[2,:] # the second row is back-scattering in the pi-direction
 	    if isempty(Pattern_all)
 	        global Pattern_all = Pattern_pi[:,:]'
 	    else
 	        global Pattern_all = [Pattern_all; Pattern_pi[:,:]']
 	    end
 	    println("Pattern_all", Pattern_all)
    end
    fileName = string("Pro_backscatter_a_",a,"_b_",b,"_inc_",theta_incdeg,"_kbmin_", kb_range[1], "_kbmax_", kb_range[end], ".dat");
    writedlm(fileName, Pattern_all , '\t' ) ;
	# Saving to disk
 	#fileName = string("Pro_a_",a,"_b_",b,"_k0_",k_0,"_k1_",k_1,"_rho1rho0_",rho10,"_inc_",theta_incdeg,"Dirichlet.dat");
 	#writedlm(fileName, Pattern , '\t' ) ;
	


	
	
	
	
	
	
	
	