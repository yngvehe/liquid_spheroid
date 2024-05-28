	#
	#	Script for calculation of far-field pattern in the liquid case (sphere)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Physical parameters
	rho10 = 1.10 ; # Density ratio
	f = 200000
	c_0 = 1500
	c_1 = 1540
	k_0 = 2*pi*f/c_0  # Wave number in media 0
	k_1 = 2*pi*f/c_1  # Wave number in media 1
	a = 0.05 ; # Major semiaxis spheroid
	b = a / 5 ; # Minor semiaxis spheroid
	theta_inc = pi/2 ; # Incidence angle in radians
	theta_incdeg= round(Int,180*theta_inc/pi); 
	# Software parameters
	# M = 10 ;
	M = 8
	method = 2 ;
	delta_eta = 0.00390625 ;  # safe value: 1/2^8 
	
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	eta_inc = cos( theta_inc ) ;

# 	for i = 1 : M + 1 # 'm' loop
#         for j = i : M + 1 # 'n' loop
# 	        Smn0 = IncidentWave_Pro_Smn_vs_eta(a, b, i, j, f, method, eta_inc, delta_eta, c_0)
# 	        fileNameSmn = string("CRIMAC_Smn_a_",a,"_b_",b,"_m_",i-1,"_n_",j-1,"_k0_",k_0,"_inc_",theta_incdeg,".dat");
# 	        writedlm(fileNameSmn, Smn0 , '\t' ) ;
# 	    end
# 	end

# 	Rmn = IncidentWave_Pro_Rmn(a, b, M, f, c_0)
# 	fileNameRmn = string("CRIMAC_Rmn_a_",a,"_b_",b,"_M_",M,"_k0_",k_0,"_inc_",theta_incdeg,".dat");
# 	writedlm(fileNameRmn, Rmn , '\t' ) ;

	# Coefficient evaluation
	(Amn, Bmn) = Coeff_Pro( rho10, k_0, k_1, a, b, M, method, eta_inc ) ;
 	
	# Far-Field pattern calculation
	Pattern = Pattern_LiquidPro( k_0, a, b, M, method, eta_inc, Amn, delta_eta ) ;
	
	# Saving to disk
	fileName = string("CRIMAC_Pro_M_",M,"_a_",a,"_b_",b,"_k0_",k_0,"_k1_",k_1,"_rho1rho0_",rho10,"_inc_",theta_incdeg,".dat");
	writedlm(fileName, Pattern , '\t' ) ;
	
	