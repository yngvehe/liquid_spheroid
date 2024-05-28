	#
	#	Script for calculation of far-field pattern in the liquid case (sphere)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Physical parameters
	rho10 = 1.05 ; # Density ratio
	c_0 = 1500
	c_1 = c_0 * 1.05
	a = 0.08 ; # Major semiaxis spheroid
	b = a / 4 ; # Minor semiaxis spheroid
	theta_inc = pi/2 ; # Incidence angle in radians
	theta_incdeg= round(Int,180*theta_inc/pi); 
	# Software parameters
	# M = 10 ;
	M = 40
	method = 2 ;
	delta_eta = 0.00390625 ;  # safe value: 1/2^8 
	
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	eta_inc = cos( theta_inc ) ;

    #h = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
    h = [55]

    d = 2 * sqrt(a^2 - b^2)
    dir = "Smn_files"
    mkpath(dir)

    for hh in h
        f = (c_0 / pi) * (1/d) * hh
        i = 1  # 'm' = 0
        j = M + 1
        #for j = i : M + 1 # 'n' loop
            Smn0 = IncidentWave_Pro_Smn_vs_eta(a, b, i, j, f, method, eta_inc, delta_eta, c_0)
            fileNameSmn = string(dir, "/CRIMAC_Smn_a_",a,"_b_",b,"_m_",i-1,"_n_",j-1,"_h_",hh,"_inc_",theta_incdeg,".dat");
            writedlm(fileNameSmn, Smn0 , '\t' ) ;
        #end

    end
	
	