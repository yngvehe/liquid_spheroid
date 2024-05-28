	#
	#	Script for calculation of far-field pattern in the liquid case (sphere)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Physical parameters
	a = 0.05 ; # Major semiaxis spheroid
    b = 0.01 ; # Minor semiaxis spheroid
	theta_inc = pi/2 ; # Incidence angle in radians
	theta_incdeg= round(Int,180*theta_inc/pi);
	# Software parameters
	m = 4 ;
	n = 8 ;
	method = 2 ;
	delta_eta = 0.00390625 ;  # safe value: 1/2^8
    f = 200000

	#precision = 256 # default precision
	# comment out to use default precision for BigFloat
	#precision = 128
	#setprecision(precision)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	eta_inc = cos( theta_inc ) ;

    f_kHz = Int(f / 1000)

    #size_scales = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8]
    size_scales = [7.5]

    for scale in size_scales
        a_scaled = a * scale
        b_scaled = b * scale
        Pattern = IncidentWave_Pro_Smn_vs_eta(a_scaled, b_scaled, m, n, f, method, eta_inc, delta_eta)

        # Saving to disk
        fileName = string("Test_angular_pro_f_",f_kHz,"a_",a_scaled,"_b_",b_scaled,"_m_",m,"_n_",n,"_inc_",theta_incdeg,"_", precision, ".dat");
        writedlm(fileName, Pattern , '\t' ) ;
    end

