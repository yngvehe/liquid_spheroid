	#
	#	Script for calculation of far-field pattern in the liquid case (sphere)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Physical parameters
	a = 0.08 ; # Major semiaxis spheroid
	b = a / 4 ; # Minor semiaxis spheroid
	theta_inc = pi/2 ; # Incidence angle in radians
	theta_incdeg= round(Int,180*theta_inc/pi);
	# Software parameters
	m = 1
	n = [6, 11, 21, 31]

# 	fmin = 83000;
# 	fmax = 97500;
# 	df = 100;
# 	xi = 20;

# 	f = 200000  # corresponds to c = 41 for a=0.05, b=0.01
# 	k = 2* pi * f / 1500
# 	d = 2*sqrt( a^2 - b^2 )
#     c_0 = ( d/2 )*k
# 	ximin = 1/sqrt( 1 - (b/a)^2 ) ; # Prolate
# 	ximax = 10
# 	N = 20  # number of cycles per 2pi cycles
# 	dxi = 2 * pi / (N * c_0)

    h = [20.0, 40.0]

#     c_0 = 20
#     d = 2*sqrt( a^2 - b^2 )
#     f = (c_0 * 1500) / (pi * d)
    ximin = 1/sqrt( 1 - (b/a)^2 ) ; # Prolate
    ximax = 2
    N = 200
    dxi = (ximax - ximin) / N

	method = 2 ;

	#precision = 256 # default precision
	# comment out to use default precision for BigFloat
	#precision = 512
	#setprecision(precision)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	eta_inc = cos( theta_inc ) ;

    asp_ratio = a / b
    for hh in h
        for nn in n
            #Pattern = IncidentWave_Pro(a, b, m, n, fmin, fmax, df, method, eta_inc, xi)
            local Pattern = IncidentWave_Pro_vs_xi(a, b, m, nn, hh, ximin, ximax, dxi, method, eta_inc)

            #f_kHz = Int(fmin / 1000)

            # Saving to disk
            #fileName = string("Test_pro_f_",f_kHz,"_a_",a,"_b_",b,"_m_",m,"_n_",n,"_inc_",theta_incdeg,"_", precision, ".dat");
            fileName = string("Rmn_xi_c_",hh,"_asp_ratio_",asp_ratio,"_m_",(m-1),"_n_",(nn-1),".dat");
            writedlm(fileName, Pattern , '\t' ) ;
        end
    end

