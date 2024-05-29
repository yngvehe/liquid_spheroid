#=
JUL_Crimac_test:
- Julia version: 
- Author: ynhe
- Date: 2023-04-24
=#
    using Printf ;
	using DelimitedFiles ;

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Prolate Spheroid coefficients (liquid)
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function IncidentWave_Pro( a::Any, b::Any, m::Any, n::Any, fmin::Any, fmax::Any, df::Any, method::Any, eta_inc::Any, xi::Any)

		#	Calculate the incident wave in terms of spheroidal wave functions.
		#
		#	input:      a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               m	: 'm' index
        #               n   : 'n' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#               eta_inc	: incidence angle
		#
		#	output:	Amn	: vector of coefficients
		#               Bmn	: vector of coefficients
		#

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Type conversion
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		convert( AbstractFloat, a ) ;
		convert( AbstractFloat, b ) ;
		convert( Int , m ) ;
		convert( Int , n ) ;
		convert( Int , method ) ;
		convert( AbstractFloat, eta_inc ) ;
		convert( AbstractFloat, fmin ) ;
		convert( AbstractFloat, fmax ) ;
		convert( AbstractFloat, df ) ;
		convert( AbstractFloat, xi ) ;

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		# Derived parameters

		d = 2*sqrt( a^2 - b^2 ) ;
		#xi_0 = 1/sqrt( 1 - (b/a)^2 ) ; # Prolate
		xi_0 = xi

        # todo: input parameter (or maybe f as input parameter)
		f_range = fmin:df:fmax
        Size = length(f_range)

		# Structure declaration

        P_inc = zeros(BigFloat, Size, 1)
        R1 = zeros(BigFloat, Size, 1)
		R2 = zeros(BigFloat, Size, 1)
		S1 = zeros( BigFloat, Size, 1)
		R3_0 = zeros( Complex{BigFloat}, Size, 1)
        Norm00 = zeros(BigFloat, Size, 1)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Calculation of incident wave at xi = xi_0, eta=eta_inc (Adelman-Gumerov-Duraiswami software)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        i = 1
        for f = f_range
			k = 2* pi * f / 1500
            c_0 = ( d/2 )*k
            pro_lambdamn_approx( c_0, m-1, n-1 ) ;
            run( `pro_sphwv_S.R.bat $c_0 $(m-1) $(n-1) $eta_inc $xi_0` )
            R_0 = ReadFileToArrayBF( "Out_R.dat", ',', 0 )
            S_0 = ReadFileToArrayBF( "Out_S.dat", ',', 0 )
			R3_00 = complex( R_0[6], R_0[12] )
			R_0 = R_0[6]
			S_0 = S_0[3]
            name = @sprintf("%08d_%03d_%03d", trunc(Int,c_0*1000), m-1, n-1 )
            pro_lambdamn_approx( c_0, m-1, n-1 )
            run(`pro_sphwv_dr.N.bat $c_0 $(m-1) $(n-1) $name`)
            Norm0 = ReadFileToArrayBF( "Out_S_N.dat", '\n', 1)
            P_inc[i] = S_0 * R_0 / Norm0
			R1[i] = R_0
			S1[i] = S_0
			R3_0[i] = R3_00
			Norm00[i] = Norm0
            i += 1
        end
		return [f_range P_inc R1 S1 R3_0 Norm00]
		end

		function IncidentWave_Pro_vs_xi( a::Any, b::Any, m::Any, n::Any, h::Any, ximin::Any, ximax::Any, dxi::Any, method::Any, eta_inc::Any )

            #	Calculate the incident wave in terms of spheroidal wave functions.
            #
            #	input:      a	: spheroid major semiaxis
            #               b	: spheroid minor semiaxis
            #               m	: 'm' index
            #               n   : 'n' index
            #               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
            #               eta_inc	: incidence angle
            #
            #	output:	Amn	: vector of coefficients
            #               Bmn	: vector of coefficients
            #

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # 	Type conversion
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            convert( AbstractFloat, a ) ;
            convert( AbstractFloat, b ) ;
            convert( Int , m ) ;
            convert( Int , n ) ;
            convert( Int , method ) ;
            convert( AbstractFloat, eta_inc ) ;
            convert( AbstractFloat, h ) ;
            convert( AbstractFloat, ximin ) ;
            convert( AbstractFloat, ximax ) ;
            convert( AbstractFloat, dxi ) ;

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # 	Preliminaries
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # Derived parameters

            #d = 2*sqrt( a^2 - b^2 ) ;

            # todo: input parameter (or maybe f as input parameter)
            xi_range = ximin:dxi:ximax
            Size = length(xi_range)

            # Structure declaration

            R1 = zeros(BigFloat, Size, 1)
            R2 = zeros(BigFloat, Size, 1)
            R1_p = zeros(BigFloat, Size, 1)
            R2_p = zeros(BigFloat, Size, 1)
            #
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # 	Calculation of radial spheroidal function of xi (Adelman-Gumerov-Duraiswami software)
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #k = 2* pi * f / 1500
            #c_0 = ( d/2 )*k
            c_0 = h
            i = 1
            for xi = xi_range
                pro_lambdamn_approx( c_0, m-1, n-1 ) ;
                run( `pro_sphwv_S.R.bat $c_0 $(m-1) $(n-1) $eta_inc $xi` )
                R_0 = ReadFileToArrayBF( "Out_R.dat", ',', 0 )
                R1p_0 = R_0[7]
                R2p_0 = R_0[13]
                R_2 = R_0[12]
                R_0 = R_0[6]
                name = @sprintf("%08d_%03d_%03d", trunc(Int,c_0*1000), m-1, n-1 )
                println("xi = ", xi)
                R1[i] = R_0
                R2[i] = R_2
                R1_p[i] = R1p_0
                R2_p[i] = R2p_0
#                 pro_lambdamn_approx( c_0, m-1, n-1 )
#                 run(`pro_sphwv_dr.N.bat $c_0 $(m-1) $(n-1) $name`)
#                 Norm0 = ReadFileToArrayBF( "Out_S_N.dat", '\n', 1)
#                 P_inc[i] = S_0 * R_0 / Norm0
#                 R1[i] = R_0
#                 S1[i] = S_0
#                 R3_0[i] = R3_00
#                 Norm00[i] = Norm0
                i += 1
            end
            return [xi_range R1 R1_p R2 R2_p]
            end

		function IncidentWave_Pro_Smn_vs_eta( a::Any, b::Any, m::Any, n::Any, f::Any, method::Any, eta_inc::Any, deta::Any, c::Any )

		#	Calculate the incident wave in terms of spheroidal wave functions.
		#
		#	input:      a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               m	: 'm' index
        #               n   : 'n' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#               eta_inc	: incidence angle
		#
		#	output:	Amn	: vector of coefficients
		#               Bmn	: vector of coefficients
		#

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Type conversion
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		convert( AbstractFloat, a ) ;
		convert( AbstractFloat, b ) ;
		convert( Int , m ) ;
		convert( Int , n ) ;
		convert( Int , method ) ;
		convert( AbstractFloat, eta_inc ) ;
		convert( AbstractFloat, deta) ;
		convert( AbstractFloat, f) ;
		convert( AbstractFloat, c ) ;

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		# Derived parameters

		d = 2*sqrt( a^2 - b^2 ) ;
		xi_0 = 1/sqrt( 1 - (b/a)^2 ) ; # Prolate


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Calculation of incident wave at xi = xi_0, eta=eta_inc (Adelman-Gumerov-Duraiswami software)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		k = 2* pi * f / c
		c_0 = ( d/2 )*k

		run(`pro_sphwv_Sgrid.Sinc.bat $c_0 $(m-1) $(n-1) $eta_inc $deta`);
		ETA = reverse( readdlm( "Out_S.dat",',','\n')[:,2] ) ;
		Smn = reverse( ReadFileToArrayBF2( "Out_S.dat", ',', :, 5 ) ) ;
		Smn_inc= ReadFileToArrayBF("Out_Sinc.dat", ',', 5 ) ;

		return [ETA Smn]
		end


# 		function call_pro_sphwv_S_R(c_0::Float64, m::Int64, n::Int64, eta_inc::Float64, xi_0::Float64)
# 			max_mem=4000
# 			prec=256
# 			n_dr=20
# 			dr_min=1.0e-100
# 			n_dr_neg=20
# 			dr_neg_min=1.0e-100
# 			n_c2k=20
# 			c2k_min=1.0e-100
#
# 			run( `pro_sphwv.exe -max_memory $max_mem -precision $prec -verbose n -c $c_0 -m $m -n $ -w lambda`)
# 		end

        function IncidentWave_Pro_Rmn( a::Any, b::Any, M::Any, f::Any, c::Any )

        convert( AbstractFloat, a ) ;
        convert( AbstractFloat, b ) ;
        convert( Int , M ) ;
        convert( AbstractFloat, f) ;
        convert( AbstractFloat, c ) ;

        Size = round( Int64, (M+1)*(M+2)/2 ) ;
        d = 2*sqrt( a^2 - b^2 ) ;
        xi_0 = 1/sqrt( 1 - (b/a)^2 ) ; # Prolate
        k = 2* pi * f / c
        c_0 = ( d/2 )*k

        R1_0 = zeros( BigFloat, Size, 1) ;
        R1p_0 = zeros( BigFloat, Size, 1) ;
        R3_0 = zeros( Complex{BigFloat}, Size, 1) ;
        R3p_0 = zeros( Complex{BigFloat}, Size, 1) ;

        for i = 1 : M + 1 # 'm' loop
            for j = i : M + 1 # 'n' loop
                indice = Index( M, i, j) ;
                pro_lambdamn_approx( c_0, i-1, j-1 ) ;
                run( `pro_sphwv_S.R.bat $c_0 $(i-1) $(j-1) $eta_inc $xi_0` );
                R_0 = ReadFileToArrayBF( "Out_R.dat", ',', 0 )  ;
                # method 2
                R1_0[indice] = R_0[6] ;
                R1p_0[indice] = R_0[7] ;
                R3_0[indice] = complex( R_0[6], R_0[12] ) ;
                R3p_0[indice] = complex( R_0[7], R_0[13] ) ;
            end
        end

        return [R1_0 R1p_0 R3_0 R3p_0]
        end