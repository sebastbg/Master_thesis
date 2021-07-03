program yield_surface_full

    use globals !, only:  Nslips !, factor
    use utils
    use crystal_plasticity !use: yldf, gmdot

    implicit none

    integer     ::  i, np, ns,  cs, list_niter(10000,49), k, counter, a1, a2, n, m, INFO, teller
    real(KIND =8) :: stress(6), dstrain(6), dt, tol, strainrate(3,3), T_list(1000,6), T11, T22, T33, &
                    strial(6), f, grad(5), hess(5,5), C11, C12, C44, Cdiag(6), a, gm_sum(96),        &
                    algmod(6,6), factors(49), lamdot, Wp(3,3), Q0(3,3), vm_strain, deformation,      &
                    R(3,3), time, Q(3,3), L(3,3), phi1, phi2, PHI, L_crystal(3,3), D_plast(5),       &
                    D(3,3), W(3,3), Dp(3,3), De(3,3), Wc(3,3), stop_criteria, stress_array(10),      &
                    scale1, scale2, array(100), gm_out(96), dphi(5), ddphi(5,5), initial,            &
                    cs1, cs2, cs3, tmp_out, T12, T13,  T23, T(3,3), input(6), out(100), gm_temp, &
                    cs10, cs20, cs30, crit, slip_out, dp_numb(1000), final_temp, final_out
    
    real(KIND=8), allocatable ::  params(:), sdv(:), crss(:), gmdot(:), crss_string(:), tmp(:)

    C11 = 233.30d3  ! MPa
    C12 = 135.50d3  ! MPa
    C44 = 118.00d3  ! MPa
    
    cs = 123
    call init_slipsystems(cs) ! this defines symmetric Schmid P vectors and number of slip systems NSLIP - both global variables
    call init_crss(cs, crss, scale1, scale2) ! Weighting for CRSS for non-octahedral slip. Documentation inside
    ! defines number of parameters in PARAMS array
    np = 5 + Nslips
    ! initialize parameters
    allocate(params(np))
    params(1) = C11         ! Elastic constant
    params(2) = C12         ! Elastic constant
    params(3) = C44         ! Elastic constant
    params(4) = 20         ! Exponent in the RSCYF - Regularized Single Crystal Yield Function
    params(5) = Nslips      ! Number of slip systems used
    params(6:np) = crss     ! Critical resolved shear stress (CRSS)
    ! initialize state-dependent variables
    ns = 4
    gm_out = 0.d0
    initial = 10.d0
    allocate(tmp(Nslips))
    allocate(gmdot(Nslips))
    allocate(sdv(ns), source=0.d0)

    L = 0.d0
    L(1,1) = 1.0d0
    L(2,2) = -0.0d0
    L(3,3) = -1.0d0
    dt = 0.01d0 / mnorm(getsym(L))

    !Reading in the 5D files.
    open (83, file='1000_5D.dat', status='unknown', action='read')

    ! writting file to array for inner loop
    do i=1, 1000
        read(83,*) input
        T_list(i,:) = input
    end do
    close(83)

    write(*,*) 'n= ', params(4)
    write(*,*) '####### start ########'

    !Getting a vector for scaling for CRSS
    call linspace(0.010d0, 1.50d0,array)
    !Getting a vector for scaling of Sigma 12 and 31
    !call linspace(0.10d0,5.0d0,stress_array)

    ! File for writing results
    open (81, file = 'Numb_slip_n20_max.dat', status = 'REPLACE')
    !open (82, file = 'Testy_raw.dat', status = 'REPLACE')
    
    do a1=1,100! !Alpha 1 loop
        scale1 = array(a1)
        crss(13:24)   = initial * scale1
        crss(61:72)   = initial * scale1
        do a2=1,100! !alpha2 loop
            scale2 = array(a2)
            crss(25:48)   = initial * scale2
            crss(73:96)   = initial * scale2
            params(6:np)  = crss 

            gm_sum = 0.d0
            gm_temp = 0.d0
            tmp_out = 0.d0
            out(a2) = 0
            
            cs1 = 0.d0
            cs2 = 0.d0
            cs3 = 0.d0
            do n=1,100!0 
                !Reading in the strain data for the step
                T11 = 1.d0/6.d0*(sqrt(3.d0)+3.d0)*T_list(n,1) + 1.d0/6.d0*(sqrt(3.d0)-3.d0)*T_list(n,2)
                T22 = 1.d0/6.d0*(sqrt(3.d0)-3.d0)*T_list(n,1) + 1.d0/6.d0*(sqrt(3.d0)+3.d0)*T_list(n,2)
                T33 = -T11-T22
                T12 = sqrt(2.d0)/2.d0*T_list(n,5)
                T13 = sqrt(2.d0)/2.d0*T_list(n,4)
                T23 = sqrt(2.d0)/2.d0*T_list(n,3)
                T   = reshape( (/ T11, T12, T13, T12, T22, T23, T13, T23, T33 /), (/3,3/) )
                
                ! Update the returnmarp parameters for the step
                stress = 0.d0
                dstrain = matrix2vec(T)
                
                call returnmap(stress, dstrain, dt, params, np, sdv, ns, algmod, INFO, .FALSE., lamdot, D_plast)
                call calc_gmdot(stress(2:), crss, params(4), gmdot, lamdot)
                
                ! store the largest set of slip systems
                !write(*,*) gm_temp
                cs10 = 0.d0
                cs20 = 0.d0
                cs30 = 0.d0

                teller = 0
                gm_temp = 0.d0
                do i =1,96
                  gm_temp  = gm_temp  + gmdot(i)
                end do

                crit = 0.0d0
                crit = gm_temp * 0.01d0 
                do i =1 ,96
                  if (gmdot(i) .GT. crit) then
                     teller = teller +1;
                  end if
                end do
                !write(*,*) teller
                dp_numb(n) = teller;
                final_temp = final_temp + teller
            end do 
            write(*,*) maxval(dp_numb)
            final_out = final_temp / 1000.d0 
            final_temp = 0.d0     
            out(a2) = final_out
            final_out = 0.d0
        end do
        !Write to files 
        write(81,'(100F14.2)') out
        !write(82,*) 
        write(*,*) a1, '%'
    end do

    close (81)
    !close (82)
    write(*,*) '########### FINISHED ##############'

    deallocate(tmp)
    deallocate(gmdot)
    deallocate(sdv)

end program yield_surface_full
