program yield_surface

    use globals !, only:  Nslips !, factor
    use utils
    use crystal_plasticity !use: yldf, gmdot

    implicit none

    integer     ::  i, np, ns,  cs, list_niter(10000,49), k, counter, a1, a2, n, m, INFO
    real(8)     ::  stress(6), dstrain(6), dt, tol, strainrate(3,3), T_list(10000,5), T11, T22, T33,  &
                    strial(6), f, grad(5), hess(5,5), C11, C12, C44, Cdiag(6), a, gm_sum,        &
                    algmod(6,6), factors(49), lamdot, Wp(3,3), Q0(3,3), vm_strain, deformation,       &
                    R(3,3), time, Q(3,3), L(3,3), phi1, phi2, PHI, L_crystal(3,3), D_plast(5), D(3,3),&
                    W(3,3), Dp(3,3), De(3,3), Wc(3,3), stop_criteria, stress_array(10), scale, xi(96),   &
                    scale1, scale2, array(200), summer, tot(96), phit, dphi(5), ddphi(5,5), initial, &
                    cs1, cs2, cs3, tmp_out(200), test(6)
    
    real(8), allocatable    ::  params(:), sdv(:), crss(:), gmdot(:), crss_string(:), yield_temp(:), tmp(:)

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
    params(4) = 10000     ! Exponent in the RSCYF - Regularized Single Crystal Yield Function
    params(5) = Nslips      ! Number of slip systems used
    params(6:np) = crss     ! Critical resolved shear stress (CRSS)
    ! initialize state-dependent variables
    ns = 4
    tot = 0.d0
    initial = 10.d0
    stress = 0.d0
    allocate(tmp(Nslips))
    allocate(gmdot(Nslips))

    !Getting a vector for scaling for CRSS
    call linspace(0.010d0,1.5d0,array)
    !Getting a vector for scaling of Sigma 12 and 31
    call linspace(-1.0d0,1.0d0,stress_array)

    open (81, file = 'Testytest.dat', status = 'REPLACE')
    write(*,*) test

    do a1=1,200 !Alpha 1 loop
        scale1 = array(a1)

        do a2=1,200 !alpha2 loop
            scale2 = array(a2)

            crss(13:24)   = initial * scale1
            crss(25:48)   = initial * scale2
            crss(61:72)   = initial * scale1
            crss(73:96)   = initial * scale2

            open (83, file='10000_5D.dat', status='unknown', action='read')
            !read(83,*)
            do n=1,10!10 !sigma13
                read(83,*) stress
                !stress(5) = stress_array(n)
                !do m=1,11 !sigma12
                stress(6) = stress(1)
                ! Yield functions 
                a = params(4)
                xi = 1.d0
                phi = 0.d0
                scale = 0.d0
                ! scaling by the maximum element
                tmp = 0.d0

                do i=1, Nslips
                    tmp(i) = max(0.d0, vdot(stress(2:), P(:,i))/CRSS(i))
                end do
                scale = maxval(tmp, dim=1)
                tmp = tmp/scale

                do i=1, Nslips
                    phi = phi + xi(i)*(tmp(i)**a)
                end do
                phi = phi**(1.d0/a)
                phi = scale*phi - 1.d0

                xi = 1.0d0
                gmdot = 0.d0
                lamdot = 1.0d0
                gm_sum = 0.d0
                
                !! Gamma prime functions
                do i=1, Nslips
                    gmdot(i) = lamdot*(xi(i)/CRSS(i)) * (tmp(i)**(a-1.d0))
                    if (isnan(gmdot(i))) gmdot(i) = 0.d0
                    gm_sum = gm_sum + gmdot(i)
                end do
                gm_sum = gm_sum/100.d0
                gmdot = gmdot/gm_sum
                !end do
            end do
            close (83)

            !Summation of the yield for eache point in the CRSS matrix
            cs1 = 0.d0
            cs2 = 0.d0
            cs3 = 0.d0
            
            !Summation for the different slip modes
            do i =1,12
                cs1 = cs1 + gmdot(i) + gmdot(i+48)
            end do
            
            do i =13,24
                cs2 = cs2 + gmdot(i) + gmdot(i+48)
            end do
    
            do i =25,48
                cs3 = cs3 + gmdot(i) + gmdot(i+48)
            end do

            ! Setting info for matrix for easy plotting
            if (cs1 .GT. 0.000000001d0) cs1 = 1.d0
            if (cs2 .GT. 0.000000001d0) cs2 = 3.d0
            if (cs3 .GT. 0.000000001d0) cs3 = 5.d0 

            tmp_out(a2) = (cs1 + cs2 + cs3)
        end do
        write(81,'(200F14.2)') tmp_out
        tmp_out = 0.d0
    end do
    close (81)

    write(*,*) '#####################################'
    !write(*,*)  gmdot

    write(*,*) '#####################################'
    !write(*,*)  tot


end program yield_surface 
