module TaylorLin

    use globals, only:  Nslips !, factor
    use utils
    use crystal_plasticity !use: yldf, gmdot

    implicit none

    contains

    subroutine test(deformation, scale1, scale2, cs, counter, exp_n)
          integer     ::  i, j, np, ns, INFO, cs, list_niter(10000,49), k, counter
          real(8)     ::  stress(6), dstrain(6), dt, tol, strainrate(3,3), T_list(10000,5), T11, T22, T33,  &
                          T12, T23, T13, strial(6), f, grad(5), hess(5,5), C11, C12, C44, Cdiag(6),         &
                          algmod(6,6), factors(49), lamdot, Wp(3,3), Q0(3,3), vm_strain, deformation,       &
                          R(3,3), time, Q(3,3), L(3,3), phi1, phi2, PHI, L_crystal(3,3), D_plast(5), D(3,3),&
                          W(3,3), Dp(3,3), De(3,3), Wc(3,3), stop_criteria, defor_change, change_criteria,  &
                          scale1, scale2, exp_n, tmp_gm, gm_out(96)
          
          intent(in)  ::  deformation, scale1, scale2, cs, counter
          logical     ::  doprint, dochange, singel_crystal

          character(len=70) :: fn, filenum, fmt, file
          CHARACTER(len=10) :: file_id
          CHARACTER(len=30) :: file_name, file_name2
          CHARACTER(len=100):: file_info
          character(len=30) :: def_print, crss_print, scale_print
          integer, parameter :: numfiles=4
          integer, parameter :: outunit=4
          integer :: z !,filenum 

          real(8), allocatable    ::  params(:), sdv(:), crss(:), a(:), b(:), c(:), gmdot(:)
          ! Elastic constants
          !FCC C11 = 106.75d3 C12 = 60.41d3  C44 = 28.34d3
          !BCC C11 = 233.30d3 C12 = 135.50d3 C44 = 118.00d3
          C11 = 233.30d3  ! MPa
          C12 = 135.50d3  ! MPa
          C44 = 118.00d3  ! MPa
     
          ! Crystallographic slip systems
          ! Octahedral slip is given by cs = 1
          ! Non-octrahedral slip is given by FCC: cs=1234 and BCC: cs=123
          ! More documentations in crystal_fcc and crystal_bcc
          !cs = 123  ! activates 24 face-centered cubic slips systems
          call init_slipsystems(cs) ! this defines symmetric Schmid P vectors and number of slip systems NSLIP - both global variables

          call init_crss(cs, crss, scale1, scale2) ! Weighting for CRSS for non-octahedral slip. Documentation inside
          ! defines number of parameters in PARAMS array
          np = 5 + Nslips
          ! initialize parameters
          allocate(params(np))
          params(1) = C11         ! Elastic constant
          params(2) = C12         ! Elastic constant
          params(3) = C44         ! Elastic constant
          params(4) = exp_n       ! Exponent in the RSCYF - Regularized Single Crystal Yield Function
          params(5) = Nslips      ! Number of slip systems used
          params(6:np) = crss     ! Critical resolved shear stress (CRSS)
          ! initialize state-dependent variables
          ns = 4
          allocate(sdv(ns), source=0.d0)
          ! initialize critical resolved shear stress
          !allocate(crss(Nslips), source=10.d0)
          allocate(gmdot(Nslips))
          !allocate(gm_out(Nslips))


          !Velocity Deformation Gradient in Global refrence system
          L = 0.d0
          L(1,1) = 1.0d0
          L(2,2) = -0.00d0
          L(3,3) = -1.0d0
          !L(1,3) = 2.d0

          !Two step deformation.
          !The point of change is given in fraction deformation and the converted to strain
          dochange = .FALSE.
          !defor_change = 2.0

          !open (83, file = 'Result_Full.dat', status = 'REPLACE')
          
          !Outer loop for more efficent parameter search
          ! build filename -- i.dat

          !write(x1, 'I5.5') counter
          !fmt = '(I5.5)'
          !write(filenum, fmt) counter
          
          !write(fn,fmt='(i5.5,a)') filenum,'.dat'
          write(file_id, '(i5.5)') counter
          file_name = 'file'//trim(adjustl(file_id))//'.dat'
          file_name2 = 'file'//trim(adjustl(file_id))//'_mode.dat'
          !file_info = 'def:'//deformation//'_s1:'//scale1//'_s:'//scale2//'_cs:'//cs


          ! open it with a fixed unit number
          open(unit=11,file=trim(file_name), form='formatted')
          open(unit=12,file=trim(file_name2), form='formatted')
          write(11, *) 'load: uniaxial ', 'def:', deformation, 'n:', exp_n
          write(11, *) 's1:', scale1,'s2:', scale2, 'cs:', cs

          !Input Euler Angels for external file
          open (85, file= 'RAND1000.ori', status = 'unknown')
          !open (85, file= 'gamma.ori', status = 'unknown')
          read(85,*)
          read(85,*)
          read(85,*) k

          write(*,*) 'Filename:', file_name
          write(*,*) 'CS:', cs,'scale 112:',scale1,'scale 123', scale2
          write(*,*) 'Exponent:', exp_n, 'Number of grains: ', k

          !converting frac Deformation (0<deformation<1) to von Mises Strain for Rolling
          !deformation = 0.99999997d0
          stop_criteria = (2.d0 / sqrt(3.d0))*LOG(1.d0/(1.d0-deformation))
          !stop_criteria = 20.d0
          !
          write(*,*) 'Defor.:', deformation*100, '%, strain(vM):', stop_criteria

          do i=1,k
            !Euler angels in degrees from external file.
            read(85,*) phi1, PHI, phi2

            !calculating Transformation matrix(Q) and Rotation matrix(R) from Euler angeles
            Q = euler2rotm(phi1, PHI, phi2)
            R = transp3(Q)

            !Itteration variables for singel crystal itterations
            time = 0.d0
            dt = 0.01d0 / mnorm(getsym(L))
            vm_strain = 0.d0

            do while (vm_strain < stop_criteria)
              !If we want to change Deformation process:
              if (dochange .AND. vm_strain > change_criteria) then
                !Sett the new L matrix
                L = 0.d0
                L(1,1) = 1.0d0
                L(2,2) = -0.50d0
                L(3,3) = -0.50d0

                dt = 0.01d0 / mnorm(getsym(L))
                dochange = .FALSE.
              end if


              !Rotating the deformation gradient into crystall coordinates.
              L_crystal = transform(L, Q)
              D = getsym(L_crystal)
              W = getskw(L_crystal)

              ! strainrate vector
              dstrain = matrix2vec(D) * dt
              ! return-map subroutine, documentation inside
              call returnmap(stress, dstrain, dt, params, np, sdv, ns, algmod, INFO, .FALSE., lamdot, D_plast)
              ! Calculating the slip rates (eq. 11)
              call calc_gmdot(stress(2:), params(6:np), params(4), gmdot, lamdot)

              gm_out = gm_out + gmdot
              ! calculating deformation tensors (eq. 3)
              call calc_deformation(D, D_plast, gmdot,  Dp, De)
              ! Calculating plastic strain tensor (eq. 12)
              call calc_spinn(W, gmdot, Wp, Wc)
              ! Calculating the new Rotation Matrix (eq. )
              call update_R(R, Wc, dt)
              !write(*, '(96E14.3)') gm_out
              !write(12, '(96E14.3)') gm_out
              ! Calculating von Mises StrainRate for
              vm_strain = vm_strain + vonMises(Dp, 'strain')*dt
              ! Updating the itteration variables
              Q = transp3(R)
              time = time + dt
            end do
            gm_out = (gm_out*100) / (sum(gm_out))
            write(12, '(96E14.3)') gm_out
            !Output data into file
            write(11, '(4F14.6)') rotm2euler(Q)
            ! write(83, '(4F14.6)') rotm2euler(Q)
            !write(*,*) phi1, PHI, phi2, rotm2euler(Q)
            ! write(*,'(3F14.5)') phi1, PHI, phi2
            ! write(*,*) 'Q:'
          end do
          !Closing the files with results
          close(85)
          close(11)
          close(12)
          write(*,*) 

    end subroutine test


end module TaylorLin
