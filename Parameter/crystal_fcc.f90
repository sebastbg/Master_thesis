module crystal_plasticity

    use globals
    use utils

    implicit none

    contains

    subroutine init_slipsystems(cs)
    ! Subroutine defines slip systems, build Schmid matrix,its symmetric part "P"
    ! and skew part "Omega", both global variables defined in module "globals"
    !
    ! Inputs: cs    -   specifies the set of slip systems used for the calculations as:
    !                   dictionary of possible crystalographic slip systems specified by "cs"
    !                              1              2               3               4
    !                   1,   '{111} <110>'
    !                   12,  '{111} <110>'  '{011} <011>'
    !                   13,  '{111} <110>'                  '{100} <011>'
    !                   14,  '{111} <110>'                                  '{112} <011>'
    !                   123, '{111} <110>'  '{011} <011>'   '{100} <011>'
    !                   124, '{111} <110>'  '{011} <011>'                   '{112} <011>'
    !                   134, '{111} <110>'                  '{100} <011>'   '{112} <011>'
    !                   1234,'{111} <110>'  '{011} <011>'   '{100} <011>'   '{112} <011>'


        real(kind=8)        :: n(36,3), b(36,3), schmid(3,3,36), schmid_sym(3,3,36), schmid_skw(3,3,36), &
                               schmid_sym_vec(5,36), schmid_skw_vec(3,36), tmp(6)
        integer             :: i, cs, s = 0
        real(kind=8), parameter     ::  inv3 = 1.d0/sqrt(3.d0),  &
                                        inv2 = 1.d0/sqrt(2.d0),  &
                                        inv6 = 1.d0/sqrt(6.d0)

        intent(in)          :: cs

        ! slip systems described by slip plane normal "n" and slip direction "b"
        ! slip plane normals "n"
        ! cs = 1 (FCC octahedral slips)
        n(1,:) = (/ 1.d0,  1.d0, -1.d0/)
        n(2,:) = (/ 1.d0,  1.d0, -1.d0/)
        n(3,:) = (/ 1.d0,  1.d0, -1.d0/)
        n(4,:) = (/ 1.d0, -1.d0, -1.d0/)
        n(5,:) = (/ 1.d0, -1.d0, -1.d0/)
        n(6,:) = (/ 1.d0, -1.d0, -1.d0/)
        n(7,:) = (/ 1.d0, -1.d0,  1.d0/)
        n(8,:) = (/ 1.d0, -1.d0,  1.d0/)
        n(9,:) = (/ 1.d0, -1.d0,  1.d0/)
        n(10,:)= (/ 1.d0,  1.d0,  1.d0/)
        n(11,:)= (/ 1.d0,  1.d0,  1.d0/)
        n(12,:)= (/ 1.d0,  1.d0,  1.d0/)
        ! cs = 2
        n(13,:)= (/ 1.d0,  1.d0,  0.d0/)
        n(14,:)= (/ 1.d0, -1.d0,  0.d0/)
        n(15,:)= (/ 1.d0,  0.d0,  1.d0/)
        n(16,:)= (/ 1.d0,  0.d0, -1.d0/)
        n(17,:)= (/ 0.d0,  1.d0,  1.d0/)
        n(18,:)= (/ 0.d0,  1.d0, -1.d0/)
        ! cs = 3
        n(19,:)= (/ 1.d0,  0.d0,  0.d0/)
        n(20,:)= (/ 1.d0,  0.d0,  0.d0/)
        n(21,:)= (/ 0.d0,  1.d0,  0.d0/)
        n(22,:)= (/ 0.d0,  1.d0,  0.d0/)
        n(23,:)= (/ 0.d0,  0.d0,  1.d0/)
        n(24,:)= (/ 0.d0,  0.d0,  1.d0/)
        ! cs = 4
        n(25,:)= (/-1.d0,  1.d0,  2.d0/)
        n(26,:)= (/-1.d0,  1.d0, -2.d0/)
        n(27,:)= (/ 1.d0,  1.d0,  2.d0/)
        n(28,:)= (/ 1.d0,  1.d0, -2.d0/)
        n(29,:)= (/ 1.d0,  2.d0, -1.d0/)
        n(30,:)= (/ 1.d0, -2.d0, -1.d0/)
        n(31,:)= (/ 1.d0,  2.d0,  1.d0/)
        n(32,:)= (/ 1.d0, -2.d0,  1.d0/)
        n(33,:)= (/ 2.d0,  1.d0, -1.d0/)
        n(34,:)= (/-2.d0,  1.d0, -1.d0/)
        n(35,:)= (/ 2.d0,  1.d0,  1.d0/)
        n(36,:)= (/-2.d0,  1.d0,  1.d0/)

        ! slip directions "b"
        ! cs = 1
        b(1,:) = (/ 0.d0,  1.d0,  1.d0/)
        b(2,:) = (/ 1.d0,  0.d0,  1.d0/)
        b(3,:) = (/ 1.d0, -1.d0,  0.d0/)
        b(4,:) = (/ 0.d0,  1.d0, -1.d0/)
        b(5,:) = (/ 1.d0,  0.d0,  1.d0/)
        b(6,:) = (/ 1.d0,  1.d0,  0.d0/)
        b(7,:) = (/ 0.d0,  1.d0,  1.d0/)
        b(8,:) = (/ 1.d0,  0.d0, -1.d0/)
        b(9,:) = (/ 1.d0,  1.d0,  0.d0/)
        b(10,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(11,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(12,:)= (/ 1.d0, -1.d0,  0.d0/)
        ! cs = 2
        b(13,:)= (/ 1.d0, -1.d0,  0.d0/)
        b(14,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(15,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(16,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(17,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(18,:)= (/ 0.d0,  1.d0,  1.d0/)
        ! cs = 3
        b(19,:)= (/ 0.d0,  1.d0,  1.d0/)
        b(20,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(21,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(22,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(23,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(24,:)= (/ 1.d0, -1.d0,  0.d0/)
        ! cas = 4
        b(25,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(26,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(27,:)= (/ 1.d0, -1.d0,  0.d0/)
        b(28,:)= (/ 1.d0, -1.d0,  0.d0/)
        b(29,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(30,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(31,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(32,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(33,:)= (/ 0.d0,  1.d0,  1.d0/)
        b(34,:)= (/ 0.d0,  1.d0,  1.d0/)
        b(35,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(36,:)= (/ 0.d0,  1.d0, -1.d0/)


        n(1:12,:)  = n(1:12,:) *inv3
        n(13:18,:) = n(13:18,:)*inv2
        n(25:36,:) = n(25:36,:)*inv6
        b = b*inv2

        ! 
        do i=1, 36
            schmid(:,:,i)       = outer3(b(i,:), n(i,:)) !finner 3x3 matrise
            schmid_sym(:,:,i)   = getsym(schmid(:,:,i))   !finner symetri
            tmp                 = matrix2vec(schmid_sym(:,:,i)) !
            schmid_sym_vec(:,i) = tmp(2:)
            schmid_skw(:,:,i)   = getskw(schmid(:,:,i))
            schmid_skw_vec(:,i) = (/ schmid_skw(2,3,i),                 &
                                     schmid_skw(1,3,i),                 &
                                     schmid_skw(1,2,i) /)
        end do

        if (allocated(P)) deallocate(P, Omega, STAT=s)
        if (s /= 0) write(*,*) 'Deallocation of P and Omega unsuccsessful.'

        select case (cs)
        case (1)
            Nslips = 24
            allocate(P(5,Nslips), Omega(3,Nslips), STAT=s)
            P(:,:12)     =  schmid_sym_vec(:,:12)
            P(:,13:)     = -schmid_sym_vec(:,:12)
            Omega(:,:12) =  schmid_skw_vec(:,:12)
            Omega(:,13:) = -schmid_skw_vec(:,:12)

        case (12)
            Nslips = 36
            allocate(P(5,Nslips), Omega(3,Nslips), STAT=s)
            P(:,:12)       =  schmid_sym_vec(:,:12)
            P(:,13:18)     =  schmid_sym_vec(:,13:18)
            P(:,19:30)     = -schmid_sym_vec(:,:12)
            P(:,31:36)     = -schmid_sym_vec(:,13:18)
            Omega(:,:12)   =  schmid_skw_vec(:,:12)
            Omega(:,13:18) =  schmid_skw_vec(:,13:18)
            Omega(:,19:30) = -schmid_skw_vec(:,:12)
            Omega(:,31:36) = -schmid_skw_vec(:,13:18)

        case (13)
            Nslips = 36
            allocate(P(5,Nslips), Omega(3,Nslips), STAT=s)
            P(:,:12)       =  schmid_sym_vec(:,:12)
            P(:,13:18)     =  schmid_sym_vec(:,19:24)
            P(:,19:30)     = -schmid_sym_vec(:,:12)
            P(:,31:36)     = -schmid_sym_vec(:,19:24)
            Omega(:,:12)   =  schmid_skw_vec(:,:12)
            Omega(:,13:18) =  schmid_skw_vec(:,19:24)
            Omega(:,19:30) = -schmid_skw_vec(:,:12)
            Omega(:,31:36) = -schmid_skw_vec(:,19:24)

        case (14)
            Nslips = 48
            allocate(P(5,Nslips), Omega(3,Nslips), STAT=s)
            P(:,:12)       =  schmid_sym_vec(:,:12)
            P(:,13:24)     =  schmid_sym_vec(:,25:36)
            P(:,25:36)     = -schmid_sym_vec(:,:12)
            P(:,37:48)     = -schmid_sym_vec(:,25:36)
            Omega(:,:12)   =  schmid_skw_vec(:,:12)
            Omega(:,13:24) =  schmid_skw_vec(:,25:36)
            Omega(:,25:36) = -schmid_skw_vec(:,:12)
            Omega(:,37:48) = -schmid_skw_vec(:,25:36)

        case (1234)
            Nslips = 72
            allocate(P(5,Nslips), Omega(3,Nslips), STAT=s)
            P(:,:36)       =  schmid_sym_vec
            P(:,37:)       = -schmid_sym_vec
            Omega(:,:36)   =  schmid_skw_vec
            Omega(:,37:72) = -schmid_skw_vec

        end select
        if (s /= 0) write(*,*) 'Allocation of P and Omega unsuccessful.'

    end subroutine init_slipsystems

    subroutine init_crss(cs, crss)
    ! Scaling of the CRSS for different slips systems
    ! Data collected form: https://doi.org/10.1007/s11661-020-05743-y
    ! Table 3: Threshold stresses of slip modes
    !   Alloy   Application   '{111} <110>'  '{011} <011>'   '{100} <011>'   '{112} <011>'
    !   6082    Rolling           1             x0.98           x1.5            x.1.4
    !   6082    Extrusion         1             x0.98           x1.5            x.1.4
    !   3104    Rolling           1             x0.75           x1.3            x.1.3
    !   3104    Plain strain      1             x0.75           x1.3            x.1.3

      use globals, only: Nslips
      !input
      integer         :: cs
      real(8)         :: initial, scale1, scale2, scale3
      real(8), allocatable  ::  crss(:)

      intent (in)     :: cs
      !Output
      intent (inout)  :: crss

      initial = 10.d0
      scale1 = 0.75d0
      scale2 = 1.30d0
      scale3 = 1.30d0

      select case (cs)
      case (1)
          allocate(crss(Nslips))
          crss(:12)  = initial
          crss(13:)  = initial

      case (12)
          allocate(crss(Nslips))
          crss(:24)     = initial
          crss(25:36)   = initial *scale1

      case (13)
          allocate(crss(Nslips))
          crss(:24)     = initial
          crss(25:36)   = initial * scale2

      case (14)
          allocate(crss(Nslips))
          crss(:24)     = initial
          crss(25:48)   = initial * scale3

      case (1234)
          allocate(crss(Nslips))
          crss(:24)     = initial
          crss(25:36)   = initial * scale1
          crss(37:48)   = initial * scale2
          crss(49:72)   = initial * scale3

      end select

    end subroutine init_crss

    subroutine yldf(S, CRSS, a, mode, phi, dphi, ddphi)

      use globals, only: Nslips, P
      ! input
      integer         ::  mode
      real(kind=8)    ::  S(5), CRSS(:), a
      ! output
      real(kind=8)    ::  phi, dphi(5), ddphi(5,5)
      !
      integer         ::  i
      real(kind=8)    ::  xi(size(CRSS)), tmp(size(CRSS)), scale

      intent(in)      ::  S, CRSS, mode, a
      intent(out)     ::  phi, dphi, ddphi

      xi = 1.d0
      phi = 0.d0
      ! scaling by the maximum element
      tmp = 0.d0
      do i=1, Nslips
          tmp(i) = max(0.d0, vdot(S, P(:,i))/CRSS(i))
      end do
      scale = maxval(tmp, dim=1)
      tmp = tmp/scale

      do i=1, Nslips
          phi = phi + xi(i)*(tmp(i)**a)
      end do
      phi = phi**(1.d0/a)
      phi = scale*phi - 1.d0

      if (mode .EQ. 0) return

      ! caluclate gradient
      tmp = 0.d0
      do i=1, Nslips
          tmp(i) = max(0.d0, vdot(S, P(:,i))/((phi+1.d0)*CRSS(i)))
      end do
      scale = maxval(tmp, dim=1)
      tmp = tmp/scale

      dphi = 0.d0
      do i=1, Nslips
          dphi = dphi + xi(i)* (tmp(i)**(a-1.d0)) *P(:,i)/CRSS(i)
      end do
      dphi = dphi*(scale**(a-1.d0))

      if (mode .EQ. 1) return

      ! calculate hessian
      ddphi = 0.d0
      do i=1, Nslips
          ddphi = ddphi + xi(i)*(tmp(i)**(a-2.d0)) * ((1.d0/CRSS(i))**2) * symouter5(P(:,i), P(:,i))
      end do
      ddphi = ddphi*(scale**(a-2.d0))
      ddphi = (a-1.d0)/(phi+1.d0) * (ddphi - symouter5(dphi, dphi))

      return
    end subroutine yldf

  !Calcualting slip rates. Gets values from returnmap function and inital constants
   subroutine calc_gmdot(S, CRSS, a, gmdot, lamdot)

      use globals, only: Nslips, P
      ! input
      real(kind=8)    ::  S(5), a, lamdot, xi(Nslips), CRSS(:), tmp(Nslips)
      ! output
      real(kind=8)    ::  gmdot(Nslips)
      !
      integer         ::  i

      intent(in)      ::  S, CRSS , a, lamdot
      intent(out)     ::  gmdot

      ! initial
      xi = 1.0d0
      gmdot = 0.d0
      tmp = 0.d0

      !eq. 11 in singel crystal document
      do i=1, Nslips
          tmp(i) = max(0.d0, vdot(S, P(:,i)) / CRSS(i))
      end do

      do i=1, Nslips
          gmdot(i) = lamdot*(xi(i)/CRSS(i)) * (tmp(i)**(a-1.d0))
      end do


    end subroutine calc_gmdot

    subroutine calc_deformation(D, D_plast, gmdot,  Dp, De)
      use globals, only: Nslips, P
      integer       :: i, k
      ! input
      real(kind=8)  :: D(3,3), D_plast(5), D_temp(6), gmdot(Nslips)
      !output
      real(kind=8)  ::  Dp(3,3), De(3,3)

      intent(in)    :: D, D_plast
      intent(out)   :: Dp, De

      D_temp = 0.d0

      !Converting the plastic deformation vector into the tensor
      D_temp(2:) = D_plast
      Dp = vec2mat(D_temp)

      !Calculating the Elastic deformation tensor
      De = D - Dp

    end subroutine calc_deformation

    subroutine calc_spinn(W, gmdot, Wp, Wc)
      !Using omega defined in globals
      use globals, only: Nslips, Omega
      integer       :: i, k
      ! input
      !real(kind=8)  ::  gmdot(24), W(3,3), Wp_vec(6), Wp_calc(3)
      real(kind=8)  ::  gmdot(Nslips), W(3,3), Wp_vec(6), Wp_calc(3)
      !output
      real(kind=8)  ::  Wp(3,3), Wc(3,3)

      intent(in)    ::  W, gmdot
      intent(out)   ::  Wp, Wc
      Wp_calc = 0.d0
      Wp_vec = 0.d0

      do i=1, Nslips
        Wp_calc = Wp_calc + (gmdot(i) * Omega(:,i))
      end do

      do i = 1, 3
        Wp_vec(i+3) = Wp_calc(i)
      end do

      Wp = vec2skewmat(Wp_vec)
      Wc = W - Wp

    end subroutine calc_spinn

    subroutine update_R(R, W_c, dt)
      use globals, only: Nslips
      integer       :: i, j

      !input
      real(kind=8)  ::  R_temp(3,3), K(3,3), W_c(3,3), spinn(3), A, dt, Z
      !input/output
      real(kind=8)  ::  R(3,3)

      intent(in)    ::  dt, W_c
      intent(inout) ::  R

      K = 0.d0
      Z = 0.d0

      do i=1, 3
        do j=1, 3
          Z = Z + (W_c(i,j))**2.d0
        end do
      end do

      R_temp = (eye(3)+(dt/(1.d0 +((dt**2.d0)/4.d0)*Z))* (W_c+(dt/2.d0)*(MATMUl(W_c, W_c))))
      R = MATMUL(R, R_temp)

     end subroutine update_R

     ! subroutine calc_vonMises(Dp, dt, vm_strain)
     !   use globals
     !   integer       :: i, j
     !   ! input
     !   real(kind=8)  :: Dp(3,3), dt, Z, vm_temp
     !   !output
     !   real(kind=8)  ::  vm_strain
     !
     !   intent(in)    :: Dp, dt
     !   intent(inout) :: vm_strain
     !   vm_temp = 0.d0
     !   Z = 0.d0
     !
     !   do i=1, 3
     !     do j=1, 3
     !       Z = Z + (Dp(i,j))**2.d0
     !     end do
     !   end do
     !
     !   vm_temp = (sqrt(Z * (2.d0/3.d0))* dt) + vm_strain
     !   vm_strain = vm_temp
     !
     ! end subroutine calc_vonMises

     subroutine returnmap(stress, dstrain, dt, params, np, sdv, ns,   &
                            algmod, INFO, doprint, lamdot, Dp)
!**********************************************************************
!      This is the return mapping predictor-corrector algorithm for
!      solving the closest-point projection in an elastic-plastic
!      problem. It is based on a Newton-Raphson method with
!      line-search. This algorithm is based on the paper by
!      Scherzinger, W. M.: A return mapping algorithm for isotropic
!      and anisotropic plasticity models using
!      a line search method, Computer Methods in Applied Mechanics
!      and Engineering, 2017, DOI: 10.1016/j.cma.2016.11.026
!      This particular implementation is for isotropic elasticity
!      and perfect plasticity.
!
!**********************************************************************
! (IN/OUT) STRESS is REAL*8 array, dimension (6)
!          It is Cauchy stress at the beginning of the time
!          increment expressed in new notation. At return it is
!          updated by the subroutine to be the stress at the end of the
!          time increment.
! (IN)     DSTRAIN is REAL*8 array, dimension (6)
!          It is the total strain increment expressed in new notation.
! (IN)     DT is REAL*8
!          It is a time increment
! (IN)     PARAMS is REAL*8 array, dimension (NP)
!          It contains user-defined necessary material parameters
!          as e.g. elastic constants, plastic anisotropy coefficients...
! (IN)     NP is INTEGER
!          It defines the length of the PARAMS array.
! (IN)     SDV is REAL*8 array, dimension (NS)
!          It contains the solution-dependent state variables as e.g.
!          equivalent plastic strain, hardening model variables...
! (IN)     NS is INTEGER
!          It defines the length of the SDV array.
! (OUT)    ALGMOD is REAL*8 array, dimension (6,6)
!          It is the algorithmic modulus defined as dstress/dstrain.
! (OUT)    INFO is INTEGER
!          INFO = 0 means that the algorithm finished successfully with
!          use of the Newton-Raphson steps only, i.e. no need for
!          line-seach
!          INFO = 1 means that the algorithm finished successfully and
!          the line-search was activated at least once
!          INFO = 2 means that the algorithm finished successfully
!          within IMAX iterations but there was at least one
!          iteration in which line-search was terminated due to
!          reaching JMAX line-search iterations
!          INFO = 3 means that there was no need for plastic corrector
!          since the yield function at STRIAL <= 0.d0
!          INFO = -1 means unsuccesfull finish, i.e. the residual PSI
!          did not get below TOL within IMAX iterations
! (IN)     DOPRINT is LOGICAL
!          DOPRINT = .TRUE., enables controlled printing of the
!          iteration process
!          DOPRINT = .FALSE., any info printing is suppressed
      use globals, only: factor
!
      implicit none
!     ! in/out
      integer   ::  info, np, ns
      real(8)   ::  stress(6), dstrain(6), dt, params(np), sdv(ns),     &
                    algmod(6,6), lamdot, Dp(5)
      logical   ::  doprint
!     ! internal
      integer   ::  i, j, k, iter, imax, jmax, Nss
      real(8)   ::  DpI(5), res(5), strial(5), f, grad(5),       &
                  hess(5,5), alpha, psi0, psi, Linv(5,5),       &
                  C11, C12, C44, a, L(5,5), beta, tol,                  &
                  dlamdot, eta, dsdev(5), lamdotI, sdevI(5),            &
                  y(5), y2(5), aparams(18), sdev(5),                    &
                  H, choldiag(5), denom, Cdiag(6)
      logical   ::  linesearch
      real(8), allocatable  ::  tau0(:), crss(:)

!     Newton-Raphson and line-search parameters
      beta = 1.d-4
      eta  = 0.1d0
      jmax = 10
      imax = 1000
      tol = 1.d-18
      linesearch = .false.
!
! --- PARAMETERS READING -------
!
!     Elastic constants for cubic symmetry
      C11 = params(1)
      C12 = params(2)
      C44 = params(3)
      Cdiag = [C11 + 2.d0*C12, C11 - C12, C11 - C12, 2.d0*C44, 2.d0*C44, 2.d0*C44]
!
!     exponent of the regularized yield function
      a = params(4)
!     number of slip systems
      Nss = int(params(5))
!     initial critical resolved shear stresses
      allocate(tau0(Nss), source=params(6:5+Nss))
      allocate(crss, source=tau0)
!
! --- ELASTIC PREDICTOR -------
!
!     hydrostatic pressure update
      stress(1) = stress(1) + Cdiag(1)*dstrain(1)
!
!     compute deviatoric trial stress
      strial = stress(2:) + Cdiag(2:)*dstrain(2:)
!
!     critical resolved shear stress (TODO)

!     yield function at STRIAL)
      call yldf(strial, crss, a, 0, f, grad, hess)
      !write(*,*) '########## f: ###############'
      !write(*,*) 'f:', f
!     if no need for the plastic corrector
      if (f .LT. 0.d0) then
!         update the stress
          stress(2:) = strial
!         update the state variables after pure elastic step

!         caluclate the algorithmic modulus
          algmod = 0.d0
          do i=1, 6
              algmod(i,i) = Cdiag(i)
            end do
!         solution flag
          INFO = 3
          goto 101
      end if
!
! --- PLASTIC CORRECTOR -------
!
!     if plastic corrector needed
!     projected trial stress onto yield surface
!     as improved initial guess for deviatoric stress
      sdev = strial/(f+1.d0)
!     initialize plastic strain rate
      do i=1, 5
          Dp(i) = -(sdev(i) - strial(i))/Cdiag(i+1)/dt
      end do
!
!     initial guess for plastic multiplier
      lamdot = vdot(Dp, sdev)
      !write(*,*) '############## Debug: ###################'
      !write(*,*) 'Debug:', lamdot
!     compute residuals
      call yldf(sdev, crss, a, 1, f, grad, hess)
      res = -Dp + lamdot*grad
      !write(*,*) '############## Debug: ###################'
      !write(*,*) 'Debug:', res
!     calculate the merit function
      psi0 = 0.5d0*(vdot(res,res)*dt**2 + f**2)
!
      if (doprint) write(*, *) 'Psi initial ', psi0

!
!     start of iterations
      do iter=1, imax
!
          if (doprint) write(*, *) 'Iter no. ', iter
!
!         start with the Newton step, i.e. alpha = 1
          alpha = 1.d0
!         calculate yield func, its gradient and hessian
          call yldf(sdev, crss, a, 2, f, grad, hess)
!         due to major symmetry of both elastic modulus and
!         the Hessian, only the upper triangle is calculated
          Linv = lamdot*hess
          do i=1, 5
              Linv(i,i) = Linv(i,i) + 1.d0/Cdiag(i+1)/dt
          end do
!         solve Linv*y = grad by Cholesky decomposition, i.e.
!         y = L*grad
!         (lower triangle of Linv will be modified but the upper
!          triangle incl. the diagonal will not)
          call chol_decomp(Linv, 5, choldiag)
          call chol_solve(Linv, 5, grad, choldiag, y)
!         plastic modulus H
          ! H = TODO
          H = 0.d0
!         calculate plastic multiplier increment
          denom = vdot(grad,y) + H
          dlamdot = (f - vdot(res,y))/denom
!
          y2 = res + dlamdot*grad
!         solve Linv*dsdev = y2 by Cholesky decomposition
          call chol_solve(Linv, 5, y2, choldiag, dsdev)
!
          dsdev = -dsdev
          sdevI = sdev + dsdev
!         update plastic strain rate
          do i=1, 5
              DpI(i) = -(sdevI(i) - strial(i))/Cdiag(i+1)/dt
          end do
!         update plastic multiplier
          lamdotI = lamdot + dlamdot
!         update the crss
          !TODO
!         yield function
          call yldf(sdevI, crss, a, 1, f, grad, hess)
!         calculate residuals
          res = -DpI + lamdotI*grad
          psi = 0.5d0*(vdot(res,res)*dt**2 + f**2)
!         Newton-Raphson step end
!
          if (doprint) write(*, *) 'Psi after Newton step ', psi
!
!         start line-search, if not good enough improvement
!         by the Newton-Raphson step
          j = 0
          do while ((psi .GT. (1.d0-2.d0*beta*alpha)*psi0) .and. (j .LE. jmax))
!             flag indicating activated line-search
              linesearch = .TRUE.
              j = j + 1
!             find new alpha as
              alpha = max(eta*alpha, (psi0*alpha**2)/(psi-(1.d0-2.d0*alpha)*psi0))
!
              if (doprint) write(*, *) 'Alpha ', alpha
!
!             update deviatoric stress
              sdevI = sdev + alpha*dsdev
!             update plastic multiplier
              lamdotI = lamdot + alpha*dlamdot
!             update plastic strain rate
              do i=1, 5
                  DpI(i) = -(sdevI(i) - strial(i))/Cdiag(i+1)/dt
              end do
!             update the crss
              !TODO
!             yield function
              call yldf(sdevI, crss, a, 1, f, grad, hess)

!             residuals
              res = -DpI + lamdotI*grad
              psi = 0.5d0*(vdot(res,res)*dt**2 + f**2)
!
              if (doprint) write(*, *) 'Psi in Linesearch ', psi
!
          end do
          ! make candidate solution being the solution of this iteration
          sdev   = sdevI
          Dp     = DpI
          lamdot = lamdotI
          psi0   = psi
          ! crss = TODO hardening stress update
          if (j .GT. jmax) INFO = 2
!         leave the loop if converged
          if (psi0 .LT. tol) exit
      end do
!
!     assign INFO
      if (psi .LT. tol) then
          if (linesearch) then
!             succesful return-map by use of line-search and
!             without reaching JMAX line-search iterations
              if (INFO .NE. 2) INFO = 1
          else
!             succesful return-map by pure Newton-Raphson steps
              INFO = 0
          end if
      else
!         not success
          INFO = -1
      end if
!
! --- STRESS UPDATE (DEVIATORIC PART) -------
!
      do i=1,5
          stress(i+1) = sdev(i)
      end do
!
! --- SOLUTION-DEPENDENT VARIABLES UPDATE -------
!
!     plastic work update
      sdv(1) = sdv(1) + lamdot*dt
!     plastic work rate
      sdv(2) = lamdot
!     number of iterations used
      sdv(3) = real(iter)
!     final value of PSI residual
      sdv(4) = psi
!
! --- ALGORITHMIC MODULUS UPDATE -------
!     (due to major symmetry only upper triangle is calculated)
      call yldf(sdev, crss, a, 2, f, grad, hess)
      Linv = lamdot*hess
      do i=1, 5
          Linv(i,i) = Linv(i,i) + 1.d0/Cdiag(i+1)/dt
      end do
      call chol_decomp(Linv, 5, choldiag)
      call chol_inverse(Linv, 5, choldiag, L)
      do i=1, 5
          y(i) = 0.d0
          do j=1, 5
              y(i) = y(i) + L(i,j)*grad(j)
          end do
      end do
      ! H = TODO
      denom = vdot(grad,y) + H
!
      algmod(1,1) = Cdiag(1)
      do j=2, 6
          algmod(1,j) = 0.d0
      end do
      do i=2, 6
          do j=i, 6
              algmod(i,j) = (L(i-1,j-1) - y(i-1)*y(j-1)/denom)/dt
          end do
      end do
!
101   return
      end subroutine returnmap



end module crystal_plasticity
