module utils

    use globals, only: pi
    implicit none

    contains
! ###################
!
!
   function euler2rotm(ANG1d, ANG2d, ANG3d) result(Q0)
    ! ANG three Euler angles in degrees
        real(kind=8)    ::  ANG(3), Q0(3,3), ANG1, ANG2, ANG3, ANG1d, ANG2d, ANG3d

        ANG1 = deg2rad(ANG1d)
        ANG2 = deg2rad(ANG2d)
        ANG3 = deg2rad(ANG3d)

        Q0(1,1) = COS(ANG1)*COS(ANG3)-SIN(ANG1)*COS(ANG2)*SIN(ANG3)
        Q0(1,2) = SIN(ANG1)*COS(ANG3)+COS(ANG1)*COS(ANG2)*SIN(ANG3)
        Q0(1,3) = SIN(ANG2)*SIN(ANG3)
        Q0(2,1) =-COS(ANG1)*SIN(ANG3)-SIN(ANG1)*COS(ANG2)*COS(ANG3)
        Q0(2,2) =-SIN(ANG1)*SIN(ANG3)+COS(ANG1)*COS(ANG2)*COS(ANG3)
        Q0(2,3) = SIN(ANG2)*COS(ANG3)
        Q0(3,1) = SIN(ANG1)*SIN(ANG2)
        Q0(3,2) =-COS(ANG1)*SIN(ANG2)
        Q0(3,3) = COS(ANG2)

    end function euler2rotm


!C**********************************************************************
!C                         FUNCTION rotm2euler                         *
!C**********************************************************************
!C Computes the Euler angles associated with the orientation matrix    *
!C in terms of the three Euler angles: phi1 (rotation about Z1), PHI   *
!C (rotation about X2) and phi2 (rotation about Z3) - Bunge notation   *
!C**********************************************************************
      function rotm2euler(Q) result(ANG)

      real(kind=8)      ::  ANG(3), Q(3,3), ANG1, ANG2, ANG3, STH

      if (abs(Q(3,3)) < 1.D0) then
        ANG2 = acos(Q(3,3))
        STH  = sin(ANG2)
        ANG1 = atan2(Q(3,1)/STH,-Q(3,2)/STH)
        ANG3 = atan2(Q(1,3)/STH,Q(2,3)/STH)
      else
        ANG1 = atan2(Q(1,2),Q(1,1))
        ANG2 = 0.D0
        ANG3 = 0.D0
      end if
      ANG(1) = rad2deg(ANG1)
      ANG(2) = rad2deg(ANG2)
      ANG(3) = rad2deg(ANG3)

      end function rotm2euler



    function vonMises(A, m) result(VM)
    ! von Mises norm of stress or strain
        real(kind=8)   ::  A(3,3), VM, fVM, trace
        character(6)   ::  m

        trace = A(1,1) + A(2,2) + A(3,3)
        A = A - trace*eye(3)/3.d0

        if (m == 'strain') then
            fVM = sqrt(2.d0/3.d0)
        elseif (m == 'stress') then
            fVM = sqrt(3.d0/2.d0)
        end if

        VM = fVM*sqrt(A(1,1)**2 + A(2,2)**2 + A(3,3)**2 + 2.d0*(A(1,2)**2 + A(1,3)**2 + A(2,3)**2))

    end function vonMises


    subroutine vonMises_dvec(A, m, MA, VM, An)
    ! von Mises norm of deviatoric stress or dev strain given in Voigt notation
    ! masked array MA
        real(kind=8)    ::  A(5), An(5), VM, fVM
        integer         ::  i, MA(5)
        character(6)    ::  m

        intent(in)      ::  A, m, MA
        intent(out)     ::  VM, An

        if (m == 'strain') then
            fVM = sqrt(2.d0/3.d0)
        elseif (m == 'stress') then
            fVM = sqrt(3.d0/2.d0)
        end if

        An = 0.d0
        forall (i=1:5, MA(i)==1) An(i) = A(i)

        VM = fVM*sqrt(An(1)**2 + An(2)**2 + (-An(1)-An(2))**2 + 2.d0*An(3)**2 + 2.d0*An(4)**2 + 2.d0*An(5)**2)

        An = A
        forall (i=1:5, MA(i)==1) An(i) = An(i)/VM

    end subroutine vonMises_dvec


    function m2voigt_dev(M) result(v)

        real(kind=8)    ::  M(3,3), v(5)
        v = (/ M(1,1), M(2,2), M(2,3), M(1,3), M(1,2) /)

    end function m2voigt_dev


    function voigt_dev2m(x) result(M)

        real(kind=8)    ::  M(3,3), x(5)
        M = reshape( (/ x(1), x(5), x(4), x(5), x(2), x(3), x(4), x(3), -x(1)-x(2) /), (/3,3/) )

    end function voigt_dev2m


    function symouter5(x, y) result(M)

	integer     ::  i, j
	real(8)     ::  x(5), y(5), M(5,5)
!
    do i=1, 5
        do j=i, 5
            M(i,j) = x(i)*y(j)
        end do
    end do
!
    end function symouter5


    function outer2vec(e1, e2) result(v)
!
!   Uses new notation
!
	integer     ::  i, j
	real(8)     ::  e1(3), e2(3), v(6), tmp(6)
!
	tmp(1) = e1(1)*e2(1)
	tmp(2) = e1(2)*e2(2)
	tmp(3) = e1(3)*e2(3)
	tmp(4) = e1(2)*e2(3)
	tmp(5) = e1(1)*e2(3)
	tmp(6) = e1(1)*e2(2)
!
	v(1) = 1.d0/sqrt(3.d0)*(tmp(1)+tmp(2)+tmp(3))
	v(2) = 1.d0/sqrt(6.d0)*(2.d0*tmp(3)-tmp(1)-tmp(2))
	v(3) = 1.d0/sqrt(2.d0)*(tmp(2)-tmp(1))
	v(4) = sqrt(2.d0)*tmp(4)
	v(5) = sqrt(2.d0)*tmp(5)
	v(6) = sqrt(2.d0)*tmp(6)
!
    end function outer2vec
!
!
    function matrix2vec(A) result(v)
!
!   Transforms 3x3 matrix into new notation
    real*8      A(3,3), v(6)
!
	v(1) = 1.d0/sqrt(3.d0)*(A(1,1)+A(2,2)+A(3,3))
	v(2) = 1.d0/sqrt(6.d0)*(2.d0*A(3,3)-A(1,1)-A(2,2))
	v(3) = 1.d0/sqrt(2.d0)*(A(2,2)-A(1,1))
	v(4) = sqrt(2.d0)*A(2,3)
	v(5) = sqrt(2.d0)*A(1,3)
	v(6) = sqrt(2.d0)*A(1,2)
!
    end function matrix2vec
!
!
    function vec2matrix(v) result(A)
!
!   Transforms 6x1 vector v written in new notation
!   into a symmetric matrix A
!
    real*8      v(6), A(3,3)
!
    A(1,1) = 1.d0/sqrt(3.d0)*v(1) - 1.d0/sqrt(6.d0)*v(2) - 1.d0/sqrt(2.d0)*v(3)
    A(2,2) = 1.d0/sqrt(3.d0)*v(1) - 1.d0/sqrt(6.d0)*v(2) + 1.d0/sqrt(2.d0)*v(3)
    A(3,3) = 1.d0/sqrt(3.d0)*v(1) + sqrt(2.d0/3.d0)*v(2)
    A(2,3) = 1.d0/sqrt(2.d0)*v(4)
    A(1,3) = 1.d0/sqrt(2.d0)*v(5)
    A(1,2) = 1.d0/sqrt(2.d0)*v(6)
    A(3,2) = A(2,3)
    A(3,1) = A(1,3)
    A(2,1) = A(1,2)
!
    end function vec2matrix
!
!
    subroutine chol_decomp(A, n, choldiag)
!**********************************************************************
!   The CHOL_DECOMP function calculates the Cholesky decomposition
!   A = L*LT of matrix A, where L is lower-triangular matrix and is
!   stored in the lower triangle of A except for the diagonal, which
!   is stored in array CHOLDIAG
!
!**********************************************************************
! (IN/OUT) A is REAL*8 array, dimension (n,n)
!          Matrix A must be a positive-definite symmetric matrix and
!          at input it is read from the upper triangle of A. At return,
!          the factorized matrix L will be stored in the lower triangle
!          of A.
! (IN)     N is INTEGER,
!          It is the dimension of A.
! (OUT)    CHOLDIAG is REAL*8 array, dimension (n)
!          It contains diagonal Cholesky factors of matrix A.
!
      implicit none
      integer     i, j, k, n
      real*8      A(n,n), choldiag(n), sum
!
!     Cholesky factorization of A
      do i=1, n
         do j=i, n
            sum = A(i,j)
            do k=i-1, 1, -1
               sum = sum - A(i,k)*A(j,k)
            end do
            if (i .EQ. j) then
               choldiag(i) = sqrt(sum)
            else
               A(j,i) = sum/choldiag(i)
            end if
         end do
      end do
!
      return
      end subroutine chol_decomp
!
!
      subroutine chol_solve(A, n, b, choldiag, x)
!**********************************************************************
!      The CHOL_SOLVE function returns an n-element vector X containing
!      the solution to the set of linear equations Ax = b. The Cholesky
!      factorization of A must be stored in the lower-diagonal of A and
!      the diagonal factors in CHOLDIAG.
!
!**********************************************************************
! (IN)     A is REAL*8 array, dimension (n,n)
!          Array A contains the original matrix A stored in the upper
!          triangle including the diagonal. The Cholesky factorization
!          of A must at input be stored in the lower triangle, and
!          the diagonal factors stored in CHOLDIAG.
! (IN)     N is INTEGER,
!          It is the dimension of the linear system Ax = b
! (IN)     B is REAL*8 array, dimension (n)
!          It is the right-hand-side vector of the linear system Ax = b
! (IN)     CHOLDIAG is REAL*8 array, dimension (n)
!          It contains diagonal Cholesky factors of matrix A.
! (OUT)    X is REAL*8, dimension (n)
!          It is the solution of the linear system Ax = b
!
      implicit none
      integer     i, j, k, n
      real*8      A(n,n), b(n), choldiag(n), x(n), sum
!
!     backsubstitution
      do i=1, n
         sum = b(i)
         do k=i-1, 1, -1
            sum = sum - A(i,k)*x(k)
         end do
         x(i) = sum/choldiag(i)
      end do
!
      do i=n, 1, -1
         sum = x(i)
         do k=i+1, n
            sum = sum - A(k,i)*x(k)
         end do
         x(i) = sum/choldiag(i)
      end do

      return
      end subroutine chol_solve
!
!
      subroutine chol_inverse(A, n, choldiag, C)
!**********************************************************************
!      The CHOL_INVERSE function returns the matrix C which is inverse of
!      a positive-definite symmetric matrix A.
!
!**********************************************************************
! (IN)     A is REAL*8 array, dimension (n,n)
!          Array A contains the original matrix A stored in the upper
!          triangle including the diagonal. The Cholesky factorization
!          of A must at input be stored in the lower triangle, and
!          the diagonal factors stored in CHOLDIAG.
! (IN)     N is INTEGER,
!          It is the dimension of the linear system Ax = b
! (IN)     B is REAL*8 array, dimension (n)
!          It is the right-hand-side vector of the linear system Ax = b
! (IN)     CHOLDIAG is REAL*8 array, dimension (n)
!          It contains diagonal Cholesky factors of matrix A.
! (OUT)    C is REAL*8, dimension (n,n)
!          It contains inverse of A
!
      implicit none
      integer     i, j, k, n, m
      real*8      A(n,n), b(n), choldiag(n), sum, C(n,n)
!
!     inverse of A by factorized A stored in lower
!     triangle of A and in CHOLDIAG
      do m=1, n
!         create the cartesian basis vectors b
          do i=1, n
              if (m .EQ. i) then
                  b(i) = 1.d0
              else
                  b(i) = 0.d0
              end if
          end do
!         fill the columns of C with the solutions x=A-1*b
          do i=1, n
             sum = b(i)
             do k=i-1, 1, -1
                sum = sum - A(i,k)*C(k,m)
             end do
             C(i,m) = sum/choldiag(i)
          end do
!
          do i=n, 1, -1
             sum = C(i,m)
             do k=i+1, n
                sum = sum - A(k,i)*C(k,m)
             end do
             C(i,m) = sum/choldiag(i)
          end do
      end do
!
      return
      end subroutine chol_inverse



    function deg2rad(angd) result(angr)

    real(kind=8)    ::  angd, angr

    angr = pi*angd/180.d0

    end function deg2rad



    function rad2deg(angr) result(angd)

    real(kind=8)    ::  angd, angr

    angd = 180.d0*angr/pi

    end function rad2deg




    function mdot(A,x) result(y)

    real(kind=8) 	            ::  A(:,:), x(:)
    real(kind=8), allocatable   ::  y(:)
    integer			            ::  i, j, tmp1(2), tmp2, statu

    intent(in)      ::  A, x

    tmp1 = shape(A)
    tmp2 = size(x)
    if (tmp1(2) /= tmp2) then
        write(*,*) "Non-compatible dimensions for dot product."
        read(*,*)
        stop
    else
        allocate( y(tmp1(1)))
        if (allocated(y)) then
            do i = 1, tmp1(1)
                y(i) = 0.d0
                do j = 1, size(x)
                    y(i) = y(i) + A(i,j)*x(j)
                end do
            end do
        end if
    end if

    end function mdot
! ###################


!C**********************************************************************
!C                         FUNCTION TRANSFORM                          *
!C**********************************************************************
!C Rotates a matrix from one frame to another                          *
!C B = P.A.P^T                                                         *
!C**********************************************************************
    function transform(A, P) result(B)

    real(kind=8) A(3,3), B(3,3), P(3,3), PT(3,3), MAT(3,3)

    intent(in)  ::  A, P

    PT = transp3(P)
    MAT = mmult(A, PT)
    B = mmult(P, MAT)

    end function transform


!C**********************************************************************
!C                         FUNCTION MAT2VEC                            *
!C**********************************************************************
!C Transforms a (3x3) symmetric matrix into a (6x1) vector			   *
!C Van Houtte's specification of 3*3 matrix into vector, 6 dimensional *
!C space.                                                              *
!C**********************************************************************
    function mat2vec(MAT) result(VEC)

    real(kind=8) MAT(3,3), VEC(6)

    intent(in)  :: MAT

    VEC(1) = MAT(1,1)
    VEC(2) = MAT(2,2)
    VEC(3) = MAT(3,3)
    VEC(4) = MAT(2,3)
    VEC(5) = MAT(1,3)
    VEC(6) = MAT(1,2)

    end function mat2vec


!C**********************************************************************
!C                         FUNCTION VEC2MAT                            *
!C**********************************************************************
!C Transforms a (6x1) vector (3x3) into a symmetric matrix             *
!C**********************************************************************
    function vec2mat(VEC) result(MAT)

    real(kind=8) MAT(3,3), VEC(6)
    intent(in)  :: VEC

    MAT(1,1) = VEC(1)
    MAT(2,2) = VEC(2)
    MAT(3,3) = VEC(3)
    MAT(2,3) = VEC(4)
    MAT(3,2) = VEC(4)
    MAT(1,3) = VEC(5)
    MAT(3,1) = VEC(5)
    MAT(1,2) = VEC(6)
    MAT(2,1) = VEC(6)

    end function vec2mat

  !C**********************************************************************
  !C                         FUNCTION VEC2SKEWMAT                            *
  !C**********************************************************************
  !C Transforms a (6x1) vector (3x3) into a symmetric matrix             *
  !C**********************************************************************
      function vec2skewmat(VEC) result(MAT)

      real(kind=8) MAT(3,3), VEC(6)
      intent(in)  :: VEC

      MAT(1,1) = VEC(1)
      MAT(2,2) = VEC(2)
      MAT(3,3) = VEC(3)
      MAT(2,3) = VEC(4)
      MAT(3,2) = -VEC(4)
      MAT(1,3) = VEC(5)
      MAT(3,1) = -VEC(5)
      MAT(1,2) = VEC(6)
      MAT(2,1) = -VEC(6)

    end function vec2skewmat


!C**********************************************************************
!C                         FUNCTION MINV                               *
!C**********************************************************************
!C Computes the inverse of a (3x3) matrix AINV = A^-1                  *
!C**********************************************************************
    function minv3(A) result(AINV)

    real(kind=8) ::  A(3,3), AINV(3,3), ADJ(3,3), DET
    integer      ::  I,J

    intent(in)   ::  A

    !C Compute the determinant of A
    DET = determ(A)
    !C Compute the adjoint matrix of A
    ADJ(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
    ADJ(1,2)=-(A(2,1)*A(3,3)-A(3,1)*A(2,3))
    ADJ(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
    ADJ(2,1)=-(A(1,2)*A(3,3)-A(3,2)*A(1,3))
    ADJ(2,2)=A(1,1)*A(3,3)-A(3,1)*A(1,3)
    ADJ(2,3)=-(A(1,1)*A(3,2)-A(3,1)*A(1,2))
    ADJ(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
    ADJ(3,2)=-(A(1,1)*A(2,3)-A(2,1)*A(1,3))
    ADJ(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
    !C Compute the transpose of the adjoint matrix of A = inverse of A
    DO I=1,3
        DO J=1,3
            AINV(I,J)=ADJ(J,I)/DET
        END DO
    END DO

    end function minv3


!C**********************************************************************
!C                         FUNCTION MMULT                              *
!C**********************************************************************
!C Computes the product of two (3x3) matrices AB = A.B                 *
!C**********************************************************************
    function mmult(A,B) result(AB)

    real(kind=8)    ::  A(3,3), B(3,3), AB(3,3), P
    integer         ::  I,J,K

    intent(in)   ::  A, B

    DO I=1,3
        DO J=1,3
            P=0.D0
            DO K=1,3
                P=P+A(I,K)*B(K,J)
            END DO
            AB(I,J)=P
        END DO
    END DO

    end function mmult


!C**********************************************************************
!C                         FUNCTION TRANSP3                            *
!C**********************************************************************
!C Computes the transpose of a (3x3) matrix AT = A^T                   *
!C**********************************************************************
    function transp3(A) result(AT)

    real(kind=8)    ::  A(3,3), AT(3,3)
    integer         ::  I,J

    intent(in)      ::  A

    DO I=1,3
        DO J=1,3
            AT(I,J)=A(J,I)
        END DO
    END DO

    end function transp3

!C**********************************************************************
!C                         FUNCTION TRANSP                            *
!C**********************************************************************
!C Computes the transpose of a (mxn) matrix AT = A^T                   *
!C**********************************************************************
    function transp(A) result(AT)

    real(kind=8)    ::  A(:,:)
    real(kind=8), allocatable ::  AT(:,:)
    integer         ::  I,J,shpA(2)

    intent(in)      ::  A

    shpA = shape(A)
    allocate(AT(shpA(2),shpA(1)))
    if (allocated(AT)) then
        DO I=1,shpA(2)
            DO J=1,shpA(1)
                AT(I,J)=A(J,I)
            END DO
        END DO
    end if

    end function transp


!C**********************************************************************
!C                         FUNCTION GETSYM                             *
!C**********************************************************************
!C Computes the symmetric part of a (3x3) matrix A                     *
!C**********************************************************************
    function getsym(A) result(ASYM)

    real(kind=8)    ::  A(3,3), ASYM(3,3)

    intent(in)      ::  A

    ASYM = 0.5d0*(A+TRANSP3(A))

    end function getsym


!C**********************************************************************
!C                         FUNCTION GETSKW                             *
!C**********************************************************************
!C Computes the symmetric part of a (3x3) matrix A                     *
!C**********************************************************************
    function getskw(A) result(ASKW)

    real(kind=8)    ::  A(3,3), ASKW(3,3)

    intent(in)      ::  A

    ASKW = 0.5d0*(A-TRANSP3(A))

    end function getskw


!C**********************************************************************
!C                         FUNCTION FLATTEN                            *
!C**********************************************************************
!C To flatten the array A (N,M) into the vector y (M*N) row-vise       *
!C**********************************************************************
    function flatten(A) result(y)

    implicit none

    REAL(kind=8)    ::  A(:,:)
    REAL(kind=8), allocatable :: y(:)
    INTEGER ::  shpA(2), i, j, k

    intent(in) :: A

    shpA = shape(A)
    allocate(y(shpA(2)*shpA(1)))
    k = 0
    if (allocated(y)) then
        do i=1,shpA(1)
            do j=1,shpA(2)
                k = k + 1
                y(k) = A(i,j)
            end do
        end do
    end if

    end function flatten

!C**********************************************************************
!C                         FUNCTION reshape_rows                          *
!C**********************************************************************
!C Reshape the vector y (M*N) into the array A (N,M) row-vise          *
!C**********************************************************************
    function reshape_rows(y, shp) result(A)

    real(kind=8)    ::  y(:)
    integer(kind=8), allocatable :: A(:,:)
    integer ::  shp(2), i, j, k

    intent(in) :: y, shp

    allocate(A(shp(1),shp(2)))
    k = 0
    if (allocated(A)) then
        do i=1,shp(1)*shp(2),shp(2)
            k = k + 1
            A(k,:) = y(i:i+shp(2)-1)
        end do
    end if

    end function reshape_rows


!C**********************************************************************
!C                         FUNCTION DETERM                             *
!C**********************************************************************
!C Computes the determinant of a (3x3) matrix                          *
!C**********************************************************************
    function determ(A) result(DET)

    real(kind=8)    ::  A(3,3), DET

    intent(in)      ::  A

    DET = A(1,1)*A(2,2)*A(3,3) + &
            A(2,1)*A(3,2)*A(1,3) + &
            A(3,1)*A(1,2)*A(2,3) - &
            A(3,1)*A(2,2)*A(1,3) - &
            A(1,1)*A(3,2)*A(2,3) - &
            A(2,1)*A(1,2)*A(3,3)

    IF (ABS(DET) < 1.D-12) THEN
        WRITE(*,*) 'Determinant is null!'
        READ(*,*)
        STOP
    END IF

    end function determ


!C**********************************************************************
!C                         FUNCTION POLAR                              *
!C**********************************************************************
!C Polar decomposition of F = R.U using the Cayley-Hamilton theorem    *
!C see Nemat-Nasser's book p.55                                        *
!C Returns R                                                           *
!C**********************************************************************
    function polar(F) result(R)

    real(kind=8)    :: F(3,3),R(3,3),C(3,3),CS(3,3),U(3,3),UI(3,3),C1,C3,P, &
                        CD11,CD22,CD33,CD2,CD3,U1,U2,U3,A,B,PHI,L1,D,E,A3,B2
    integer         :: I,J,K

    intent(in)      :: F

    !C Compute stretch tensor: C = F^T.F
    DO J=1,3
        DO I=1,3
            C(I,J)=0.D0
            DO K=1,3
                C(I,J)=C(I,J)+F(K,I)*F(K,J)
            END DO
        END DO
    END DO
    !C Compute C^2
    CS = MMULT(C, C)
    !C Compute invariants
    C1=C(1,1)+C(2,2)+C(3,3)

    C3=C(1,1)*(C(2,2)*C(3,3)-C(2,3)*C(3,2))+ &
        C(1,2)*(C(2,3)*C(3,1)-C(2,1)*C(3,3))+ &
        C(1,3)*(C(2,1)*C(3,2)-C(2,2)*C(3,1))

    !C Invariants of the deviatoric part CD of tensor C
    P=(C(1,1)+C(2,2)+C(3,3))/3.D0

    CD11=C(1,1)-P
    CD22=C(2,2)-P
    CD33=C(3,3)-P

    CD2=CD11*CD22+CD11*CD33+CD22*CD33- &
        (C(1,2)*C(2,1)+C(1,3)*C(3,1)+C(2,3)*C(3,2))
    CD3=CD11*(CD22*CD33-C(2,3)*C(3,2))+ &
        C(1,2)*(C(2,3)*C(3,1)-C(2,1)*CD33)+ &
        C(1,3)*(C(2,1)*C(3,2)-CD22*C(3,1))
    !C Invariants of U
    U3=sqrt(C3)
    A=-CD2/3.D0
    B=CD3/2.D0
    A3=A**3.D0
    B2=B*B

    IF (ABS(A3-B2).GT.1.D-12) THEN
        PHI=ACOS(B/A**(3.D0/2.D0))
        L1=SQRT(C1/3.D0+2.D0*SQRT(A)*COS(PHI/3.D0))
    ELSE
        L1=SQRT(C1/3.D0)
    END IF

    U1=L1+SQRT(C1-L1*L1+2.D0*U3/L1)
    U2=0.5D0*(U1*U1-C1)
    !C Computes U
    D=U3-U1*U2
    E=U1**2.D0-U2
    DO I=1,3
        DO J=1,3
            U(I,J)=(CS(I,J)-E*C(I,J))/D
            IF (I.EQ.J) THEN
                U(I,J)=U(I,J)-U1*U3/D
            END IF
        END DO
    END DO
    UI = minv3(U)
    R = mmult(F,UI)

    end function polar


!C**********************************************************************
!C                         function MNORM                              *
!C**********************************************************************
!C Computes the norm of a (RxC) matrix using double contraction        *
!C Norm = (Tr(A.A^T))^1/2 = (A:A)^1/2                                  *
!C**********************************************************************
    function mnorm(A) result(norma)

    integer         :: sha(2)
    real(kind=8)    :: A(:,:), norma
    integer         :: I, J

    intent(in)      :: A

    shA = shape(A)
    norma = 0.D0
    DO I=1,shA(1)
        DO J=1,shA(2)
            norma = norma + A(I,J)**2
        END DO
    END DO

    norma = SQRT(norma)

    end function mnorm


! ###################
    pure function vdot(x,y) result(z)

    real(kind=8)	::  x(:), y(:), z
    integer			::  i

    intent(in)      ::  x, y

	z = 0.d0
    do i = 1, size(x)
        z = z + x(i)*y(i)
    end do

    end function vdot
! ###################

! ###################
	pure function cross(x,y) result(z)

	real(kind=8)	::  x(3), y(3), z(3)

    intent(in)      ::  x, y

	z(1) = x(2)*y(3)-y(2)*x(3)
	z(2) = x(3)*y(1)-y(3)*x(1)
	z(3) = x(1)*y(2)-y(1)*x(2)

    end function cross
! ###################

! ###################
    pure function norm(x) result(y)

    real(kind=8)	::  x(:), y
    integer			::  i

    intent(in)      ::  x

	y = 0.d0
    do i = 1, size(x)
        y = y + x(i)**2
    end do
	y = sqrt(y)

    end function norm
! ###################

! ###################
    function eye(n) result(A)

    integer			::  n, i
    real(kind=8)	::  A(n,n)

    intent(in)      ::  n

	A = 0.d0
    do i = 1, n
        A(i,i) = 1.d0
    end do

    end function eye
! ###################

! ###################
    function outer3(x,y) result(A)

    integer			::  i, j
    real(kind=8)	::  x(3), y(3), A(3,3)

    intent(in)      ::  x, y

	A = 0.d0
    do i = 1,3
        do j = 1,3
            A(i,j) = x(i)*y(j)
        end do
    end do

    end function outer3
! ###################

! ###################
    function inner(A,B) result(c)

    integer			::  i, j
    real(kind=8)	::  A(3,3), B(3,3), c

    intent(in)      ::  A, B

	c = 0.d0
    do i = 1,3
        do j = 1,3
            c = c + A(i,j)*B(i,j)
        end do
    end do

    end function inner
! ###################

! ###################
    function mean(x) result(a)

    real(kind=8)	::  x(:), a

    intent(in)      ::  x

    a = sum(x)/size(x)

    end function mean
! ###################


! ###################
    function schmidt(e2, e1) result(e2a)
    ! Schmidt ortho-normalize e2 against e1

    real(kind=8)	::  e2(3), e1(3), e2a(3), e3(3), l3

    intent(in)      ::  e2, e1

    e3 = cross(e1, e2)
    l3 = norm(e3)
    if (l3 < 1.d-10) then
        e2a = e2
        e2a(1) = e2a(1) + 1.d0
        e3 = cross(e1,e2a)
        l3 = norm(e3)
        if (l3==0) then
            e2a = e2
            e2a(2) = e2a(2) + 1.d0
            e3 = cross(e1, e2a)
            l3 = norm(e3)
        end if
    end if
    e3 = e3/l3

    e2a = cross(e3, e1)
    e2a = e2a/norm(e2a)
    end function schmidt
! ###################


! ###################
    subroutine minverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
    integer        ::  n
    real(kind=8)   ::  a(n,n), c(n,n)
    real(kind=8)   ::  L(n,n), U(n,n), b(n), d(n), x(n)
    real(kind=8)   ::  coeff
    integer        ::  i, j, k

    intent(in)     ::  n
    intent(out)    ::  c

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 allows such operations on matrices
    L = 0.d0
    U = 0.d0
    b = 0.d0

    ! step 1: forward elimination
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
        L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
        do i=1,j
        U(i,j) = a(i,j)
        end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
        b(k)=1.0
        d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
        d(i)=b(i)
        do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
        end do
        end do
    ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
        end do
    ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
        c(i,k) = x(i)
        end do
        b(k)=0.0
    end do

    end subroutine minverse


	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

!  Purpose: Computes all eigenvalues and eigenvectors of a real
!     symmetric matrix A, which is of size N by N, stored in a
!     physical NP by NP array.  On output, elements of A above the
!     diagonal are destroyed.  D returns the eigenvalues of A in
!     its first N elements.  V is a matrix with the same logical and
!     physical dimensions as A whose columns contain, on output, the
!     normalized eigenvectors of A.  NROT returns the number of Jacobi
!     rotations which were required.
!
!  Source: W. H. Press et al., "Numerical Recipes", 1989, p. 346.
!
!  Modifications:
!
!     1. Double precision version
!
!  Prepared by J. Applequist, 10/23/91

    IMPLICIT REAL*8(A-H,O-Z)
    INTEGER, PARAMETER :: NMAX=100
    INTEGER :: I, IP, IQ, J, N, NP, NROT
    DIMENSION :: A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

!     Initialize the identity matrix.

    DO 12 IP=1,N
    DO 11 IQ=1,N
    V(IP,IQ)=0.D0
11   CONTINUE
    V(IP,IP)=1.D0
12   CONTINUE

!     Initialize B and D to the diagonal of A.

    DO 13 IP=1,N
    B(IP)=A(IP,IP)
    D(IP)=B(IP)
    Z(IP)=0.D0
13   CONTINUE
    NROT=0
    DO 24 I=1,50
    SM=0.D0

!     Sum off-diagonal elements.

    DO 15 IP=1,N-1
    DO 14 IQ=IP+1,N
    SM=SM+DABS(A(IP,IQ))
14   CONTINUE
15   CONTINUE
    IF (SM.EQ.0.D0) RETURN
    IF (I.LT.4) THEN
    TRESH=0.2D0*SM/N**2
    ELSE
    TRESH=0.D0
    ENDIF
    DO 22 IP=1,N-1
    DO 21 IQ=IP+1,N
    G=100.D0*DABS(A(IP,IQ))

!     After four sweeps, skip the rotation if the off-diagonal
!     element is small.

    IF ((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP))) &
    .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ)))) THEN
    A(IP,IQ)=0.D0
    ELSE IF (DABS(A(IP,IQ)).GT.TRESH) THEN
    H=D(IQ)-D(IP)
    IF (DABS(H)+G.EQ.DABS(H)) THEN
    T=A(IP,IQ)/H
    ELSE
    THETA=0.5D0*H/A(IP,IQ)
    T=1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
    IF (THETA.LT.0.D0) T=-T
    ENDIF
    C=1.D0/DSQRT(1.D0+T**2)
    S=T*C
    TAU=S/(1.D0+C)
    H=T*A(IP,IQ)
    Z(IP)=Z(IP)-H
    Z(IQ)=Z(IQ)+H
    D(IP)=D(IP)-H
    D(IQ)=D(IQ)+H
    A(IP,IQ)=0.D0
    DO 16 J=1,IP-1
    G=A(J,IP)
    H=A(J,IQ)
    A(J,IP)=G-S*(H+G*TAU)
    A(J,IQ)=H+S*(G-H*TAU)
16   CONTINUE
    DO 17 J=IP+1,IQ-1
    G=A(IP,J)
    H=A(J,IQ)
    A(IP,J)=G-S*(H+G*TAU)
    A(J,IQ)=H+S*(G-H*TAU)
17   CONTINUE
    DO 18 J=IQ+1,N
    G=A(IP,J)
    H=A(IQ,J)
    A(IP,J)=G-S*(H+G*TAU)
    A(IQ,J)=H+S*(G-H*TAU)
18   CONTINUE
    DO 19 J=1,N
    G=V(J,IP)
    H=V(J,IQ)
    V(J,IP)=G-S*(H+G*TAU)
    V(J,IQ)=H+S*(G-H*TAU)
19   CONTINUE
    NROT=NROT+1
    ENDIF
21   CONTINUE
22   CONTINUE
    DO 23 IP=1,N
    B(IP)=B(IP)+Z(IP)
    D(IP)=B(IP)
    Z(IP)=0.D0
23   CONTINUE
24   CONTINUE
    WRITE (6,600)
600  FORMAT(/'50 ITERATIONS OCCURRED IN SUBROUTINE JACOBI.')
    RETURN
    END SUBROUTINE JACOBI




! ###################


    subroutine linspace(min, max, array)
        real(kind=8)	::  step, dx
        integer         ::  n_steps, i

        real(kind=8), intent(in)      ::  min, max
        real(kind=8), intent(out)     ::  array(:)


        n_steps = size(array)

        if (n_steps == 0) return

        if (n_steps == 1) then
            array(1) = min
            return
        end if

        dx = (max - min) / real(n_steps - 1)

        do i=1,n_steps
            array(i) = min + dx*(i-1)
        end do

    end subroutine linspace




end module utils
