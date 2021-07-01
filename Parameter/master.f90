program Master
    use globals
    use utils
    use crystal_plasticity
    use TaylorLin

    implicit none

    INTEGER :: i, z, zz, zzz, counter, cs
    real(8) :: deformation(1), scale1(11), scale2(4), exp_n(4), start, finished

    call cpu_time(start)
    !deformation = (/ 0.50d0, 0.70d0, 0.8d0, 0.9d0/)
    deformation = (/ 0.90d0/)
    exp_n = (/20.d0, 50.d0, 100.d0, 1000.d0/)
    !exp_n = (/1000.d0/)
    
    ! Scaling the CRSS values 
    ! Scale 1 = {211}<111> /{110}<111>. The range of effect is Cos(30) to 1/cos(30) for CO-slip {211} & {110}
    ! Scale 2 = {321}<111> /{110}<111>. The range of effect is Cos(19.1) to cos(10.9)/cos(30) for CO-slip {321} & {110}
    ! For CO-slip on {321}{211}{110} the limit is given by scale1 = scale2/cos(10.9) to scale1 = scale2 *cos(30) /cos(10.9)

    ! scale2 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.d0/)
    ! scale1 = (/0.90d0, 0.925d0, 0.95d0, 0.975d0, 1.0d0, 1.025d0, 1.05d0, 1.075d0, 1.10d0/)
    ! counter = 0
    ! cs = 123
    ! scale1 = (/0.85d0, 0.90d0, 0.925d0, 0.95d0, 0.975d0, 1.0d0, 1.025d0, 1.05d0, 1.075d0, 1.10d0, 1.20d0/)
    ! scale2 = (/1.0d0, 0.99d0, 1.00d0, 1.01d0/)
    ! do zzz=6, 6
    !     do i=1, 1
    !         do z=1, 4
    !             counter = counter + 1
    !             call test(deformation(1), scale1(zzz), scale2(i), cs, counter, exp_n(z))
    !         end do
    !     end do
    !     !counter = zzz
    ! end do
    counter = 0
    deformation = (/ 0.99d0/)
    do z=1, 3
        counter = counter + 1
        call test(deformation(1), 0.945d0, 0.985d0, 123, counter, exp_n(z))
    end do

    ! counter = 1000
    ! cs = 123
    ! scale1 = (/0.85d0, 0.90d0, 0.925d0, 0.95d0, 0.975d0, 1.0d0, 1.025d0, 1.05d0, 1.075d0, 1.10d0, 1.20d0/)
    ! scale2 = (/1.0d0, 0.99d0, 1.00d0, 1.01d0/)
    ! do zzz=1, 1
    !     do i=1, 1
    !         do z=1, 4
    !             counter = counter + 1
    !             call test(deformation(1), scale1(zzz), scale2(i), cs, counter, exp_n(z))
    !         end do
    !     end do
    !     !counter = zzz
    ! end do
    ! ! CS 1 simulations 
    ! cs = 1
    ! scale2 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
    ! scale1 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
    ! do zzz=1,1
    !     do i=1,4
    !         do z=1, 1
    !             counter = counter + 1
    !             call test(deformation(zzz), scale1(z), scale2(z), cs, counter, exp_n(i))
    !         end do
    !     end do
    ! end do

    ! CS 2 simulations 
!     counter = 100
!     cs = 12
!     scale2 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
!     scale1 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
!     do zzz=1,1
!         do i=1,4
!             do z=1, 1
!                 counter = counter + 1
!                 call test(deformation(zzz), scale1(z), scale2(z), cs, counter, exp_n(i))
!             end do
!         end do
!     end do
! ! CS 3 simulations 
!     counter = 200
!     cs = 13
!     scale2 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
!     scale1 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
!     do zzz=1,1
!         do i=1,4
!             do z=1, 1
!                 counter = counter + 1
!                 call test(deformation(zzz), scale1(z), scale2(z), cs, counter, exp_n(i))
!             end do
!         end do
!     end do
    
    ! counter = 1000
    ! ! CS 12 simulations 
    ! cs = 12
    ! scale2 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.d0/)
    ! scale1 = (/0.90d0, 0.925d0, 0.95d0, 0.975d0, 1.0d0, 1.025d0, 1.05d0, 1.075d0, 1.10d0/)
    ! do zzz=1,1
    !     do i=1,4
    !         do z=1, 9
    !             counter = counter + 1
    !             call test(deformation(zzz), scale1(z), scale2(z), cs, counter, exp_n(i))
    !         end do
    !     end do
    ! end do

    ! ! CS 13 simulations 
    ! counter = 2000
    ! cs = 13
    ! scale1 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.d0/)
    ! scale2 = (/0.90d0, 0.925d0, 0.95d0, 0.975d0, 1.0d0, 1.025d0, 1.05d0, 1.075d0, 1.10d0/)
    ! do zzz=1,1
    !     do i=1,4
    !         do z=1, 9
    !             counter = counter + 1
    !             call test(deformation(zzz), scale1(z), scale2(z), cs, counter, exp_n(i))
    !         end do
    !     end do
    ! end do


    !CS 123 simulations 
    ! counter = 10000
    ! cs = 123
    ! scale1 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.d0/)
    ! scale2 = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.d0/)
    ! exp_n = (/5.d0, 10.d0, 100.d0/)
    ! do zzz=1, 1
    !     do i=1, 3
    !         do z=1, 1
    !             counter = counter + 1
    !             call test(deformation(zzz), scale1(z), scale2(z), cs, counter, exp_n(i))
    !         end do
    !     end do
    ! end do

    call cpu_time(finished)
    print '("Time = ",f6.3," seconds.")',finished-start
    print '(''Ferdig'')'
end program Master
