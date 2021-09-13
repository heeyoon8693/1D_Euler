!===========================================================================
!                     1D Euler Sod Shock tube problem
!                       common parameter for this problem
!                           2021.09.07 
!                       made by Hee Yoon / Republic of Korea
! ==========================================================================
! delcare of parameter
!---------------------------------------------------------------------------
 module parameter_mod
 implicit real(a-h,o-z)
    integer, parameter :: im = 101, im1 = im-1 ! number of computing cell at x direction  
    real, parameter  :: gm = 1.4d0, gm1 = 0.4d0
    real, parameter  :: dt = 0.001d0 
    real, parameter  :: leng_x = 1.0d0, x_diap = leng_x/2
    real, parameter  :: prel_ini = 1.0d0, prer_ini = 0.1d0 ! initial condition of pressure
    real, parameter  :: vell_ini = 0.0d0, velr_ini = 0.0d0 ! initial condition of velocity
    real, parameter  :: rhol_ini = 1.0d0, rhor_ini = 0.125d0  ! initial condition of density
    real :: dx
    integer, parameter  :: iflux =5, iend =140
    integer, parameter :: iorder=3 
   ! integer :: i
    real, dimension(im) :: x, pre, vel, rho, e
    real, dimension(im, 3) :: q
 end module


