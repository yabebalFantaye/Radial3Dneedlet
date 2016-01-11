!==============================================================================
! FILE:  fftw_m.f90 - module contains fftw3 1D and 2D wrappers:
!    Functions:   fftw, dfftw           - 1D subroutines (single, double)
!    Subroutines: fftw2D, dfftw2D       - 2D subroutines (single, double)
!                 dfftw_r2c, dfftw_c2r  - real->complex, complex->real
!==============================================================================

module fftw_m

contains

!*******************************************************************************
! function fftw(in[, direction]) - Returns the output of fftw's FFT of 'in'.
!      -Input and output are real complex arrays, which are converted to
!         (and then back from) double precision for the actual FFT.
!      -Direction is an optional argument:  if present, and specified as 1,
!         direction is backward; otherwise, direction is forward.
! ELW-20090820
!*******************************************************************************
function fftw(in, direction)
  complex, dimension(0:),intent(in)      :: in
  integer(4), intent(in), optional       :: direction
  complex, dimension(0:size(in)-1)       :: fftw
  double complex,dimension(0:size(in)-1) :: dblin, dblout
  integer(4) :: dir
  integer(8) :: N, i, plan

#include "fftw3.f"

  N = size(in)
  dblin = dcmplx(in)            ! Convert cmplx real 'in' array to cmplx double

  if (present(direction)) then
        dir = direction         ! Take direction from parameter
  else
      dir = 0                   ! Default forward direction
  endif

  if (dir .ne. 1) then
      call dfftw_plan_dft_1d(plan,N,dblin, dblout, FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute( plan )
      dblout = dblout/dble(N)   ! Normalize
  else
      call dfftw_plan_dft_1d(plan,N,dblin, dblout, FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute( plan )
  endif
  call dfftw_destroy_plan( plan )
  fftw = cmplx(dblout)          ! Convert cmplx double array to cmplx real
end function fftw

!*******************************************************************************
! fftw2D() - subroutine to compute 2D FFT in the given direction, as with fftw,
!       except that output is via the "out" variable, for which memory is
!       assumed to have been allocated by the caller.
!       -in and out are single precision complex 2D arrays
!       -in/out are converted to double for use w/ dfft -- should change the
!        subroutine to use a single precision fftw call instead, if possible.
!*******************************************************************************
subroutine fftw2D(in, out, direction)
  complex,dimension(0:,0:),intent(in)      :: in
  complex,dimension(0:size(in(:,0))-1,0:size(in(0,:))-1),intent(inout) :: out
  double complex,dimension(0:size(in(:,0))-1,0:size(in(0,:))-1) :: dblin,dblout
  integer(4), intent(in), optional :: direction
  integer(4) :: dir
  integer(8) :: Nx, Ny, i, plan

#include "fftw3.f"
  Nx = size(in(:,0)) ; Ny = size(in(0,:))

  dblin = dcmplx(in)            ! Convert cmplx real 'in' array to cmplx double

  if (present(direction)) then
        dir = direction         ! Take direction from parameter
  else
      dir = 0                   ! Default forward direction
  endif

  if (dir .ne. 1) then
      call dfftw_plan_dft_2d( plan, Nx,Ny, dblin, dblout, FFTW_FORWARD, FFTW_ESTIMATE )
      call dfftw_execute( plan )
      dblout = dblout/dble(Nx*Ny)   ! Normalize
  else
      call dfftw_plan_dft_2d( plan, Nx,Ny, dblin, dblout, FFTW_BACKWARD, FFTW_ESTIMATE )
      call dfftw_execute( plan )
  endif
  call dfftw_destroy_plan( plan )
  out = cmplx(dblout)          ! Convert cmplx double array to cmplx real
end subroutine fftw2D


!*******************************************************************************
! function dfftw(in[, direction]) - Returns the output of fftw's FFT of 'in'.
!      -Input and output are double-precision complex arrays.
!      -Direction is an optional argument:  if present, and specified as 1,
!         direction is backward; otherwise, direction is forward.
! ELW-20090819
!*******************************************************************************
function dfftw(in, direction)
  double complex,dimension(0:),intent(in) :: in
  double complex,dimension(0:size(in)-1)  :: dfftw
  integer(4), intent(in), optional :: direction
  integer(4) :: dir
  integer(8) :: N, i, plan

#include "fftw3.f"

  N = size(in)

  if (present(direction)) then
      dir = direction           ! Take direction from parameter
  else
      dir = 0                   ! Default forward direction
  endif

  if (dir .ne. 1) then
      call dfftw_plan_dft_1d( plan, N, in, dfftw, FFTW_FORWARD, FFTW_ESTIMATE )
      call dfftw_execute( plan )
      dfftw = dfftw/dble(N)             ! Normalize
  else
      call dfftw_plan_dft_1d( plan, N, in, dfftw, FFTW_BACKWARD, FFTW_ESTIMATE )
      call dfftw_execute( plan )
  endif
  call dfftw_destroy_plan( plan )
end function dfftw

!*******************************************************************************
! subroutine dfftw_r2c(in, out) - Returns the FFT of 'in' in 'out'.
!      -Forward direction implied
!      -N and Plan are computed on first call, and never ever again.
!          Therefore, use this function ONLY when you never expect to see a
!          different sized or shaped array in the program.
! ELW-20090921
!*******************************************************************************
subroutine dfftw_r2c(in, out)
  real(8),dimension(0:),intent(in)            :: in
  complex(8),dimension(0:size(in)-1),intent(inout) :: out
  real(8), save :: Ndouble
  integer(8), save :: N, plan_r2c
  logical, save :: reset = .true.
#include "fftw3.f"

  if (reset .eq. .true.) then
      N = size(in)
      Ndouble = dble(N)
      call dfftw_plan_dft_r2c_1d( plan_r2c, N, in, out, FFTW_ESTIMATE )
      call dfftw_execute( plan_r2c )
      reset = .false.
  else
      call dfftw_execute_dft_r2c( plan_r2c, in, out )
  endif
  out = out/Ndouble     ! Normalize
! call dfftw_destroy_plan( plan_r2c )
end subroutine dfftw_r2c

!*******************************************************************************
! subroutine dfftw_c2r(in, out) - Returns the FFT of 'in' in 'out'.
!      -Implied reverse direction
!      -Q&D: N and Plan are computed on first call, and never ever again.
!          Therefore, use this function ONLY when you never expect to see a
!          different sized or shaped array in the program.
! ELW-20090921
!*******************************************************************************
subroutine dfftw_c2r(in, out)
  complex(8),dimension(0:),intent(in)           :: in
  real(8),dimension(0:size(in)-1),intent(inout) :: out
  integer(8), save :: N, plan_c2r
  logical, save :: reset = .true.
#include "fftw3.f"

  if (reset .eq. .true.) then
      N = size(in)
      call dfftw_plan_dft_c2r_1d( plan_c2r, N, in, out, FFTW_ESTIMATE )
      call dfftw_execute( plan_c2r )
      reset = .false.
  else
      call dfftw_execute_dft_c2r( plan_c2r, in, out )
  endif

! call dfftw_destroy_plan( plan_c2r )
end subroutine dfftw_c2r


!*******************************************************************************
! dfftw2D() - subroutine to compute 2D FFT in the given direction, as with fftw,
!       except that output is via the "out" variable, for which memory is
!       assumed to have been allocated by the caller.
!       -in and out are double precision complex 2D arrays.
!*******************************************************************************
subroutine dfftw2D(in, out, direction)
  double complex,dimension(0:,0:),intent(in) :: in
  double complex,dimension(0:size(in(:,0))-1,0:size(in(0,:))-1),intent(inout) :: out
  integer(4), intent(in), optional :: direction
  integer(4) :: dir
  integer(8) :: Nx, Ny, i, plan

#include "fftw3.f"

  Nx = size(in(:,0)) ; Ny = size(in(0,:))

  if (present(direction)) then
      dir = direction           ! Take direction from parameter
  else
      dir = 0                   ! Default forward direction
  endif

  if (dir .ne. 1) then
      call dfftw_plan_dft_2d( plan, Nx,Ny, in, out, FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute( plan )
      out = out/dble(Nx*Ny)             ! Normalize
  else
      call dfftw_plan_dft_2d( plan, Nx,Ny, in, out, FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute( plan )
  endif
  call dfftw_destroy_plan( plan )
end subroutine dfftw2D

end module fftw_m
