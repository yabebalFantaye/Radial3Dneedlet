module fftw_wrapper

  use healpix_types
  use healpix_fft
  implicit none  

  include "fftw3.f"

  integer*8 plan, plan1, plan2
  
contains

  integer*8  function make_plan(in,out,backward,ncpu)
    COMPLEX(DPC), dimension(:) :: in, out
    integer*8 :: p
    integer :: N, direction
    integer :: iret, id, nthreads
    logical :: backward,use_fftw_omp
    integer, optional :: ncpu
    INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS


    N = size(in)
    if (backward) then
       direction=FFTW_BACKWARD
    else
       direction=FFTW_FORWARD
    endif
    
    !----------------------------------------------------
    ! here, we determine how many thready we have available
    !----------------------------------------------------  
    use_fftw_omp=.false.
    if (present(ncpu)) then
       use_fftw_omp=.true.

       !$omp parallel private(id)
       id = omp_get_thread_num()
       !$omp barrier
       if ( id == 0 ) then
          nthreads = omp_get_num_threads()
          if (ncpu==0) ncpu=nthreads !if 0 use all threads
          write (*,*) "fftw wrapper: using nthreads of total OMP_threads=", ncpu,nthreads
       end if
       !$omp end parallel  
  
       call dfftw_init_threads(iret)
       if (iret==0) call dfftw_plan_with_nthreads(ncpu)
    endif
    print*, 'fftw plan: direction,FFTW_FORWARD, N,use_fftw_omp',direction,FFTW_FORWARD,use_fftw_omp
    call dfftw_plan_dft_1d(p,N,in,out,direction,FFTW_ESTIMATE)

   make_plan=p
   return 

  end function make_plan
!=============================

subroutine fftw_complex(arr,backward, pp)
  COMPLEX(DPC), dimension(1:) :: arr
  logical :: backward
  integer*8 :: p
  integer*8, optional::pp

  if (.not. present(pp)) then
     p = make_plan(arr, arr, backward)
  else
     p=pp
  endif
  call dfftw_execute_dft(p, arr, arr)

end subroutine fftw_complex


!======================
subroutine fftw_complex_lm(arr, backward, fout)
  COMPLEX(DPC), dimension(0:,0:,0:) :: arr
  COMPLEX(DPC), dimension(:), allocatable :: temp
  REAL(DP) :: norm
  logical :: backward,flag
  integer*8 :: p
  integer :: i, j, k,l,m,id, lmax, nshell, unit
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  character(len=*), optional :: fout
 
  !applying the forward and then the backward transform will multiply the input by n
  lmax=size(arr(0,:,0))-1
  nshell=size(arr(:,0,0))

  !omp does not like input args
  flag=backward

  if (backward) then
     !norm=1d0
     norm=1d0/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
     allocate(temp(0:nshell-1))
  else
     norm=1d0
     !allocate(temp(-3:nshell-1+3)) !with padding
     allocate(temp(0:nshell-1))
  endif

  temp=0
  temp(0:nshell-1) = arr(0:nshell-1,lmax,lmax)

  !create fftw plan
  !p = make_plan(temp, temp, backward) 
  
  !--------- excute fftw in OMP parellel ----
  !$OMP PARALLEL DEFAULT(PRIVATE), &     
  !$OMP& SHARED(arr,p,lmax, nshell, norm), private(id)
  id = omp_get_thread_num()
  if (id==0) print*, 'fftw_complex_lm, lmax, nshell,nthreads', lmax, nshell,OMP_GET_NUM_THREADS()

  !$OMP DO PRIVATE(l,m,temp,flag)
  do m=0,lmax
     do l=m, lmax
        temp=0
        temp(0:nshell-1)=arr(0:nshell-1,l,m) 
        !excute fft
        !call dfftw_execute_dft(p, temp, temp)
        !if to use healpix fft
        call complex_fft(temp,flag)
        arr(0:nshell-1,l,m)=temp(0:nshell-1)*norm
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL     
  !-----------------------------------------------
  deallocate(temp)

  print*, 'fftw done: sum(almn):',sum(arr)  
  !call dfftw_destroy_plan(p)

  !save almn if fout is given
  if (present(fout)) then
     unit=12
     write(*,*), 'save almn to file: ',trim(fout)
     open(unit,file=trim(fout),status='unknown',form='unformatted')
     write(unit) arr
     close(unit)
  endif

end subroutine fftw_complex_lm

!======================
subroutine fftw_filter_complex_lm(arr, gln,backward, fout)
  REAL(DP), DIMENSION(0:,0:), INTENT(IN) :: gln
  COMPLEX(DPC), dimension(0:,0:,0:) :: arr
  COMPLEX(DPC), dimension(:), allocatable :: temp
  REAL(DP) :: norm
  logical :: backward,ifft_flag
  integer*8 :: p
  integer :: i, j, k,l,m,id, lmax, nshell, unit, nnmax
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  character(len=*), optional :: fout
  character(len=*), parameter :: code='FFT_ALMN_GLN'

  !applying the forward and then the backward transform will multiply the input by n
  lmax=size(arr(0,:,0))-1
  nshell=size(arr(:,0,0))
  nnmax=size(gln(0,:))
  if (nnmax < nshell) stop code//' nnmax /= nnshell' !,nnmax, nshell

  !omp does not like input args
  ifft_flag=backward

  if (backward) then
     !norm=1d0
     norm=1d0/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
     allocate(temp(0:nshell-1))
  else
     norm=1d0
     !allocate(temp(-3:nshell-1+3)) !with padding
     allocate(temp(0:nshell-1))
  endif

  temp=0
  temp(0:nshell-1) = arr(0:nshell-1,lmax,lmax)

  !create fftw plan
  !p = make_plan(temp, temp, backward) 

  !--------- excute fftw in OMP parellel ----
  !$OMP PARALLEL DEFAULT(PRIVATE), &     
  !$OMP& SHARED(arr,gln,p,lmax, nshell, norm), private(id,ifft_flag)
  id = omp_get_thread_num()
  if (id==0) print*, 'fftw_complex_lm, lmax, nshell,nthreads', lmax, nshell,OMP_GET_NUM_THREADS()

  if (ifft_flag) then !backward transform

     !$OMP DO PRIVATE(l,m,temp)
     do m=0,lmax
        do l=m, lmax
           temp=0
           temp(0:nshell-1)=arr(0:nshell-1,l,m)*gln(l,0:nshell-1) 
           !excute fft
           !call dfftw_execute_dft(p, temp, temp)
           !if to use healpix fft
           call complex_fft(temp,ifft_flag)
           arr(0:nshell-1,l,m)=temp(0:nshell-1)*norm
        enddo
     enddo
     !$OMP END DO

  else !forward transform

     !$OMP DO PRIVATE(l,m,temp)
     do m=0,lmax
        do l=m, lmax
           temp=0
           temp(0:nshell-1)=arr(0:nshell-1,l,m) 
           !excute fft
           !call dfftw_execute_dft(p, temp, temp)
           !if to use healpix fft
           call complex_fft(temp,ifft_flag)
           arr(0:nshell-1,l,m)=temp(0:nshell-1)*gln(l,0:nshell-1)*norm
        enddo
     enddo
     !$OMP END DO

  endif

  !$OMP END PARALLEL     
  !-----------------------------------------------
  deallocate(temp)

  print*, 'fftw done: sum(almn):',sum(arr)  
  !call dfftw_destroy_plan(p)

  !save almn if fout is given
  if (present(fout)) then
     unit=12
     write(*,*), 'save almn to file: ',trim(fout)
     open(unit,file=trim(fout),status='unknown',form='unformatted')
     write(unit) arr
     close(unit)
  endif

end subroutine fftw_filter_complex_lm

!==============================

!======================
subroutine ifft_almn_gln(arr, arr_out, gln,fout)
  COMPLEX(DPC), dimension(0:,0:,0:),intent(in) :: arr
  COMPLEX(DPC), dimension(0:,0:,0:,0:),intent(out) :: arr_out
  REAL(DP), dimension(0:,0:,0:),intent(in) :: gln

  !
  COMPLEX(DPC), dimension(:), allocatable :: temp
  REAL(DP) :: norm
  logical :: backward
  integer*8 :: p
  integer ::lmax, nshell,nj, nnmax
  integer :: i, j, k,l,m,id, unit
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

  character(len=*), optional :: fout
  character(len=*), parameter :: code='IFFT_ALMN_GLN'

  !Note:  it seems omp does not like input args

  !applying the forward and then the backward transform will multiply the input by n
  lmax=size(arr(0,:,0))-1
  nshell=size(arr(:,0,0))
  nnmax=size(gln(0,:,0))
  nj=size(gln(0,0,:))


  if (nnmax < nshell) stop code//' nnmax /= nnshell' !,nnmax, nshell

  !only backward transform. 
  backward=.true.

  if (backward) then
     !norm=1d0
     norm=1d0/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
     allocate(temp(0:nshell-1))
  else
     norm=1d0
     !allocate(temp(-3:nshell-1+3)) !with padding
     allocate(temp(0:nshell-1))
  endif

  temp=0
  temp(0:nshell-1) = arr(0:nshell-1,lmax,lmax)

  !create fftw plan
  !p = make_plan(temp, temp, backward) 

  !--------- excute fftw in OMP parellel ----
  !$OMP PARALLEL DEFAULT(PRIVATE), &     
  !$OMP& SHARED(arr,arr_out, gln,p,lmax,nj, nshell, norm), private(id)
  id = omp_get_thread_num()
  if (id==0) print*, code,' lmax, nshell,norm,nthreads', lmax, nshell,norm,OMP_GET_NUM_THREADS()

  !$OMP DO PRIVATE(l,m,temp,j)
  do m=0,lmax
     do l=m, lmax
        do j=0,nj-1
           temp=0
           temp(0:nshell-1)=arr(0:nshell-1,l,m)*gln(l,0:nshell-1,j) !filter
           !if to use healpix fft
           call complex_fft(temp,.true.) !only backward transform
           arr_out(0:nshell-1,j,l,m)=temp(0:nshell-1)*norm
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL     
  !-----------------------------------------------
  print*, 'last temp.sum',sum(temp)
  deallocate(temp)

  print*, code//': done:,sum(almn), sum(almnj),gln.sum:',sum(arr),sum(arr_out),sum(gln)

  !save almn if fout is given
  if (present(fout)) then
     unit=12
     write(*,*), 'save almn to file: ',trim(fout)
     open(unit,file=trim(fout),status='unknown',form='unformatted')
     write(unit) arr_out
     close(unit)
  endif

end subroutine ifft_almn_gln


!================

end module fftw_wrapper
