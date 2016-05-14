module almnutil

  use healpix_types
  use alm_tools
  use fftw_wrapper
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp
  implicit none  


contains

!===============================
subroutine almn2cln(arr,cln, fname)
  COMPLEX(DPC), dimension(0:,0:,0:) :: arr
  REAL(DP), dimension(0:,0:) :: cln
  REAL(DP), dimension(:,:), allocatable :: cl
  integer :: i, j, k,l,m,id, lmax, nshell,unit
  character(len=*), optional :: fname
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

  lmax=size(arr(0,:,0))-1
  nshell=size(arr(:,0,0))

  allocate(cl(0:lmax, 1))
  cl=0d0
  cln=0d0

!  !$OMP PARALLEL DEFAULT(PRIVATE), &     
!  !$OMP& SHARED(arr, cln, lmax, nshell), private(id)
!  id = omp_get_thread_num()
!  if (id==0) print*, 'almn2cln, lmax, nshell', lmax, nshell

!  !$OMP DO PRIVATE(i, cl)
  do i=0,nshell-1
     !compute power spectrum from mean alm
     call alm2cl(lmax, lmax, arr(i:i, 0:lmax, 0:lmax),cl)
     cln(i,0:lmax) = cln(i,0:lmax) + cl(0:lmax,1)
  enddo
!  !$OMP END DO
!  !$OMP END PARALLEL     
  deallocate(cl)

  if (present(fname)) then
     unit=12
     print*, '******* saving 2d power spectrum cln **'
     print*, 'cln(3:10,3:10)',cln(3:10,3:10)
     print*, trim(fname)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
     write(unit) cln
     close(unit)
  endif

end subroutine almn2cln
!==============================

subroutine almr2cl(arr,meancl,fname)
  COMPLEX(DPC), dimension(0:,0:,0:) :: arr
  REAL(DP), dimension(0:,1:) :: meancl
  REAL(DP), dimension(:,:), allocatable :: tempcl
  COMPLEX(DPC), dimension(:,:,:), allocatable :: alm
  integer :: i, j, k,l,m,id, lmax, nshell,unit
  character(len=*), optional :: fname
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

  lmax=size(arr(0,:,0))-1
  nshell=size(arr(:,0,0))

  allocate(tempcl(0:lmax, 1))
  allocate(alm(1,0:lmax,0:lmax))

  alm=0
  tempcl=0d0
  meancl(:,1:2)=0d0

  !compute mean cl and mean alm 
  do i=0,nshell-1
     alm(1,0:lmax,0:lmax) = alm(1,0:lmax,0:lmax)+arr(i, 0:lmax, 0:lmax)
     !compute power spectrum at each shell
     call alm2cl(lmax, lmax, arr(i:i, 0:lmax, 0:lmax),tempcl)
     meancl(0:lmax,1) = meancl(0:lmax,1) + tempcl(0:lmax,1)
  enddo

  !compute power spectrum from mean alm
  print*, 'almr2cl: sum(alm), max(alm_r), min(alm_r), max(alm_i),min(alm_i)',sum(alm),maxval(real(alm)),minval(real(alm)),maxval(aimag(alm)),minval(aimag(alm))

  call alm2cl(lmax, lmax, alm,tempcl) 
  meancl(0:lmax,2)=tempcl(0:lmax,1)

  if (present(fname)) then
     unit=12
     print*, '******* saving mean cl **'
     print*, 'cl_mean_clr(3:10)',meancl(3:10,1)
     print*, 'cl_mean_alm_cl(3:10)',meancl(3:10,2)
     print*, trim(fname)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
     write(unit) meancl(:,1:2)
     close(unit)
  endif

  deallocate(tempcl, alm)

end subroutine almr2cl

!==========================
subroutine write_almx(arr, froot,pio)
  COMPLEX(DPC), dimension(0:,0:,0:) :: arr
  COMPLEX(DPC), dimension(:,:), allocatable :: alm
  integer :: i, j, k,l,m,id, lmax, nshell, unit
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  character(len=30) :: str
  character(len=2000) :: fname,fout
  character(len=*) :: froot
  logical, optional ::pio

  lmax=size(arr(0,:,0))-1
  nshell=size(arr(:,0,0))

  !omp does not like input/optional params
  fout=trim(froot)

  if (present(pio)) then 
     !write in parellal
     
     allocate(alm(0:lmax,0:lmax))
     
     !$OMP PARALLEL DEFAULT(PRIVATE), &
     !$OMP& SHARED(arr,nshell,lmax), private(id, unit)
     id = omp_get_thread_num()
     unit=12+id
     
     !$OMP DO PRIVATE(i,fname,fout,alm,str)
     do i=0,nshell-1
        write(str,*)i
        fname = trim(adjustl(fout))//'_'//trim(adjustl(str))
        write(*,*), 'save almn to file: ',trim(fname)
        open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
        alm(0:lmax,0:lmax)=arr(i,0:lmax,0:lmax)
        write(unit) alm
        close(unit)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     
     deallocate(alm)
  else 
     !serial all write
     unit=12
     fname = trim(adjustl(fout))
     write(*,*), 'save almn to file: ',trim(fname)
     open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
     write(unit) arr
     close(unit)
  endif

end subroutine write_almx

!================
subroutine read_almx(arr, froot,pio)
  COMPLEX(DPC), dimension(0:,0:,0:) :: arr
  COMPLEX(DPC), dimension(:,:), allocatable :: alm
  integer :: i, j, k,l,m,id, lmax, nshell, unit
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  character(len=30) :: str
  character(len=2000) :: fname,fout
  character(len=*) :: froot
  logical, optional ::pio

  lmax=size(arr(0,:,0))-1
  nshell=size(arr(:,0,0))

  !omp does not like input/optional params
  fout=trim(froot)

  if (present(pio)) then 
     !write in parellal
     
     allocate(alm(0:lmax,0:lmax))
     
     !$OMP PARALLEL DEFAULT(PRIVATE), &
     !$OMP& SHARED(arr,nshell,lmax), private(id, unit)
     id = omp_get_thread_num()
     unit=12+id
     
     !$OMP DO PRIVATE(i,fname,fout,alm,str)
     do i=0,nshell-1
        write(str,*)i
        fname = trim(adjustl(fout))//'_'//trim(adjustl(str))
        write(*,*), 'read almn/r from file: ',trim(fname)
        open(unit,file=trim(fname),status='old',form='unformatted', ACCESS="STREAM")
        read(unit) alm
        arr(i,0:lmax,0:lmax)=alm(0:lmax,0:lmax)
        close(unit)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     
     deallocate(alm)
  else 
     !serial all write
     fname = trim(adjustl(fout))
     write(*,*), 'read almn/r from file: ',trim(fname)
     open(unit,file=trim(fname),status='old',form='unformatted', ACCESS="STREAM")
     read(unit) arr
     close(unit)
  endif

end subroutine read_almx




!------------
end module almnutil
