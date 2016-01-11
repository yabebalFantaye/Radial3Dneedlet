  
  
!!!IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
  use healpix_types
  use alm_tools
  use pix_tools
  use extension

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,nj,large,first_j,last_j
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: map_TQU
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(SP) :: fwhm,sigma
  CHARACTER(LEN=128) :: healpixdir,dir,filenum
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)



  call getEnvironment("HEALPIX",healpixdir)
  dir=''!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  open(12,file=trim(adjustl(dir))//'nside_beta.unf',status='old',form='unformatted')
  read(12) nside
  close(12)
  open(12,file=trim(adjustl(dir))//'lmax_beta.unf',status='old',form='unformatted')
  read(12) lmax
  close(12)
  open(12,file=trim(adjustl(dir))//'nj.unf',status='old',form='unformatted')
  read(12) nj
  close(12)
  open(12,file=trim(adjustl(dir))//'large.unf',status='old',form='unformatted')
  read(12) large
  close(12)

  polar=0

  npix=nside**2*12

  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=0
  if (me.ne.0) first_j=sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1

  !print*,'proc:',me,nj_pp(me),first_j,last_j

  ALLOCATE(map_TQU(0:npix-1,1:1+2*polar),alm_TGC(1:1+2*polar,0:lmax,0:lmax))
  ALLOCATE(w8ring_TQU(1:2*nside,1:1+2*polar))


  w8ring_TQU=1d0
  open(12,file=trim(adjustl(dir))//'w8_beta.unf',status='old',form='unformatted')
  read(12) w8ring_TQU
  close(12)
  zbounds = (/ -1.0_dp, 1.0_dp /)

  do j=0,nj_pp(me)-1
     !print*,me,first_j-j0+j

     WRITE(filenum,'(I4)') first_j+j+1000
     if (large.eq.0) open(12,file=trim(adjustl(dir))//'beta'//filenum(2:4)//'.unf',status='old',form='unformatted')
     if (large.eq.1) open(12,file=trim(adjustl(dir))//'large/beta'//filenum(2:4)//'.unf',status='old',form='unformatted')
     read(12) map_TQU(:,1)
     close(12)

     !print*,'map2alm:',me,j
     if (polar.eq.0) call map2alm(nside, lmax, lmax, map_TQU(:,1), alm_TGC, zbounds, w8ring_TQU)
     !if (polar.eq.1) call map2alm(nside, lmax, lmax, map_TQU, alm_TGC, zbounds, w8ring_TQU)

     if (large.eq.0) open(12,file=trim(adjustl(dir))//'betaalm'//filenum(2:4)//'.unf',status='unknown',form='unformatted')
     if (large.eq.1) open(12,file=trim(adjustl(dir))//'large/betaalm'//filenum(2:4)//'.unf',status='unknown',form='unformatted')
     write(12) alm_TGC
     close(12)

  enddo


  CALL MPI_FINALIZE(ierr)


end program

