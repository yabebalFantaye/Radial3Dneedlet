program idl_an2ball

  !yabebal fantaye 12 Jan 2013

  use healpix_types
  use pix_tools
  use extension
  use alm_tools
  use wav_ball_mod

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large,itype
  REAL(DP), DIMENSION(:), ALLOCATABLE :: funcr
  CHARACTER(LEN=128) :: healpixdir,dir,str,stype
  CHARACTER(LEN=200) :: fname

  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src, s0, sn,stat

  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,an
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,kk
  real(dp) ::rmin,rmax,dr,rmid,y
  real(dp):: s2n, s2nm1,xxx


  ! CALL MPI_INIT(ierr)

  ! CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  me = 0
  ntasks=1

  unit=15+me


  call getEnvironment("HEALPIX",healpixdir)
  dir='large/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  !read parameters
  open(12,file=trim(adjustl(dir))//'params_an2ball.unf',status='old',form='unformatted')
  read(12) large
  read(12) nnmax
  read(12) nshell
  read(12) itype
  read(12) rmin,rmax
  close(12)

  write(stype,*)itype

  !check correctness of input parameters
  if (me == 0) then
     print*, '************** almn2ball.f90 *****'
     print*, '---------------------------'
     print*, 'large = ', large 
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, 'itype = ',itype
     print*, 'rmin, rmax = ',real(rmin),real(rmax)
     print*, '---------------------------'


     if ( rmin == rmax ) then
        stop 'idl_ball2almn : rmin=rmax'
     end if

     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be even and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif


  !make r array
  allocate(rvec(0:nshell))
  call r8vec_even ( nshell, rmin, rmax, rvec )

  ALLOCATE(funcr(0:nshell))
  allocate(an(0:nnmax))



  !load almn array
  open(unit,file=trim(adjustl(dir))//'an_ball2an_stype'//trim(adjustl(stype))//'.unf',status='old',form='unformatted')
  read(unit) an
  close(unit)



  funcr = 0.

  do ii=0,nshell !nmap_pp(me)-1

     i = ii

     !central r value for the ith radial pixel
     !rmid = (rvec(i+1) + rvec(i))/2d0

     rmid=rvec(ii)

     do jj=0,nnmax

        y=jj*argpi
        !print*, 'an2ball: y, argpi',y,argpi

        !print*, 'n, rmid = , sin(pi*n*rmid)/rmid= ',jj, rmid,sinc
        if (itype==1) funcr(ii) = funcr(ii) + (an(jj)*y/pifac*sinc(y*rmid))
        if (itype==2) funcr(ii) = funcr(ii) + (an(jj)*sin(y*rmid))
        if (itype==3) funcr(ii) = funcr(ii) + (an(jj)*cos(y*rmid))

     enddo


  enddo

  !print*, funcr


  fname=trim(adjustl(dir))//'ballmap_an2ball_stype'//trim(adjustl(stype))//'.unf'
  
  print*, '************* saving ballmap array to **'
  print*, trim(fname)
  print*, '---------------------------'
  
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) funcr
  close(unit)
  



end program idl_an2ball

