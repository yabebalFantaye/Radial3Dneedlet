program test_filon_quad

  !yabebal fantaye 12 Jan 2013

  use healpix_types


  IMPLICIT NONE

  !  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large
  REAL(DP), DIMENSION(:), ALLOCATABLE :: f,ff,fff,rf,fn

  CHARACTER(LEN=128) :: healpixdir,dir,str
  CHARACTER(LEN=200) :: fname,fout


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp

  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix,rpaw
  real(dp) ::rmin,rmax,dr,rmid
  real(dp):: s2n, s2nm1

  !   CALL MPI_INIT(ierr)

  !   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  !   CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  me=0
  ntasks=1
  unit=12+me


!  call getEnvironment("HEALPIX",healpixdir)
  dir=''!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  if (iargc() > 0) then
     call getarg(1,fname)
  else
     stop 'usage test_filon_quad paramffile.txt'
  endif
  open(12,file=trim(adjustl(fname)),status='old')
  read(12,*) fout
  read(12,*) rpaw
  read(12,*) nnmax
  read(12,*) nshell
  read(12,*) rmin,rmax
  read(12,*) nrsample
  close(12)


  !check correctness of input parameters
  if (me == 0) then
     print*, '**************** ball2almn.f90 *****'
     print*, '---------------------------'
     print*, 'large = ', large 
     print*, 'nside = ',nside 
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, 'rmin, rmax = ',real(rmin),real(rmax)
     print*, 'nrsample = ',nrsample
     print*, '---------------------------'


     if ( rmin == rmax ) then
        stop 'idl_ball2almn : rmin=rmax'
     end if

     if ( nrsample < 0 .or. mod ( nrsample, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'WARNING NRSAMPLE (per shell) must be even and >= 0. Assuming max(0,nrsample+1)'
        nrsample = max(0,nrsample+1)
        write ( *, '(a,i8)' ) '  new nrsample  = ',nrsample
        stop
     end if

     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be evec and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif


  nrpix = nshell*(nrsample+1)

  if  (me == 0) then
     print*, 'nnmax, nrpix = ',nnmax, nrpix
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell
     print*, 'calling map2alm in MPI distributed loop'
  endif


  ALLOCATE(f(0:nshell), rf(0:nshell), fff(0:nshell), fn(nnmax))

  !get radial values (min and max of each shell) at nrpix points
  allocate(rvec(0:nrpix),ff(0:nrpix))

  print*, 'starting filon_sin integration ff --> fn'


  dr = 1d0 / real( nrpix-1, dp )


  do jj = 1,nnmax

     call fil_simp_points(rmin,rmax,jj*PI,nshell,1,rvec,ff) !sin(xy)/xy

     f = rvec*rvec**rpaw !assume f = r*f where f=r
     fn(jj) = sum(ff*jj*PI*f*rvec)

  enddo

  print*, 'rvec =', rvec(0:10)

  fn = sqrt2*fn
  print*, 'starting inverse filon_sin  fn --> ff'

  fff = 0d0
  do ii=0,nshell
     do jj=1,nnmax
        fff(ii) = fff(ii) + fn(jj)*sqrt2*sin(PI*jj*rvec(ii)) 
     enddo
  enddo

!print*, 'fn = ',fn
  print*, '----------------------------'
  print*, '----------------------------'
  print*, 'orginal func (f): ',f
  print*, '----------------------------'
  print*, '----------------------------'
  print*, 'inv sync expansion func (fff): ',fff
  print*, '----------------------------'
  print*, '----------------------------'
  print*, 'abs(f - fff): ',abs(f - fff)/max(f,1e-10)
  print*, '----------------------------'

  print*, '----------------------------'
  


end program test_filon_quad


