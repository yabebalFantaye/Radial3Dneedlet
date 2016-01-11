program idl_ball2an

  !yabebal fantaye 12 Jan 2013

  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod

  IMPLICIT NONE

!  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large,itype

  REAL(DP), DIMENSION(:), ALLOCATABLE :: an
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: anm

  REAL(DP), DIMENSION(:), ALLOCATABLE :: funcr

  CHARACTER(LEN=128) :: healpixdir,dir,str,stype
  CHARACTER(LEN=200) :: fname


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src,sn,s0,stat

  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili,wfili2
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix
  real(dp) ::rmin,rmax,dr,rmid,y,y2
  real(dp):: s2n, s2nm1,xxx,anorm

  ! CALL MPI_INIT(ierr)

  ! CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  me = 0
  ntasks=1

  unit=51+me


  call getEnvironment("HEALPIX",healpixdir)
  dir='large/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  open(12,file=trim(adjustl(dir))//'params_ball2an.unf',status='old',form='unformatted')
  read(12) large
  read(12) nnmax
  read(12) nshell
  read(12) rmin,rmax
  read(12) itype
  close(12)

  write(stype,*)itype

  !check correctness of input parameters
  if (me == 0) then
     print*, '**************** ball2an.f90 *****'
     print*, '---------------------------'
     print*, 'large = ', large 
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, 'rmin, rmax = ',real(rmin),real(rmax)
     print*, 'itype = ',itype
     print*, '---------------------------'


     if ( rmin == rmax ) then
        stop 'idl_ball2almn : rmin=rmax'
     end if

     ! if ( nrsample < 1  ) then
     !    write ( *, '(a)' ) 'WARNING NRSAMPLE (per shell) must be >= 1. Assuming max(0,nrsample+1)'
     !    nrsample = 1
     ! end if

     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be eve and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif


  nrpix = nshell+1 !*(nrsample)


  ALLOCATE(funcr(0:nshell))

  open(unit,file=trim(adjustl(dir))//'ballmap.unf',status='old',form='unformatted')
  read(unit) funcr
  close(unit)


  allocate(an(0:nnmax))



  print*, 'nnmax, nrpix = ',nnmax, nrpix
  print*, 'ntasks, nshell = ',nshell


  !allocate Filon-simpson quadrature nodes and weights
  allocate(rfili(0:nshell),wfili(0:nshell),wfili2(0:nshell))


  print*, 'starting filon_sin integration, itype=',itype


  an = 0
  do jj = 0,nnmax

     y = jj*argpi

     if (itype==1) y = y/pifac
     anorm=1d0
     if (jj==0) anorm=2d0

     !print*, 'y, argpi',y,argpi
     call fil_simp_points(rmin,rmax,y,nshell,1,rfili,wfili,itype=itype) !sin(xy)/xy 


     do ii=0,nshell-1
        if (itype==1) then
           an(jj) =  an(jj) + rfili(ii)**2*y*wfili(ii)*funcr(ii) 
        else
           an(jj) =  an(jj) + wfili(ii)*funcr(ii)/anorm
        endif
     enddo


  enddo

  an = 2d0*an
  !print*, 'ball2an: an = ',an

  fname = trim(adjustl(dir))//'an_ball2an_stype'//trim(adjustl(stype))//'.unf' 
  print*, '******* saving an array to **'
  print*, trim(fname)
  print*, '---------------------------'

  !print*, size(an)
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) an
  close(unit)

  !!=========== check identitiy ==========

  allocate(anm(0:nnmax,0:nnmax))
  anm = 0


  do kk = 0,nnmax
     y2 = kk*argpi
     call fil_simp_points(rmin,rmax,y2,nshell,1,rfili,wfili2,itype=itype) !sin(xy)/xy 
     
     do jj = 0,nnmax
        
        y = jj*argpi
        call fil_simp_points(rmin,rmax,y,nshell,1,rfili,wfili,itype=itype) !sin(xy)/xy 
        
        
        do ii=0,nshell-1
           if (itype==1) anm(jj,kk) =  anm(jj,kk) + 2d0*wfili(ii)*sinc(y2*rfili(ii))*y*y2*rfili(ii)**2 !identitiy test !*funcr(ii) !
           if (itype==2) anm(jj,kk) =  anm(jj,kk) + 2d0*wfili(ii)*sin(y2*rfili(ii)) 
           if (itype==3) anm(jj,kk) =  anm(jj,kk) + 2d0*wfili(ii)*cos(y2*rfili(ii)) 
        enddo
        
        
     enddo
  enddo



  fname = trim(adjustl(dir))//'anm_check_stype'//trim(adjustl(stype))//'.unf' 
  print*, '******* saving an array to **'
  print*, trim(fname)
  print*, '---------------------------'
  print*,'*** load it as'
  print*, 'anm=dblarr(nnmax+1,nnmax+1) '
  print*, 'runf, anm, '//trim(fname)
  print*, '---------------------------'

  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) anm
  close(unit)

  print*, '******* saving an array done! **'

end program idl_ball2an

