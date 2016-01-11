program idl_ball2almn

  !yabebal fantaye 12 Jan 2013

  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE :: map_TQU
  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(SP) :: fwhm,sigma
  CHARACTER(LEN=128) :: healpixdir,dir,str
  CHARACTER(LEN=200) :: fname


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src,sn,s0
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,cnt_pp

  COMPLEX(SPC), DIMENSION(:,:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC,alm_TGC_rpp
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix
  real(dp) ::rmin,rmax,dr,rmid
  real(dp):: s2n, s2nm1

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me


  call getEnvironment("HEALPIX",healpixdir)
  dir=''!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  open(12,file=trim(adjustl(dir))//'params_ball2almn.unf',status='old',form='unformatted')
  read(12) large
  read(12) nside
  read(12) polar
  read(12) lmax
  read(12) nnmax
  read(12) nshell
  read(12) rmin,rmax
  read(12) nrsample
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
     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'

     if ( rmin == rmax ) then
        stop 'idl_ball2almn : rmin=rmax'
     end if

     if ( nrsample < 1  ) then
        write ( *, '(a)' ) 'WARNING NRSAMPLE (per shell) must be >= 1. Assuming max(0,nrsample+1)'
        nrsample = 1
     end if

     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be eve and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif

  polar=0  
  npix=nside2npix(nside) !**2*12
  nrpix = nshell*(nrsample+1)

  ALLOCATE(nmap_pp(0:ntasks-1),cnt_pp(0:ntasks-1))
  nmap_pp=nshell/ntasks
  if(mod(nshell,ntasks).ne.0) nmap_pp(0:mod(nshell,ntasks)-1)=nmap_pp(0:mod(nshell,ntasks)-1)+1
  

  ALLOCATE(map_TQU(0:npix-1,1:1+2*polar,0:nshell-1))
  ALLOCATE(w8ring_TQU(1:2*nside,1:1+2*polar))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))

  w8ring_TQU=1d0
  open(12,file=trim(adjustl(dir))//'w8.unf',status='old',form='unformatted')
  read(12) w8ring_TQU
  close(12)
  zbounds = (/ -1.0_dp, 1.0_dp /)

  open(12,file=trim(adjustl(dir))//'ballmap.unf',status='old',form='unformatted')
  read(12) map_TQU
  close(12)

  sn=nmap_pp(me)
  allocate(alm_TGC_rpp(0:sn-1,1:1+2*polar,0:lmax,0:lmax))
  allocate(alm_TGC_r(0:nrpix-1,1:1+2*polar,0:lmax,0:lmax))
  allocate(almn_TGC(1:nnmax,1:1+2*polar,0:lmax,0:lmax))
  almn_TGC = 0
  alm_TGC=0


  if  (me == 0) then
     print*, 'nnmax, nrpix = ',nnmax, nrpix
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell, nmap_pp
     print*, 'calling map2alm in MPI distributed loop'
  endif

  do ii=0,nmap_pp(me)-1

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif

     !print*, 'nshell = , max(map) =  ', i, maxval(map_TQU(:,1,i))

     if (polar.eq.0) call map2alm(nside, lmax, lmax, map_TQU(:,1,i), alm_TGC, zbounds, w8ring_TQU)
     if (polar.eq.1) call map2alm(nside, lmax, lmax, map_TQU(:,:,i), alm_TGC, zbounds, w8ring_TQU)

     alm_TGC_rpp(ii,:,:,:) = alm_TGC

     !print*, 'done with nshell =  from jj:kk= ',i,(nrsample+1)*i,(nrsample+1)*(i+1)-1
  enddo
  
  !deallocate(map_TQU, alm_TGC)

  !get radial values (min and max of each shell) at nrpix points
  allocate(rvec(0:nrpix-1))
  call r8vec_even (nrpix, rmin, rmax, rvec )
  !if  (me == 0) print*, 'rvec: ',REAL(rvec)

  !r spacing
  dr = ( rmax - rmin ) / real( nrpix-1, dp )

  !allocate coeficinets for Filon type quadrature 
  allocate(alpha(nnmax),beta(nnmax),gamma(nnmax))
  call filon_params(nnmax,dr,alpha,beta,gamma)

  if  (me == 0) then 
     print*, 'starting filon_sin integration'
  endif

  !multiply alm(r) by r so that the integrand would be f(r)*sin(n*pi*r)*dr
  do ii=0,nrpix-1
     alm_TGC_r(ii,:,:,:) = alm_TGC_r(ii,:,:,:)*rvec(ii)
  enddo

  do jj = 1,nnmax

     do ii=0,(nrpix-1)/2

        oddInd = 2*ii !zero based and knowning nrpix is odd
        evenInd = 2*(ii+1)-1 !zero based 
        

        if (oddInd <=nrpix-1) almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + dr * ( beta(jj)*alm_TGC_r(oddInd,:,:,:)*sin(PI*jj*rvec(oddInd)) ) !odd sum
        if (evenInd <=nrpix-1) almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + dr * ( gamma(jj)*alm_TGC_r(evenInd,:,:,:)*sin(PI*jj*rvec(evenInd)) )  !even sum

!        if (jj==1 .and. me==0 .and. oddInd <=nrpix-1) print*, 'odd ind, r, sin(pi*n*r)', oddInd, real(rvec(oddInd)), real(sin(PI*jj*rvec(oddInd))),   dr * ( beta(jj)*alm_TGC_r(oddInd,1,1,1)*sin(PI*jj*rvec(oddInd)) )
!        if (jj==1 .and. me==0 .and. evenInd <=nrpix-1) print*, 'even ind, r, sin(pi*n*r)', evenInd,real(rvec(evenInd)), real(sin(PI*jj*rvec(evenInd))),  alm_TGC_r(evenInd,:,:,:)
        
     enddo

      almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) -  dr * ( &
           beta(jj)/2d0*( alm_TGC_r(0,:,:,:)*sin(PI*jj*rvec(0)) + alm_TGC_r(nrpix-1,:,:,:)*sin(PI*jj*rvec(nrpix-1)) ) &
          + alpha(jj)*( alm_TGC_r(0,:,:,:)*cos(PI*jj*rvec(0)) - alm_TGC_r(nrpix-1,:,:,:)*cos(PI*jj*rvec(nrpix-1)) ) )

  enddo


  !write(str,*)ii
  !'//trim(adjustl(str))//'
  if (me==0) then
     if (large.eq.1)  fname = trim(adjustl(dir))//'large/almn.unf'
     if (large.eq.0)  fname = trim(adjustl(dir))//'almn.unf' 
     print*, '******* saving almn array to **'
     print*, trim(fname)
     !print*, 'size(almn), shape(almn) = ',size(almn_TGC), shape(almn_TGC)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) almn_TGC
     close(unit)
  endif

  CALL MPI_FINALIZE(ierr)

end program idl_ball2almn



!central r value for the ith radial pixel
!i = ii-1
!rmid = (rvec(i+1) + rvec(i))/2d0
!almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + alm_TGC_r(i,:,:,:)*sin(PI*jj*rmid)*rmid*dr
