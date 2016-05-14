program idl_ball2almn

  !yabebal fantaye 12 Jan 2013

  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large,itype
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: map_TQU
  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(DP) :: fwhm,sigma
  CHARACTER(LEN=128) :: healpixdir,dir,str,stype
  CHARACTER(LEN=200) :: fname


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src,sn,s0,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp

  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix
  real(dp) ::rmin,rmax,dr,rmid,y
  real(dp):: s2n, s2nm1,anorm

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me


  call getEnvironment("HEALPIX",healpixdir)
  dir='large/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  open(12,file=trim(adjustl(dir))//'params_ball2almn.unf',status='old',form='unformatted')
  read(12) large 
  read(12) nside
  read(12) polar
  read(12) lmax
  read(12) nnmax
  read(12) nrpix
  read(12) rmin,rmax
  read(12) nrsample
  close(12)

  nshell=nrpix-1

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

     if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be eve and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif

  polar=0  
  npix=nside2npix(nside) 


  !task distribution across processors: each processor recieves nmap_pp(me) shells
  ALLOCATE(nmap_pp(0:ntasks-1),nct_pp(0:ntasks-1))
  nmap_pp=nrpix/ntasks
  if(mod(nrpix,ntasks).ne.0) nmap_pp(1:mod(nrpix,ntasks)-1)=nmap_pp(1:mod(nrpix,ntasks)-1)+1
  nct_pp = nmap_pp*(1+2*polar)*(lmax+1)*(lmax+1) !data count per each processor


  ALLOCATE(w8ring_TQU(1:2*nside,1:1+2*polar))


  w8ring_TQU=1d0
  open(12,file=trim(adjustl(dir))//'w8.unf',status='old',form='unformatted')
  read(12) w8ring_TQU
  close(12)
  zbounds = (/ -1.0_dp, 1.0_dp /)


  sn=nmap_pp(me)
  allocate(alm_TGC_rpp(0:sn-1,1:1+2*polar,0:lmax,0:lmax))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))
  alm_TGC=0


  if  (me == 0) then
     print*, 'nnmax, lmax, nrpix = ',nnmax, lmax, nrpix
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell, nmap_pp
     print*, 'calling map2alm in MPI distributed loop'
  endif

!   !All shells in one file case
!   ALLOCATE(map_TQU(0:npix-1,1:1+2*polar,0:nshell))
!   open(12,file=trim(adjustl(dir))//'ballmap.unf',status='old',form='unformatted')
!   read(12) map_TQU
!   close(12)

  !shell by shell read case
  ALLOCATE(map_TQU(0:npix-1,1:1+2*polar,0:nmap_pp(me)-1))


  !in each processor loop over radial shells, and do harmonic expansion
  do ii=0,nmap_pp(me)-1

     open(12,file=trim(adjustl(dir))//'ballmap.unf',status='old',form='unformatted')
     read(12) map_TQU
     close(12)

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif

     !print*, 'nshell = , max(map) =  ', i, maxval(map_TQU(:,1,i))
     if (i == 0) then !the 0th layer is an empty map
        alm_TGC=0.
     else
        if (polar.eq.0) call map2alm(nside, lmax, lmax, map_TQU(:,1,ii), alm_TGC, zbounds, w8ring_TQU)
        if (polar.eq.1) call map2alm(nside, lmax, lmax, map_TQU(:,:,ii), alm_TGC, zbounds, w8ring_TQU)
     endif
     alm_TGC_rpp(ii,:,:,:) = alm_TGC(:,:,:)

     !print*, 'done with nshell =  from jj:kk= ',i,(nrsample+1)*i,(nrsample+1)*(i+1)-1
  enddo

  deallocate(map_TQU)


  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (me==0) then
     !allocate big arrays
     allocate(alm_TGC_r(0:nshell,1:1+2*polar,0:lmax,0:lmax))
     alm_TGC_r = 0

     print*, 'MPI send revciev alm_TQU_rpp from work processors'
     print*, ntasks,me
  endif


  if (ntasks>1) then

     dest=0
     if (me .ne. 0) then
        tag = me
        cnt=nct_pp(me)
        print*, 'sending proc=',me
        call  MPI_SSEND(alm_TGC_rpp,cnt,MPI_DOUBLE_COMPLEX,dest,tag,MPI_COMM_WORLD,ierr)
     else

        do i=0,ntasks-1

           if (i==0) then
              s0 = 0
              sn = nmap_pp(0)-1
              alm_TGC_r(s0:sn,:,:,:)=alm_TGC_rpp
           else

              src = i
              tag = i
              cnt=nct_pp(i)
              s0 = sum(nmap_pp(0:i-1))
              sn = sum(nmap_pp(0:i))-1
              print*, 'recieving from proc=',i
              call MPI_RECV(alm_TGC_r(s0:sn,:,:,:),cnt,MPI_DOUBLE_COMPLEX,src,tag,MPI_COMM_WORLD,stat,ierr)
           endif
        enddo
        print*, 'MPI send revciev alm_TQU_rpp done!'        
     endif

  else
     alm_TGC_r = alm_TGC_rpp
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !print*,'deallocating alm_TGC_rpp ',me
  deallocate(alm_TGC_rpp)

  if  (me == 0) then

     !damp alm(r) to a file

     fname = trim(adjustl(dir))//'almr.unf' 
     print*, '******* saving alm_r array to **'
     print*, trim(fname)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) alm_TGC_r
     close(unit)

     print*, 'lmax, nrpix, nshell = ',lmax,nrpix,nshell
     !print*, 'size(almn), shape(almn) = ',size(alm_TGC), shape(alm_TGC)

     !allocate big array 
     allocate(almn_TGC(1:1+2*polar,0:lmax,0:lmax,0:nnmax))
     almn_TGC = 0

     !allocate Filon-simpson quadrature nodes and weights
     allocate(rfili(0:nshell),wfili(0:nshell))


     print*, 'starting filon_sin integration'

     do itype=2,3

        intflag=itype
        write(stype,*) itype
        if (me==0) print*, 'computing almn_ball for itype=',itype

        almn_TGC = 0

        do jj = 0,nnmax

           anorm=1d0
           if (jj==0) anorm=2d0

           y = argpi*jj

           call fil_simp_points(rmin,rmax,y,nshell,1,rfili,wfili,itype=itype) !sin or cos
           !alm_TGC = 0
           do ii=0,nshell-1
              almn_TGC(1,0:lmax,0:lmax,jj) = almn_TGC(1,0:lmax,0:lmax,jj) + wfili(ii)*alm_TGC_r(ii,1,0:lmax,0:lmax)/anorm
           enddo

        enddo

        almn_TGC = 2d0*almn_TGC  !2d0 normalization 


        !if (large.eq.1)  
        fname = trim(adjustl(dir))//'almn_stype'//trim(adjustl(stype))//'.unf'
        !if (large.eq.0)  fname = trim(adjustl(dir))//'almn.unf' 
        print*, '******* saving almn array to **'
        print*, trim(fname)
        print*, 'size(almn), shape(almn) = ',size(almn_TGC), shape(almn_TGC)
        print*, '---------------------------'
        open(unit,file=trim(fname),status='unknown',form='unformatted')
        write(unit) almn_TGC
        close(unit)

     enddo
  endif



  CALL MPI_FINALIZE(ierr)

end program idl_ball2almn



!central r value for the ith radial pixel
!i = ii-1
!rmid = (rvec(i+1) + rvec(i))/2d0
!almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + alm_TGC_r(i,:,:,:)*sin(PI*jj*rmid)*rmid*dr
