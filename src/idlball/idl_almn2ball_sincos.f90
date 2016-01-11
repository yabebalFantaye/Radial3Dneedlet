program idl_almn2ball

  !yabebal fantaye 12 Jan 2013

  use healpix_types
  use pix_tools
  use extension
  use alm_tools
  use wav_ball_mod
  use sincm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large,itype
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: map_TQU_pp,map_TQU
  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(DP) :: fwhm,sigma
  CHARACTER(LEN=128) :: healpixdir,dir,str,stype
  CHARACTER(LEN=200) :: fname

  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src, s0, sn,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp

  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almn_TGC
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrpix
  real(dp) ::rmin,rmax,dr,rmid
  real(dp):: s2n, s2nm1,y


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me


  call getEnvironment("HEALPIX",healpixdir)
  dir='large/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  !read parameters
  open(12,file=trim(adjustl(dir))//'params_almn2ball.unf',status='old',form='unformatted')
  read(12) large
  read(12) nside
  read(12) polar
  read(12) lmax
  read(12) nnmax
  read(12) nrpix
  read(12) rmin,rmax
  close(12)

  nshell=nrpix-1

  !check correctness of input parameters
  if (me == 0) then
     print*, '************** almn2ball.f90 *****'
     print*, '---------------------------'
     print*, 'large = ', large 
     print*, 'nside = ',nside 
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, 'rmin, rmax = ',real(rmin),real(rmax)
     print*, '---------------------------'

     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'

     if ( rmin == rmax ) then
        stop 'idl_ball2almn : rmin=rmax'
     end if

     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be even and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif
  polar=0  


  npix = nside2npix(nside)

  !define per processor task array
  ALLOCATE(nmap_pp(0:ntasks-1), nct_pp(0:ntasks-1))
  nmap_pp=nrpix/ntasks
  if(mod(nrpix,ntasks).ne.0) nmap_pp(1:mod(nrpix,ntasks)-1)=nmap_pp(1:mod(nrpix,ntasks)-1)+1
  nct_pp = nmap_pp*npix*(1+2*polar)

  !make r array
  allocate(rvec(0:nshell))
  !call r2dr_vec (nshell, rmin, rmax, rvec)
  call r8vec_even (nshell, rmin, rmax, rvec)


  sn=nmap_pp(me)

  ALLOCATE(map_TQU_pp(0:npix-1,1:1+2*polar,0:sn-1))
  allocate(almn_TGC(1:1+2*polar,0:lmax,0:lmax,0:nnmax))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))
  !allocate big array to collect all maps at the master node 
  if (me==0) ALLOCATE(map_TQU(0:npix-1,1:1+2*polar,0:nshell))

  !allocate(alm_TGC_r(0:nshell-1,1:1+2*polar,0:lmax,0:lmax))

  do itype=2,3

     intflag=itype
     write(stype,*) itype
     if (me==0) print*, 'reading beta and reconstructing ball for itype=',itype

     !load almn array
     open(unit,file=trim(adjustl(dir))//'almn_stype'//trim(adjustl(stype))//'.unf',status='old',form='unformatted')
     read(unit) almn_TGC
     close(unit)


     if  (me == 0) print*, 'calling alm2map in MPI distributed loop.'


     do ii=0,nmap_pp(me)-1

        if (me.ne.0) then
           i=ii+sum(nmap_pp(0:me-1))
        else
           i=ii
        endif
        !r value for the ith radial pixel
        rmid = rvec(i)
        if  (me == 0) print*, 'almn2ball: me, i, rmid = ',me,i, rmid


        alm_TGC = 0
        !
        if (i >= nshell) then !set by hand the map at the 0th shell to be 0
           map_TQU_pp(:,1,ii)=0

        else

           if (itype==2) then 
              do jj=0,nnmax
                 y=jj*argpi
                 alm_TGC = alm_TGC + almn_TGC(:,:,:,jj)*sin(y*rmid)
              enddo
           else
              do jj=0,nnmax
                 y=jj*argpi
                 alm_TGC = alm_TGC + almn_TGC(:,:,:,jj)*cos(y*rmid)
              enddo
           endif

           if (polar.eq.0) call alm2map(nside,lmax,lmax,alm_TGC,map_TQU_pp(:,1,ii)) 
           if (polar.eq.1) call alm2map(nside,lmax,lmax,alm_TGC,map_TQU_pp(:,:,ii)) 
        endif

        write(str,*)i
        fname = trim(adjustl(dir))//'almr_stype'//trim(adjustl(stype))//'.unf_'//trim(adjustl(str)) 
        print*, '******* almn2ball: saving alm_r array to **'
        print*, trim(fname)
        print*, '---------------------------'

        open(unit,file=trim(fname),status='unknown',form='unformatted')
        write(unit) alm_TGC
        close(unit)

     enddo


     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     if (me==0) print*, 'MPI send revcieve maps from work processors'

     if (ntasks>1) then

        dest=0
        if (me .ne. 0) then
           tag = me
           cnt=nct_pp(me)
           call  MPI_SSEND(map_TQU_pp,cnt,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr)
        else
           do i=0,ntasks-1

              if (i==0) then
                 s0 = 0
                 sn = nmap_pp(0)-1
                 map_TQU(:,:,s0:sn)=map_TQU_pp
              else

                 src = i
                 tag = i
                 cnt=nct_pp(i)
                 s0 = sum(nmap_pp(0:i-1))
                 sn = sum(nmap_pp(0:i))-1

                 call MPI_RECV(map_TQU(:,:,s0:sn),cnt,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,stat,ierr)
              endif
           enddo
        endif

     else
        map_TQU = map_TQU_pp
     endif

     if (me==0) then

        !if (large.eq.1) 
        fname=trim(adjustl(dir))//'ballmap_stype'//trim(adjustl(stype))//'.unf'
        !if (large.eq.0) fname=trim(adjustl(dir))//'ballmap.unf'

        print*, '************* saving ballmap array to **'
        print*, trim(fname)
        print*, '---------------------------'

        open(unit,file=trim(fname),status='unknown',form='unformatted')
        write(12) npix,nrpix,4,1
        write(12) map_TQU
        close(12)
     endif


  enddo

  !free alm memory
  deallocate(alm_TGC,almn_TGC)

  CALL MPI_FINALIZE(ierr)

end program idl_almn2ball

