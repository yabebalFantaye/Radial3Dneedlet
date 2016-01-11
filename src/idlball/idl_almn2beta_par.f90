  
!!!IDL extension to HEALPix [TESTVERSION] 2D needlet: Written by Frode K. Hansen, Nov. 2006
!! Modified for 3D needlet by Y. Fantaye Sep 2013 
  
  use healpix_types
  use pix_tools
  use extension
  use alm_tools
  use wav_ball_mod

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,nmax,i,j,l,m,polar,p,wav,large,m1,m2,multialm,lmaxx
  INTEGER(I4B) :: j0,nj,first_j,last_j, v0,nv,first_v, last_v, in
  REAL(SP), DIMENSION(:,:,:,:), ALLOCATABLE :: map_TQU
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln
  COMPLEX(SPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almn_TGC
  REAL(SP) :: fwhm,sigma
  REAL(DP) :: bb,norm
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)


  call getEnvironment("HEALPIX",healpixdir)
  dir=''!trim(adjustl(healpixdir))//'/src/idl/almtools/'
  open(12,file=trim(adjustl(dir))//'bb.unf',status='old',form='unformatted')
  read(12) bb
  close(12)
  open(12,file=trim(adjustl(dir))//'j0.unf',status='old',form='unformatted')
  read(12) j0
  close(12)
  open(12,file=trim(adjustl(dir))//'nj.unf',status='old',form='unformatted')
  read(12) nj
  open(12,file=trim(adjustl(dir))//'j0.unf',status='old',form='unformatted')
  read(12) v0
  close(12)
  open(12,file=trim(adjustl(dir))//'nj.unf',status='old',form='unformatted')
  read(12) nv
  close(12)
  open(12,file=trim(adjustl(dir))//'nside.unf',status='old',form='unformatted')
  read(12) nside
  close(12)
  open(12,file=trim(adjustl(dir))//'lmax.unf',status='old',form='unformatted')
  read(12) lmax
  close(12)
  open(12,file=trim(adjustl(dir))//'nmax.unf',status='old',form='unformatted')
  read(12) nmax
  close(12)

  ! open(12,file=trim(adjustl(dir))//'p.unf',status='old',form='unformatted')
  ! read(12) p
  ! close(12)
  ! open(12,file=trim(adjustl(dir))//'wav.unf',status='old',form='unformatted')
  ! read(12) wav
  ! close(12)
  ! open(12,file=trim(adjustl(dir))//'fwhm_beta.unf',status='old',form='unformatted')
  ! read(12) fwhm
  ! close(12)
  ! open(12,file=trim(adjustl(dir))//'large.unf',status='old',form='unformatted')
  ! read(12) large
  ! close(12)
  ! open(12,file=trim(adjustl(dir))//'m1.unf',status='old',form='unformatted')
  ! read(12) m1
  ! close(12)
  ! open(12,file=trim(adjustl(dir))//'m2.unf',status='old',form='unformatted')
  ! read(12) m2
  ! close(12)
  ! open(12,file=trim(adjustl(dir))//'multialm.unf',status='old',form='unformatted')
  ! read(12) multialm
  ! close(12)
  !open(12,file=trim(adjustl(dir))//'polar.unf',status='old',form='unformatted')
  !read(12) polar
  !close(12)

  !!necessary to be defined params
  polar=0
  wav=0
  p=0
  m1=1
  m2=0

  npix=nside**2*12
  nv=int(bb**(j0+nj-1) - bb**(j0) + 1)

  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=j0
  if (me.ne.0) first_j=j0+sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1

  first_v=v0
  last_v=first_v+nv

  !print*,'ntasks:',ntasks,nj_pp(me),first_j,last_j

  ALLOCATE(map_TQU(0:npix-1,1:1+polar*2,0:nv-1, 0:nj_pp(me)-1))
  ALLOCATE(gln(0:lmax,0:nmax,j0:j0+nj-1))
  ALLOCATE(almn_TGC(0:nmax, 1:1+polar*2,0:lmax,0:lmax))

  !sigma=fwhm/60./180.*pi/sqrt(8.*log(2.))
  !allocate(beam(0:lmax,1:1+polar*2))
  ! open(12,file=trim(adjustl(dir))//'beam_beta.unf',status='old',form='unformatted')
  ! read(12) beam
  ! close(12)

  open(12,file=trim(adjustl(dir))//'almn_beta.unf',status='old',form='unformatted')
  read(12) almn_TGC
  close(12)

  ! do l=0,lmax
  !    alm_TGC(1,l,0:l,:)=alm_TGC(1,l,0:l,:)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,1)
  !    if (polar.eq.1) then
  !       alm_TGC(2,l,0:l,:)=alm_TGC(2,l,0:l,:)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,2)
  !       alm_TGC(3,l,0:l,:)=alm_TGC(3,l,0:l,:)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,3)
  !    endif
  ! enddo

  if (wav.eq.0) wavtyp='standard'
  if (wav.eq.1) wavtyp='mexican'

  CALL calc_f2


  !print*,'alm2beta1:',j0,nj,lmax,bb,wavtyp,p
  CALL calc_gln(j0,nj,lmax,nmax,gln,bb,wavtyp,p)

  if ((m1.eq.1).and.(m2.eq.0)) then
     do l=0,lmax
        do in=0,nmax
           norm=sum(gln(l,in,:)**3)
           if (norm.gt.0d0) gln(l,in,:)=gln(l,in,:)/sqrt(norm)
        enddo
     enddo
  endif

  ! if (m2.eq.1) then
  !    norm=sum(gl(2:lmax/8,:)**2)/(REAL(lmax/8-2+1,DP))
  !    gl=gl/sqrt(norm)
  ! endif


  CALL synmap_multi_ball(almn_TGC,map_TQU,nside,nmax, lmax, bb,j0,gln(:,:, first_j:last_j),nj_pp(me),nv)

  do j=0,nj_pp(me)-1
     WRITE(filenum,'(I4)') first_j-j0+j+1000
     !nv =  bb**(first_j-j0+j)
     open(12,file=trim(adjustl(dir))//'beta_ball_j'//filenum(2:4)//'.unf',status='unknown',form='unformatted')
     write(12) map_TQU(:,1,0:nv-1,j)
     close(12)
  enddo

  !print*,'done',me

  CALL MPI_FINALIZE(ierr)

end program

