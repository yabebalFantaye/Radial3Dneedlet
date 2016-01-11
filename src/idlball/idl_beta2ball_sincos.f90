  
  
!!!IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod
  use sincm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,j0,nj,p,wav,nside_out
  INTEGER(I4B) :: npix_out,large,m1,m2,first_j,last_j,savealm
  INTEGER(I4B) :: nv, v0, nmax,nvj,nshell,ir, in,iv,nrpix, itype
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: map_slet
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: map_ball
  REAL(DP), DIMENSION(:), ALLOCATABLE :: map,map1,map2,map3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam,rvec

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE ::gln
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: arlm,arlm_all
  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: arlmj,almn,alm_r
  REAL(DP), DIMENSION(1:2) :: zbounds

  REAL(DP),DIMENSION(:,:),allocatable :: w8ring_map

  REAL(DP) :: fwhm,sigma,y
  REAL(DP) :: bb,norm,rmin,rmax,rmid,vval,blnj_sn,temp
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn,stype
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)



  call getEnvironment("HEALPIX",healpixdir)
  dir='temp/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  open(12,file=trim(adjustl(dir))//'params_beta2ball.unf',status='old',form='unformatted')
  read(12) bb
  read(12) j0
  read(12) nj
  read(12) nv
  read(12) nside
  read(12) nside_out
  read(12) polar
  read(12) lmax
  read(12) nmax
  read(12) nrpix
  read(12) rmin,rmax
  read(12) feedback
  close(12)


  nshell=nrpix-1
  if (nv < nshell) nv=nshell

  polar=0
  wav=0
  p=0
  m1=0
  m2=0

  npix=nside**2*12
  npix_out=nside_out**2*12


  if (me == 0) then
     print*, '**************** ball2almn.f90 *****'
     print*, '---------------------------'
     print*, 'itype (1-sinc, 2-sin, 3-cos, 4-sin+cos)= ',2,3
     print*, 'nside, nside_out = ',nside,nside_out
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nmax
     print*, 'nshell = ',nshell
     print*, 'rmin, rmax = ',real(rmin),real(rmax)
     print*, 'nv, nj, j0 = ',nv,nj,j0
     print*, '---------------------------'
     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'

     if ( rmin == rmax ) then
        stop 'idl_ball2beta : rmin=rmax'
     end if

     if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be eve and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif


  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=j0
  if (me.ne.0) first_j=j0+sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1

  if (feedback>2) print*,'proc:',me,nj_pp(me),first_j,last_j

  ALLOCATE( w8ring_map(2*nside,1+2*polar) )
  w8ring_map=1d0
  open(12,file=trim(adjustl(dir))//'w8.unf',status='old',form='unformatted')
  read(12) w8ring_map
  close(12)
  zbounds = (/ -1.0_dp, 1.0_dp /)


  if (feedback>2) print*,'beta map read!'
  if (wav.eq.0) wavtyp='standard'
  if (wav.eq.1) wavtyp='mexican'

  CALL calc_f2


  ALLOCATE(gln(0:lmax,0:nmax,j0:j0+nj-1))
  CALL calc_gln(j0,nj,lmax,nmax,gln,bb,wavtyp,p)

  if (feedback>2) print*,'gln computed!'


  
  do itype=2,3
     
     intflag=itype
     write(stype,*) itype
     if (me==0) print*, 'reading beta and reconstructing ball for itype=',itype
     

     ALLOCATE(map(0:npix-1),map2(0:npix-1), map_slet(0:npix-1,1:1+2*polar,0:nv,0:nj_pp(me)-1))
     ALLOCATE(map_ball(0:npix_out-1,1:1+2*polar,0:nshell))

     ALLOCATE(alm(1:1+2*polar,0:lmax,0:lmax))
     ALLOCATE(arlm(1:1+2*polar,0:lmax,0:lmax))
     allocate(almn(1:1+2*polar,0:lmax,0:lmax,0:nmax))
     allocate(arlm_all(1:1+2*polar,0:lmax,0:lmax))

     arlm=0.
     almn=0.
     alm=0.


     do j=0,nj_pp(me)-1
        !print*,me,first_j-j0+j
        WRITE(filenum,'(I4)') first_j-j0+j !+1000
        open(12,file=trim(adjustl(dir))//'beta_ball_stype'//trim(adjustl(stype))//'_j'//trim(adjustl(filenum))//'.unf',status='old',form='unformatted')

        !nvj=bb**(first_j-j0+j)
        read(12) map_slet(:,1,0:nv,j)
        close(12)
     enddo
     
     do in=0,nmax
        arlm=0.        
        
        do j=first_j,last_j
           
           map=0.
           do iv=0,nv-1
              vval= dble(iv)
              map = map + map_slet(:,1,iv, j-first_j)*snrv_sincos(itype,in, vval,j,bb,nv)
           enddo
           
           !map = map_slet(:,1,iv, j-first_j)
           call map2alm(nside, lmax, lmax, map, alm, zbounds, w8ring_map)
           
           do l=0,lmax
              arlm(1,l,0:l) = arlm(1,l,0:l) + alm(1,l,0:l)*gln(l,in,j)
              !print*, 'j,n,l,glnj,sn,sum_alm, sum_arlm: ',j,in,l, gln(l,in,j),temp,sum(abs(alm)),sum(abs(arlm))
           enddo
           
        enddo !end j
        
        cnt=(lmax+1)**2
        !if ( me.eq.0) arlm_all=arlm
        CALL MPI_ALLREDUCE(arlm,arlm_all,cnt,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
        
        almn(1,0:lmax,0:lmax,in) = arlm_all(1,0:lmax,0:lmax)
        
     enddo !end in
     
     deALLOCATE(map,map2, map_slet)
     deallocate(arlm_all)
     
  !=== for now compute sum over lmn in main processor ===

     !make r array
     allocate(rvec(0:nshell))
     !  call r2dr_vec (nshell, rmin, rmax, rvec)
     call r8vec_even (nshell, rmin, rmax, rvec)
     
     if ( me.eq.0) then 
        print*, 'me, allj sum_almn: ',sum(abs(almn))
        
        open(12,file=trim(adjustl(dir))//'almn_from_beta_stype'//trim(adjustl(stype))//'.unf',status='unknown',form='unformatted')
        write(12) almn
        close(12)
        !endif

        allocate(alm_r(0:nshell,1:1+2*polar,0:lmax,0:lmax))
        alm_r = 0
        
        do ir=0,nshell-1  !the zero r is set to zero by hand
           rmid = rvec(ir)
           
           !almn2ball
           arlm = 0
           do in=0,nmax
              y = in*argpi
              if (itype==1) temp = y*sinc(y*rmid)
              if (itype==2) temp = sin(y*rmid)
              if (itype==3) temp = cos(y*rmid)
              
              arlm(1,0:lmax,0:lmax) = arlm(1,0:lmax,0:lmax) + almn(1,0:lmax,0:lmax,in)*temp
              
              !if (feedback>4) print*,'ir, rmid, in, arlm, almn, Sn(r): ', ir,rmid,in, sum(abs(arlm)), sum(abs(almn(1,:,:,in))),temp
              
           enddo !end in
           
           alm_r(ir,1,0:lmax,0:lmax)=arlm(1,0:lmax,0:lmax)
           
           !if (feedback>3) print*,'calling alm2map: ir, r, sum_abs_alm: ',ir,rmid
           
           call alm2map(nside_out,lmax,lmax,arlm,map_ball(:,1,ir))
           if (feedback>3) print*,'calling alm2map done!'
           
        enddo !end ir
        
        !print*,'writing file to :',trim(adjustl(dir))//'map_ball_beta.unf'
        
        
        open(unit,file=trim(adjustl(dir))//'almr_from_beta_stype'//trim(adjustl(stype))//'.unf' ,status='unknown',form='unformatted')
        write(unit) alm_r
        close(unit)
        
        !if (me.eq.0) then
        open(12,file=trim(adjustl(dir))//'map_ball_beta_stype'//trim(adjustl(stype))//'.unf',status='unknown',form='unformatted')
        write(12) npix,nrpix,4,1
        write(12) map_ball
        close(12)
     endif
     
     deallocate(rvec)
     deALLOCATE(arlm,alm,alm_r)
     deALLOCATE(almn)
     deALLOCATE(map_ball)

  enddo

  deALLOCATE(nj_pp)
  deALLOCATE(gln)


     
  print*, 'Calling MPI finalize'
  CALL MPI_FINALIZE(ierr)


end program

