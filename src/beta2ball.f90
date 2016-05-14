  
  
!!!IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use fitstools

  use parsers
  use needlet_kernel
  use mapio
  use almnutil
  use fftw_wrapper
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp


  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,nside_jmap, npix_jmap
  INTEGER(I4B) :: lmax,i,j,k,l,m,polar,balltype
  INTEGER(I4B) :: ifirstshell,ilastshell,iwidth
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd
  INTEGER(I4B) :: nv,in,iv,ir
  integer(i4b) :: feedback=1

  !needlet
  CHARACTER(LEN=30) :: wavtyp
  INTEGER(I4B) :: p,j0,nj,glnpow,wav, first_j,last_j,first_v, last_v
  REAL(DP) :: fwhm,sigma,norm, bb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl, cl, cln
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln


  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: map_slet
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map_ball,map
  REAL(DP), DIMENSION(:), ALLOCATABLE :: map1,map2,map3

  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almv,anlm,anlm_all
  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: anlmj
  REAL(DP), DIMENSION(1:2) :: zbounds

  REAL(DP),DIMENSION(:,:),allocatable :: w8ring_map

  !strings
  CHARACTER(LEN=*),parameter :: code='BETA2BALL_PROGRAM'
  CHARACTER(LEN=30) :: filestr,str,strn,filenum
  CHARACTER(LEN=FILELEN) :: fname,parfile,inputfile,outputfile,healpixdir,tempdir,finext

  !MPI params
  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src, s0, sn,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp,nj_pp

  !timing
  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1


  if (iargc() > 0) then
     call getarg(1,parfile)
  else
     stop 'usage: '//code//' params_filename.dat'
  endif


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)


  call wall_clock_time(time0)
  call cpu_time(ptime0)

  call getEnvironment("HEALPIX",healpixdir)

  polar=0
  wav=0
  p=0


  !read parameter file
  call almnparser(parfile,nside, lmax, nnmax, nshell,ifirstshell,ilastshell,iwidth,norm,balltype,inputfile,outputfile,finext,tempdir)
  call betaparser(parfile, bb, j0, nj, glnpow, wav, nside_jmap)

  if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
     if (me == 0) write ( *, '(a)' ) 'Fatal error! NSHELL must be eve and > 1.'
     if (me == 0) write ( *, '(a,i8)' ) '  NTAB = ',nshell
     stop
  end if

  if ( 2*nnmax>nshell ) then
     if (me == 0) print *,  'WARNNING! setting nnmax=NSHELL/2. Given nnmax=',nnmax
     nnmax=nshell/2
  end if

  npix=12*nside*nside
  npix_jmap=12*nside_jmap*nside_jmap

  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=j0
  if (me.ne.0) first_j=j0+sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1

  if (feedback>2) print*,'proc:',me,nj_pp(me),first_j,last_j


  !if jmaps are squared, glnpow>1, do simple sum. Otherwise, do the following
  if (glnpow>0) then


     ALLOCATE( w8ring_map(2*nside,1+2*polar) )
     w8ring_map=1d0
     call read_ringweights(nside, 1 ,healpixdir, w8ring_map )
     zbounds = (/ -1.0_dp, 1.0_dp /)


     allocate(map_slet(0:npix_jmap-1,0:nshell-1,0:nj_pp(me)-1))
     ALLOCATE(map(0:npix-1,1))
     ALLOCATE(alm(1:1,0:lmax,0:lmax))

     !read map based on balltype format
     if (balltype == slice_balltype) then
        do j=0,nj_pp(me)-1
           do k=0,nshell-1
              WRITE(filenum,*) first_j-j0+j !+1000
              WRITE(strn,*) k
              fname=trim(adjustl(inputfile))//'_j'//trim(adjustl(filenum))//'_r'//trim(adjustl(strn))//'.unf'
              if (glnpow>1) fname=trim(adjustl(inputfile))//'_gln2_j'//trim(adjustl(filenum))//'_r'//trim(adjustl(strn))//'.unf'

              open(12,file=fname,status='old',form='unformatted', ACCESS="STREAM")
              read(12) map_slet(0:npix_jmap-1,k,j)
              close(12)
           enddo
        enddo
     else !full ball map
        do j=0,nj_pp(me)-1
           WRITE(filenum,*) first_j-j0+j !+1000
           fname=trim(adjustl(inputfile))//'_j'//trim(adjustl(filenum))//'.unf'
           if (glnpow>1) fname=trim(adjustl(inputfile))//'_gln2_j'//trim(adjustl(filenum))//'.unf'

           open(12,file=fname,status='old',form='unformatted', ACCESS="STREAM")
           read(12) map_slet(0:npix-1,0:nshell-1,j)
           close(12)
        enddo
     endif



     if (feedback>2) print*,'beta map read!'
     if (wav.eq.0) wavtyp='standard'
     if (wav.eq.1) wavtyp='mexican'

     ALLOCATE(gln(0:lmax,0:nshell-1,j0:j0+nj-1))
     CALL calc_f2
     !print*,'beta2map1:',j0,nj,lmax,bb,wavtyp,p
     CALL calc_gln(j0,nj,lmax,nnmax,gln,bb,wavtyp,p,gln4fft=1)
     if (feedback>2) print*,'gln computed!'


     ALLOCATE(almv(0:nshell-1,0:lmax,0:lmax))
     ALLOCATE(anlm(0:nshell-1,0:lmax,0:lmax))

     anlm=CMPLX(0.0_dp, 0.0_dp, KIND=DP)

     do j=0,nj_pp(me)-1 

        map=0.
        do iv=0,nshell-1
           map(0:npix-1,1) = map_slet(0:npix-1,iv, j)
           if (count( map(:,1)>1e-13 )>10) call map2alm_iterative(nside, lmax, lmax,3, map, alm, zbounds, w8ring_map) !sum over k
           almv(iv,0:lmax,0:lmax) = alm(1,0:lmax,0:lmax)
        enddo

        !do parallel fft forward
        call fftw_filter_complex_lm(almv, gln(0:lmax,0:nshell-1,first_j+j),.false.)


        anlm = anlm + almv !sum over local j

     enddo !end local j

     deALLOCATE(w8ring_map)
     deALLOCATE(map_slet)
     deALLOCATE(gln)
     deALLOCATE(almv)


     !sum over globalj    - used only at MASTER 0
     allocate(anlm_all(0:nshell-1,0:lmax,0:lmax))


     if (ntasks >1 ) then
        cnt=(nshell)*(lmax+1)**2
        if (me<2) print*,'calling MPI_ALLREDUCE, cnt=',cnt

        !take the sum per node and sum over all nodes of alm
        call MPI_REDUCE(anlm,anlm_all,cnt, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD,ierr)  
        !CALL MPI_ALLREDUCE(anlm,anlm_all,cnt,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        print*,'No MPI_ALLREDUCE'
        anlm_all=anlm
     endif

     !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     if (me<2) print*,'calling MPI_ALLREDUCE done, iproc=',me

     deallocate(anlm)


     !=== for now compute sum over lmn in main processor ===
     if ( me.eq.0) then 

        !-- compute and save C_ln from a_nlm --
        allocate(cln(0:nshell-1,0:lmax))
        fname = trim(adjustl(tempdir))//'cln.unf'
        call almn2cln(anlm_all,cln,fname)     
        deallocate(cln)


        !=== for now compute alm_r in all processor ===
        !do parallel fft backward
        call fftw_complex_lm(anlm_all, .true.) 
        
        !-- compute and save mean C_l frm a_lm(r)
        allocate(cl(0:lmax,2))
        fname = trim(adjustl(tempdir))//'meancl.unf'
        call almr2cl(anlm_all,cl,fname)
        deallocate(cl)


        !ALLOCATE(map_ball(0:npix-1,0:nshell-1))

        print*, 'doing the final alm2map and saving map ..'
   

        do ir=0,nshell-1
           alm(1,0:lmax,0:lmax) =  anlm_all(ir,0:lmax,0:lmax)
           call alm2map(nside,lmax,lmax,alm,map(:,1)) !sum over lm

           write(strn,*)ir
           fname=trim(adjustl(outputfile))//'_'//trim(adjustl(strn))

           open(12,file=fname,status='unknown',form='unformatted', ACCESS="STREAM")
           write(12) map 
           close(12)
           if (mod(ir, 10)==0) print*, 'saved map at ishell=',ir
        enddo

     endif

     if (me<10) print*,'calling alm2map done!'

     deALLOCATE(anlm_all)

  else !glnpow<=0


     if (me==0) then
        allocate(map1(0:npix_jmap-1))
        ALLOCATE(map2(0:npix_jmap-1))

        do ir=0,nshell-1
           write(strn,*)ir
           map2=0d0
           do j=0,nj-1
              WRITE(filenum,*) j
              fname=trim(adjustl(inputfile))//'_gln2_j'//trim(adjustl(filenum))//'_r'//trim(adjustl(strn))//'.unf'
              open(12,file=fname,status='old',form='unformatted', ACCESS="STREAM")
              read(12) map1 
              close(12)
              map2=map2 + map1 !sum over j
           enddo

           fname=trim(adjustl(outputfile))//'_'//trim(adjustl(strn))
           open(12,file=fname,status='unknown',form='unformatted', ACCESS="STREAM")
           write(12) map2 
           close(12)
        enddo
     endif

  endif !glnpow condition




  deALLOCATE(nj_pp)



  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) " Clock and CPU time [s]: ", clock_time, ptime


  print*, 'Calling MPI finalize'
  CALL MPI_FINALIZE(ierr)


end program

