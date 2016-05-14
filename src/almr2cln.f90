program idl_ball2almn

  !yabebal fantaye 12 Jan 2013

  use parsers
  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod

  use fitstools
  use head_fits
  use mapio

  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  USE paramfile_io, only: paramfile_handle, parse_init, parse_int, &
       parse_string, parse_double, parse_lgt, parse_summarize, parse_check_unused, &
       parse_finish, concatnl, scan_directories, get_healpix_data_dir, get_healpix_test_dir, &
       get_healpix_pixel_window_file

  IMPLICIT NONE


  INCLUDE 'mpif.h'

  type(paramfile_handle) :: handle

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,balltype
  INTEGER(I4B) :: ifirstshell,ilastshell,iwidth
  INTEGER(I4B) :: ifstr,ffstr
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map_TQU

  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(DP) :: fwhm,sigma,norm
  CHARACTER(LEN=30) :: filestr,str
  CHARACTER(LEN=FILELEN) :: fname,parfile,inputfile,outputfile,healpixdir,tempdir,finext


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src,sn,s0,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp

  COMPLEX(DPC) :: fr_mean
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cln
  REAL(DP), DIMENSION(:), ALLOCATABLE :: cl_proj
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cl
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix
  real(dp) ::rmin,rmax,dr,rmid,y
  real(dp):: s2n, s2nm1

  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1
  !type(planck_fft2_plan) :: plan_ball2beta
  !call make_fft2_plan (plan_ball2beta,length,fft2_forward)

  if (iargc() > 0) then
     call getarg(1,parfile)
  else
     !stop 'usage: box2ball box_filename.dat'
     parfile=''
  endif

  handle = parse_init(parfile)

  inputfile = parse_string(handle, 'inputfile', default='almr.unf', descr='input almr filename. Default almr.unf', filestatus='old')
  outputfile = parse_string(handle, 'outputfile', default="ps", descr='basename Power spectra. Default is ./ps', filestatus='new')
  lmax = parse_int(handle, 'lmax', default=2000, descr='angular max multipole. Default is lmax=2000.')
  nshell = parse_int(handle, 'nshell', default=256, descr='number of radial shells. Default is nshell=256.')
  nnmax=nshell/2

  call parse_check_unused(handle, code='ALMr2CL')
  call parse_summarize(handle,code='ALMr2CL')
  call parse_finish(handle)

  unit=12+me



  call wall_clock_time(time0)
  call cpu_time(ptime0)


  call getEnvironment("HEALPIX",healpixdir)


  !check correctness of input parameters
  if (me == 0) then
     print*, '******* Input parameter ******'
     print*, '---------------------------'
     print*, 'nside = ',nside 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, '---------------------------'

     if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'Fatal error! NSHELL must be even and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if

     if ( nnmax==0 .or. nshell <= 2*nnmax ) then
        print *,  'WARNNING! setting nnmax=NSHELL/2. Given nnmax=',nnmax
        nnmax=nshell/2
     end if

  endif


  !allocate big arrays
  allocate(alm_TGC_r(0:nshell-1,0:lmax,0:lmax))
  alm_TGC_r = 0

  !allocate signle r alm array 
  allocate(alm_TGC(1:1,0:lmax,0:lmax))
  alm_TGC = 0

  !read almr 
  fname = trim(adjustl(inputfile))
  print*, '******* reading almn_r array from **'
  print*, trim(fname)
  print*, '---------------------------'
  open(unit,file=trim(fname),status='old',form='unformatted')
  read(unit) alm_TGC_r
  close(unit)


  print*, 'starting fftw integration'
  allocate(fr(-3:nshell-1+3)) !pad -3:0, nshell:nshell+2
  allocate(cln(0:nshell-1,0:lmax))
  allocate(cl_proj(0:lmax),cl(0:lmax,1:1))
  cln=0d0
  fr=0d0


  !compute projected power spectrum
  cl_proj=0
  do ii=0,nshell-1 
     alm_TGC(1,0:lmax,0:lmax) =alm_TGC(1,0:lmax,0:lmax) + alm_TGC_r(ii,0:lmax,0:lmax)
     call alm2cl(lmax, lmax, alm_TGC_r(ii:ii,0:lmax,0:lmax),cl)
     cl_proj=cl_proj + cl(:,1) !sum over Cl(r)
  enddo
  call alm2cl(lmax, lmax, alm_TGC,cl) !cl from mean alm 


  do l=0,lmax
     do m=0,l
        fr=0d0
 
        fr(0:nshell-1) = alm_TGC_r(0:nshell-1,l,m)


        !Healpix forward fft
        call complex_fft(fr,backward=.false.) 

        !      if (l<10 .and. l>5) then
        !         print*, 'l,m',l,m
        !         print*, 'f_n',fr(nnmax/2:nnmax-1)
        !         !print*, '-----------'
        !         !print*, 'conj(f_n)',fr(nshell-1:nshell-nnmax+1:-1)
        !      endif

        !DC mode n=0
        !Nyquist mode n=n_nyquist=nshell/2
        !positive frequency 0:nshell/2 -  n=0,1,2..nshell/2
        !negative frequency nshell/2+1:nshell-1 - n=[1,2,3...nshell/2-1]-nshell/2,
        cln(0:nshell-1,l) = cln(0:nshell-1,l) + real(fr(0:nshell-1)*conjg(fr(0:nshell-1)),kind=dp )/(l+1d0)
     enddo
  enddo
  deallocate(fr)



  fname = trim(adjustl(outputfile))//'.cln'
  print*, '******* saving 2d power spectrum cln **'
  print*, 'sum_n(cln(:,lmax/2:(lmax/2+10))',sum(cln(:,lmax/2:(lmax/2+10)),1)
  print*, trim(fname)
  print*, 'size(almn), shape(almn) = ',size(cln), shape(cln)
  print*, '---------------------------'
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) cln
  close(unit)

  fname = trim(adjustl(outputfile))//'.cl_mean_clr'
  print*, '******* saving sum cl(r) or sum over n of cln **'
  print*, 'cl_proj(lmax/2:(lmax/2+10))',cl_proj(lmax/2:(lmax/2+10))
  print*, trim(fname)
  !print*, 'size(almn), shape(almn) = ',size(cln), shape(cln)
  print*, '---------------------------'
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) cl_proj
  close(unit)

  fname = trim(adjustl(outputfile))//'.cl_mean_almr'
  print*, '******* saving cl from mean alm **'
  print*, 'cl_mean_alm(lmax/2:(lmax/2+10))',cl(lmax/2:(lmax/2+10),1)
  print*, trim(fname)
  !print*, 'size(almn), shape(almn) = ',size(cln), shape(cln)
  print*, '---------------------------'
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) cl(:,1)
  close(unit)




  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  if (me<=10) write(*,*) " Clock and CPU time [s]: ", clock_time, ptime


end program idl_ball2almn



!central r value for the ith radial pixel
!i = ii-1
!rmid = (rvec(i+1) + rvec(i))/2d0
!almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + alm_TGC_r(i,:,:,:)*sin(PI*jj*rmid)*rmid*dr
