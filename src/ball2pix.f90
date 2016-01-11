program ball2pix

  use alm_tools
  use pix_tools
  USE healpix_types
  use parsers
  use common_module

  IMPLICIT NONE
  character(LEN=FILELEN) :: fparamin,fname_box,fname_ball,fmapout,ordering
  INTEGER(I4B) :: ii,jj,count,ndata,balltype
  INTEGER(I4B) :: nrint,nside,nshell,npix,ipix
  REAL(DP) :: Rbox,maxcount,rnmin,rnmax
  REAL(DP) :: x,y,z,r,theta, phi,rmax,h
  REAL(DP), DIMENSION(:,:), ALLOCATABLE ::ball
  REAL(SP), DIMENSION(:,:), ALLOCATABLE ::map
  REAL(DP), DIMENSION(:), ALLOCATABLE ::rtab


  if (iargc() > 0) then
     call getarg(1,fparamin)
  else
     stop 'usage: box2ball box_filename.dat'
  endif

  rnmin=0d0
  rnmax=1d0
  ordering='RING'
  nrint=32
  call ball2pixparser(fparamin,fname_box,fmapout,nrint,nside,ordering,rnmin,rnmax)
  print *, 'fname_box: ',trim(fname_box)
  print *, 'fmapout: ', trim(fmapout)


  !! load data
  open(13,file=trim(fname_box)//'.unf',status='old',form='unformatted')
  read(13) ndata,rmax
  print*, 'ndata, rmax = ',ndata,rmax
  allocate(ball(ndata,3))
  read(13) ball
  close(13)

  !!prepare arrays to pixelize loaded data
  nshell = 2*nrint
  npix = 12*nside*nside

  print *, 'nrint,nside,ordering, nshell: ',nrint,nside,trim(ordering),nshell
  print*, 'rnmin, rnmax = ',rnmin,rnmax

  allocate(map(0:npix-1,0:nshell))
  allocate(rtab(0:nshell))
  !call r2dr_vec ( nshell, rnmin, rnmax, rtab)
  call r8vec_even ( nshell, rnmin, rnmax, rtab)


  map=0
  maxcount = 0.0
  do ii=1,ndata
     r = ball(ii,1)
     if (r<0 .or. r>1) stop 'r is outside [0,1]'

     theta = ball(ii,2)
     phi = ball(ii,3)


     if (ordering=='NESTED') then 
        call ang2pix_nest(nside,theta,phi,ipix)
     else
        call ang2pix_ring(nside,theta,phi,ipix)
     endif

     exist = .false.
     do jj=0,nshell
        if (jj==0) then 
           if ((rtab(jj)-r)>=0) then
              print*, 'zero r case: gridr, datar:',rtab(jj),r
              exist=.true.
              exit
           endif
        else
           if ((r-rtab(jj-1))>=0 .and. (rtab(jj)-r)>=0) then
              exist=.true.
              exit
           endif
        endif
     enddo
     if (.not. exist) then
        print*, 'rtab = ',rtab
        print*, 'r= ',r
        stop 'ball2pix: r could not be assigned an index in rtab'
     endif

     if (jj==0) print* ,'jj=0 found: rdata=',r

     !if (ipix > 3070 .and. jj<2) print*,'ipix, jj',ipix, jj
     
     map(ipix,jj) = map(ipix,jj)+1
     !get maximum count
     if (map(ipix,jj)>maxcount) maxcount = map(ipix,jj)

  enddo

  !do jj=0,nshell
  !   print*, 'map count per shell: ishell, count',jj, sum(map(:,jj))
  !enddo

  print*, 'pixelized ball will be writting in: ',trim(fmapout)
  print*,'nside, nrint: ',nside,nrint
  open(13,file=trim(fmapout),status='unknown',form='unformatted')
  write(13) npix,nshell+1,4,int(maxcount)
  map=map/real(maxcount)
  write(13) map
  close(13)



end program ball2pix
