program ball2almn

  use alm_tools
  use pix_tools
  USE healpix_types
  use parsers

  IMPLICIT NONE
  character(LEN=FILELEN) :: fparamin,fname_box,fname_ball,fmapout,ordering
  INTEGER(I4B) :: ii,jj,count,ndata,balltype
  INTEGER(I4B) :: nrpaw,nside,ntab,npix,ipix
  REAL(DP) :: Rbox,maxcount
  REAL(DP) :: x,y,z,r,theta, phi,rmax,h
  REAL(DP), DIMENSION(:,:), ALLOCATABLE ::ball,map
  REAL(DP), DIMENSION(:), ALLOCATABLE ::rtab


  if (iargc() > 0) then
     call getarg(1,fparamin)
  else
     stop 'usage: box2ball box_filename.dat'
  endif

  call ball2pixparser(fparamin,fname_box,fmapout,nrpaw,nside,ordering)
  print *, 'fname_box: ',trim(fname_box)
  print *, 'fmapout: ', trim(fmapout)
  print *, 'nrpaw,nside,ordering: ',nrpaw,nside,trim(ordering)

  !! load data
  open(13,file=trim(fname_box)//'.unf',status='old',form='unformatted')
  read(13) ndata,rmax
  print*, 'ndata, rmax = ',ndata,rmax
  allocate(ball(ndata,3))
  read(13) ball
  close(13)

  !!prepare arrays to pixelize loaded data
  ntab = 2**nrpaw
  npix = 12*nside*nside

  allocate(map(npix,ntab))
  allocate(rtab(ntab+1))
  call r8vec_even ( ntab+1, 0d0, 1d0, rtab)


  map=0.0
  maxcount = 0.0
  do ii=1,ndata
     r = ball(ii,1)
     if (r<0 .or. r>1) stop 'r is outside [0,1]'

     theta = ball(ii,2)
     phi = ball(ii,3)

     call ang2pix_ring(nside,theta,phi,ipix)
     if (ordering=='NESTED') call ang2pix_ring(nside,theta,phi,ipix)

     exist = .false.
     do jj=1,ntab
        if ((rtab(jj+1)-r)>=0 .and. (r-rtab(jj))>=0) then
           exist=.true.
           exit
        endif
     enddo
     if (.not. exist) then
        print*, 'rtab = ',rtab
        print*, 'r= ',r
        stop 'ball2pix: r could not be assigned an index in rtab'
     endif

     map(ipix,jj) = map(ipix,jj)+1
     !get maximum count
     if (map(ipix,jj)>maxcount) maxcount = map(ipix,jj)

  enddo

  print*, 'reading pixelized ball from: ',trim(fmapin)
  open(13,file=trim(fmapin),status='old',form='unformatted')
  read(13) npix,ntab,8,int(maxcount)
  read(13) map
  map=map*maxcount
  close(13)



end program ball2almn
