program box2ball

  USE healpix_types
  use parsers

  IMPLICIT NONE
  character(LEN=FILELEN) :: fparamin,fname_box
  INTEGER(I4B) :: ii,jj,count,ndata,balltype
  REAL(DP) :: Rbox
  REAL(DP) :: x,y,z,r,theta, phi,rmax
  REAL(DP), DIMENSION(:,:), ALLOCATABLE ::ball,ballf


  if (iargc() > 0) then
     call getarg(1,fparamin)
  else
     stop 'usage: box2ball paramfile.ini'
  endif

  call box2ballparser(fparamin,fname_box,Rbox,ndata,balltype)

  print *,'filename  to read: ', trim(fname_box)
  print*, 'ndata, box size: ',ndata, Rbox

  open(unit=12,file=trim(fname_box),status='old')


  !   rmin=1e-32
  !   rmax = -1e32

  !   xmin=1e-32
  !   xmax = -1e32

  !   ymin=1e-32
  !   ymax = -1e32

  !   zmin=1e-32
  !   zmax = -1e32

  allocate(ball(ndata,3))

  !outscribe ball with have radius equal to the distance to the corner of the box
  !while inscribe box will have radius equal to the distance to the face of the box
  rmax = sqrt(3.0*Rbox*Rbox/4.0) 
  if (balltype==inscribe) rmax=Rbox/2.0

  count=0

100 CONTINUE 
  read(12,*,end=200)x,y,z

  !! move the center from the corner of the box to the center of the box
  x=  mod(Rbox+x,Rbox)- Rbox/2.0
  y=  mod(Rbox+y,Rbox)- Rbox/2.0
  z=  mod(Rbox+z,Rbox)- Rbox/2.0

  !! r is dimensionless and is in between [0, 1]
  r = sqrt(x*x+y*y+z*z)

  if (r<=rmax) then 
     count = count+1

     theta = acos(z/r)
     phi = atan2(y,x)
     if (phi<0) phi=TWOPI + phi !now phi is in [0,2pi]

     ball(count,1) = r/rmax
     ball(count,2) = theta
     ball(count,3) = phi
  endif

  GO TO 100
200 CONTINUE
  close(12)

  allocate(ballf(count, 3))
  do ii=1,count
     ballf(ii,1:3) = ball(ii,1:3)
  enddo

  print*, 'count, rmax = ',count, rmax

  !!write r, theta, phi data
  open(13,file=trim(fname_box)//'.unf',status='unknown',form='unformatted')
  write(13) count,rmax
  write(13) ballf
  close(13)

end program box2ball
