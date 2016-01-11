program test_filon_quad

  !yabebal fantaye 12 Jan 2013

  use healpix_types


  IMPLICIT NONE

  !  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large
  REAL(DP), DIMENSION(:), ALLOCATABLE :: f,ff,fff,rf,fn

  CHARACTER(LEN=128) :: healpixdir,dir,str
  CHARACTER(LEN=200) :: fname


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp

  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix
  real(dp) ::rmin,rmax,dr,rmid
  real(dp):: s2n, s2nm1

  !   CALL MPI_INIT(ierr)

  !   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  !   CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  me=0
  ntasks=1
  unit=12+me


  call getEnvironment("HEALPIX",healpixdir)
  dir=''!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  if (iargc() > 0) then
     call getarg(1,fname)
  else
     stop 'usage test_filon_quad paramffile.txt'
  endif
  open(12,file=trim(adjustl(fname)),status='old')
  read(12,*) large
  read(12,*) nside
  read(12,*) polar
  read(12,*) lmax
  read(12,*) nnmax
  read(12,*) nshell
  read(12,*) rmin,rmax
  read(12,*) nrsample
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


     if ( rmin == rmax ) then
        stop 'idl_ball2almn : rmin=rmax'
     end if

     if ( nrsample < 0 .or. mod ( nrsample, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'WARNING NRSAMPLE (per shell) must be even and >= 0. Assuming max(0,nrsample+1)'
        nrsample = max(0,nrsample+1)
        write ( *, '(a,i8)' ) '  new nrsample  = ',nrsample
        stop
     end if
     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 1 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be odd and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif


  nrpix = nshell*(nrsample+1)

  if  (me == 0) then
     print*, 'nnmax, nrpix = ',nnmax, nrpix
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell
     print*, 'calling map2alm in MPI distributed loop'
  endif


  ALLOCATE(f(0:nshell-1), rf(0:nshell), fff(0:nshell-1), fn(nnmax))
  call r8vec_even (nshell, dble(0.0001), dble(0.999), f )
  f = f**(polar+1)



  !get radial values (min and max of each shell) at nrpix points
  allocate(rvec(0:nrpix-1),ff(0:nrpix-1))
  call r8vec_even (nrpix, rmin, rmax, rvec )

  do i=0,nshell-1
     oddInd = (nrsample+1)*i
     evenInd = ((nrsample+1)*(i+1)-1)
     do jj=oddInd, evenInd
        ff(jj) = f(i)
     enddo
  enddo



  print*, 'starting filon_sin integration ff --> fn'


  dr = ( rmax - rmin ) / real( nrpix-1, dp )
  allocate(alpha(nnmax),beta(nnmax),gamma(nnmax))
  call filon_params(nnmax,dr,alpha,beta,gamma)


  do jj = 1,nnmax

     do ii=0,(nrpix-1)/2

        oddInd = 2*ii !zero based and knowning nrpix is odd
        evenInd = 2*(ii+1)-1 !zero based 


        if (oddInd <=nrpix-1) fn(jj) =  fn(jj) + dr * ( beta(jj)*ff(oddInd)*sin(PI*jj*rvec(oddInd)) ) !odd sum
        if (evenInd <=nrpix-1) fn(jj) =  fn(jj) + dr * ( gamma(jj)*ff(evenInd)*sin(PI*jj*rvec(evenInd)) )  !even sum

        if (jj==1 .and. me==0 .and. oddInd <=nrpix-1) print*, 'odd ind, r, sin(pi*n*r)', oddInd, real(rvec(oddInd)), real(sin(PI*jj*rvec(oddInd))),  ff(oddInd)
        if (jj==1 .and. me==0 .and. evenInd <=nrpix-1) print*, 'even ind, r, sin(pi*n*r)', evenInd,real(rvec(evenInd)), real(sin(PI*jj*rvec(evenInd))),  ff(evenInd)

     enddo

     fn(jj) =  fn(jj) -  dr * ( &
          beta(jj)/2d0*( ff(0)*sin(PI*jj*rvec(0)) + ff(nrpix-1)*sin(PI*jj*rvec(nrpix-1)) ) &
          + alpha(jj)*( ff(0)*cos(PI*jj*rvec(0)) - ff(nrpix-1)*cos(PI*jj*rvec(nrpix-1)) ) )

  enddo

  print*, 'starting inverse filon_sin  fn --> ff'

  call r8vec_even (nshell+1, rmin, rmax, rf )
  do ii=1,nshell
     rmid = (rf(ii+1) + rf(ii))/2d0
     do jj=1,nnmax
        fff(ii) = fff(ii) + fn(jj)*sin(PI*jj*rmid)/rmid
     enddo
  enddo

  print*, '----------------------------'
  print*, '----------------------------'
  print*, 'orginal func (f): ',f
  print*, '----------------------------'
  print*, '----------------------------'
  print*, 'inv sync expansion func (ff): ',fff
  print*, '----------------------------'
  print*, '----------------------------'
  print*, 'abs(f - ff): ',abs(f - ff)
  print*, '----------------------------'
  print*, '----------------------------'
  
  
  if (me==0) then
     if (large.eq.1)  fname = trim(adjustl(dir))//'large/ff.unf'
     if (large.eq.0)  fname = trim(adjustl(dir))//'f_finv.unf' 
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) f
     write(unit) ff
     close(unit)

  endif

!  CALL MPI_FINALIZE(ierr)

end program test_filon_quad



subroutine filon_params(nnmax,h,alpha,beta,gamma)

  implicit none

  integer ( i4b ), intent(in) :: nnmax
  real(dp), intent(in)::  h
  real(dp), dimension(nnmax), intent(out) :: alpha,beta, gamma
  real(dp), dimension(nnmax) :: cost,sint,theta
  integer ::ii


do ii=1,nnmax
  theta(ii) = PI*ii*h
enddo

  sint = sin ( theta )
  cost = cos ( theta )

  where ( 6.0D+00 * abs ( theta ) <= 1.0D+00 ) 

     alpha = 2.0D+00 * theta**3 /   45.0D+00 &
          - 2.0D+00 * theta**5 /  315.0D+00 &
          + 2.0D+00 * theta**7 / 4725.0D+00

     beta =  2.0D+00            /     3.0D+00 &
          + 2.0D+00 * theta**2 /    15.0D+00 &
          - 4.0D+00 * theta**4 /   105.0D+00 &
          + 2.0D+00 * theta**6 /   567.0D+00 &
          - 4.0D+00 * theta**8 / 22275.0D+00

     gamma = 4.0D+00            /      3.0D+00 &
          - 2.0D+00 * theta**2 /     15.0D+00 &
          +           theta**4 /    210.0D+00 &
          -           theta**6 /  11340.0D+00

  elsewhere

     alpha = ( theta**2 + theta * sint * cost &
          - 2.0D+00 * sint**2 ) / theta**3

     beta = ( 2.0D+00 * theta + 2.0D+00 * theta * cost**2 &
          - 4.0D+00 * sint * cost ) / theta**3

     gamma = 4.0D+00 * ( sint - theta * cost ) / theta**3

     
  endwhere



! s2n = sum ( ftab(1:ntab:2) * sin ( t * xtab(1:ntab:2) ) ) &
!      - 0.5D+00 * ( ftab(ntab) * sin ( t * xtab(ntab) ) &
!      + ftab(1) * sin ( t * xtab(1) ) )

! s2nm1 = sum ( ftab(2:ntab-1:2) * sin ( t * xtab(2:ntab-1:2) ) )

! result = h * ( &
!      alpha * ( ftab(1) * cos ( t * xtab(1) ) &
!      - ftab(ntab) * cos ( t * xtab(ntab) ) ) &
!      + beta * s2n &
!      + gamma * s2nm1 )


end subroutine filon_params

!------------------------------

subroutine r8vec_even( ntab, alo, ahi, a )

  implicit none

  integer(i4b) :: i,ntab
  real(dp):: a(ntab), ahi, alo

  if ( ntab == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, ntab
      a(1) = ( real(ntab-i, dp) * alo + real(i-1, dp) * ahi )/real(ntab-1, dp)
    end do

  end if

  return

end subroutine r8vec_even



!central r value for the ith radial pixel
!i = ii-1
!rmid = (rvec(i+1) + rvec(i))/2d0
!almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + alm_TGC_r(i,:,:,:)*sin(PI*jj*rmid)*rmid*dr
