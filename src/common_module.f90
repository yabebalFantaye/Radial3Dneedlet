module common_module


  use healpix_types

  IMPLICIT NONE



contains


  !=================================================

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
    real(dp):: a(0:ntab), ahi, alo,h


    h = (ahi-alo)/real(ntab, dp)

    do i = 0, ntab
       a(i) = real(i, dp)*h !( real(ntab-i, dp) * alo + real(i-1, dp) * ahi )/real(ntab-1, dp)
    end do

    ! if ( ntab == 1 ) then

    !    a(0) = alo
    !    a(1) = ahi !0.5D+00 * ( alo + ahi )

    ! else

    !    do i = 0, ntab
    !       a(i) = ( real(ntab-i, dp) * alo + real(i-1, dp) * ahi )/real(ntab-1, dp)
    !    end do

    ! end if

    return

  end subroutine r8vec_even


  !============================

  subroutine r2dr_vec( ntab, alo, ahi, rvec,iedge )

    implicit none
    integer(i4b), optional::iedge !if passed the return is the edge of the radial bins 
    integer(i4b) :: i,ntab,nhalf,k,kk
    real(dp):: rvec(0:ntab), ahi, alo,dr,rr,dr0,dv,ak,bk,h


    nhalf = ntab/2
    if(2*nhalf .ne. ntab) stop 'filon_simp_r2drpoints: nn must be an even number'

    dv=(ahi**3-alo**3)/dble(ntab) !equal volume division

    !print*, 'dv = ',dv
    dr=dv !first bin dr would be just dv

    !special condition for the first bin

    do k = 0,ntab
       dr=0d0
       !if (.not.present(iedge) .and. k.ne.0) dr=(dv)**(1./3.)/(dble(k))**(2./3.)
       rvec(k)=(dble(k)*dv)**(1./3.) - dr/2d0  
       !print*,'r2dr: dr=',k, dr,(dble(k)*dv)**(1./3.),rvec(k)
    enddo

    !print*, 'r2dr: rvec=',rvec

    return

  end subroutine r2dr_vec


  !================
  subroutine r2dr_vec_ir( ntab, alo, ahi, ir, rval)

    implicit none
    integer(i4b) :: i,ntab,nhalf,k,kk,ir
    real(dp):: rval, ahi, alo,dr,rr,dr0,dv,ak,bk,h


    nhalf = ntab/2
    if(2*nhalf .ne. ntab) stop 'filon_simp_r2drpoints: nn must be an even number'

    dv=(ahi**3-alo**3)/dble(ntab) !equal volume division

    dr=dv !first bin dr would be just dv

    !special condition for the first bin

    dr=0d0
    rval=(dble(ir)*dv)**(1./3.) - dr/2d0  


    return

  end subroutine r2dr_vec_ir


  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !     
  subroutine fil_simp_r2drpoints(a,b,y,nn,j,rvec,wvec)
    !     
    !     calculates nn (even) + 1 Filon-Simpson points "xfisi" 
    !     and weights "wfisi" for numerical integration of
    !     
    !     I_1[f] = int_a^b dx f(x) sin(x*y)/(x*y)              j = 1
    !     
    !     I_2[f] = int_a^b dx f(x) 4 sin^2(x*y/2)/(x^2*y^2)    j = 2
    !     
    !     at fixed y
    !     
    !     input:    a -- lower limit of integral
    !     b -- upper limit of integral
    !     y -- frequency parameter 
    !     nn -- number of intervals (even, < 1000)
    !     j  -- type of integral (see above)
    !     
    !     output:   xfisi(n) -- equidistant integration points (n = 0,1, ... nn)
    !     wfisi(n) -- weights for Filon-Simpson integration 
    !     (n = 0,1, ... nn)
    !     
    !     
    !     I_j[f] approx \sum_{n=0}^{nn} wfisi(n)*f(xfisi(n))
    !     
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     
    !     author: R. Rosenfelder
    !     Paul-Scherrer-Institut (PSI)
    !     CH-5232 Villigen PSI
    !     Switzerland
    !     
    !     reference: R. Rosenfelder, arXiv: hep-ph/0603161
    !     
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     
    !     21 March, 2006   
    !     
    !     
    integer, parameter :: nmax=30000
    integer :: nn,nhalf,k,kk,j
    real*8:: a,b,h,ak,bk,y,x0,x1,x2
    real*8, dimension(0:nn) ::xfisi,wfisi,weight
    real*8, dimension(0:nn) ::rvec,wvec


    ! if(nn .gt. nmax) then
    !    write(6,*) 'nn > nmax'
    !    stop
    ! endif
    !     
    nhalf = nn/2
    if(2*nhalf .ne. nn) stop 'filon_simp_r2drpoints: nn must be an even number'

    !     
    call r2dr_vec( nn, a, b, rvec ) !initialize radial nodes

    wfisi(0) = 0.d0
    do k = 1,nhalf
       kk = 2*(k - 1)

       x0 = rvec(kk)
       x1 = rvec(kk+1)
       x2 = rvec(kk+2)
       !     
       call fil_simp_r2dr(x0,x1,x2,y,j,weight)
       !     
       wfisi(kk) =  wfisi(kk) + weight(0)
       !     
       wfisi(kk+1) = weight(1)
       !     
       wfisi(kk+2) = weight(2)
    enddo

    wfisi(nn) = 2*wfisi(nn) !make the end to similar order as previous nodes

    wvec = wfisi(0:nn)

    return

  end subroutine fil_simp_r2drpoints

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  subroutine fil_simp_r2dr(x0,x1,x2,y,j,wfisi)
    !  
    !  calculates weights wfisi(i), i = 0, 1, 2 so that
    !
    !  int_a^{b=a+2*h} dx f(x) O_j(x*y)
    !
    !               = wfisi(0)*f(a) + wfisi(1)*f(a+h) + wfisi(2)*f(b)
    !                                        
    !  is exact for monomials x^0, x^1, x^2
    !
    !  oscillating weight function:   O_1(x*y) = sin(x*y)/(x*y)
    !                                 O_2(x*y) = 4*sin^2(x*y/2)/(x*y)^2
    !
    !  for y = 0 : Simpson weights 
    !
    !  input:    a -- lower limit of integral
    !            b -- upper limit of integral
    !            y -- frequency parameter 
    !            j  -- type of integral (see above)
    !
    !  output:   wfisi(i) -- Filon-Simpson weights (i = 0,1,2)
    !
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !  author: R. Rosenfelder
    !          Paul-Scherrer-Institut (PSI)
    !          CH-5232 Villigen PSI
    !          Switzerland
    !
    !  reference: R. Rosenfelder, arXiv: hep-ph/0603161
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !  21 March 2006
    !
    !
    integer::  j
    real*8 :: a,b,wfisi(0:2),y
    real*8 :: h, h2,fak,aint0,aint1,aint2,x0,x1,x2

    !
    !    moments
    !
    aint0 = aux(j,0,x2,y) - aux(j,0,x0,y)
    aint1 = aux(j,1,x2,y) - aux(j,1,x0,y)
    aint2 = aux(j,2,x2,y) - aux(j,2,x0,y)
    !
    !     weights
    !

    wfisi(0) = (x1*x2*aint0 - (x1+x2)*aint1 + aint2)/((x2-x0)*(x1-x0))
    wfisi(1) = (-x0*x2*aint0 + (x2+x0)*aint1 - aint2)/((x2-x1)*(x1-x0))
    wfisi(2) = (x0*x1*aint0 - (x0+x1)*aint1 + aint2)/((x2-x1)*(x2-x0))
    !
    return
  end subroutine fil_simp_r2dr
  !

  !=====================
  subroutine fil_simp_points(a,b,y,nn,j,rvec,wvec,itype)
    !     
    !     calculates nn (even) + 1 Filon-Simpson points "xfisi" 
    !     and weights "wfisi" for numerical integration of
    !     
    !     I_1[f] = int_a^b dx f(x) sin(x*y)/(x*y)              j = 1
    !     
    !     I_2[f] = int_a^b dx f(x) 4 sin^2(x*y/2)/(x^2*y^2)    j = 2
    !     
    !     at fixed y
    !     
    !     input:    a -- lower limit of integral
    !     b -- upper limit of integral
    !     y -- frequency parameter 
    !     nn -- number of intervals (even, < 1000)
    !     j  -- type of integral (see above)
    !     
    !     output:   xfisi(n) -- equidistant integration points (n = 0,1, ... nn)
    !     wfisi(n) -- weights for Filon-Simpson integration 
    !     (n = 0,1, ... nn)
    !     
    !     
    !     I_j[f] approx \sum_{n=0}^{nn} wfisi(n)*f(xfisi(n))
    !     
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     
    !     author: R. Rosenfelder
    !     Paul-Scherrer-Institut (PSI)
    !     CH-5232 Villigen PSI
    !     Switzerland
    !     
    !     reference: R. Rosenfelder, arXiv: hep-ph/0603161
    !     
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     
    !     21 March, 2006   
    !     
    !     
    integer, parameter :: nmax=30000
    integer :: nn,nhalf,k,kk,j,iflag
    real*8:: a,b,h,ak,bk,y
    real*8, dimension(0:nmax) ::xfisi,wfisi,weight
    real*8, dimension(0:nn) ::rvec,wvec
    integer,optional ::itype


    !itype=1 (sinx/x), 2 (sinx) 3(cosx)
    iflag=2
    if (present(itype)) iflag=itype

    if(nn .gt. nmax) then
       write(6,*) 'nn > nmax'
       stop
    endif
    !     
    nhalf = nn/2
    if(2*nhalf .ne. nn) then
       write(6,*) 'nn odd'
       stop
    endif

    !call r8vec_even( nn, a, b, rvec ) !initialize radial nodes
    !     
    h = (b-a)/nn
    wfisi(0) = 0.d0
    do k = 1,nhalf
       kk = 2*(k - 1)
       ak = a + kk*h
       bk = ak + h + h
       !     
       if (iflag==1) call fil_simp(ak,bk,y,j,weight)
       if (iflag==2) call fil_simp_sin(ak,bk,y,j,weight)
       if (iflag==3) call fil_simp_cos(ak,bk,y,j,weight)
       !     
       if(kk .eq. 0) xfisi(kk) = ak
       wfisi(kk) =  wfisi(kk) + weight(0)
       !     
       xfisi(kk+1) = ak + h
       wfisi(kk+1) = weight(1)
       !     
       xfisi(kk+2) = bk
       wfisi(kk+2) = weight(2)
    enddo

    rvec = xfisi(0:nn)
    wvec = wfisi(0:nn)

    return

  end subroutine fil_simp_points
  !
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !
  subroutine fil_simp(a,b,y,j,wfisi)
    !  
    !  calculates weights wfisi(i), i = 0, 1, 2 so that
    !
    !  int_a^{b=a+2*h} dx f(x) O_j(x*y)
    !
    !               = wfisi(0)*f(a) + wfisi(1)*f(a+h) + wfisi(2)*f(b)
    !                                        
    !  is exact for monomials x^0, x^1, x^2
    !
    !  oscillating weight function:   O_1(x*y) = sin(x*y)/(x*y)
    !                                 O_2(x*y) = 4*sin^2(x*y/2)/(x*y)^2
    !
    !  for y = 0 : Simpson weights 
    !
    !  input:    a -- lower limit of integral
    !            b -- upper limit of integral
    !            y -- frequency parameter 
    !            j  -- type of integral (see above)
    !
    !  output:   wfisi(i) -- Filon-Simpson weights (i = 0,1,2)
    !
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !  author: R. Rosenfelder
    !          Paul-Scherrer-Institut (PSI)
    !          CH-5232 Villigen PSI
    !          Switzerland
    !
    !  reference: R. Rosenfelder, arXiv: hep-ph/0603161
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !  21 March 2006
    !
    !
    integer::  j
    real*8 :: a,b,wfisi(0:2),y
    real*8 :: h, h2,fak,aint0,aint1,aint2,x0,x1,x2

    !
    h = 0.5d0*(b-a)
    h2 = h*h
    fak = 0.5d0/h2
    !
    !    moments
    !
    aint0 = aux(j,0,b,y) - aux(j,0,a,y)
    aint1 = aux(j,1,b,y) - aux(j,1,a,y)
    aint2 = aux(j,2,b,y) - aux(j,2,a,y)
    !
    !     weights
    !
    x0 = a
    x1 = a + h
    x2 = b

    wfisi(0) = fak*(x1*x2*aint0 - (x1+x2)*aint1 + aint2)
    wfisi(1) = 2.d0*fak*(-x0*x2*aint0 + 2.d0*x1*aint1 - aint2)
    wfisi(2) = fak*(x0*x1*aint0 - (x0+x1)*aint1 + aint2)
    !
    return
  end subroutine fil_simp
  !
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  subroutine fil_simp_sin(a,b,y,j,wfisi)
    !
    integer::  j
    real*8 :: a,b,wfisi(0:2),y
    real*8 :: h, h2,fak,aint0,aint1,aint2,x0,x1,x2

    !
    h = 0.5d0*(b-a)
    h2 = h*h
    fak = 0.5d0/h2
    !
    !    moments
    !
    aint0 = aux_sin(j,0,b,y) - aux_sin(j,0,a,y)
    aint1 = aux_sin(j,1,b,y) - aux_sin(j,1,a,y)
    aint2 = aux_sin(j,2,b,y) - aux_sin(j,2,a,y)
    !
    !     weights
    !
    x0 = a
    x1 = a + h
    x2 = b

    wfisi(0) = fak*(x1*x2*aint0 - (x1+x2)*aint1 + aint2)
    wfisi(1) = 2.d0*fak*(-x0*x2*aint0 + 2.d0*x1*aint1 - aint2)
    wfisi(2) = fak*(x0*x1*aint0 - (x0+x1)*aint1 + aint2)
    !
    return
  end subroutine fil_simp_sin

  !
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  subroutine fil_simp_cos(a,b,y,j,wfisi)
    !
    integer::  j
    real*8 :: a,b,wfisi(0:2),y
    real*8 :: h, h2,fak,aint0,aint1,aint2,x0,x1,x2

    !
    h = 0.5d0*(b-a)
    h2 = h*h
    fak = 0.5d0/h2
    !
    !    moments
    !
    aint0 = aux_cos(j,0,b,y) - aux_cos(j,0,a,y)
    aint1 = aux_cos(j,1,b,y) - aux_cos(j,1,a,y)
    aint2 = aux_cos(j,2,b,y) - aux_cos(j,2,a,y)
    !
    !     weights
    !
    x0 = a
    x1 = a + h
    x2 = b

    wfisi(0) = fak*(x1*x2*aint0 - (x1+x2)*aint1 + aint2)
    wfisi(1) = 2.d0*fak*(-x0*x2*aint0 + 2.d0*x1*aint1 - aint2)
    wfisi(2) = fak*(x0*x1*aint0 - (x0+x1)*aint1 + aint2)
    !
    return
  end subroutine fil_simp_cos


  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !
  function  aux(j,k,z,y) result(ddd) 
    !
    !  auxiliary function  int_0^z dx x^k O_j(x*y) ; k = 0, 1, 2
    !
    !  Feb 28, 2006


    integer :: j, k
    integer ::n,k1,n2k,n2
    real*8, parameter :: gamma = 0.577215664901533d0
    real*8 :: ddd,z,y,zy,sinzy,coszy,fak,cizy,sizy
    REAL*8 ::c0,zy2,c1,faz
    !
    k1 = k + 1
    zy = z*y
    if (zy .lt. 1.d-3) go to 100
    !
    sinzy = dsin(zy)
    coszy = dcos(zy)
    fak = j/(y**k1)

    call cisiluke(zy,cizy,sizy)
    !
    if(j .eq. 1) then
       if(k .eq. 0) then
          ddd = fak*sizy
          return
       endif
       if(k .eq. 1) then
          ddd = fak*(1.d0 - coszy)
          return
       endif
       if(k .eq. 2) then
          ddd = fak*(sinzy - zy*coszy)
          return
       endif
    endif
    !
    if(j .eq. 2) then
       if(k .eq. 0) then
          ddd = fak*(sizy - (1.d0-coszy)/zy)
          return
       endif
       if(k .eq. 1) then
          ddd = fak*(gamma + dlog(zy) - cizy)
          return
       endif
       if(k .eq. 2) then
          ddd = fak*(zy- sinzy)
          return
       endif
    endif
    !
    !  expansion for small arguments
    !
100 continue
    !

    c0 = 1.d0/k1
    ddd = c0
    n = 0
    zy2 = zy*zy
    !
10  n2 = 2*n
    n2k = n2 + k1
    fak = -zy2*n2k/((n2+j+1.d0)*(n2+j+2.d0)*(n2k+2.d0))
    c1 = c0*fak
    ddd = ddd + c1
    if(dabs(c1) .lt. 1.d-12) go to 20
    c0 = c1
    n = n + 1
    go to 10
    !
20  faz = 0.d0
    if(z .gt. 1.d-12) faz = z**k1
    ddd = ddd*faz

    return
  end  function aux

  !
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !

  function  aux_sin(j,k,z,y) result(ddd) 
    !
    !  auxiliary function  int_0^z dx x^k O_j(x*y) ; k = 0, 1, 2
    !
    !  Feb 28, 2006


    integer :: j, k
    integer ::n,k1,n2k,n2,k2
    real*8, parameter :: gamma = 0.577215664901533d0
    real*8 :: ddd,z,y,zy,sinzy,coszy,fak,cizy,sizy
    REAL*8 ::c0,zy2,c1,faz
    !
    k1 = k + 1
    zy = z*y
    if (zy .lt. 1.d-3) go to 100
    !
    sinzy = dsin(zy)
    coszy = dcos(zy)
    fak = 1.d0/(y**k1)
    !
    !if(j .eq. 1 or j.eq.2) then
    if(k .eq. 0) then
       ddd = fak*(1.d0-coszy)
       return
    endif
    if(k .eq. 1) then
       ddd = fak*(sinzy - zy*coszy)
       return
    endif
    if(k .eq. 2) then
       ddd = fak*(2.d0*zy*sinzy -(zy**2-2.d0)*coszy - 2.d0)
       return
    endif
    !
    !  expansion for small arguments
    !  
100 continue
    !
    k2=(k1+1.d0)  !since x^ksin(x) instead of x^ksin(x)/x
    c0 = zy/k2
    ddd = c0
    n = 0
    zy2 = zy*zy
    !
10  n2 = 2*n
    n2k = n2 + k2
    !genius way to take care of (2n+1)! n0 here --> n1 in series
    fak = -zy2*n2k/((n2+1.d0)*(n2+2.d0)*(n2k+2.d0))
    c1 = c0*fak
    ddd = ddd + c1
    if(dabs(c1) .lt. 1.d-12) go to 20
    c0 = c1
    n = n + 1
    go to 10
    !
20  faz = 0.d0
    if(z .gt. 1.d-12) faz = z**k1
    ddd = ddd*faz

    return
  end  function aux_sin
  !
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !
  function  aux_cos(j,k,z,y) result(ddd) 
    !
    !  auxiliary function  int_0^z dx x^k cos(x*y) ; k = 0, 1, 2
    !
    !  Feb 28, 2006


    integer :: j, k
    integer ::n,k1,n2k,n2,k2
    real*8, parameter :: gamma = 0.577215664901533d0
    real*8 :: ddd,z,y,zy,sinzy,coszy,fak,cizy,sizy
    REAL*8 ::c0,zy2,c1,faz
    !
    k1 = k + 1
    zy = z*y
    if (zy .lt. 1.d-3) go to 100
    !
    sinzy = dsin(zy)
    coszy = dcos(zy)
    fak = 1.d0/(y**k1)
    !
    !if(j .eq. 1 or j.eq.2) then
    if(k .eq. 0) then
       ddd = fak*sinzy
       return
    endif
    if(k .eq. 1) then
       ddd = fak*(zy*sinzy + coszy - 1.d0)
       return
    endif
    if(k .eq. 2) then
       ddd = fak*((zy**2-2.d0)*sinzy + 2.d0*zy*coszy)
       return
    endif
    !
    !  expansion for small arguments
    !  
100 continue
    !
    
    c0 = 1.d0/k1
    ddd = c0
    n = 0
    zy2 = zy*zy
    !
10  n2 = 2*n
    n2k = n2 + k1
    !genius way to take care of (2n+1)! n0 here --> n1 in series
    fak = zy2*n2k/((n2+1.d0)*(n2+2.d0)*(n2k+2.d0))
    c1 = c0*fak
    ddd = ddd + c1
    if(dabs(c1) .lt. 1.d-12) go to 20
    c0 = c1
    n = n + 1
    go to 10
    !
20  faz = 0.d0
    if(z .gt. 1.d-12) faz = z**k1
    ddd = ddd*faz

    return
  end  function aux_cos
  !
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !
  subroutine cisiluke(x,ci,si)
    !
    !   Chebyshev expansion of cosine and sine integral
    !   according to Luke, vol. 2, p. 325, 327; 
    !   coefficients retained down to 10^(-16)
    !
    !   Si(x) =  int_0^x dt sin(t)/t
    !   Ci(x) = gamma_Euler + ln x +  int_0^x dt [ cos(t) - 1]/t
    !
    !   definition of "handbook" p. 231 : Si(x)  --> pi/2 for large x; 
    !                                     Ci(x)  --> 0    for large x
    !
    !
    !   input:   x > 0 argument
    !   
    !   output:  ci -- Ci(x) cosine integral
    !            si -- Si(x) sine integral
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !  author: R. Rosenfelder
    !          Paul-Scherrer-Institut (PSI)
    !          CH-5232 Villigen PSI
    !          Switzerland
    !
    !  reference: R. Rosenfelder, arXiv: hep-ph/0603161
    !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !   21 March 2006
    !
    real*8 :: x,ci,si
    real*8 :: gamma, a(0:15), b(0:15),e(0:30),f(0:30)
    real*8 :: t(0:61)
    integer ::n,n2
    REAL*8 :: y,y2,fc,fs,p,q,help,hc,hs,x2

    data gamma /0.577215664901533d0/
    data a /1.9405491464835549d0,9.413409132865213d-1, &
         -5.798450342929928d-1,3.091572011159271d-1, &
         -9.16101792207713d-2,1.644374075154625d-2, &
         -1.9713091952164d-3,1.692538850835d-4,-1.09393295731d-5, &
         5.522385748d-7,-2.23994933d-8,7.465332d-10,-2.08183d-11, &
         4.931d-13,-1.005d-14,1.8d-16/
    data b /1.9522209759530711d0,-6.884042321257154d-1, &
         4.551855132255848d-1,-1.804571236838779d-1, &
         4.10422133758592d-2,-5.9586169555889d-3, &
         6.001427414144d-4,-4.44708329107d-5,2.53007823075d-6, &
         -1.141307593d-7,4.1857839d-9,-1.273471d-10,3.2674d-12, &
         -7.17d-14,1.36d-15,2.d-17/
    data e /9.6074783975203596d-1,-3.71138962123981d-2, &
         1.9414398889919d-3,-1.716598842515d-4,2.11263775323d-5, &
         -3.2716325671d-6,6.006921161d-7,-1.258679440d-7, &
         2.93256346d-8,-7.4569592d-9,2.0410548d-9,-5.950223d-10, &
         1.832297d-10,-5.92051d-11,1.99652d-11,-6.9951d-12, &
         2.5369d-12,-9.493d-13,3.655d-13,-1.445d-13,5.85d-14, &
         -2.42d-14,1.02d-14,-4.4d-15,1.9d-15,-8.7d-16,3.9d-16, &
         -1.8d-16,8.d-17,-4.d-17,2.d-17/
    data f /9.860406569623826d-1,-1.34717382082952d-2, &
         4.532928411652d-4,-3.06728865166d-5,3.1319919760d-6,  &
         -4.2110196496d-7,6.90724483d-8,-1.31832129d-8, &
         2.8369743d-9,-6.732923d-10,1.733969d-10,-4.78694d-11, &
         1.40323d-11,-4.3350d-12,1.4027d-12,-4.731d-13, &
         1.656d-13,-5.99d-14,2.24d-14,-8.6d-15,3.4d-15, &
         -1.4d-15,5.6d-16,-2.4d-16,1.d-16,-4.d-17,2.d-17, &
         0.,0.,0.,0./
    !data pi /3.1415926535897932d0/
    !

    if(x .le. 0.d0) then
       write(6,*) 'x < 0'
       stop
    endif
    !
    if(x .lt. 8.d0) then
       y = 0.125d0*x
       y2 = 2.d0*y
       !
       t(0) = 1.d0
       t(1) = y
       do  n = 1,30
          t(n+1) = y2*t(n) - t(n-1)
       enddo
       !
       fc = 0.d0
       fs = 0.d0
       do  n = 0,15
          n2 = 2*n
          fc = fc + a(n)*t(n2)
          fs = fs + b(n)*t(n2+1)
       enddo
       ci = gamma + dlog(x) - fc
       si = fs
    endif
    !

    if(x .ge. 8.d0) then
       y = 8.d0/x
       y2 = 2.d0*y
       !
       t(0) = 1.d0
       t(1) = y
       do  n = 1,60
          t(n+1) = y2*t(n) - t(n-1)
       enddo
       !
       p = 0.d0
       q = 0.d0
       do  n = 0,30
          help = t(2*n)
          p = p + e(n)*help
          q = q + f(n)*help
       enddo
       hc = dcos(x)
       hs = dsin(x)
       x2 = x*x
       ci = hs*q/x - hc*p/x2
       si = 0.5d0*pi - hs*p/x2 - hc*q/x
    endif
    return
  end subroutine cisiluke
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



end module common_module



!!========================================
!!    SINC(X) function with array/scalar single/double precision 
!!========================================
module sincm
  interface sinc
     module procedure sincas, sincss,sincad, sincsd
  end interface
contains

!!------ single precision array and scalar sinc ---
  function sincas(x) result(z) ! array
    implicit none
    real, dimension(:) :: x
    real, dimension(size(x)) :: z
    where(abs(x).lt.1e-10)
       z=1.
    elsewhere
       z=sin(x)/x
    endwhere
  end function sincas
  function sincss(x) result(z) ! scalar
    implicit none
    real :: x,z
    if(abs(x).lt.1e-10) then
       z=1.
    else
       z=sin(x)/x
    endif
  end function sincss

!!------ double precision array and scalar sinc ---
  function sincad(x) result(z) ! array
    implicit none
    real(kind=8), dimension(:) :: x
    real(kind=8), dimension(size(x)) :: z
    where(abs(x).lt.1e-10)
       z=1.
    elsewhere
       z=sin(x)/x
    endwhere
  end function sincad
  function sincsd(x) result(z) ! scalar
    implicit none
    real(kind=8) :: x,z
    if(abs(x).lt.1e-10) then
       z=1.
    else
       z=sin(x)/x
    endif
  end function sincsd
end module sincm






