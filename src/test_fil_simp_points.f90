!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!         **************************************
!         *                                    *
!         *   program test_fil_simp_points.f   *
!         *                                    *
!         **************************************
!
!
!  tests Filon-Simpson integration program fil_simp_points.f      
!  for ( Dwight p. 234, 235 )
!
!  j = 1:
!
!  I1_f0 = int_0^infinity dx exp(-x) 1/x sin(y*x)/y = arctan(y)/y 
!
!  I1_f1 = int_0^infinity dx exp(-x) sin(y*x)/y = 1/(1 + y^2)
!
!*********************************************************************
!
!  j = 2: 
!
!  I2_f0 = int_0^infinity dx exp(-x) 1/x^2 4*sin^2(y*x/2)/y^2 = 
!                                       2*arctan(y)/y  - ln(1+y^2/)/y^2
!
!  I2_f1 = int_0^infinity dx exp(-x) 1/x 4*sin^2(y*x/2)/y^2=ln(1+y^2)/y^2
!
!
!***********************************************************************
!
!
!  input:    b -- upper limit of integral to be evaluated numerically
!            nn -- number of intervals (even, < 1000)
!            dy -- increment of y-values
!            ny -- number of y-values
! 
!  leading order asymptoti! contributions for integrals from b...infty
!  added to numerically evaluated integrals
!
!  output:   y, relative deviations from exact value
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  author: R. Rosenfelder
!          Paul-Scherrer-Institut (PSI)
!          CH-5232 Villigen PSI
!          Switzerland
!
!  reference: R. Rosenfelder, arXiv: hep-th/0603...
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  March 17, 2006: 

program test_fil_simp_points!
  implicit real*8 (a-h,o-z)
  parameter(nmax=1000)
  dimension xfisi(0:nmax),wfisi(0:nmax)
  dimension dev0(2),dev1(2)

  pi =  3.14159265358979324E0;

  !write(6,*) 'enter: b,nn(even),dy,ny'
  !read(5,*) b,nn,dy,ny
!  b=1
!  nn=10
!  dy=pi

b = 20
nn = 288
dy = 0.5
ny=10
  write(6,10)  b,nn




10 format(//'   Filon-Simpson integration with xmax = ',f5.1,5x, 'nn = ',i4//)

  !write(6,11)
11 format('   relative deviations:')
  write(6,12)
12 format(/'       y',8x,'I1_f0',8x,'I1_f1', 8x,'I2_f0',8x,'I2_f1'/)

  xmax = b
  h = xmax/nn

  do  n = 0,ny
     y = dy*n
     y2 = y*y

     do  j = 1,2
        s0 = 0.d0
        s1 = 0.d0
        !
        call fil_simp_points(0.d0,xmax,y,nn,j,xfisi,wfisi)   
        !
        do  i = 0,nn
           x = xfisi(i)
           x2 = x*x
           hilf = dexp(-x)
           s0 = s0 + wfisi(i)*hilf
           s1 = s1 + wfisi(i)*hilf*x
        enddo
        !
        !  exact values
        !
        hilf = 1.d0 + y2
        helf = 1.d0
        if(y2 .gt. 1.d-20) helf = datan(y)/y
        sinxy = dsin(xmax*y)
        cosxy = dcos(xmax*y)
        !
        if(j .eq. 1) then
           s1ex = 1.d0/hilf
           s0ex = helf
           !
           !  asymptotic contribution to I1
           !
           fak = xmax
           if(n .ne. 0) fak = sinxy/y
           fak = cosxy + fak
        endif

        if(j .eq. 2) then
           if(n .eq. 0) then
              s1ex = 1.d0
              s0ex = 1.d0
           else
              s1ex = dlog(hilf)/y2              
              s0ex = 2.d0*helf - dlog(hilf)/y2
           endif
           !
           !  asymptotic contribution to I2
           !
           fak = 0.5d0*xmax*xmax + xmax + 1.d0
           if(n .ne. 0) fak = (hilf - cosxy + y*sinxy)/y2
           fak = fak/xmax
        endif
        !
        asy = j*dexp(-xmax)*fak/hilf
        s0 = s0 + asy/xmax 
        s1 = s1 + asy
        !
        !  relative deviation
        !
        dev0(j) = (s0 - s0ex)/s0ex
        dev1(j) = (s1 - s1ex)/s1ex
     enddo
     !
     write(6,20) y,dev0(1),dev1(1), dev0(2),dev1(2)
20   format(3x,f6.1,4(3x,d10.3))
  enddo

  stop

!  include'fil_simp_points.f'
  !
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
end program test_fil_simp_points!


