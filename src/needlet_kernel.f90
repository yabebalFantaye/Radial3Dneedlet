module needlet_kernel

  use healpix_types

  IMPLICIT NONE

  INTEGER(I4B) :: nmax_f2=10000000
  REAL(DP), DIMENSION(:), ALLOCATABLE :: f2


contains


  subroutine calc_f2

    INTEGER(I4B) :: i
    REAL(DP) :: x

    ALLOCATE(f2(0:nmax_f2-1))

    x=-1.d0
    !f2(0)=exp(-1.d0/(1.d0-x**2))
    f2(0)=0.d0
    do i=1,nmax_f2-1
       x=2.d0/REAL(nmax_f2,DP)*i-1.d0
       f2(i)=f2(i-1)+exp(-1.d0/(1.d0-x**2))
    enddo

    f2=f2/f2(nmax_f2-1)

  end subroutine calc_f2

!--------------------------------

  subroutine calc_phiphi(phiphi,x,aa)

    REAL(DP), INTENT(OUT) :: phiphi
    REAL(DP), INTENT(IN) :: x,aa
    INTEGER(I4B) :: ii
    phiphi=0.d0

    if ((x.gt.-1.d0).and.(x.lt.-1d0/aa)) then
       !ii=NINT(REAL(n,DP)/2.d0*(1.d0+3.d0+4.d0*x))
       stop 'STOP, x lt -1/aa, need formula for this case!'
       if (ii.gt.(nmax_f2-1)) ii=nmax_f2-1
       phiphi=f2(ii)
    endif

    if ((x.ge.-1d0/aa).and.(x.le.1d0/aa)) phiphi=1.d0

    if ((x.le.1.d0).and.(x.gt.1d0/aa)) then
       !ii=NINT(REAL(nmax_f2,DP)/2.d0*(1.d0+4.d0*x-3.d0))
       ii=NINT(REAL(nmax_f2,DP)*(1d0-(aa/(aa-1d0)*(x-1d0/aa))))  !!!the figure 2 in here has disappeared because we try to find an index which runs from 0 to n-1 instead of a number between -1 and 1.
       if (ii.gt.(nmax_f2-1)) ii=nmax_f2-1
       phiphi=f2(ii)
    endif

  end subroutine calc_phiphi

!--------------------------------

  subroutine calc_phi(phi,x,aa)

    REAL(DP), INTENT(OUT) :: phi
    REAL(DP), INTENT(IN) :: x,aa
    REAL(DP) :: xx,phi1,phi2

    xx=x/aa
    call calc_phiphi(phi1,xx,aa)
    xx=x
    call calc_phiphi(phi2,xx,aa)
    phi=phi1-phi2
    !print*,phi1,phi2
  end subroutine calc_phi

!--------------------------------

  subroutine calc_gl_stand(j0,nj,lmax,gl,aa)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:lmax,j0:j0+nj-1), INTENT(OUT) :: gl
    INTEGER(I4B) :: j,l0,l1,l
    REAL(DP) :: xx,phi


    gl=0.d0

    do j=j0,j0+nj-1
       l0=NINT(aa**(j-1))
       l1=NINT(aa**(j+1))
       if (l0.gt.lmax) stop 'first l for j greater than lmax!!!'
       if (l1.gt.lmax) l1=lmax
       do l=l0,l1
          xx=REAL(l,DP)/aa**REAL(j,DP)
!!!xx=(3.d0*aa*xx+aa**2-4.d0)/(2.d0*(aa**2-1.d0))
          !!if (j.eq.1) print*,'aa',l,xx
          call calc_phi(phi,xx,aa)
          gl(l,j)=phi
       enddo
    enddo

    gl=sqrt(gl)

    if ((aa.lt.2d0).and.(j0.eq.0)) gl(0,0)=1d0

  end subroutine calc_gl_stand

!--------------------------------

  subroutine calc_gl_mex(j0,nj,lmax,gl,aa,p)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,p
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:lmax,j0:j0+nj-1), INTENT(OUT) :: gl
    INTEGER(I4B) :: j,l0,l1,l


    gl=0.d0

    do j=j0,j0+nj-1
       do l=0,lmax
          gl(l,j)=(REAL(l,DP)/aa**j)**(2d0*p)*exp(-REAL(l,DP)**2/aa**(2d0*j))
       enddo
    enddo

  end subroutine calc_gl_mex

!--------------------------------

  subroutine calc_gl(j0,nj,lmax,gl,aa,wavtyp,p,nside,smhw)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,p
    INTEGER(I4B), INTENT(IN), OPTIONAL :: nside
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:lmax,j0:j0+nj-1), INTENT(OUT) :: gl
    REAL(DP),  DIMENSION(j0:j0+nj-1), INTENT(IN), OPTIONAL :: smhw
    CHARACTER(LEN=128), INTENT(IN) :: wavtyp

    if (trim(adjustl(wavtyp)).eq.'standard') call calc_gl_stand(j0,nj,lmax,gl,aa)
    if (trim(adjustl(wavtyp)).eq.'mexican')  call calc_gl_mex(j0,nj,lmax,gl,aa,p)
    if ((trim(adjustl(wavtyp)).eq.'smhw').and.(present(nside)).and.(present(smhw))) &
         & call calc_smhw(nside*4,lmax,nj,smhw,gl,.false.,.false.,wavtyp)

    if ((trim(adjustl(wavtyp)).eq.'smhw').and.((.not.present(nside)).or.&
         &(.not.present(smhw))))  stop 'SMHW, but not all necessary input to calc_gl!'

    if ((trim(adjustl(wavtyp)).ne.'standard').and.(trim(adjustl(wavtyp)).ne.'mexican').and.&
         &(trim(adjustl(wavtyp)).ne.'smhw')) stop 'wavelet type not specified in calc_gl'

  end subroutine calc_gl

!--------------------------------

  subroutine calc_smhw(nside,lmax,nscales,scales,g_real_tot,loadwav,savewav,filename)

    use pix_tools
    use alm_tools

    IMPLICIT NONE


    INTEGER(I4B), INTENT(IN) :: nside,lmax,nscales
    REAL(DP),DIMENSION(0:lmax,0:nscales-1), INTENT(OUT) :: g_real_tot
    CHARACTER(LEN=128), INTENT(IN) :: filename
    LOGICAL(LGT), INTENT(IN) :: loadwav,savewav
    REAL(DP),DIMENSION(0:nscales-1), INTENT(IN) :: scales
    COMPLEX(SPC), ALLOCATABLE, DIMENSION(:,:,:) :: almw
    REAL(SP),ALLOCATABLE,DIMENSION(:) ::wav_pix
    REAL(DP) :: r,nr,y,theta,phi
    INTEGER(I4B) :: n_pix_T,i_pix,irisol,i, n
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring
    REAL(DP), DIMENSION(1:2) :: z = 0.0d0


    if (loadwav) then

       open(12,file=filename,status='old',form='unformatted')
       read(12) g_real_tot
       close(12)

    else


       n_pix_T=nside**2*12

       ALLOCATE(w8ring(1:2*nside,1:1))
       w8ring=1d0

       ALLOCATE(wav_pix(0:n_pix_T-1))
       ALLOCATE(almw(1:1,0:lmax,0:lmax))

       DO irisol=0,nscales-1

          !!     r = (irisol+1)*r0

          if (scales(irisol).ne.0.) then

             r=scales(irisol)/60.d0/180.d0*pi
             nr=r*SQRT(1.+r**2./2.+r**4./4.) 

             DO i_pix = 0,n_pix_T-1 

                CALL pix2ang_ring(nside,i_pix,theta,phi)

                y =2.*SIN(theta/2.)/COS(theta/2.)

                if ((-y**2./2./r**2.).gt.(-700.d0)) then
                   wav_pix(i_pix)=1./nr/SQRT(2.*pi)*(1.+(y/2.)**2.)**2.*(2.- (y/r)**2. )*EXP(-y**2./2./r**2.)
                else
                   wav_pix(i_pix)=0.d0
                endif


             END DO     !i_pix


             CALL map2alm(nside, lmax, lmax, wav_pix, almw, z, w8ring)  

             g_real_tot(:,irisol) = REAL(almw(1,:,0),DP)


             DO i=0,lmax
                g_real_tot(i,irisol) = g_real_tot(i,irisol)/SQRT(2d0*REAL(i,DP)+1)*SQRT(4d0*pi)
             END DO

          else
             g_real_tot(0:lmax,irisol)=1.
          endif


       END DO  !irisol

       if (savewav) then
          open(12,file=filename,status='unknown',form='unformatted')
          write(12) g_real_tot
          close(12)
       endif

       DEALLOCATE(wav_pix,almw,w8ring)

    endif


  END subroutine calc_smhw

  !==============================================
  !----- 3D Radial Needlet window functions -----
  !==============================================

  subroutine calc_gln_stand(j0,nj,lmax,nmax,gln,aa,gln4fft)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,nmax
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:,0:,j0:), INTENT(OUT) :: gln
    INTEGER(I4B) :: j,eln0,eln1,l,n,isfft,two_nmax
    REAL(DP) :: xx,phi,eln,elnmax,norm
    INTEGER(I4B),optional::gln4fft


    isfft=0
    if (present(gln4fft)) isfft=1

    !if isfft then we should store gln such that the radial frequency ordering is 
    !0:nngrid = 1,2 .. N/2, -N/2-1,..-2, -1
    two_nmax=2*nmax

    gln=0.d0

    elnmax = real(nmax**2,DP)+REAL(lmax,DP)*(REAL(lmax,DP)+1)

    do j=j0,j0+nj-1
       !if (feedback>2) print*,'gln j=',j

       eln0=NINT(aa**(j-1))
       eln1=NINT(aa**(j+1))

       if (eln0.gt.elnmax) print*, 'WARNNING: first eln for j greater than elnmax!!!, eln0, elnmax: ', eln0, elnmax
       if (eln1.gt.elnmax) print*, 'WARNNING: last eln for j greater than elnmax!!!, eln1m elnmax: ', eln0, elnmax

       do n=0,nmax

          do l=0,lmax

             eln = real(n**2,DP)+REAL(l,DP)*(REAL(l,DP)+1)


             xx=sqrt(real(n**2,DP)+REAL(l,DP)*(REAL(l,DP)+1))/aa**REAL(j,DP)

             call calc_phi(phi,xx,aa)
             gln(l,n,j)=phi

             !print*,'j, n, l, aa^j, sqrt(n2+l(l+1)),xx:',j,n,l,aa**REAL(j,DP),xx*aa**REAL(j,DP)
          enddo
       enddo
    enddo

    do n=0,nmax
       do l=0,lmax
          norm=sum(gln(l,n,:))
          if (norm.gt.0d0) gln(l,n,:)=gln(l,n,:)/norm
       enddo
    enddo
    gln=sqrt(gln)

    !if isfft then we should store gln such that the radial frequency ordering is 
    !0:nngrid = 1,2 .. N/2, -N/2-1,..-2, -1
    if (isfft==1) then
       do n=nmax+1,two_nmax-1
          gln(0:lmax,n,j0:j0+nj-1) = gln(0:lmax,two_nmax-n,j0:j0+nj-1)
       enddo
    endif

    if ((aa.lt.2d0).and.(j0.eq.0)) gln(0,:,0)=1d0

  end subroutine calc_gln_stand

  !==============================================

  subroutine calc_gln_mex(j0,nj,lmax,nmax,gln,aa,p,gln4fft)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,p,nmax
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:,0:,j0:), INTENT(OUT) :: gln
    REAL(DP) :: xx,phi,eln,elnmax,norm
    INTEGER(I4B) :: j,l0,l1,l,n
    INTEGER(I4B) :: isfft,two_nmax
    INTEGER(I4B),optional::gln4fft

    isfft=0
    if (present(gln4fft)) isfft=1

    !if isfft then we should store gln such that the radial frequency ordering is 
    !0:nngrid = 1,2 .. N/2, -N/2-1,..-2, -1
    two_nmax=2*nmax


    gln=0.d0

    do j=j0,j0+nj-1
       do n=0,nmax
          do l=0,lmax
             gln(l,n,j)=(sqrt(real(n**2,DP)+REAL(l,DP)*(REAL(l,DP)+1))/aa**j)**(2d0*p)*exp(-REAL(l,DP)**2/aa**(2d0*j))
          enddo
       enddo
    enddo

    do n=0,nmax
       do l=0,lmax
          norm=sum(gln(l,n,:))
          if (norm.gt.0d0) gln(l,n,:)=gln(l,n,:)/norm
       enddo
    enddo
    gln=sqrt(gln)

    !if isfft then we should store gln such that the radial frequency ordering is 
    !0:nngrid = 1,2 .. N/2, -N/2-1,..-2, -1
    if (isfft==1) then
       do n=nmax+1,two_nmax-1
          gln(0:lmax,n,j0:j0+nj-1) = gln(0:lmax,two_nmax-n,j0:j0+nj-1)
       enddo
    endif


  end subroutine calc_gln_mex

  !==============================================

  subroutine calc_gln(j0,nj,lmax,nmax,gln,aa,wavtyp,p,nside,smhw,gln4fft)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,p,nmax
    INTEGER(I4B), INTENT(IN), OPTIONAL :: nside,gln4fft
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:,0:,j0:), INTENT(OUT) :: gln
    REAL(DP),  DIMENSION(j0:j0+nj-1), INTENT(IN), OPTIONAL :: smhw
    CHARACTER(LEN=*), INTENT(IN) :: wavtyp

    if (present(gln4fft)) then
       if (trim(adjustl(wavtyp)).eq.'standard') call calc_gln_stand(j0,nj,lmax,nmax,gln,aa,gln4fft)
    else
       if (trim(adjustl(wavtyp)).eq.'standard') call calc_gln_stand(j0,nj,lmax,nmax,gln,aa)
    endif

    if (trim(adjustl(wavtyp)).eq.'mexican')  call calc_gln_mex(j0,nj,lmax,nmax,gln,aa,p)

    if ((trim(adjustl(wavtyp)).eq.'smhw').and.(present(nside)).and.(present(smhw))) &
         &  stop 'calc_swhm not implemented for needlet in a ball' 

    if ((trim(adjustl(wavtyp)).eq.'smhw').and.((.not.present(nside)).or.(.not.present(smhw)))) &
         & stop 'SMHW, but not all necessary input to calc_gl!'

    if ((trim(adjustl(wavtyp)).ne.'standard').and.(trim(adjustl(wavtyp)).ne.'mexican').and.&
         &(trim(adjustl(wavtyp)).ne.'smhw')) stop 'wavelet type not specified in calc_gl'

  end subroutine calc_gln



end module needlet_kernel
