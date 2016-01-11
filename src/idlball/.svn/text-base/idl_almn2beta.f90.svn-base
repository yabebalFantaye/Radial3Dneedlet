
!!!IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006

use healpix_types
use pix_tools
use extension
use alm_tools
use wav_mod

IMPLICIT NONE

INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,j0,nj,p,wav,large,m1,m2
REAL(SP), DIMENSION(:,:,:), ALLOCATABLE :: map_TQU
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
REAL(SP) :: fwhm,sigma
REAL(DP) :: bb,norm
CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum

print*, 'idl_alm2beta code.f90 starting ... ' 

call getEnvironment("HEALPIX",healpixdir)
dir=''!trim(adjustl(healpixdir))//'/src/idl/almtools/'
open(12,file=trim(adjustl(dir))//'bb.unf',status='old',form='unformatted')
read(12) bb
close(12)
open(12,file=trim(adjustl(dir))//'j0.unf',status='old',form='unformatted')
read(12) j0
close(12)
open(12,file=trim(adjustl(dir))//'nj.unf',status='old',form='unformatted')
read(12) nj
close(12)
open(12,file=trim(adjustl(dir))//'p.unf',status='old',form='unformatted')
read(12) p
close(12)
open(12,file=trim(adjustl(dir))//'wav.unf',status='old',form='unformatted')
read(12) wav
close(12)
open(12,file=trim(adjustl(dir))//'nside.unf',status='old',form='unformatted')
read(12) nside
close(12)
open(12,file=trim(adjustl(dir))//'lmax.unf',status='old',form='unformatted')
read(12) lmax
close(12)
open(12,file=trim(adjustl(dir))//'fwhm_beta.unf',status='old',form='unformatted')
read(12) fwhm
close(12)
open(12,file=trim(adjustl(dir))//'large.unf',status='old',form='unformatted')
read(12) large
close(12)
open(12,file=trim(adjustl(dir))//'m1.unf',status='old',form='unformatted')
read(12) m1
close(12)
open(12,file=trim(adjustl(dir))//'m2.unf',status='old',form='unformatted')
read(12) m2
close(12)
!open(12,file=trim(adjustl(dir))//'polar.unf',status='old',form='unformatted')
!read(12) polar
!close(12)
polar=0

print*, 'files read ..'

sigma=fwhm/60./180.*pi/sqrt(8.*log(2.))
npix=nside**2*12
 
ALLOCATE(map_TQU(0:npix-1,1:1+polar*2,0:nj-1),alm_TGC(1:1+polar*2,0:lmax,0:lmax),beam(0:lmax,1:1+polar*2))
ALLOCATE(gl(0:lmax,j0:j0+nj-1))

print*, 'map and gl arrays allocated!'

open(12,file=trim(adjustl(dir))//'beam_beta.unf',status='old',form='unformatted')
read(12) beam
close(12)
if (large.eq.0) open(12,file=trim(adjustl(dir))//'alm_beta.unf',status='old',form='unformatted')
if (large.eq.1) open(12,file=trim(adjustl(dir))//'large/alm_beta.unf',status='old',form='unformatted')
read(12) alm_TGC
close(12)
 
print*, 'alm arrays read!'

do l=0,lmax
alm_TGC(1,l,0:l)=alm_TGC(1,l,0:l)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,1)
if (polar.eq.1) then
alm_TGC(2,l,0:l)=alm_TGC(2,l,0:l)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,2)
alm_TGC(3,l,0:l)=alm_TGC(3,l,0:l)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,3)
endif
enddo

print*, 'alm_TGC loop done!'

if (wav.eq.0) wavtyp='standard'
if (wav.eq.1) wavtyp='mexican'

CALL calc_f2


!print*,'alm2beta1:',j0,nj,lmax,bb,wavtyp,p

CALL calc_gl(j0,nj,lmax,gl,bb,wavtyp,p)

if ((m1.eq.1).and.(m2.eq.0)) then
do l=0,lmax
norm=sum(gl(l,:)**2)
if (norm.gt.0d0) gl(l,:)=gl(l,:)/sqrt(norm)
enddo
endif

if (m2.eq.1) then
norm=sum(gl(2:lmax/8,:)**2)/(REAL(lmax/8-2+1,DP))
gl=gl/sqrt(norm)
endif

print*, 'calling synmap_multi ...'

CALL synmap_multi(alm_TGC(1,:,:),map_TQU,nside,lmax,gl,nj)

print*, 'synmap_multi done'
!print*,'alm2beta2:',sum(gl),sum(alm_TGC),sum(map_TQU)
!call synmap(alm,map,nside,lmax)
!if (polar.eq.0) call alm2map(nside,lmax,lmax,alm_TGC,map_TQU(:,1))
!if (polar.eq.1) call alm2map(nside,lmax,lmax,alm_TGC,map_TQU)

do j=0,nj-1
WRITE(filenum,'(I4)') j + 1000
if (large.eq.1) open(12,file=trim(adjustl(dir))//'large/beta'//filenum(2:4)//'.unf',status='unknown',form='unformatted')
if (large.eq.0) open(12,file=trim(adjustl(dir))//'beta'//filenum(2:4)//'.unf',status='unknown',form='unformatted')
write(12) map_TQU(:,1,j)
close(12)
enddo

end

