
!!!IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006

use healpix_types
use pix_tools
use extension
use alm_tools
use wav_mod

IMPLICIT NONE

INTEGER(I4B) :: lmax,i,j,l,j0,nj,p,wav,nside,nmax
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln
REAL(DP), DIMENSION(:), ALLOCATABLE :: smhw
REAL(DP) :: bb
CHARACTER(LEN=128) :: healpixdir,dir,wavtyp

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
open(12,file=trim(adjustl(dir))//'lmax.unf',status='old',form='unformatted')
read(12) lmax
close(12)

open(12,file=trim(adjustl(dir))//'nmax.unf',status='old',form='unformatted')
read(12) nmax
close(12)


if (wav.eq.2) j0=0

ALLOCATE(gln(0:lmax,0:nmax,j0:j0+nj-1))

if (wav.eq.0) wavtyp='standard'
if (wav.eq.1) wavtyp='mexican'
if (wav.eq.2) wavtyp='smhw'

if (wav.eq.2) then

   print *,'idl_get_needlet_ball: swhm not implemented for a ball'
   stop
   ! open(12,file=trim(adjustl(dir))//'nside.unf',status='old',form='unformatted')
   ! read(12) nside
   ! close(12)
   ! ALLOCATE(smhw(0:nj-1))
   ! open(12,file=trim(adjustl(dir))//'smhw.unf',status='old',form='unformatted')
   ! read(12) smhw
   ! close(12)
   ! CALL calc_gln(j0,nj,lmax,nmax,gln,bb,wavtyp,p,nside,smhw)


else


   CALL calc_f2


   !print*,'alm2beta1:',j0,nj,lmax,bb,wavtyp,p
   CALL calc_gln(j0,nj,lmax,nmax,gln,bb,wavtyp,p)

endif

open(12,file=trim(adjustl(dir))//'gln.unf',status='unknown',form='unformatted')
write(12) gln
close(12)


end

