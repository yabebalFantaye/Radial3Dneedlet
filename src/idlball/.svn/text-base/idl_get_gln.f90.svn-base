
!!!IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006

use healpix_types
use pix_tools
use extension
use alm_tools
use wav_mod

IMPLICIT NONE

INTEGER(I4B) :: n0,nn,nmax,nr,ii,jj
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sn
REAL(DP), DIMENSION(:), ALLOCATABLE :: rr
REAL(DP) ::sqrt2_to_3
CHARACTER(LEN=128) :: healpixdir,dir

call getEnvironment("HEALPIX",healpixdir)
dir=''!trim(adjustl(healpixdir))//'/src/idl/almtools/'

open(12,file=trim(adjustl(dir))//'n0.unf',status='old',form='unformatted')
read(12) n0
close(12)
open(12,file=trim(adjustl(dir))//'nn.unf',status='old',form='unformatted')
read(12) nn
close(12)

open(12,file=trim(adjustl(dir))//'nmax.unf',status='old',form='unformatted')
read(12) nmax
close(12)

open(12,file=trim(adjustl(dir))//'nr.unf',status='old',form='unformatted')
read(12) nr
close(12)

ALLOCATE(rr(1:nr))
open(12,file=trim(adjustl(dir))//'rr.unf',status='old',form='unformatted')
read(12) rr
close(12)

ALLOCATE(sn(nr,n0:n0+nn-1))
sqrt2_to_3 = sqrt(2d0/3d0)


do jj=n0,n0+nn-1
   do ii=1,nr
      sn(ii,jj) = sqrt2_to_3*sin(pi*n*rr(ii))/rr(ii)
   enddo
enddo


open(12,file=trim(adjustl(dir))//'sn.unf',status='unknown',form='unformatted')
write(12) sn
close(12)


end

