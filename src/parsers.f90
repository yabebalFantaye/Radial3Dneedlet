
MODULE parsers

use healpix_types

IMPLICIT NONE
integer(i4b), parameter :: FILELEN=128,inscribe=1,outscribe=2
logical :: exist

CONTAINS


SUBROUTINE box2ballparser(filename,fname_box,Rbox,ndata,balltype)
  character(len=FILELEN) :: line,name,value
  character(len=FILELEN) :: filename,fname_box
  real(dp) :: Rbox
  integer(i4b) :: ndata,i,balltype

 
    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if


    ! --- set defaults:

    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('ndata')
          read (value,*) ndata
       case ('fname_box')
          fname_box = trim(value)
       case ('Rbox')
          read (value,*) Rbox
       case ('balltype')
          read (value,*) balltype
       end select
    end do
11  close (1)

  end SUBROUTINE box2ballparser

!--------------------------------

SUBROUTINE ball2pixparser(filename,fname_box,fmapout,nrint,nside,ordering,rnmin,rnmax)
  character(len=FILELEN) :: line,name,value
  character(len=FILELEN) :: filename,fname_box,fmapout,ordering
  integer(i4b) :: i, nrint,nside
  real(dp) ::rnmin,rnmax
 
    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if


    ! --- set defaults:

    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('fname_box')
          fname_box = trim(value)
       case ('fmapout')
          fmapout = trim(value)
       case ('ordering')
          ordering = trim(value)
       case ('nrint')
          read (value,*) nrint
       case ('nside')
          read (value,*) nside
       case ('rnmin')
          read (value,*) rnmin
       case ('rnmax')
          read (value,*) rnmax
       end select
    end do
11  close (1)

  end SUBROUTINE ball2pixparser



end MODULE parsers
