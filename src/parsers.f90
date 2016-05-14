
MODULE parsers

use healpix_types

IMPLICIT NONE
integer(i4b), parameter :: FILELEN=2000,inscribe=1,outscribe=2
integer(i4b), parameter :: full_balltype=2,slice_balltype=1

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

!-------------------------------
SUBROUTINE almnparser(filename,nside, lmax, nnmax, nshell,ifirstshell,ilastshell,iwidth,norm,balltype,inputfile,outputfile,finext,tempdir)
  character(len=FILELEN) :: line,name,value,finext,tempdir
  character(len=FILELEN) :: filename,inputfile,outputfile
  real(dp) :: norm
  integer(i4b) :: balltype,i
  integer(i4b) :: nside, lmax, nnmax, nshell,ifirstshell,ilastshell,iwidth
 
    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if


    ! --- set defaults:
    nside=0
    lmax=0
    nnmax=0
    nshell=0
    ifirstshell=0
    ilastshell=0
    iwidth=4
    norm=1d0
    balltype=1
    inputfile=''
    outputfile=''
    finext='unf'
    tempdir=''

    !now read parfile. Consider checking for sanity
    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('nside')
          read (value,*) nside
       case ('lmax')
          read (value,*) lmax
       case ('nnmax')
          read (value,*) nnmax
       case ('nshell')
          read (value,*) nshell
       case ('ifirstshell')
          read (value,*) ifirstshell
       case ('ilastshell')
          read (value,*) ilastshell
       case ('iwidth')
          read (value,*) iwidth
       case ('norm')
          read (value,*) norm
       case ('balltype')  !1 - shells per file; not 1 - all shells in one file
          read (value,*) balltype
       case ('inputfile')
          inputfile=trim(adjustl(value))
       case ('finext')
          finext=trim(adjustl(value))
       case ('outputfile')
          outputfile=trim(adjustl(value))
       case ('tempdir')
          tempdir=trim(adjustl(value))
       end select

    end do
11  close (1)

    nnmax=nshell

  end SUBROUTINE almnparser

!--------------------------------

SUBROUTINE betaparser(filename, bb, j0, nj, glnpow, wav,nside_jmap)
  character(len=FILELEN) :: line,name,value,filename
  integer(i4b) :: i
  !params
  character(len=FILELEN) :: finext,tempdir
  character(len=FILELEN) :: inputfile,outputfile
  real(dp) :: bb
  integer(i4b) :: j0, nj, glnpow, wav,nside_jmap
 
    inquire(file=filename, exist=exist)
    if (.not. exist) then
       print 1, trim(filename)
1      format(/,&
            & " Error in parser:",/,&
            & " File ",a," does not exist.")
       stop 1
    end if


    ! --- set defaults:
    j0=0
    bb=0
    nj=0
    glnpow=1
    wav=0
    nside_jmap=0

    !now read parfile. Consider checking for sanity
    open  (1,file=filename)
    do
       read (1,'(a)',end=11) line
       i=scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))

       select case (trim(name))
       case ('nside_jmap')
          read (value,*) nside_jmap
       case ('j0')
          read (value,*) j0
       case ('nj')
          read (value,*) nj
       case ('wav')
          read (value,*) wav
       case ('glnpow')
          read (value,*) glnpow
       case ('bb')
          read (value,*) bb
       end select
    end do
11  close (1)


  end SUBROUTINE betaparser

!--------------------------------

end MODULE parsers
