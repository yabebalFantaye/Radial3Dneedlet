

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006

pro alm2beta_ball,alm,map,bb,j0,nj,wavtyp=wavtyp_in,mexican_p=p_in,nside=nside_in,beam=beam_in,fwhm=fwhm_in,silent=silent,large=large_in,modif=modif,softmod=softmod,nproc=nproc_in

dir='$HEALPIX/src/idl/almtools/'

dir0=dir
ldir=''
if (keyword_set(large_in)) then begin
    if (size(large_in,/tn) eq 'INT') then begin
        if (large_in eq 1) then ldir='large/'
    endif
    if (ldir eq '') then ldir=trim(string(large_in),2)+'/'
    if (ldir ne 'large/') then dir=dir+ldir
endif
direxe=dir


;direxe='$HEALPIX/bin/'
unit=99
status = FSTAT(unit)
if (status.open ne 0) then begin
    WHILE(status.open ne 0) do begin
        unit=unit-1
        status = FSTAT(unit)
    ENDWHILE
endif

routine = 'alm2beta'
uroutine = strupcase(routine)
if (n_params() ne 5) then begin
    PRINT, 'Wrong number of arguments in ' + uroutine
    print,'Syntax : '
    print,uroutine+', alm, map, B, j0,nj'
    print,'              [wavtyp=, mexican_p=, NSIDE=, BEAM=, FWHM=,SILENT] '
;    print,' Type '+uroutine+', /help '
;    print,'   for an extended help'
    return
endif

if (n_elements(alm) lt 1) then begin
    print,'alm MUST BE of the form alm=complexarr(1+2*pol,lmax+1,lmax+1) or (1+2*pol,lmax+1,lmax+1,nj)'
    return
endif

salm=size(alm)
salm=LONG(salm[0])
if ((salm lt 2) and (salm gt 4)) then begin
    print,'alm MUST BE of the form alm=complexarr(1+2*pol,lmax+1,lmax+1) or (1+2*pol,lmax+1,lmax+1,nj)'
    return
endif

multialm=0l
almx=alm
if (salm eq 4) then begin
    if (not KEYWORD_SET(nproc_in)) then begin
        print,'must set nproc if using multialm! Only par.version can be used!'
        return
    endif
    multialm=1l
    alm00=alm
    almx=(alm(*,*,*,0))
    salm=3
endif

if (salm eq 2) then lmax=n_elements(reform(almx(*,0)))
if (salm eq 3) then lmax=n_elements(reform(almx(0,*,0)))

if (((salm eq 2) and (lmax ne n_elements(reform(almx(0,*))))) or ((salm eq 3) and (lmax ne n_elements(reform(almx(0,0,*)))))) then begin
    print,'alm MUST BE of the form alm=complexarr(1+2*pol,lmax+1,lmax+1) or (1+2*pol,lmax+1,lmax+1,nj)'
    return
endif


polar=0l
;if (KEYWORD_SET(pol)) then polar=1l
if (KEYWORD_SET(pol)) then print,'polarization is not implemented yet!'


if ((polar eq 1) and (salm ne 3)) then begin
    print,'alm MUST BE of the form alm=complexarr(1+2*pol,lmax+1,lmax+1) or (1+2*pol,lmax+1,lmax+1,nj)'
    return
endif


if (salm eq 3) then begin
    if  ((n_elements(reform(almx(*,0,0))) ne 3l) and (n_elements(reform(almx(*,0,0))) ne 1l)) then begin
        print,'alm MUST BE of the form alm=complexarr(1+2*pol,lmax+1,lmax+1) or (1+2*pol,lmax+1,lmax+1,nj)'
        return
    endif
endif

if ((polar eq 1) and ((n_elements(reform(almx(*,0,0))) ne 3l))) then begin
    print,'alm MUST BE of the form alm=complexarr(1+2*pol,lmax+1,lmax+1) or (1+2*pol,lmax+1,lmax+1,nj)'
    return
endif

if (multialm eq 0l) then begin
    alm00=complexarr(1l+polar*2l,lmax,lmax)
    if ((salm eq 2) and (polar eq 0)) then alm00(0,*,*)=alm
    if ((salm eq 3) and (polar eq 0)) then alm00(0,*,*)=alm(0,*,*)
    if ((salm eq 3) and (polar eq 1)) then alm00=alm
endif


lmax=lmax-1
if (lmax lt 1) then begin
    print,'lmax MUST BE GREATER THAN 0!'
    return
endif
if (KEYWORD_SET(nside_in)) then begin
    nside=LONG(nside_in)
    map=fltarr(nside^2*12l)
endif else begin
    nside=sqrt(DOUBLE(n_elements(map))/12d0)
endelse


if (KEYWORD_SET(nproc_in)) then begin
    if (nproc_in gt nj) then begin
        print,'nproc must be smaller or equal to nj!'
        return
    endif
endif

if (KEYWORD_SET(wavtyp_in)) then begin
    if ((wavtyp_in ne 'standard') and (wavtyp_in ne 'mexican')) then begin
        print,'wavtyp must be standard or mexican'
        return
    endif
    wav=0l
    if (wavtyp_in eq 'mexican') then wav=1l
endif else begin
    wav=0l
endelse

if ((KEYWORD_SET(p_in)) and (wav eq 0)) then begin
    print,'mexican_p can only be set when wavelet type is mexican!'
    return
endif

p=1l
if (KEYWORD_SET(p_in) and (wav eq 1)) then p=LONG(p_in)

if (j0 lt 0) then begin
    print,'j0 must be larger than 0!'
    return
endif
j0=LONG(j0)

if (nj lt 1) then begin
    print,'nj must be larger than 0!'
    return
endif
nj=LONG(nj)

if (bb lt 0d0) then begin
    print,'B must be larger than 0!'
    return
endif
bb=DOUBLE(bb)


cns=alog(DOUBLE(nside))/alog(2d0)
if ((abs(cns-NINT(cns)) gt 1d-6) or (nside lt 1)) then begin
    print,'nside MUST BE of the form 2^N'
    return
endif
if (lmax gt 3*nside) then begin
    if (NOT (KEYWORD_SET(silent))) then print,'WARNING:lmax>3*nside will use lmax=3*nside',lmax, nside
    lmax=3*nside
endif
openw,unit,dir+'nside.unf',/f77
writeu,unit,LONG(nside)
close,unit
openw,unit,dir+'lmax.unf',/f77
writeu,unit,LONG(lmax)
close,unit
fwhm0=0.
if (KEYWORD_SET(fwhm_in)) then fwhm0=FLOAT(fwhm_in)
beam0=dblarr(lmax+1,1l+polar*2l)
beam0(*,*)=1d0
if (KEYWORD_SET(beam_in)) then begin
    beam=beam_in
    sbeam=SIZE(beam0)
    sbeam=LONG(sbeam[0])
    if ((sbeam lt 1) or (sbeam gt 2)) then begin
        print,'BEAM must be of the form beam=dblarr(lmax+1,1+pol)'
        return
    endif
    if (sbeam eq 1) then beam0(*,0)=beam(0:lmax)
    if (sbeam eq 2) then beam0(*,0)=beam(0:lmax,0)
    if (polar eq 1l) then begin
        if (sbeam eq 1) then beam0(*,1)=beam(0:lmax)
        if (sbeam eq 1) then beam0(*,2)=beam(0:lmax)
        if (sbeam eq 2) then beam0(*,1)=beam(0:lmax,1)
        if (sbeam eq 2) then beam0(*,2)=beam(0:lmax,1)
    endif
endif

openw,unit,dir+'multialm.unf',/f77
writeu,unit,multialm
close,unit
openw,unit,dir+'beam_beta.unf',/f77
writeu,unit,beam0
close,unit
openw,unit,dir+'fwhm_beta.unf',/f77
writeu,unit,fwhm0
close,unit
;openw,unit,dir+'polar.unf',/f77
;writeu,unit,polar
;close,unit
if (ldir eq 'large/') then openw,unit,dir+'large/alm_beta.unf',/f77 else openw,unit,dir+'alm_beta.unf',/f77
;writeu,unit,COMPLEX(alm00(0:polar*2l,0:lmax,0:lmax))
writeu,unit,COMPLEX(alm00(0,0:lmax,0:lmax,*))
close,unit
openw,unit,dir+'bb.unf',/f77
writeu,unit,bb
close,unit
openw,unit,dir+'j0.unf',/f77
writeu,unit,j0
close,unit
openw,unit,dir+'nj.unf',/f77
writeu,unit,nj
close,unit
openw,unit,dir+'p.unf',/f77
writeu,unit,p
close,unit
openw,unit,dir+'wav.unf',/f77
writeu,unit,wav
close,unit
m1=0l
if (KEYWORD_SET(modif)) then m1=1l
m2=0l
if (KEYWORD_SET(softmod)) then m2=1l
openw,unit,dir+'m1.unf',/f77
writeu,unit,m1
close,unit
openw,unit,dir+'m2.unf',/f77
writeu,unit,m2
close,unit
openw,unit,dir+'large.unf',/f77
if (ldir eq 'large/') then writeu,unit,1l else writeu,unit,0l
close,unit
if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_alm2beta '+direxe
if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_alm2beta_par '+direxe

if (not KEYWORD_SET(nproc_in)) then SPAWN,'cd '+direxe+' ; ./idl_alm2beta'
if (KEYWORD_SET(nproc_in)) then begin
    nodelist=direxe+'nodelist'
    str='cd '+direxe+' ; mpiexec -machinefile '+nodelist+'  -n '+trim(string(nproc_in),2)+' '+direxe+'idl_alm2beta_par'
    SPAWN,str
endif

;map=fltarr(nside^2*12l,1+polar*2l,nj)
if (ldir eq 'large/') then dir=dir+'large/'
map=fltarr(nside^2*12l,nj)
mp=fltarr(nside^2*12l)
for j=0l,nj-1l do begin
    nm='00'+trim(string(j),2)
    lnm=strlen(nm)
    openr,unit,dir+'beta'+strmid(nm,lnm-3,3)+'.unf',/f77
    readu,unit,mp
    close,unit
    map(*,j)=mp
endfor
spawn, 'rm '+dir+'beta'+strmid(nm,lnm-3,3)+'.unf'
end









