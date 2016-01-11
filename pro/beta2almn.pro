

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006

pro beta2alm_ball,beta,alm,lmax=lmax_in,pol=pol,all=all,w8=w8,large=large_in,nproc=nproc_in,mask=mask_in

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

routine = 'beta2alm'
uroutine = strupcase(routine)
if (n_params() ne 2) then begin
    PRINT, 'Wrong number of arguments in ' + uroutine
    print,'Syntax : '
    print,uroutine+', beta, alm,'
    print,'              [large=, LMAX=, POL= ,ALL=,W8=,nproc=,mask= ] '
;    print,' Type '+uroutine+', /help '
;    print,'   for an extended help'
    return
endif

polar=0l
;if (KEYWORD_SET(pol)) then polar=1l
if (KEYWORD_SET(pol)) then print,'polarization not implemented!'


smap=size(beta)
smap=LONG(smap[0])
if ((smap lt 2) or (smap gt 3)) then begin
print,'beta MUST BE of the form (npix,pol,nj) OR (npix,nj)'
return
endif

if (smap eq 2) then nj=LONG(n_elements(reform(beta(0,*))))
if (smap eq 3) then nj=LONG(n_elements(reform(beta(0,0,*))))

npix=n_elements(beta(*,0))
nside=npix2nside(npix)

mask=fltarr(npix)
mask(*)=1.
if (KEYWORD_SET(mask_in)) then begin
npix2=n_elements(mask_in)
npix2=npix2[0]
if (npix2 ne npix) then begin
print,'mask must have same nside as beta!'
return
endif
mask=mask_in
endif

if (KEYWORD_SET(lmax_in)) then begin
lmax=LONG(lmax_in)
endif else begin
lmax=3l*nside
endelse

if (lmax lt 1) then begin
print,'lmax MUST BE GREATER THAN 0!'
return
endif
if (lmax gt 3*nside) then begin
print,'WARNING:lmax>3*nside will use lmax=3*nside'
lmax=3*nside
endif

if (KEYWORD_SET(nproc_in)) then begin
if (nproc_in gt nj) then begin
print,'nproc must be smaller or equal to nj!'
return
endif
endif


if (KEYWORD_SET(w8)) then begin
wfile='weight_ring_n'
wnumfile=trim('00000'+trim(string(LONG(nside)),2),2)
wsize=strlen(wnumfile)
wfile='$HEALPIX/data/'+wfile+strmid(wnumfile,wsize-5,wsize)+'.fits'
read_fits_map,wfile,wl
wl=wl+1d0
endif else begin
wl=dblarr(2*nside,3)
wl(*,*)=1d0
endelse
wunf,wl,dir+'w8_beta.unf'


openw,unit,dir+'nside_beta.unf',/f77
writeu,unit,LONG(nside)
close,unit
openw,unit,dir+'lmax_beta.unf',/f77
writeu,unit,LONG(lmax)
close,unit
for j=0l,nj-1l do begin
nm='00'+trim(string(j),2)
lnm=strlen(nm)
if (ldir eq 'large/') then openw,unit,dir+'large/beta'+strmid(nm,lnm-3,3)+'.unf',/f77 else  openw,unit,dir+'beta'+strmid(nm,lnm-3,3)+'.unf',/f77
if (smap eq 2) then map=reform(float(beta(*,j)))
if (smap eq 3) then map=float(reform(beta(*,0,j)))
map=map*mask
writeu,unit,map
close,unit
endfor
;openw,unit,dir+'polar.unf',/f77
;writeu,unit,polar
;close,unit
openw,unit,dir+'large.unf',/f77
if (ldir eq 'large/') then writeu,unit,1l else writeu,unit,0l
close,unit
openw,unit,dir+'nj.unf',/f77
writeu,unit,nj
close,unit

if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_beta2alm_par '+direxe

if (not KEYWORD_SET(nproc_in)) then nproc_in=1
nodelist=direxe+'nodelist'
str='cd '+direxe+' ; mpiexec -machinefile '+nodelist+' -n '+trim(string(nproc_in),2)+' '+direxe+'idl_beta2alm_par'
SPAWN,str

alm=complexarr(lmax+1,lmax+1,nj)
alm0=complexarr(lmax+1,lmax+1)
for j=0l,nj-1l do begin
nm='00'+trim(string(j),2)
lnm=strlen(nm)
if (ldir eq 'large/') then openr,unit,dir+'large/betaalm'+strmid(nm,lnm-3,3)+'.unf',/f77 else openr,unit,dir+'betaalm'+strmid(nm,lnm-3,3)+'.unf',/f77
readu,unit,alm0
alm(*,*,j)=alm0
close,unit
endfor

end




