

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006

pro get_3dneedlet,bb,j0,nj,lmax,nmax,gln,wavtyp=wavtyp_in,mexican_p=p_in,modif=modif,softmod=softmod,nside=nside_in,smhw=smhw_in

dir= './temp/' ;;'$HEALPIX/src/idl/almtools/'
spawn, 'mkdir -p '+dir


unit=99
status = FSTAT(unit)
if (status.open ne 0) then begin
WHILE(status.open ne 0) do begin
unit=unit-1
status = FSTAT(unit)
ENDWHILE
endif

routine = 'get_needlet'
uroutine = strupcase(routine)
if (n_params() ne 6) then begin
    PRINT, 'Wrong number of arguments in ' + uroutine
    print,'Syntax : '
    print,uroutine+',B,j0,nj,lmax,nmax,gl'
    print,'              [wavtyp=, mexican_p=,nside=] '
;    print,' Type '+uroutine+', /help '
;    print,'   for an extended help'
    return
endif

if (KEYWORD_SET(wavtyp_in)) then begin
if ((wavtyp_in ne 'standard') and (wavtyp_in ne 'mexican') and (wavtyp_in ne 'smhw')) then begin
print,'wavtyp must be standard or mexican'
return
endif
wav=0l
if (wavtyp_in eq 'mexican') then wav=1l
if (wavtyp_in eq 'smhw') then wav=2l
endif else begin
wav=0l
endelse

if ((KEYWORD_SET(p_in)) and (wav ne 1)) then begin
print,'mexican_p can only be set when wavelet type is mexican!'
return
endif

p=1l
if (KEYWORD_SET(p_in) and (wav eq 1)) then p=LONG(p_in)

;; if (j0 lt 0) then begin
;; print,'j0 must be larger than 0!'
;; return
;; endif
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

if (((KEYWORD_SET(nside_in)) or (KEYWORD_SET(smhw_in))) and (wav ne 2)) then begin
print,'nside and smhw should be set only when using SMHW!'
return
endif


if (wav eq 2) then begin
if ((not KEYWORD_SET(nside_in)) or (not KEYWORD_SET(smhw_in))) then begin
print,'nside and smhw must be set when using SMHW!'
return
endif

openw,unit,dir+'smhw.unf',/f77
writeu,unit,DOUBLE(smhw_in)
close,unit
openw,unit,dir+'nside.unf',/f77
writeu,unit,LONG(nside_in)
close,unit

endif


openw,unit,dir+'lmax.unf',/f77
writeu,unit,LONG(lmax)
close,unit

openw,unit,dir+'nmax.unf',/f77
writeu,unit,LONG(nmax)
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


SPAWN,'./idl_get_3dneedlet'
gln=dblarr(lmax+1l,nmax+1l,nj)
openr,unit,dir+'gln.unf',/f77
readu,unit,gln
close,unit

;; if (KEYWORD_SET(modif) and (NOT KEYWORD_SET(softmod))) then begin
;; for l=0,lmax do begin
;; norm=total(gl(l,*)^2)
;; if (norm gt 0d0) then gln(l,*)=gl(l,*)/sqrt(norm)
;; endfor
;; endif

;; if (KEYWORD_SET(softmod)) then begin
;; norm=total(gl(2:lmax/8.,*)^2)/(DOUBLE(lmax)/8d0-2d0+1d0)
;; gl=gl/sqrt(norm)
;; endif


end









