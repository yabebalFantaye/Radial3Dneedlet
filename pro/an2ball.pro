;;yabebal fantaye 12 Jan 2013

pro an2ball,an,map,nside=nside_in,beam=beam_in,fwhm=fwhm_in,nopixwin=nopixwin,pol=pol,$
            silent=silent,d1=d1,d2=d2,der1=der1_in,der2=der2_in,large=large_in,rmin=rmin,rmax=rmax,$
            lmax=lmax_in,nmax=nmax_in,nshell=nshell_in,nproc=nproc_in,itype=itype

                                ;dir='$HEALPIX/src/idl/almball/'
  dir='./'
  large_in=1
  spawn,'mkdir -p large'


  if not keyword_set(itype) then itype=2

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
  if (ldir eq 'large/') then dir=dir+ldir

;direxe='$HEALPIX/bin/'
  unit=99
  status = FSTAT(unit)
  if (status.open ne 0) then begin
     WHILE(status.open ne 0) do begin
        unit=unit-1
        status = FSTAT(unit)
     ENDWHILE
  endif

  nproc = 4
  if keyword_set(nproc_in) then nproc=nproc_in


  routine = 'an2ball'
  uroutine = strupcase(routine)
  print, '--------- '+uroutine+' ----------'

  if (n_params() ne 2) then begin
     PRINT, 'Wrong number of arguments in ' + uroutine
     print,'Syntax : '
     print,uroutine+', an, map,'
     print,'              [NSIDE=, BEAM=, FWHM=, NOPIXWIN=, POL= , SILENT] '
;    print,' Type '+uroutine+', /help '
;    print,'   for an extended help'
     return
  endif


  if (KEYWORD_SET(nmax_in)) then nmax=LONG(nmax_in)



  if (nmax lt 1) then begin
     print,'nmax MUST BE GREATER THAN 0!'
     return
  endif

  if not keyword_set(rmin) then rmin=0d0
  if not keyword_set(rmax) then rmax=1d0



  ;;check for nshell 
  nshell=32
  if size(map,/type) ne 7 then nshell=n_elements(map)-1
  if keyword_set(nshell_in) then nshell=LONG(nshell_in)
  print, 'nshell = ',nshell

  ;; load an if a file name is passed
  if size(an,/type) eq 7 then begin

     if not KEYWORD_SET(nmax_in) then begin
        print, 'nmax must be set when passing filename to an.'
        return
     endif    
     fname = an
     an=dblarr(nmax+1)

     openr,unit,fname,/f77
     readu,unit, an
     close,unit
  endif



  ;;write parameter file 
  openw,unit,dir+'params_an2ball.unf',/f77
  writeu,unit,LONG((ldir eq 'large/'))
  writeu,unit,LONG(nmax)
  writeu,unit,LONG(nshell)
  writeu,unit,LONG(itype)  
  writeu,unit,double(rmin),double(rmax)
  close,unit

  openw,unit,dir+'an_ball2an_stype'+strn(itype)+'.unf',/f77
  writeu,unit,double(an)
  close,unit


  if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_an2ball '+direxe

  str='cd '+direxe+'; time ./idl_an2ball'
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'
  SPAWN,str


  ;; if size(map, /type) eq 7 then begin
  ;;    spawn, 'cp '+dir+'ballmap.unf '+map
  ;; endif else begin
  
  openr,unit,dir+'ballmap_an2ball_stype'+strn(itype)+'.unf',/f77
  map=dblarr(nshell+1)
  readu,unit,map
  close,unit
;;  endelse



end





