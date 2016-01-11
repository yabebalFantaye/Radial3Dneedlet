;;yabebal fantaye 12 Jan 2013
pro ball2an,map,an,lmax=lmax_in,nmax=nmax_in,cl=cl_in,nshell=nshell_in,$
              mask=mask_in,pol=pol,all=all,w8=w8,large=large_in,alpha=alpha,$
              nproc=nproc_in,rmin=rmin,rmax=rmax,nrsample=nrsample_in,itype=itype,amat=amat

  ;dir='$HEALPIX/src/idl/almball/'
  dir='./'
  large_in=1
  spawn,'mkdir -p large'

  if not keyword_set(alpha) then alpha=3

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

  routine = 'ball2an'
  uroutine = strupcase(routine)
  print, '--------- '+uroutine+' ----------'
  if (n_params() ne 2) then begin
     PRINT, 'Wrong number of arguments in ' + uroutine
     print,'Syntax : '
     print,uroutine+', map, alm, '
     print,'              [LMAX=, NMAX=, CL=, MASK= ,POL= ,ALL=,W8= ] '
;    print,' Type '+uroutine+', /help '
;    print,'   for an extended help'
     return
  endif

  ;;r sampling per shel
  nrsample=0l
  nshell=32
  if (KEYWORD_SET(nrsample_in)) then nrsample=nrsample_in
  if (KEYWORD_SET(nshell_in)) then nshell=nshell_in
  if (KEYWORD_SET(nmax_in)) then nmax=LONG(nmax_in)


  if (nmax lt 1) then begin
     print,'nmax MUST BE GREATER THAN 0!'
     return
  endif

  if not keyword_set(rmin) then rmin=0d0
  if not keyword_set(rmax) then rmax=1d0

  ;;map = dblarr(nshell)
  ;;for ii=0,nshell-1 do begin
  ;;   map[ii] = (ii*1d0/(nshell-1d0))^alpha    ;float(max([ii,0.001]))/float(nshell-1)
  ;;endfor


  openw,unit,dir+'params_ball2an.unf',/f77
  writeu,unit,LONG((ldir eq 'large/'))
  writeu,unit,LONG(nmax)
  writeu,unit,LONG(nshell)
  writeu,unit,double(rmin),double(rmax)
  writeu,unit,LONG(itype)
  close,unit

  openw,unit,dir+'ballmap.unf',/f77
  writeu,unit,double(map)
  close,unit


  if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_ball2an '+direxe

  str='time ./idl_ball2an'
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'
  spawn,'time ./idl_ball2an'
  ;print, stderr


  if size(an,/type) eq 7 then begin
     ;;spawn, 'ls -l '+dir+'an.unf '
     print,'copying an_ball2an_stype'+strn(itype)+'.unf to ',an
     spawn, 'cp '+dir+'an_ball2an_stype'+strn(itype)+'.unf '+ an 
  endif else begin
     print,'reading an_ball2an_stype'+strn(itype)+'.unf'
     an=dblarr(nmax+1)
     openr,unit,dir+'an_ball2an_stype'+strn(itype)+'.unf',/f77
     readu,unit,an
     close,unit
     ;;cl_in=an2cl(an,all=all)
  endelse


     print,'reading anm_check_stype'+strn(itype)+'.unf'
     amat=dblarr(nmax+1,nmax+1)
     openr,unit,dir+'anm_check_stype'+strn(itype)+'.unf',/f77
     readu,unit,amat
     close,unit

end


