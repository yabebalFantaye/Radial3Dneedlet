pro browse_ballpix,arr,dataout=dataout,rmax=rmax,norm=norm,title=title

if size(arr,/type) eq 7 then begin
    read_pixball,arr,map,maxcount=maxcount
    arr = map
endif 

if keyword_set(norm) then arr= norm*arr

npix = n_elements(arr[*,0])
ntab = n_elements(arr[0,*])

if not keyword_set(rmax) then rmax = 1.
rvec = rmax*findgen(ntab)/float(ntab)

; ----------------------------------
;; define grid for plotting
;----------------------
defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

!P.BACKGROUND = 1               ; white background
!P.COLOR = 0                    ; black foreground
bad_data= !healpix.bad_value
if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)

hist_equal=1 & log=0 & asinh=0
mode_col = keyword_set(hist_equal)
mode_col = mode_col + 2*keyword_set(log) + 4*keyword_set(asinh)
obs_npix = npix
nside = npix2nside(npix)

do_true = keyword_set(truecolors) ;;defines RGB color per pixel
truetype = do_true ? truecolors : 0
du_dv = 2.                      ; aspect ratio
fudge = 1.02               ; spare some space around the Mollweide egg
xsize = 800L
ysize = xsize/2L
size = (do_true) ? 3 : 1
n_uv = xsize*ysize
xll= 0 & xur =  xsize-1
yll= 0 & yur =  ysize-1
xc = 0.5*(xll+xur) & dx = (xur - xc)
yc = 0.5*(yll+yur) & dy = (yur - yc)


; -------------------------------------------------
; converts the position on the sphere into pixel number
; -------------------------------------------------

planmap_ind = indgen(xsize,ysize,/long)
planmap_ind[*,*] = -1
planmap_val =  planmap_ind

count = 0
yband = LONG(5.e5 / FLOAT(xsize))
for ystart = 0, ysize - 1, yband do begin
    yend   = (ystart + yband - 1) < (ysize - 1)
    nband = yend - ystart + 1
    u = FINDGEN(xsize)     # REPLICATE(1,nband)
    v = REPLICATE(1,xsize) # (FINDGEN(nband) + ystart)
    u =  du_dv*(u - xc)/(dx/fudge) ; in [-2,2]*fudge
    v =        (v - yc)/(dy/fudge) ; in [-1,1] * fudge 

    ellipse  = WHERE( (u^2/4. + v^2) LE 1. , nellipse)

    if (nellipse gt 0) then begin
        u1 =  u(ellipse)
        v1 =  v(ellipse)
        u = 0 & v = 0
        s1 =  SQRT( (1-v1)*(1+v1) )
        a1 =  ASIN(v1)

        z = 2./!PI * ( a1 + v1*s1)
        phi = (flipconv *!Pi/2.) * u1/s1 ; lon in [-pi,pi], the minus sign is here to fit astro convention
        sz = SQRT( (1. - z)*(1. + z) )
        vector = [[sz * COS(phi)], [sz * SIN(phi)], [z]]
        u1 = 0 & v1 = 0 & s1 = 0 & a1 = 0 & z = 0 & phi = 0 & sz = 0

        VEC2PIX_RING, nside, vector, id_pix
        planmap_val[count] = id_pix
        planmap_ind[count] = ystart*xsize+ellipse
        count = count + 1
    endif 
endfor

; -------------------------------------------------
; generate the (u,v) position on the mollweide map 
; project the corresponding data value on the map
; ------------------------------------------------

dataout = MAKE_ARRAY(/BYTE, ntab,xsize, ysize,Value = !P.BACKGROUND) ; white
temp = MAKE_ARRAY(/BYTE, xsize, ysize,Value = !P.BACKGROUND) ; white

for nn = 0,ntab-1 do begin
    planmap = temp

    find_min_max_valid, arr[*,nn], mindata, maxdata, valid=Obs, bad_data= 0.9 * bad_data
    print, 'ir, min, max',nn,mindata,maxdata
    data = COLOR_MAP(arr[*,nn], mindata, maxdata, Obs, $
                     color_bar = color_bar, mode=mode_col, $
                     minset = min_set, maxset = max_set, /silent )

    Tmin = mindata & Tmax = maxdata
    
    planmap[planmap_ind[0:count-2]] = data[planmap_val[0:count-2]]
    dataout[nn,*,*] = planmap
endfor

upos = lindgen(xsize)
vpos = lindgen(ysize)

(scope_varfetch('dataset', level=-1, /enter)) = dataout
(scope_varfetch('r', level=-1, /enter)) = rvec
(scope_varfetch('x', level=-1, /enter)) = upos
(scope_varfetch('y', level=-1, /enter)) = vpos

browse, dataout, rvec, upos,vpos,ttitle='r',xtitle='phi',ytitle='theta',data_title=title

end
