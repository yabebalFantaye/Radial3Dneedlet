lmax=200
nmax=200
j0=0
bb=2d0
nj=10

get_3dneedlet,bb,j0,nj,lmax,nmax,gln

wunf, double(gln), '../python/gln_data.unf'


ax=['xtt',' ytt','zttttt']
ay=['$\ell$','$n$','$b(e_{\ell n}/B^j)$']

 LoadCT, 13                     ;25, /Brewer, /Reverse
.r
for j=7,7 do begin

   data=gln[*,*,j]

   fname = '../figures/gln_3d_lmax'+strn(lmax)+'_nmax'+strn(nmax)+'_j'+strn(j)+'.ps'
   ps_start_planck,file=fname,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xs,ysize=ys,font=18 ;,medium=medium,large=large ;, /large ;,/medium ;large

   surface, data,xtitle=ax[0],ytitle=ax[1],ztitle=ax[2],$
            charsize=0.6,xmargin=[20,16],ymargin=[16,16],zmargin=[16,16], POSITION = [0.2, 0.2, 0.9, 0.9] ;,Shades=BytScl(data)

;;,xcharsize=0.6,ycharsize=0.7,zcharsize=0.7

   ps_end_planck, /png,feps=fname,/latex,xsize=xs,ysize=ys ,axisx=ax,axisy=ay ;,/nops

   ;; fname = '../figures/gln_3d_lmax'+strn(lmax)+'_nmax'+strn(nmax)+'_j'+strn(j)+'.ps'
   

   ;; cgLoadCT, 25, /Brewer, /Reverse
   ;; cgSurf, data, Shades=BytScl(data), /NoErase, xtitle=ayy[0],ytitle=ayy[1],ztitle=ayy[2],output=fname

   ;;  ps_to_any, /png,feps=fname,/latex,axisx=ax,axisy=ay ;,/nops
endfor
end
