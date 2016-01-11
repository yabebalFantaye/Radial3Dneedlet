maxcount_mocbox=1 ;519

nrint=64
nside=16
nshell=2*nrint ;;number of radial shells 
nnmax=30 ;;maximum radial multipole
lmax=5 ;;angular multipole
nproc=1 
nrsample=1 ;;r sampling factor for ball2almn
figdir='../figures_nshell'+strn(nshell)+'_nnmax'+strn(nnmax)+'_nside'+strn(nside)+'_lmax'+strn(lmax)+'/'
spawn, 'mkdir -p '+figdir

fball_orig = '../output/ballmap_mockbox_almn2ball_nrint'+strn(nrint)+'_2.unf'
fball = '../output/ballmap_mockbox_almn2ball_nrint'+strn(nrint)+'_3.unf'
falmn = '../output/almn_mockbox_2.unf'

;; fball_orig='../output/ballmap_mockbox_nrint'+strn(nrint)+'.unf'
;; fball = '../output/ballmap_mockbox_almn2ball_nrint'+strn(nrint)+'.unf'
;; falmn = '../output/almn_mockbox.unf'

fkey = '_nshell'+strn(nshell)+'_nmax'+strn(nnmax)

almr=1
ralmr=1
outalmn=1
ball2almn,fball_orig,falmn,lmax=lmax,nmax=nnmax,/large,nproc=nproc,nrsample=nrsample,almr=almr,outalmn=outalmn,rmin=0d0,rmax=1d0
almn2ball,falmn,fball,nside=nside,nshell=nshell,lmax=lmax,nmax=nnmax,/large,nproc=nproc,almr=ralmr,rmin=0d0,rmax=1d0

;return

x = indgen(nside2npix(nside))+1
y = indgen(nshell+1)

;;save 2D projected density maps

axisx = ['xlbl','ylbl','zlbl']
axisy = ['Healpix pixel number', 'Radial pixel number','$\Delta N/N [\%]$']
CTFILE = 'Planck_CT.tbl'

;fball_orig[*,0]=0.01
;fball[*,0]=0.01

;; ;;--------------------------
;; fname=figdir+'test_original_ball'+fkey+'.ps'
;; ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

;; hfi_ct, CTDIR='./', CTFILE=CTFILE, /LOAD

;; density, fball_orig, x, y, /colorbar, title='Original 3D ball: from VLS simulation Ngal = 1841599',$
;;              bar_charsize=2, xtitle='xlbl',ytitle='ylbl',xra=[1,nside2npix(nside)],yra=[1,nshell],/xs,/ys

;; ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy

;; ;;--------------------------
;; fname=figdir+'test_reconstructed_ball'+fkey+'.ps'
;; ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

;; hfi_ct, CTDIR='./', CTFILE=CTFILE, /LOAD

;; density, fball, x, y, /colorbar, title='Reconstructed 3D ball',$
;;              bar_charsize=2, xtitle='xlbl',ytitle='ylbl',xra=[1,nside2npix(nside)],yra=[1,nshell],/xs,/ys

;; ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy

;; ;;--------------------------
;; fname=figdir+'test_difference_ball'+fkey+'.ps'
;; ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

;; hfi_ct, CTDIR='./', CTFILE=CTFILE, /LOAD

;; density, (fball_orig - fball), x, y, zmin=-2, zmax=2, /colorbar, title='Difference  3D ball',$
;;              bar_charsize=2, xtitle='xlbl',ytitle='ylbl',xra=[1,nside2npix(nside)],yra=[1,nshell],/xs,/ys

;; ;*100./(fball_orig)

;; ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy
;; ;;---------------------------

fname=figdir+'frac_diff_raidal_orig_vs_reconstruct_ball'+fkey+'.ps'
ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

tek_color

plot, y,(fball_orig[10,*]-fball[10,*])*100./fball_orig[10,*], title=' ',color=0,$
             xtitle='ylbl',ytitle='zlbl',xra=[0,nshell],yra=[-1,1],/xs,/ys

oplot, y,(fball_orig[100,*]-fball[100,*])*100./fball_orig[100,*],color=2
oplot, y,(fball_orig[500,*]-fball[500,*])*100./fball_orig[500,*],color=3
oplot, y,(fball_orig[1000,*]-fball[1000,*])*100./fball_orig[1000,*],color=4
oplot, y,(fball_orig[2000,*]-fball[2000,*])*100./fball_orig[2000,*],color=5

;*100./(fball_orig)
;; legend,['Original','Reconstructed'],color=[0,2],textcolor=[0,2],/top,/right,box=0,$                                                                                                                                                                
;;        line=[0,0],pspacing=1,thick=4,charsize=lcharsize

legend,['pixel=10','pixel=100','pixel=500','pixel=1000','pixel=2000'],color=[0,2,3,4,5],textcolor=[0,2,3,4,5],/top,/left,box=0,$                                                                                                                                                                
       line=[0,0,0,0,0],pspacing=1,thick=4,charsize=lcharsize

ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy
;;--------------------------


;; browse_ballpix,fball
;; browse_ballpix,fball_orig

;; fball_diff = (abs(fball-fball_orig))
;; print, 'mean diff = ',mean(fball_diff)

;; browse_ballpix,fball_diff

end


;; density, z, x, y, zmin=zmin, zmax=zmax, zlog=zlog, colorbar=colorbar, $
;;              nlevels=nlevels, fill=fill, bar_charsize=bar_charsize, $
;;              bar_pad=bar_pad, bar_title=bar_title, bar_format=bar_format, scientific=scientific, $
;;              _extra=extra
