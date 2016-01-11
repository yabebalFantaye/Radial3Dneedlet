maxcount_mocbox=519


forig = '../output/func_r4.unf'
fan = '../output/an_func_r4.unf'
fout = '../output/an2func_func_r4.unf'

nnmax=20 ;;max n for sin(n*pi*r)/r
nshell = 100 ;;sampling on r for original function
nshell_out = nshell ;;radial sampling for reconstructed
fkey = '_nshell'+strn(nshell)+'_nmax'+strn(nnmax)


forig=dblarr(nshell)
ntilda=3
r=findgen(nshell)/float(nshell-1)

.r
for ii=1,nshell-1 do begin
   snr = -4.*r[ii]*(r[ii]-1);sinc_x(!pi*ntilda*r[ii])*!pi*ntilda
   forig[ii]=snr
endfor
end

ball2an,forig,fan,nmax=nnmax,nrsample=1,nshell=nshell,alpha=100
an2ball,fan,fout,nshell=nshell_out,nmax=nnmax

tek_color
axisx = ['xlbl','ylbl']
axisy = ['$r$', '$f(r)=r^3$']

r = (forig)^(1./3.)

;plot, fball              
fname='../figures/test_radial_function_alpha100_reconstruction'+fkey+'.ps'
ps_start_planck,file=fname,/medium,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

plot, r, forig,color=0,xtitle='xlbl',ytitle='ylbl',xra=[0,1],yra=[-0.1,1.1],/xs,/ys
oplot,r, fout,color=2

legend,['Original','Reconstructed'],col=[0,2],textcolor=[0,2],line=[0,0],pspacing=1,/left,/top

ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy

fname='../figures/test_radial_function_alpha100_ratio'+fkey+'.ps'
ps_start_planck,file=fname,/medium,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

plot, r,(forig-fout)/forig,color=0,yra=[-1,1],ytitle='Ratio',xtitle='xlbl',xra=[0,1],/xs,/ys
oplot, r, r*0,col=3

ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy



end
