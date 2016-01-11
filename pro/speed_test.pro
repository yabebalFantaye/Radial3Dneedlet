
mmetric=dblarr(10,5)   
metric=dblarr(10,5)   

.r
nsim=20
for ii=1,nsim do begin
   fname='../speed_measure/metric_almn2ball_ball2almn_r7_11_ns5_9_sim'+strn(ii)+'.unf'
   openr, 12, fname,/f77
   readu, 12,metric
   close, 12
   mmetric=mmetric+metric/nsim
endfor
end

mmetric4=dblarr(4,5)   
metric4=dblarr(4,5)   

.r
nsim=20
for ii=1,nsim do begin
   fname='../speed_measure/metric_almn2ball_ball2almn_r7_8_ns5_6_sim'+strn(ii)+'.unf'
   openr, 12, fname,/f77
   readu, 12,metric4
   ;;print, 'metric4 isim='+strn(ii),metric4
   close, 12
   mmetric4=mmetric4+metric4/nsim
endfor
end


mt4 = mmetric4[0:1,*]
mr4 = mmetric4[2:3,*] ;;

mt = mmetric[0:4,*]
mr = mmetric[5:9,*] ;;

mt[0,*]=mt4[0,*]

print, 'mt = ',mt
print, 'mr = ',mr

print, 'mt4 = ',mt4
print, 'mr4 = ',mr4


;; ind=sort(mt[*,0])
;; mt=mt[ind,*]

;; ind=sort(mr[*,0])
;; mr=mr[ind,*]


;; mr[0,3]=9.10
;; mr[0,3]=1.10
;; mr[1,3]=4.20

;metric=transpose(metric)
;;0-npix 1-lmax 2-nmax 3-time 4-maxdiff

wunf,mt,'../figsOct14/cpu_almn2ball_ball2almn_mt_55_20imc.unf'
wunf,mr,'../figsOct14/cpu_almn2ball_ball2almn_mr_55_20imc.unf'

LS_circle,sz=0.5,color=2
planck_plot, mr[*,0],mr[*,3],xra=[1e6,1e9],yra=[0,210],xover=mt[*,0],yover=mt[*,3],ax=['xttl','yttl'],ay=['$N_{pix}$','$t_s$[sec]'],$
             fname='../figsOct14/cpu_almn2ball_ball2almn.ps',/psend,color=[0,2],psym=[-8,-8] ,/xlog,/xstyle
planck_plot, mr[*,0],mr[*,4],xra=[1e6,1e9],yra=[1e-6,3e-5],xover=mt[*,0],yover=mt[*,4],ax=['xttl','yttl'],ay=['$N_{pix}$','$max(|a_{\ell mn}-a_{\ell mn}^{rec}|)$'],$
             fname='../figsOct14/almn_diff_almn2ball_ball2almn.ps',/psend,color=[0,2],psym=[-8,-8] ,/xlog,/xstyle



almnsim = dcomplexarr(256,129,129)
openr, 12, 'temp/almn_sim.unf',/f77
readu, 12, almnsim
close, 12

almnrec=almnsim
openr, 12, 'temp/almn_rec.unf',/f77
readu, 12, almnrec
close, 12


print, 'max sim: ',max(reform(abs(almnsim[3,*,*])))
print, 'max rec: ',max(reform(abs(almnrec[3,*,*])))
print, 'max diff n=2',max(reform(abs(almnsim[2,*,*]-almnrec[2,*,*])))
print, 'max_diff n=4',max(reform(abs(almnsim[4,*,*]-almnrec[4,*,*])))




;cl_file='../cls/Cl_model_bin5.dat'
;readcol,cl_file, ll,clmm,clmm_err
;cl=clmm*1e5

lmax=64
nside=128
cl=1d0+findgen(129)
cl[0]=0
cl[1]=0
create_alm,cl,alm,lmax=lmax,12312
alm2map,alm,map,/nopix,nside=nside
map2alm,map,alm2,lmax=lmax,/w8
print, reform(abs(alm2[0,*,0]-alm[0,*,0])) 
print, reform(abs(alm2[0,*,2]-alm[0,*,2])) 

print, 'max sim: ',max(reform(abs(alm[0,*,*])))
print, 'max rec: ',max(reform(abs(alm2[0,*,*])))

print, 'max diff n=2',max(reform(abs(alm[0,*,*]-alm2[0,*,*])))
print, 'max_diff n=4',max(reform(abs(alm[0,*,*]-alm2[0,*,*])))



end
