pro read_pixball,filename,arr,maxcount=maxcount,npix=npix,ntab=ntab

;;read ball format data
unit=get_funit()
openr,unit,filename,/f77
npix = lonarr(1)
ntab=npix
maxcount=npix

readu,unit,npix,ntab,prec,maxcount
npix = npix[0] & ntab=ntab[0] & maxcount=maxcount[0]
arr=fltarr(npix,ntab)
readu,unit,arr
close,unit

end
