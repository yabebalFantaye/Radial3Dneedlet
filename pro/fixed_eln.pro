pro fixed_eln,lmax,nmax,eln_out=eln_out,count=count, lind=lind,nind=nind

if not keyword_set(eln_out) then eln_out=long(3*3+4*(4+1))

lind=make_array(long(lmax+1)*long(nmax+1),/long)
nind=lind


count=0
for l=0,lmax do begin
   for n=0,nmax do begin
      
      eln = long(n*n+l*(l+1))
      
      if eln eq eln_out  then begin
         lind[count]=l
         nind[count]=n
         count=count+1
      endif
      
   endfor
endfor


print, 'nu,ber indices of l,n to have equal e_ln= '+strn(eln_out)+' are '+strn(count)
lind = lind[0:count-1]
nind = nind[0:count-1]

end
