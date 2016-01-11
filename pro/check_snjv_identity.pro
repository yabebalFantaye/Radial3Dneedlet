function snjvv, nn, v, j,bb ;;,rvbj=rvbj 

  if bb^(double(j)-1d0) lt nn then begin
     xx=0d0
  endif else begin
     
     rvj=double(v)/bb^double(j)
     ;;if keyword_set(rvbj) then 
     
     xx=sqrt(2d0)*!pi*nn*sinc_x(!pi*nn*rvj)
  endelse
  
  return, xx
end



function check_snjv_identity,jz=jz, jj=jj,bb=bb,nvmaps=nvmaps,nnmax=nnmax,vvec=vvec


  if not keyword_set(nnmax) then nnmax=20
  if not keyword_set(j0) then j0=1d0
  if not keyword_set(jj) then jj=1d0
  if not keyword_set(bb) then bb=2d0
  if not keyword_set(nvmaps) then nvmaps=128
  
  nvmaps=round(bb^double(jj))
;  nnmax = nvmaps-1
  mat=dblarr(nnmax, nnmax)

  for in=0,nnmax do begin
     for in2=0,nnmax do begin

        ssprod=0d0
        vvec=dblarr(nvmaps)
        count=0
        if nnmax lt nvmaps then begin
           for iv=0,nvmaps-1 do begin
              vval= double(iv)*(bb^double(jj)-1d0)/double(nvmaps-1d0)
              vvec[iv]=vval     ;/bb^double(jj)
              rv2=(double(vval)/bb^double(jj))^2
              
              arg1=snjvv(in,vval,jj,bb)
              arg2=snjvv(in2,vval,jj,bb)
              ssprod=ssprod + (rv2/bb^double(jj))*arg1*arg2
              count=count+1
           endfor
        endif
                                ;;if in eq 1 and in2 eq 1 then print, vvec
        mat(in-1, in2-1) = ssprod ;/(bb^double(jj)*double(count))

     endfor
  endfor

  return, mat

end
