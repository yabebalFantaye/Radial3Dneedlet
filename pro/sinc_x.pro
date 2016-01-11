;+
; NAME:
;
;         SINC
;
; PURPOSE:
;
;         Function to return the value of the SINC function,
;         i.e., sin(x)/x.
;
; CALLING SEQUENCE:
;
;         Result = SINC(X)
; 
; INPUTS:
;
;         X - Input value.  Scalar or array.
;
; OUTPUTS:
;
;         Result - Value of SIN(X)/X. 
;
; PROCEDURE:
;
;         Straightforward; except Result is explicitly set to
;         one when X=0.
;
; MODIFICATION HISTORY:
;
;         David L. Windt, Bell Laboratories, May 1997
;         windt@bell-labs.com
;
;-
function sinc_x,x

wh=where(x eq 0,count)
if count ne 0 then x(wh)=1.
result=sin(x)/x
if count ne 0 then result(wh)=1.
return,result
end
