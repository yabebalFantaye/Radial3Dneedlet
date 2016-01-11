function get_funit

unit=99
status = FSTAT(unit)
if (status.open ne 0) then begin
    WHILE(status.open ne 0) do begin
        unit=unit-1
        status = FSTAT(unit)
    ENDWHILE
endif

return, unit
end
