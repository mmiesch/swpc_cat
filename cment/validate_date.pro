function validate_date,year,month,day,hour,minute

;20110628 anewman - started

compile_opt idl2

if year lt 1000 or year gt 2999 then return,0
if month lt 1 or month gt 12 then return,0
mlen = MonthLen(month,year4=year)
if day lt 1 or day gt mlen then return,0
if hour lt 0 or hour gt 23 then return,0
if minute lt 0 or minute gt 59 then return,0

return,1
end