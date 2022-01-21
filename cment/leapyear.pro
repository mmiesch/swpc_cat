function LeapYear,Year
;+
; Name:        LeapYear.pro
; Usage: n = LeapYear(Year)
;
; Returns 1 if Year is a leap year, 0 if not
;
; Input:        Year as 4-digit integer (may be array)
; Output:       1 if Year is a leap year, or 0 otherwise (may be array)
;
; Calls:        None.
; Called by:    function MonthLen
; Written by:   Sue Greer
; Date: 5/10/2004
; Mods: 11/05/2004 Works with array of years (S.Greer).
;       6/28/2011 anewman: use in this program
;-
    compile_opt idl2

    n = n_elements(year)
    leap = intarr(n) & leap[*] = 0
    for i=0,n_elements(year)-1 do begin
        if Year[i] mod 400 eq 0 then begin
            leap[i] = 1
        endif else begin
            if Year[i] mod 4 eq 0 then begin
                 if Year[i] mod 100 ne 0 then leap[i] = 1
            endif
        endelse
    endfor
    return,leap

end