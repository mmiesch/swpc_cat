function MonthLen,month,year4=year4
;+
; Name:     MonthLen.pro
; Usage: n = MonthLen(month,year4=year4)
;
; Returns the number of days in the month.
;
; Input:        month       Month as int (1..12) or string.
;               year4       Optional 4-digit year as int or string.
;                           (If not supplied, current year is used).
; Output:       returns number of days in the month
;
; Calls:        SystimeX.pro,LeapYear.pro
; Called by:    function validate_date
; Written by:   Sue Greer
; Date: 5/10/2004
; Mods: 6/28/2011 anewman: modified for this program
;
; If year is omitted, current year is used.
; (Year must be 4 digits or error is returned.)
; If input(s) are supplied as strings, they
; are converted to numbers.
;-
  compile_opt idl2

  if not keyword_set(year4) then begin
    year4 = SystimeX(/year,/str)
  endif
  year4 = fix(year4)
  if size(month,/type) eq 7 then begin
    month = strmid(month,0,3)
    mstr = ["xx","Jan","Feb","Mar","Apr","May","Jun", $
        "Jul","Aug","Sep","Oct","Nov","Dec"]
    mnum = where(month eq mstr,cnt)
    month = mnum
  endif
  mlen = [0,31,28,31,30,31,30,31,31,30,31,30,31]
  if LeapYear(year4) then mlen[2] = 29
  return,mlen[month]
end