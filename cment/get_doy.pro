function get_doy, Year, Month, Day, str=str
;+
; Name:         get_doy.pro
; Usage: retval = get_doy(Year, Month, Day, str=str)
;
; Returns day of year, given Year,Month,Day
; If keyword str is set, returns I3.3 string
;
; Input:        Year as 4-digit integer
;               Month as integer
;               Day	(day of month) as integer
;               str     flag:   0 or non-existent = return an integer
;				                1 = return a 3-digit string
; Output:       Return DOY as either an integer or 3-digit string
;
; Calls:        LeapYear.pro
; Called by:    DataOps_eventcb.pro, DateTime__define.pro, SystimeX.pro,
;               AgreeWithTextBoxes.pro, AOval__define.pro,
;               DateTime__AgreeWithTextBoxes.pro, Pmap__PmapPlot.pro
; Written by:   Sue Greer
; Date:	5/10/2004
; Revisions
; 20080703 ahn rename this function get_doy.
;-
	compile_opt idl2

      leap = LeapYear(Year)
      lmo = Month + 2
      lday = (lmo*3055L)/100 + Day - 91
      if lday gt 61 then lday = lday + leap - 2
      IF Year eq 1582 and lday ge 278 then lday = lday - 10
      retval = lday
      if keyword_set(str) then retval=string(format='(I3.3)',lday)
      return, retval
end
