function SystimeX,dtstr=dtstr,year=year,month=month,doy=doy,dom=dom, $
    hour=hour,minute=minute,second=second,date=date,time=time, $
    str=str,fmt=fmt,to_fmt=to_fmt
;+
; Name: SystimeX.pro
; Usage: retval = SystimeX(dtstr=dtstr,year=year,month=month, $
;                   doy=doy,dom=dom,hour=hour,minute=minute, $
;                   second=second,date=date,time=time, $
;                   str=str,fmt=fmt,to_fmt=to_fmt)

; Date manipulation with three different modes of operation:
; conversion, extraction and substitution.
;
; --------------------------------------------------------------
; Keywords:
;   dtstr   : input date/time string in a recognized format (see below)
; The following keywords apply to Extraction mode:
;   year    : return the year only from the date string
;   month   : return the month only (numeric)
;   doy     : return the day of year
;   dom     : return the day of month
;   hour    : return the hour
;   minute  : return the minute
;   second  : return the second (and .sss if available)
;   str     : return the value as a trimmed string
;   subst   : return dtstr with substituted element, e.g. year
; ---------------------------------------------------------------
; NOTE: if date or time, AND multiple element keywords present,
;   str keyword is automatically applied; returned array will be strings.
; ---------------------------------------------------------------
;
;   date    : return the date portion (all values returned as str)
;   time    : return the time portion (all values returned as str)
; ---------------------------------------------------------------
; The following keywords apply to Conversion mode:
;   fmt     : format identifier of the input string (default is 0)
;   to_fmt  : format identifier of out string (requires keyword fmt)
; If converting *to* format 0, day of week is supplied as 'DoW'.
; ---------------------------------------------------------------
; Formats:
;   0 = IDL systime()   :'Mon Jan 05 16:35:22 2004'
;   1 = SECDateTime     :'YYYY-MM-DD hh:mm:ss.sss'
;
; ---------------------------------------------------------------
; Returns:
;   Three return-value modes of operation:
;       (a) Conversion: string (date,time or date/time) is converted
;           from fmt to to_fmt
;       (b) Extraction: element of date/time string, e.g. year, month, etc.
;           Mode (b) can return value as number or string (using str keyword)
;       (c) Substitute single element in dtstr
;
; ---------------------------------------------------------------
; Calls:        DoY.pro
; Called by:    AgreeWithTextBoxes.pro, DataOps_eventcb.pro,
;               MonthLen.pro, PopulateYearDLs.pro,
; Written by:   Sue Greer
; Date: 5/10/2004
;-

    compile_opt idl2

    mm = ['xx','Jan','Feb','Mar','Apr','May','Jun', $
         'Jul','Aug','Sep','Oct','Nov','Dec']
    fmt0='(A3," ",A3," ",I2.2," ",I2.2,":",I2.2,":",I2.2," ",I4)'
    fmt0d = '(A3," ",I2.2," ",I4)'
    fmt0t = '(I2.2,":",I2.2,":",I2.2)'
    fmt1='(I4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",F6.3)'
    fmt1d = '(I4,"-",I2.2,"-",I2.2)'
    fmt1t = '(I2.2,":",I2.2,":",F6.3)'

    ; Check supplied keywords
    kstr = ['year','doy','month','dom','hour','minute','second','date','time']
    k = lonarr(9)
    k[0] = keyword_set(year)
    k[1] = keyword_set(doy)
    k[2] = keyword_set(month)
    k[3] = keyword_set(dom)
    k[4] = keyword_set(hour)
    k[5] = keyword_set(minute)
    k[6] = keyword_set(second)
    k[7] = keyword_set(date)
    k[8] = keyword_set(time)

    kw = where(k eq 1,cnt)          ; Return array
    if cnt eq 1 then kw = kw[0]     ; Return single value

    if not keyword_set(dtstr) then dt = systime() else dt = dtstr
    if not keyword_set(fmt) then fmt = 0
    if not keyword_set(str) then str = 0
    if not keyword_set(subst) then subst = 0

    varr = strarr(9)

;------------------------------------------------------------------------
    ; Extraction of requested elements
;------------------------------------------------------------------------
    case fmt of
       0:    begin  ;'Mon Jan 05 16:35:22 2004'
          if keyword_set(dtstr) then begin
              ; If user supplies fmt=0, and dtstr other than systime,
              ; report error and exit.
              ch3 = strmid(dt,0,3)
              dowStr = ['Mon','Tue','Wed','Thu','Fri','Sat','Sun']
              in = where (ch3 eq dowStr,cnt)
              if cnt eq 0 then begin
                 msg = [dt+' is incompatible with fmt=0', $
                   'Did you mean to use "to_fmt=0"?', $
                   'together with "fmt=<number>"', $
                   'where <number> is the format ',$
                   'identifier of '+dt+'?']
                 res = dialog_message(msg,/error)
                 return,dt
              endif
          endif
          year = fix(strmid(dt,3,4,/reverse_offset))
          varr[0] = strtrim(year,2)
          month = strmid(dt,4,3)
          imon = where(mm eq month,cnt1)
          month = mm[imon[0]]
          varr[2] = string(format='(I2.2)',imon[0])
          dom = fix(strmid(dt,8,2))
          varr[3] = string(format='(I2.2)',dom)
          doy = get_doy(year,imon[0],dom)
          varr[1] = string(format='(I3.3)',doy)
          hour = fix(strmid(dt,11,2))
          varr[4] = string(format='(I2.2)',hour)
          minute = fix(strmid(dt,14,2))
          varr[5] = string(format='(I2.2)',minute)
          second = float(strmid(dt,17,6))
          varr[6] = string(format='(I2.2)',second)
                date = string(format=fmt0d,month,dom,year)
                varr[7] = date
                time = string(format=fmt0t,hour,minute,second)
                varr[8] = time
                if kw gt -1 then varr = varr[kw] else $
                if not k[7] and not k[8] and not str then $
                	varr = fix(varr[0:6])
         end
       1:    begin ;'YYYY-MM-DD hh:mm:ss.sss'
          year = fix(strmid(dt,0,4))
          varr[0] = strtrim(year,2)
          month = fix(strmid(dt,5,2))
          varr[2] = string(format='(I2.2)',month)
          dom = fix(strmid(dt,8,2))
          varr[3] = string(format='(I2.2)',dom)
          doy = get_doy(year,month,dom)
          varr[1] = string(format='(I3.3)',doy)
          hour = fix(strmid(dt,11,2))
          varr[4] = string(format='(I2.2)',hour)
          minute = fix(strmid(dt,14,2))
          varr[5] = string(format='(I2.2)',minute)
          second = float(strmid(dt,17,6))
          varr[6] = string(format='(I2.2)',second)
                date = string(format=fmt1d,year,month,dom)
                varr[7] = date
                time = string(format=fmt1t,hour,minute,second)
                varr[8] = time
                if kw gt -1 then varr = varr[kw] else $
                if not k[7] and not k[8] then varr = fix(varr[0:6])
         end
       else:
    endcase

    ; Handle requests for extraction of a single element
    if cnt eq 1 then begin
       if not str then begin
       case kw of
         ; Numeric output (except date & time)
         0:val = year
         1:val = doy
         2:val = month
         3:val = dom
         4:val = hour
         5:val = minute
         6:val = second
            7:val = date    ;string
            8:val = time    ;string
         else:
       endcase
       endif else begin
         ; String output
       case kw of
         0:val = strtrim(year,2)
         1:val = string(format='(I3.3)',doy)
         2:if fmt eq 0 then val = month else val = string(format='(I2.2)',month)
         3:val = string(format='(I2.2)',dom)
         4:val = string(format='(I2.2)',hour)
         5:val = string(format='(I2.2)',minute)
         6:begin
          if fmt eq 0 then sfmt='(I2.2)' else sfmt = '(F6.3)'
          val = string(format=sfmt,second)
           end
            7:val = date    ;string
            8:val = time    ;string

       endcase
       endelse
       if n_elements(to_fmt) eq 0 then return,val
    endif else begin
        ; Multiple elements requested
       if n_elements(to_fmt) eq 0 then return,varr
    endelse
;-----------------------------------------------------------------------
    ; Conversion of requested date, time, or date/time
;-----------------------------------------------------------------------
    if n_elements(to_fmt) gt 0 then begin
       if fmt eq to_fmt then return,dt     ; No conversion done
       case to_fmt of
         0:   begin ; Output IDL systime-like string
                       ; NOTE:  IDL systime format cannot have milliseconds

              mon = mm[month]
              val = string(format=fmt0,'DoW',mon,dom, $
                 hour,minute,second,year)
                    if not k[7] and not k[8] then return,val
              if keyword_set(date) then $
              vald = string(format=fmt0d,mon,dom,year)
              if keyword_set(time) then $
              valt = string(format=fmt0t,hour,minute,second)
                    if k[7] and not k[8] then return,vald
                    if k[8] and not k[7] then return,valt
                    return,[vald,valt]
          end
         1:   begin
              mon = where(mm eq month,cnt)
              mon = mon[0]
              val = string(format=fmt1,year,mon,dom, $
                 hour,minute,second)
                    if not k[7] and not k[8] then return,val
              if keyword_set(date) then $
              vald = string(format=fmt1d,year,mon,dom)
              if keyword_set(time) then $
              valt = string(format=fmt1t,hour,minute,second)
                    if k[7] and not k[8] then return,vald
                    if k[8] and not k[7] then return,valt
                    return,[vald,valt]
          end
       endcase
    endif
end
