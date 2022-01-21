pro check_dates,date_array,start_julian,end_julian,max_interval_in_days, $
                dates_are_fine,problem_string
;20110628 anewman - modified to assure dates are valid dates...

compile_opt idl2

problem_string = 'ok'
dates_are_fine = 1
start_julian = -1
end_julian = -1

;help, date_array

date_array_int = fix(date_array)

;print, 'date_array_int ',date_array_int

start_year = date_array_int[0]
start_month = date_array_int[1]
start_day = date_array_int[2]
start_hour = date_array_int[3]
start_minute = date_array_int[4]
end_year = date_array_int[5]
end_month = date_array_int[6]
end_day = date_array_int[7]
end_hour = date_array_int[8]
end_minute = date_array_int[9]

;validate start date:
valid = validate_date(start_year,start_month,start_day,start_hour,start_minute)
if (valid eq 0) then begin
 dates_are_fine = 0
 problem_string = 'Start date is not a valid date'
endif else begin

   ;validate end date:
   valid = validate_date(end_year,end_month,end_day,end_hour,end_minute)
   if (valid eq 0) then begin
    dates_are_fine = 0
    problem_string = 'End date is not a valid date'
   endif else begin

      start_julian = JULDAY(start_month, start_day, start_year, start_hour,start_minute,0)
      end_julian = JULDAY(end_month, end_day, end_year, end_hour,end_minute,0)

      if end_julian le start_julian then begin
       dates_are_fine = 0
       problem_string = 'End date needs to be after the Start date'
      endif
   endelse
endelse

if dates_are_fine then begin ;check if time interval is within a reasonable limit   
   if (end_julian - start_julian GT max_interval_in_days) then begin
       dates_are_fine = 0
       problem_string = 'The time between start and end cannot exceed ' $
                        + STRTRIM(STRING(max_interval_in_days),2) + ' days.'   
   endif
endif

end