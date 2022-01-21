PRO find_images_by_date, ddate,info2,input_image_list
;20110701 anewman - start
compile_opt idl2

;get the name of the folder/subdir named for the date:
;e.g. folder 20110701
year = STRTRIM(STRING(ddate[0]),2 )
month = STRTRIM(STRING(ddate[1]),2 )
if STRLEN(month) eq 1 then month = '0' + month
day = STRTRIM(STRING(ddate[2]),2 )
if STRLEN(day) eq 1 then day = '0' + day
date_folder = year + month + day

;associate the telescope name with the folder name where images are stored:
;assumes a telescope is specified in same location in 
;info2.telescope_array as in info2.image_in_folder_array
res = STRCMP(info2.which_telescope,info2.telescope_array)
ele = where(res,count)

if count eq 1 then begin
  info2.telescope_folder = info2.image_in_folder_array[ele] 
endif else begin
  info2.telescope_folder = ''
  msg1 = 'ERROR: folder for telescope ' + info2.which_telescope +' is not defined.'
  msg2 = 'There should be matching entries for each telescope '
  msg3 = 'in lines "telescopes" and "image_in_folders" in file' 
  msg4 = info2.input_file + '.'
  msg5 = ' '
  msg6 = 'Please contact the SWPC Enlil software support team...'
  ok = DIALOG_MESSAGE([msg1,msg2,msg3,msg4,msg5,msg6],/ERROR,/center)
  EXIT   
endelse  

;assemble the location to look for files:
image_location = info2.image_in_root + info2.sep + date_folder + info2.sep + info2.telescope_folder

result = file_test(image_location,/directory)
if result eq 0 then begin
  msg1 = 'Warning: The subdir ' + image_location + ' does not exist;'
  msg2 = 'this means there are no images available for telescope ' + info2.which_telescope
  msg3 = 'on date ' + date_folder  + ' [yet].'
  ok = DIALOG_MESSAGE([msg1,msg2,msg3],/information,/center)
  RETURN
  
endif

;get the files from this folder:
;  ( do not get the STEREO polarization images, for now...)
if (STRCMP( info2.which_telescope, 'STEREO', 6 )) then begin
   names_of_all_files_in_folder = file_search(image_location,'*d7*.fts')
endif else begin
   names_of_all_files_in_folder = file_search(image_location,'*.fts')
endelse
number_of_files_in_folder = size(names_of_all_files_in_folder)
number_of_files_in_folder = number_of_files_in_folder[1]
if number_of_files_in_folder lt 1 then begin
   msg1 = 'Warning...No FITS files in this folder: ' + image_location
   ok = DIALOG_MESSAGE(msg1,/information,/center)
   RETURN
   
endif else begin
   ;Keep only the files between the user-selected 
   ;start and end date interval, inclusive
   
   number_of_files = n_elements(names_of_all_files_in_folder)
   
   progressbar = Obj_New('progressbar', Color='Charcoal', Text='Checking Image files',/nocancel, background = 'White', $
                    xsize = 500 , ysize = 20 , title='Getting Images', group_leader = info2.tlb2)
   progressBar -> Start
   
;   help,names_of_all_files_in_folder
   datetime_strings = STRMID(FILE_BASENAME(names_of_all_files_in_folder,'.fts'),0,15)
   i=0
   while (i lt n_elements(names_of_all_files_in_folder)) do begin
   
      data_year = FIX(STRMID(datetime_strings[i],0,4))
      data_month = FIX(STRMID(datetime_strings[i],4,2))
      data_day = FIX(STRMID(datetime_strings[i],6,2))
      data_hour = FIX(STRMID(datetime_strings[i],9,2))
      data_minute = FIX(STRMID(datetime_strings[i],11,2))
      data_second = FIX(STRMID(datetime_strings[i],13,2))
      data_dt = JULDAY(data_month,data_day,data_year,data_hour,data_minute,data_second)
      
      progress_percent = 100.*(float(i+1)/float(number_of_files))
      
      progressBar -> Update, progress_percent, Text='Checking Image files for ' + STRMID(datetime_strings[i],0,4)+ '-' + STRMID(datetime_strings[i],4,2)+'-' + STRMID(datetime_strings[i],6,2)

      ;Add to the list of images that are between start and end dates, inclusive
      if ((data_dt GE info2.start_julian) && (data_dt LE info2.end_julian)) then begin
         good_image = check_image(names_of_all_files_in_folder[i])
         ;print,'good_image=',good_image
         if good_image then $
            input_image_list.Add,names_of_all_files_in_folder[i],/NO_COPY
      endif
      i=i+1
   endwhile
   
   progressBar -> Destroy
     
endelse

END

;****


PRO get_fits_header_data, fits_header, which_telescope, telescope_folder, scaling_factor, $
                            rotation, SunX, SunY, datetime_string, rsun, pixel_scale
  
  compile_opt idl2
  
  date_string = FXPAR(fits_header,'DATE-OBS')
  
  
  ;The FITS header for LASCO and STEREO differ; STEREO does not have a TIME-OBS item
  if (STRMID(which_telescope,0,5) eq 'LASCO') then begin
    ;use dashes, not slashes, between YYYY MM DD
    date_string = strjoin(strsplit(date_string,'/',/EXTRACT),'-')
    time_string = FXPAR(fits_header,'TIME-OBS')
    time_string = strmid(time_string,0,8) ;strip off msec
    datetime_string = date_string + 'T' + time_string ;assemble datetime string
    naxis1 = fix(sxpar(fits_header, 'NAXIS1'))
    ;thanks to Doug Biesecker for the following:
    rotation = float(sxpar(fits_header, 'CROTA1')); specified in degrees CW
    pixel_scale = float(sxpar(fits_header, 'CDELT1')) ; assumes CDELT1 and CDELT2 are identical
    binning1 = fix(sxpar(fits_header, 'SUMROW')); 0 signifies no binning
    binning2 = fix(sxpar(fits_header, 'LEBXSUM')); 1 signifies no binning
    if binning2 mod 2 eq 1 then binning2 = binning2 - 1
    binning = max([binning1, binning2]); binning isn't needed so long as CDELT1 is used
    sun_x = float(sxpar(fits_header, 'CRPIX1')); confirm this is x-value
    sun_y = float(sxpar(fits_header, 'CRPIX2')); confirm this is y-value
    rsun = 941.; solar radius from L1 in arcseconds, a default value
    
  endif else if (STRMID(which_telescope,0,6) eq 'STEREO') then begin
  
    ;time is stored in DATE-OBS, for STEREO FITS files
    ;want this format: YYYY-MM-DDTHH:MM:SS
    time_string = strmid(date_string,11,8) ;strip off msec
    date_string = strmid(date_string,0,10) ;strip off time
    datetime_string = date_string + 'T' + time_string ;reassemble datetime string
    naxis1 = fix(sxpar(fits_header, 'NAXIS1'))
    ;thanks to Doug Biesecker for the following:
    rotation = float(sxpar(fits_header, 'CROTA'))
    pixel_scale = float(sxpar(fits_header, 'CDELT1'))
    binning = fix(sxpar(fits_header, 'SUMMED')); dimension = original/(2^(SUMMED-1)); 1 signifies no binning
    ; NOTE: binning isn't needed since CDELT1 is updated properly
    sun_x = float(sxpar(fits_header, 'CRPIX1'))
    sun_y = float(sxpar(fits_header, 'CRPIX2'))
    rsun = float(sxpar(fits_header, 'RSUN')); The fits files contain the solar radius in arcsec
    
  endif
  
  ;The input FITS files we get from NOC do not have 512x512 images, but
  ; we convert each image to 512x512 (see rebin command elsewhere in this
  ; program). We need to keep track of the
  ; scaling factor for the image conversion,
  ; for use in rotation and leading edge calculations
  
  
  scaling_factor = 1.
  if (naxis1 eq 256) then scaling_factor = 2. $         ;256x256 to 512x512: multiply by 2
  else if (naxis1 eq 1024) then scaling_factor = 0.5   ;1024x1024 to 512x512: divide by 2 (i.e. multiply by 0.5)
  
  if telescope_folder eq 'LASCO_C2_MERGE' or telescope_folder eq 'LASCO_C3_MERGE' then begin
    SWPC_rebin_version1_factor = 0.5
  endif else begin
    SWPC_rebin_version1_factor = 1.0
  endelse
  

sunX = sun_x * SWPC_rebin_version1_factor
sunY = sun_y * SWPC_rebin_version1_factor
    
    ; datetime_string looks like 'yyyy-mm-ddThh:mm:ss', with a 'T' between date and time
    ; this will be the datetime string stored in the output named-CME FITS file,
    ; if the user selects this image when "Defining CME Name"
cme_datetime = datetime_string
rsun = rsun
    
pixel_scale = pixel_scale / SWPC_rebin_version1_factor
  
END





PRO reset_for_images,info2,sensitive_val,num_images
; sensitive_val    0 or 1  when new images have been obtained, make some widgets sensitive  
;                      but when a new telescope has been selected, set those widgets insensitive
; num_images       number of images obtained, based on the datetime interval
;                  This number will be 0, if a new telescope has been selected.


compile_opt idl2

info2.i_representative_image_has_been_chosen = 0
reset_CME_name, info2
reset_CME_image, info2
info2.clock_angle_degrees = 0
clock_angle_string_value = ' 0'
info2.clock_angle_string->SetProperty,strings=clock_angle_string_value
info2.clock_angle_deg_object ->Setproperty, strings = clock_angle_string_value
info2.rotate_x = 0.
info2.rotate_y = 0.

info2.leading_edge_in_rs = 0
lead_edge_string = ' 0'
info2.leading_edge_string->SetProperty,strings=lead_edge_string
info2.leading_edge_txt_object->SetProperty,strings=lead_edge_string

info2.latest_click_X = 0
info2.latest_click_Y = 0
info2.latest_release_X = 0
info2.latest_release_Y = 0
info2.click_and_drag = 0
info2.previous_click_X = 0
info2.previous_click_Y = 0

info2.the_action = 0

widget_control,info2.widget_image_time_slider, set_slider_min = 1
widget_control,info2.widget_image_time_slider, set_slider_max = num_images
widget_control,info2.widget_image_time_slider, set_value = 1
widget_control,info2.widget_image_time_slider, sensitive=sensitive_val




;widget_control,info2.widget_alpha_slider, set_value = 100
;widget_control,info2.widget_alpha_slider,sensitive=sensitive_val

widget_control,info2.widget_botSlider, set_value = 0
widget_control,info2.widget_botSlider, sensitive=sensitive_val
widget_control,info2.widget_topSlider, set_value = 100
widget_control,info2.widget_topSlider, sensitive=sensitive_val
widget_control,info2.widget_gammaSlider, set_value = 10
widget_control,info2.widget_gammaSlider, sensitive=sensitive_val
widget_control,info2.widget_saturationSlider, sensitive=sensitive_val
widget_control,info2.widget_reset_image_controls, sensitive=sensitive_val

;these should be sensitive for new images; nonsensitive for new telescope...
widget_control,info2.widget_image_menu, sensitive=sensitive_val
widget_control,info2.widget_define_CME_name, sensitive=sensitive_val
widget_control,info2.widget_full_halo_cme, sensitive=sensitive_val

info2.number_of_images = 0
info2.current_image_number = 0
info2.full_halo = 0
info2.i_images_are_loaded = sensitive_val
info2.which_diff = 0
info2.background_image_number = -1

if info2.telescope_code eq 'AC2' or info2.telescope_code eq 'BC2' then begin
     info2.image_saturation_value = 200.
     widget_control,info2.widget_saturationSlider,set_slider_max = 4000.
     widget_control,info2.widget_saturationSlider,set_value = info2.image_saturation_value
     info2.sat_initial_value = info2.image_saturation_value
endif else begin
     info2.image_saturation_value = 50.
     widget_control,info2.widget_saturationSlider,set_slider_max = 800.
     widget_control,info2.widget_saturationSlider,set_value = info2.image_saturation_value
     info2.sat_initial_value = info2.image_saturation_value
endelse


END



PRO get_images, top_id,info2,number_of_images

;20110630 anewman - modified
;This procedure works with FITS files only.
;Modified to get input images automatically after the user selects a telescope and datetime range.

compile_opt idl2

;Widget_Control, event.top, Get_UValue=info2, /No_Copy  ;if selecting get_images from a menu

;print,'At beginning of get_images'
;help,info2
info2.image_type = 'fits'   

;create an empty list to hold the full-path-names of all input images
input_image_list = list()

;get 1st image from the start date folder
find_images_by_date,info2.start_date,info2,input_image_list

;get beginning of start_date,end_date, in julian date format
start_0000 = JULDAY(info2.start_date[1],info2.start_date[2],info2.start_date[0])
end_0000 = JULDAY(info2.end_date[1],info2.end_date[2],info2.end_date[0])

;get remaining images, if the date/time interval spans multiple days
next_date = start_0000
if (start_0000 NE end_0000) then begin  ;start date is not the same as end date
    while next_date LT end_0000 do begin
      ;get next date
      next_date = next_date + 1.0D
      CALDAT,next_date,next_month,next_day,next_year,next_hour,next_minute,next_sec
      next_date_array = [next_year,next_month,next_day]
      find_images_by_date,next_date_array,info2,input_image_list
    endwhile    

endif
;print,'list=',input_image_list

number_of_images = n_elements(input_image_list)
if number_of_images eq 0 then begin
;   help,info2.start_date
;   help,info2.end_date
   msg1='No images found for ' + info2.which_telescope + ' for interval'
   msg2=info2.start_str + ' through ' + info2.end_str + '.'
   msg3='Please use "Set Dates" to select another date/time interval.'
   ok = DIALOG_MESSAGE([msg1,msg2,msg3],/information,/center)
endif else begin

*info2.pointer_to_list_of_image_files = input_image_list
info2.number_of_images = number_of_images

; Now we know the number of images, we can define the following arrays......

image_time_Julian = dblarr(number_of_images)
full_time_string = strarr(number_of_images)

;  Loop through all the images,
;  Calculate a Julian time for each of the images (from the FITS header variables).
;  Get also the datetime string that will appear in the upper right of each image

for it = 0 , number_of_images - 1 do begin

image_data = READFITS(input_image_list[it],fits_header,/silent)

image_data = rebin(image_data,512,512)
image_data = float(image_data)
get_fits_header_data, fits_header, info2.which_telescope, info2.telescope_folder, scaling_factor, $
                            rotation, SunX, SunY, datetime_string, rsun, pixel_scale

info2.pixel_scale = pixel_scale
info2.scaling_factor = scaling_factor
                            
image_data = rot(image_data, -1.*rotation, 1.0, $
                 sunX*scaling_factor,  $
                 sunY*scaling_factor, /interp)
                            
exptime_string = FXPAR(fits_header,'EXPTIME')
exposure_time = float(exptime_string)
offset_string = FXPAR(fits_header,'OFFSET')
offset = float(offset_string)

divide_by_2_factor = 1.
if info2.telescope_code eq 'AC2' or info2.telescope_code eq 'BC2' then begin
  OBS_ID = fix(FXPAR(fits_header,'OBS_ID'))
  if OBS_ID eq 1694 then divide_by_2_factor = 2.
endif
image_data = image_data/divide_by_2_factor - offset
image_data = image_data * 20. / exposure_time
info2.all_image_data_list.add, image_data



date_string = FXPAR(fits_header,'DATE-OBS')
year_string = strmid(date_string,0,4)
month_string = strmid(date_string,5,2)
day_string = strmid(date_string,8,2)

;The FITS header for LASCO and STEREO differ; STEREO does not have a TIME-OBS item
if (STRMID(info2.which_telescope,0,5) eq 'LASCO') then begin
   time_string = FXPAR(fits_header,'TIME-OBS')
   hour_string = strmid(time_string,0,2)
   min_string = strmid(time_string,3,2)
   sec_string = strmid(time_string,6)
endif else if (STRMID(info2.which_telescope,0,6) eq 'STEREO')then begin
   ;time is stored in DATE-OBS, for STEREO FITS files
   hour_string = strmid(date_string,11,2)
   min_string = strmid(date_string,14,2)
   sec_string = strmid(date_string,17)   
endif

year = fix(year_string)
month = fix(month_string)
day = fix(day_string)
hours = fix(hour_string)
mins = fix(min_string)

secs = float(sec_string)
secs = round(secs)
sec_string = strtrim(string(secs),2)
if secs lt 10 then sec_string = '0' + sec_string

image_time_Julian[it] = JULDAY(month, day, year, hours, mins, secs)

date_string = year_string + '-' + month_string + '-' + day_string

full_time_string[it] = date_string + ' ' + hour_string + ':' + min_string

endfor


*info2.pointer_to_image_time_Julian = image_time_Julian
*info2.pointer_to_full_time_string = full_time_string

; add the time string to the initial image.....

info2.ut_string_object -> setproperty, strings = full_time_string[0]

;*******************
;set initial clock angle and leading edge lines:
; code from PRO set_sun_position
thedata = fltarr(2,4)

;A small yellow 9-pixel square (3 pixel x 3 pixel) that looks
;like a yellow dot marks the center of the sun image
;Specify the square:

half_square_size = 1.5 ;half of 3 pixels

xvals = [255-half_square_size,255+half_square_size,255+half_square_size,255-half_square_size]
yvals = [255+half_square_size,255+half_square_size,255-half_square_size,255-half_square_size]

thedata[0,*] = xvals
thedata[1,*] = yvals

info2.marker_Sun_position -> SetProperty, data=thedata

;Reset for new images:
reset_for_images,info2,1,number_of_images

; unhide our next and previous icon arrows....
info2.goto_next_image_icon->SetProperty,hide = 0
info2.goto_previous_image_icon->SetProperty,hide = 0
info2.top_title_background->SetProperty,hide = 0

info2.current_image_number = 1
info2.which_diff = 1
actually_change_the_image, info2

; update the main window....
info2.Main_Window->Draw, info2.Main_View


endelse


END

;****




pro get_config,source_path,sep,telescope_ary,image_in_folder_ary, $
              image_in_root,image_out_location, $
              max_interval_in_days,input_file
;20110630 anewman - start

compile_opt idl2

 CASE STRUPCASE(!Version.OS_FAMILY)OF
   'WINDOWS' : sep = '\'
   'UNIX'    : sep = '/'   ;if running on linux, IDL returns this...
   ELSE      : sep = '/'
 ENDCASE

telescopes = ''        ;coronagraph telescopes having images to be processed
image_in_folders = ''  ;names of the folders containing the input FITS images
image_in_root = ''     ;root of the location of the input FITS image folders on the server, e.g. on nas-bes
image_out_location = ''     ;location of the output FITS image on the server, e.g. on nas-bes
max_interval_in_days = ''   ;integer; days: upper limit of the datetime interval when selecting input images
line = ''

;get configuration info from input file
input_file = source_path + sep + 'cment.in'
openr,lun,input_file,/GET_LUN
WHILE NOT EOF(lun) DO BEGIN
  READF, lun, line
  pos = STRPOS(line,'=')
  if STRCMP( line,'telescopes', pos, /FOLD_CASE ) then telescopes = STRMID(line,pos+1)
  if STRCMP( line,'image_in_folders', pos, /FOLD_CASE ) then image_in_folders = STRMID(line,pos+1)
  if STRCMP( line,'image_in_root', pos, /FOLD_CASE ) then image_in_root = STRMID(line,pos+1)
  if STRCMP( line,'image_out_location', pos, /FOLD_CASE ) then image_out_location = STRMID(line,pos+1)
  if STRCMP( line,'max_interval_in_days', pos, /FOLD_CASE ) then max_interval_in_days = STRMID(line,pos+1)
     
ENDWHILE
close,lun
free_lun, lun

if telescopes eq '' or image_in_folders eq '' or image_in_root eq '' then begin
  print,'ERROR: configuration file ',input_file,' is missing some important information.'
endif

;convert some strings to arrays:
telescope_ary = strsplit(telescopes,',',/EXTRACT)
image_in_folder_ary = strsplit(image_in_folders,',',/EXTRACT)

;convert max interval to integer
max_interval_in_days = FIX(max_interval_in_days)

END

;
;
;+
; NAME:
;       sourcepath
;
; PURPOSE:
; This procedure returns the directory path associated with
; the routine calling this function.  This is useful for
; building applications that need to bootstrap resource and
; configuration files when the installation directory may not
; be known until run time.  Use this function in conjunction
; with FILEPATH to build platform-independent file path strings
; to your resources. <br>
; For example, <pre>
;   b = WIDGET_BUTTON(tlb, /BITMAP, $
;     VALUE=FILEPATH('up.bmp', ROOT = SourcePath(), SUBDIR = ['resource'])</pre>
; This will search for a file named "up.bmp" in the subdirectory named
; "resource" below the directory in which is located the source code
; (or SAVE file) for the routine containing the above statement.
;
; @Keyword
;   Base_Name {out}{optional}{type=string}
;       Set this keyword to a named variable to retrieve the
;       base file name of the routine's source.
; @Keyword
;   Extra {in}{optional}
;       Any extra keywords are passed to the FILE_DIRNAME
;       function, for example /MARK_DIRECTORY.
;
; @Returns
;   The return value is the root directory path to
;   the calling routine's source file or SAVE file.
;
; @Examples <pre>
;   Create a file myapp.pro with the contents and run it.
;     PRO MYAPP
;     PRINT, SourcePath()
;     END
;   The printed output will be the full path to the
;   directory in which abc.pro was created, regardless of
;   IDL's current working directory.</pre>
;
; MORE DETAILS (from the ITTVIS IDL Code Library, retrieved on 3/18/2011 by anewman):
;  The SOURCEPATH function, in combination with FILEPATH, allows a program to locate 
;  other files within a routine source file's related directory tree. For example, 
;  an IDL routine file named C:\myapp\abc.pro calls SOURCEPATH as in 
;     PRO ABC
;        PRINT, SOURCEPATH() 
;     END 
;  the resulting output will be the string "C:\myapp". If data associated with the 
;  application are in C:\myapp\mydata, a data file in this directory can be located in code via 
;        datafile = FilePath('data.dat',$
;                ROOT=SourcePath(), $
;                SUBDIR=['data']) 
;   The programmer can distribute the application to another user who may install the
;   original directory tree into "D:\app". No code modifications would be required for
;   this user to successfully locate the data.dat file. If the routine ABC were compiled 
;   and saved to an IDL SAVE file and distributed, the SOURCEPATH function will return 
;   the path to the SAVE file instead. This function supercedes the SOURCEROOT function,
;   found elsewhere in the IDL codebank, as of IDL 6.2. SOURCEPATH uses a "supported" method 
;   for retrieving the file path, while SOURCEROOT parses the output from the HELP procedure,
;   a technique that is not generally recommended. 
;   October 10, 2005 - Added _EXTRA keyword to be passed to FILE_DIRNAME.


;
; @History
;   03/18/2005  JLP, RSI - Original version <br>
;   10/10/2005 JLP, RSI - On Ben Tupper's suggestion, added _EXTRA
;-

Function SourcePath, Base_Name = BaseName, _Extra = Extra
Compile_Opt StrictArr
On_Error, 2
Stack = Scope_Traceback(/Structure)
Filename = Stack[N_elements(Stack) - 2].Filename
If (Arg_Present(BaseName)) then Begin
    BaseName = File_BaseName(Filename)
EndIf
Return, File_DirName(Filename, _Extra = Extra)
End
;-------------------------------------------------------------------------
;


PRO select_telescope, event
;20110628 anewman - start

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

Widget_Control, event.id, Get_UValue=which_telescope

info2.clock_angle_model->SetProperty, hide = 1



       if which_telescope eq -1 then begin
       print,'which_telescope is undefined...'
       info2.telescope_code = 'ERR'
       endif else begin
;       print,'which_telescope=',which_telescope
       Widget_Control, event.id, Get_Value=telescope_name
       info2.which_telescope = telescope_name
       info2.telescope_string_object -> setproperty, strings = telescope_name
       info2.telescope_string_object->SetProperty,hide = 0
       info2.Main_Window->Draw, info2.Main_View
;       print,'telescope_name=',telescope_name
       CASE info2.which_telescope  OF
          'LASCO C2': info2.telescope_code = 'SC2'
          'LASCO C3': info2.telescope_code = 'SC3'
          'STEREO A COR2': info2.telescope_code = 'AC2'
          'STEREO B COR2': info2.telescope_code = 'BC2'
          ELSE: info2.telescope_code = 'ERR'
       ENDCASE

       ;Reset some items for new telescope selected:
       ;
       ;blank out any image, from a previous telescope:
       data = [0,0]
       data = rebin(data,512,512)
       info2.diff_Image_object -> SetProperty, DATA = data  
       
       ;many items that are reset for new images also must be reset for a new telescope
       reset_for_images,info2,0,0

       ;allow user to set or select date ranges now
;       widget_control, info2.widget_select_dates, sensitive=1       
;       widget_control, info2.widget_select_recent, sensitive=1
       info2.ut_string_object -> setproperty, strings = ''
       *info2.pointer_to_list_of_image_files = list()
       info2.start_julian = -1.D
       info2.end_julian = -1.D
       info2.start_str = ''
       info2.end_str = ''
       info2.cme_datetime = ''
       ; Having chosen a telescope this menu is now greyed out....       
       widget_control,info2.widget_select_telescope, sensitive=0
       
       
       
       
info2.clock_angle_model->SetProperty, hide = 1
info2.clock_angle_deg_object->SetProperty, hide = 1
info2.leading_edge_txt_object->SetProperty, hide = 1

       date_array = datebox_ghm_no_basedate(Title='Image times', Cancel=cancelled)


       ;check the validity of the dates entered
       checking = 1
       while checking eq 1 do begin

          IF cancelled THEN begin
;          print, 'Cancelled '
           checking = 0
          endif else begin
;          help, date_array
          ;print, 'the text ',date_array

;          for i = 0 , 9 do begin
;          print, i , ' ', date_array[i], '*'
;          endfor

          ;print,'ready to call check_dates'
          check_dates,date_array,start_julian,end_julian,info2.max_interval_in_days, $
                      dates_are_fine,problem_string

          if dates_are_fine eq 0 then begin  ;there's a problem with one or both dates...
          Result = DIALOG_MESSAGE(problem_string,/center)
          date_array = datebox_ghm_no_basedate(Title='Image times', Cancel=cancelled, date_array=date_array)
          endif  else begin
             checking = 0  ;have valid dates; so, done checking them
          endelse
          endelse
       endwhile
       
       IF cancelled THEN begin
;       print, 'Cancelled '
        ;user pressed the "Cancel" button
       endif else begin

; At this point we should have valid dates.......

       ;print, 'VALID DATES :'
;       print,' date_array ', date_array
       info2.date_array_int = fix(date_array)
       ;print, 'date_array_int ',info2.date_array_int

;       Result = DIALOG_MESSAGE(date_array,/center)

       info2.start_date = info2.date_array_int[0:4]
       info2.end_date = info2.date_array_int[5:9]
       info2.start_julian = start_julian
       info2.end_julian = end_julian
       
       str_month = strtrim(string(info2.start_date[1]),2)
       if strlen(str_month) eq 1 then str_month = '0' + str_month
       str_day = strtrim(string(info2.start_date[2]),2)
       if strlen(str_day) eq 1 then str_day = '0' + str_day       
       str_hour = strtrim(string(info2.start_date[3]),2)
       if strlen(str_hour) eq 1 then str_hour = '0' + str_hour
       str_min = strtrim(string(info2.start_date[4]),2)
       if strlen(str_min) eq 1 then str_min = '0' + str_min
       
       info2.start_str = str_month + '/' + str_day + '/' + $
                         strtrim(string(info2.start_date[0]),2) + ' ' + str_hour + ':' +str_min
       
       str_month = strtrim(string(info2.end_date[1]),2)
       if strlen(str_month) eq 1 then str_month = '0' + str_month
       str_day = strtrim(string(info2.end_date[2]),2)
       if strlen(str_day) eq 1 then str_day = '0' + str_day
       str_hour = strtrim(string(info2.end_date[3]),2)
       if strlen(str_hour) eq 1 then str_hour = '0' + str_hour
       str_min = strtrim(string(info2.end_date[4]),2)
       if strlen(str_min) eq 1 then str_min = '0' + str_min

       info2.end_str = str_month + '/' + str_day + '/' + $
                         strtrim(string(info2.end_date[0]),2) + ' ' + str_hour + ':' +str_min
       
          ;get input images automatically after user selects telescope and datetime range:
          get_images,event.top,info2,number_of_images
          
          info2.number_of_images = number_of_images
          
          info2.clock_angle_model->SetProperty, hide = 0
          info2.clock_angle_deg_object->SetProperty, hide = 0
          info2.leading_edge_txt_object->SetProperty, hide = 0
          info2.clock_angle_marker_line->SetProperty,hide = 0
          info2.clock_angle_string->SetProperty,hide = 0
          info2.leading_edge_marker_line->SetProperty,hide = 0
          info2.leading_edge_string->SetProperty,hide = 0
          info2.telescope_string_object->SetProperty,hide = 0
          info2.ut_string_object->SetProperty,hide = 0

       endelse  ;not cancelled
       
       
       
       
       
       
       
       
       
       endelse  ; which_telescope
      
Widget_Control, event.top, Set_UValue=info2, /No_Copy

END











PRO color_scale_top_slider, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

info2.image_color_palette -> setproperty,top_stretch = event.value

info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy
END




PRO color_scale_bottom_slider, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

info2.image_color_palette -> setproperty,bottom_stretch = event.value

info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy
END 




PRO color_scale_gamma_slider, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

info2.image_color_palette -> setproperty,gamma = float(event.value)/ 10.

info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy
END


PRO calc_clock_angle, info2, x2, y2

   compile_opt idl2

   if info2.full_halo then $
      ;When the CME is full halo,
      ; the degrees are always 360...
      set_clock_angle,info2,360 $
   else  begin
      ;Get the clock angle, based on previous location
      ;and current location of the red clock angle model
      ;line
      x1 = info2.previous_click_X
      y1 = info2.previous_click_Y
      xc = 255
      yc = 255
      theta2 = atan((y2 - yc) , (x2 - xc))
      theta1 = atan((y1 - yc) , (x1 - xc))
      rotation_angle_degrees = (theta2 - theta1) * 180. / !pi
      info2.previous_click_X = x2
      info2.previous_click_Y = y2
         
      info2.clock_angle_model -> rotate, [0,0,1], rotation_angle_degrees, /premultiply
      info2.clock_angle_model -> getProperty, transform=transform
      
      val1=round(( atan(transform[1,0], transform[0,0]) * 180. / !pi ))
      if val1 gt 0 then begin
         val2 = 360 - val1
      endif else begin
         val2 = 0 - val1
      endelse
 
      set_clock_angle,info2,val2
  
   endelse

END



PRO set_clock_angle,info2,degrees

   compile_opt idl2

   ;save the clock angle, in degrees, and display it on the image
   info2.clock_angle_degrees = degrees
   angle_string = strcompress(string(info2.clock_angle_degrees))
   if info2.clock_angle_degrees lt 100 then angle_string = ' ' + angle_string
   if info2.clock_angle_degrees lt 10 then angle_string = ' ' + angle_string
   ;Displays new value of the clock angle degrees near the outer part of the rotating line
   info2.clock_angle_string->SetProperty,strings=angle_string
   ;Display new value of the clock angle degrees in lower left corner of image:
    info2.clock_angle_deg_object->SetProperty,strings=angle_string

END


PRO set_full_halo_cme, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

if (event.select eq 1) then begin
  ;print, 'Full Halo ON\n'
  info2.full_halo = 1
  
  ;For full halo CMEs, the clock angle is always 360 degrees
  set_clock_angle,info2,360 
  
  
endif else begin
  ;print, 'Full Halo OFF\n'
  info2.full_halo = 0
  
  ;calculate and display the current clock angle in degrees
  calc_clock_angle,info2,info2.rotate_x,info2.rotate_y
  
  
endelse

;Must change CME Name if it has already been selected...
if strlen(info2.cme_name) gt 0 then begin
   image_datetime = STRMID(info2.cme_name,0,13)
   build_CME_name,info2,image_datetime
   info2.Info_Window->Draw, info2.Info_View
endif

info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy
END



PRO change_transparency, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

alpha_value = float(event.value) / 100. 
info2.clock_angle_marker_line ->SetProperty, alpha_channel = alpha_value
info2.clock_angle_string ->SetProperty, alpha_channel = alpha_value
info2.leading_edge_marker_line ->SetProperty, alpha_channel = alpha_value
info2.leading_edge_string ->SetProperty, alpha_channel = alpha_value
info2.Main_Window->Draw, info2.Main_View                                   

Widget_Control, event.top, Set_UValue=info2, /No_Copy

END





PRO diff_format, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

Widget_Control, event.id, Get_UValue=which_diff_format

CASE which_diff_format OF

   'Normal': BEGIN
   
      info2.which_diff = 0

      ENDCASE
      
   'Diff: Running': BEGIN
   
      info2.which_diff = 1

      ENDCASE
      
   'Diff: Set Current Image as Background': BEGIN
   
      info2.which_diff = 2
      info2.background_image_number = info2.current_image_number

      ENDCASE

ENDCASE

actually_change_the_image, info2

info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy

END







PRO change_the_image_saturation, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

info2.image_saturation_value = event.value

actually_change_the_image, info2

info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy

END



PRO reset_image_controls, event

  compile_opt idl2
  
  Widget_Control, event.top, Get_UValue=info2, /No_Copy
  
  info2.image_color_palette -> setproperty,bottom_stretch = info2.bot_initial_value
  info2.image_color_palette -> setproperty,top_stretch = info2.top_initial_value
  info2.image_color_palette -> setproperty,gamma = (float(info2.gamma_initial_value))/10.
  info2.image_saturation_value = info2.sat_initial_value
  
  widget_control,info2.widget_botSlider, set_value = info2.bot_initial_value
  widget_control,info2.widget_topSlider, set_value = info2.top_initial_value
  widget_control,info2.widget_gammaSlider, set_value = info2.gamma_initial_value
  widget_control,info2.widget_saturationSlider, set_value = info2.sat_initial_value
  
  actually_change_the_image, info2
  
  info2.Main_Window->Draw, info2.Main_View
  
  Widget_Control, event.top, Set_UValue=info2, /No_Copy
  
END




PRO exit_the_prog, event

;*****************************************************************
;
;  This subroutine is executed when the 'Exit' button is clicked
;  it kills everything.
;  The cleanup routine (next) is also called.
;
;****************************************************************

compile_opt idl2

   ; Exit the program. This will cause the CLEANUP
   ; routine to be called automatically.

Widget_Control, event.top, Get_UValue=info2, /No_Copy

openw,lun,info2.cment_preferences_file, /get_lun
printf,lun,info2.tlb_position
close,lun
free_lun, lun


Widget_Control, event.top, Set_UValue=info2, /No_Copy

Widget_Control, event.top, /Destroy

END





PRO stp_Cleanup2, tlb2

;*****************************************************************
;
;  cleanup - cleans up - Frees all created objects.
;
;****************************************************************

compile_opt idl2

Widget_Control, tlb2, Get_UValue=info2
IF N_Elements(info2) NE 0 THEN begin
;print, 'destroying '
;help,/mem
;heap_free,info2
;help,/mem
;heap_gc
;help,/mem
ENDIF
END




PRO Reset, event

;*****************************************************************
;
;  Reset - intended to be used to allow the user to reset everything
;  and start again without quitting the application.
;
;****************************************************************

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

info2.clock_angle_model->SetProperty, transform = info2.clock_angle_model_original_transform
info2.clock_angle_deg_object->SetProperty,hide = 1
info2.leading_edge_txt_object->SetProperty,hide = 1
info2.clock_angle_marker_line->SetProperty,hide = 1
info2.clock_angle_string->SetProperty,hide = 1
info2.leading_edge_marker_line->SetProperty,hide = 1
info2.leading_edge_string->SetProperty,hide = 1
info2.telescope_string_object->SetProperty,hide = 1
info2.ut_string_object->SetProperty,hide = 1
widget_control,info2.widget_select_telescope, sensitive=1
;widget_control, info2.widget_select_dates, sensitive=0       
;widget_control, info2.widget_select_recent, sensitive=0
image_data = bytarr(512,512)
info2.diff_image_object -> setproperty, data = image_data
info2.i_representative_image_has_been_chosen = 0
reset_for_images,info2,0,0
info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy

END





PRO name_another_cme, event
  
  compile_opt idl2
  
  Widget_Control, event.top, Get_UValue=info2, /No_Copy
  
  info2.clock_angle_model->SetProperty, transform = info2.clock_angle_model_original_transform
  info2.clock_angle_deg_object->SetProperty,hide = 0
  info2.leading_edge_txt_object->SetProperty,hide = 0
  info2.clock_angle_marker_line->SetProperty,hide = 0
  info2.clock_angle_string->SetProperty,hide = 0
  info2.leading_edge_marker_line->SetProperty,hide = 0
  info2.leading_edge_string->SetProperty,hide = 0
  info2.telescope_string_object->SetProperty,hide = 0
  info2.ut_string_object->SetProperty,hide = 0
  info2.full_halo = 0
  info2.i_representative_image_has_been_chosen = 0
  info2.goto_next_image_icon->SetProperty,hide = 0
  info2.goto_previous_image_icon->SetProperty,hide = 0
  info2.top_title_background->SetProperty,hide = 0
  widget_control,info2.widget_image_menu, sensitive=1
  
  ; reset i_cme_has_been_named to let us name another one.....
  
  info2.i_cme_has_been_named = 0
  info2.i_images_are_loaded = 1
  widget_control,info2.widget_define_CME_name, sensitive=1
  widget_control,info2.widget_image_time_slider, sensitive=1
  widget_control,info2.widget_botSlider, sensitive=1
  widget_control,info2.widget_topSlider, sensitive=1
  widget_control,info2.widget_gammaSlider, sensitive=1
  widget_control,info2.widget_saturationSlider, sensitive=1
  widget_control,info2.widget_reset_image_controls, sensitive=1
  widget_control,info2.widget_full_halo_cme, sensitive=1
  widget_control,info2.widget_name_another_cme, sensitive=0
  info2.leading_edge_in_rs = 0
  
  info2.clock_angle_string->SetProperty,strings=' 0'
  info2.clock_angle_deg_object->SetProperty,strings=' 0'
  info2.leading_edge_string->SetProperty,strings=' 0'
  info2.leading_edge_txt_object->SetProperty,strings=' 0'
  
;  info2.clock_angle_deg_object->SetProperty,strings=info2.clock_angle_string_value
  
  info2.Main_Window->Draw, info2.Main_View
  info2.Info_Window->Draw, info2.Info_View
  
  Widget_Control, event.top, Set_UValue=info2, /No_Copy
  
END








FUNCTION compare_image_datetimes, cme_naming_file, cme_representative_file
;input: 
;filename for the CME naming file
;filename for the CME representative file
;Assumes that the filename starts with the datetime, in
;this format:
;YYYYMMDD_HHMM
;where YYYYMMDD is the year month day e.g. 20110825
;and HHMM is the time hour minute, e.g. 2354
;for example 20110825_2354...
;RETURNS
;returns 1 if the images are in the correct order
;returns 0 if the CME representative image precedes the CME naming image

compile_opt idl2

;get the datetime of the CME naming file, and convert it to Julian day number
;print,'CME named FITS = ', cme_naming_file
yy = fix(strmid(cme_naming_file,0,4))
mo = fix(strmid(cme_naming_file,4,2))
da = fix(strmid(cme_naming_file,6,2))
hh = fix(strmid(cme_naming_file,9,2))
mi = fix(strmid(cme_naming_file,11,2))
cme_name_julian = JULDAY(mo,da,yy,hh,mi,0)

;get the datetime of the CME Representative Image file, and convert it to Julian day number
;print,'Representative Image filename =' ,cme_representative_file
yy = fix(strmid(cme_representative_file,0,4))
mo = fix(strmid(cme_representative_file,4,2))
da = fix(strmid(cme_representative_file,6,2))
hh = fix(strmid(cme_representative_file,9,2))
mi = fix(strmid(cme_representative_file,11,2))
rep_image_julian = JULDAY(mo,da,yy,hh,mi,0)

;datetime of CME-naming image should be earlier 
;or equal to the datetime of CME representative image

if cme_name_julian gt rep_image_julian then return,0 $
else return,1

END


PRO build_CME_name,info2,image_datetime
;Assemble the name of the FITS output image file.
    
   compile_opt idl2

   ;clock angle degrees range from 0-360 degrees; make the string 3-char. 
   angle_string = strtrim(string(info2.clock_angle_degrees),2)
   if info2.clock_angle_degrees lt 100 then angle_string = '0' + angle_string
   if info2.clock_angle_degrees lt 10 then angle_string = '0' + angle_string

   ;CME leading edge in solar radii range from 0-99; make the string 2-char.
   leading_edge_string = strtrim(string(info2.leading_edge_in_rs),2)
   if info2.leading_edge_in_rs lt 10 then leading_edge_string = '0' + leading_edge_string
   
   image_datetime_iso8601 = strmid(image_datetime,0,8) + 'T' + strmid(image_datetime,9,4)
   
   i = image_datetime_iso8601
   
   info2.cme_catalog_datetime_for_fits_header = strmid(i,0,4) + '-' + strmid(i,4,2) + '-' + strmid(i,6,5) + ':' + strmid(i,11,2) + ':00'
   
   info2.cme_name = image_datetime_iso8601 + '-' + $
                 angle_string + '-' + $
                 leading_edge_string + '-' + $
                 info2.telescope_code
                 
   ;print,'CME_name = ',info2.cme_name


   info2.cme_name_string-> SetProperty, strings = info2.cme_name

END


PRO reset_CME_name, info2
;reset the CME name to an empty string, to force the user to press
; the "Define CME Name" button.
;This procedure should be called whenever the user changes
; something that would affect the name of the CME, e.g., if the
; user changes the clock rotation angle, or the leading edge value

  compile_opt idl2
  
  info2.cme_name = ''
  info2.cme_name_string-> SetProperty, strings = ' '
  info2.i_cme_has_been_named = 0

  ;WARNING: DO NOT reset info2.cme_datetime here, as it won't
  ; get set again if the user does not select another image.
 
  info2.Info_Window->Draw, info2.Info_View


END


PRO define_CME_name, event
; A CME named file is in this format:
; YYYYMMDDTHHMM-XXX-YY-ZZZ.fts
; where YYYYMMDD is the date string
;       HHMM is the time string
;       XXX is a 3-char clock angle string ( 000 through 360 )
;       YY is a 2-character leading edge string ( 00 through 99 )
;       ZZZ is the 3-char telescope code
;       .fts is the FITS extension

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

;Get the datetime string from the input FITS filename
list_of_image_files = *info2.pointer_to_list_of_image_files
filename = list_of_image_files[info2.current_image_number]
pos = strpos(filename,info2.sep,/reverse_search)
filename = STRMID(filename,pos+1)

ordered_correctly = 1

;if the CME Representative Image has been set
if strlen(info2.cme_image_name) gt 0 then begin
  ;check to make sure the CME naming file does not follow the CME representative image file
  ordered_correctly = compare_image_datetimes(filename, info2.cme_image_name)
endif

if ~ordered_correctly then begin
    msg1 = 'The CME-naming image must not be later than the CME Representative Image.'
    msg2 = 'Please select another image.'
    ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
endif else begin

if info2.leading_edge_in_rs eq 0 then begin
    msg1 = 'A leading edge value for the CME has not been defined.'
    msg2 = 'Move the yellow line to mark the position of the leading edge,'
    msg3 = 'then press the "Define CME Name" button again.'
    ok = DIALOG_MESSAGE([msg1,msg2,msg3],/information,/center)
endif else begin


   image_datetime = STRMID(filename,0,13)

   build_CME_name,info2,image_datetime
   
   info2.Info_Window->Draw, info2.Info_View

   widget_control,info2.widget_define_CME_name, sensitive=0
   widget_control,info2.widget_full_halo_cme, sensitive=0
   widget_control,info2.widget_define_CME_representative_image, sensitive=1
   info2.i_cme_has_been_named = 1
   
   info2.clock_angle_deg_object->SetProperty,hide = 0
   info2.leading_edge_txt_object->SetProperty,hide = 0
   info2.clock_angle_marker_line->SetProperty,hide = 0
   info2.clock_angle_string->SetProperty,hide = 0
   info2.leading_edge_marker_line->SetProperty,hide = 0
   info2.leading_edge_string->SetProperty,hide = 0
   info2.telescope_string_object->SetProperty,hide = 1
;   info2.ut_string_object->SetProperty,hide = 1
   
   info2.Main_Window->Draw, info2.Main_View
   
endelse

endelse

Widget_Control, event.top, Set_UValue=info2, /No_Copy

END




PRO reset_CME_image, info2

  ;Return CME Representative Image to the state it was in
  ;before the user selected a CME Representative Image;
  ; - and make the button unselectable
  
  compile_opt idl2  
  
  info2.cme_image_name = ''
  info2.cme_image_string-> SetProperty, strings = ' '
  info2.i_representative_image_has_been_chosen = 0
  
  widget_control,info2.widget_define_CME_representative_image, sensitive=0
  
  ;Also: Do not allow the user to save the CME
  widget_control,info2.widget_save_cme, sensitive=0
  
  info2.Info_Window->Draw, info2.Info_View


END



PRO define_cme_representative_image, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

;print, ' in Define CME image '

list_of_image_files = *info2.pointer_to_list_of_image_files
filename = list_of_image_files[info2.current_image_number]
pos = strpos(filename,info2.sep,/reverse_search)

filename = STRMID(filename,pos+1)

;check to make sure the CME representative image file does not preceed the CME naming file
ordered_correctly = compare_image_datetimes(info2.cme_name, filename)

if ~ordered_correctly then begin
    msg1 = 'The Representative Image must not be earlier than the CME-naming image.'
    msg2 = 'Please select another image.'
    ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
endif else begin

   info2.cme_image_name = filename

   cme_image_name_string = strmid(filename,0,4)+'-'+strmid(filename,4,2)+'-'+strmid(filename,6,2) + ' ' + strmid(filename,9,2)+':'+strmid(filename,11,2)

   info2.cme_image_string-> SetProperty, strings = 'image time : ' + cme_image_name_string

   info2.Info_Window->Draw, info2.Info_View

   widget_control,info2.widget_define_CME_representative_image, sensitive=0
   widget_control,info2.widget_image_menu, sensitive=0
   widget_control,info2.widget_image_time_slider, sensitive=0
   widget_control,info2.widget_botSlider, sensitive=0
   widget_control,info2.widget_topSlider, sensitive=0
   widget_control,info2.widget_gammaSlider, sensitive=0
   widget_control,info2.widget_saturationSlider, sensitive=0
   widget_control,info2.widget_reset_image_controls, sensitive=0
   info2.ut_string_object->SetProperty, hide = 1
   
   widget_control,info2.widget_save_cme, sensitive=1
   info2.i_representative_image_has_been_chosen = 1
   info2.goto_next_image_icon->SetProperty,hide = 1
   info2.goto_previous_image_icon->SetProperty,hide = 1


   
endelse

   Widget_Control, event.top, Set_UValue=info2, /No_Copy
END


PRO create_cme_catalog_file, info2

;Write the FITS output image file:
;Create the FITS header;
;Save the output FITS file in its final location,
; with its "Representative Image" data, and its header
;20110718 anewman - start
;
;using FITS I/O in the IDL Astronomy Library
;http://idlastro.gsfc.nasa.gov/fitsio.html
; search for needed IDL procedures/functions in google. e.g.: IDL WRITEFITS


compile_opt idl2


if info2.i_confirm_cme_save_to_spi eq 1 then begin
msg1 = 'You are about to catalog the CME : '
msg2 = info2.cme_name
msg3 = 'Are you sure?'
dialog_title = 'CME ' + info2.cme_name
result = DIALOG_MESSAGE([msg1,msg2,msg3],/question,/center,title=dialog_title)
endif else begin
  result = 'Yes'
endelse

if result eq 'Yes' then begin

;print,'In create_cme_catalog_file'

;This is the "Representative Image" the user selected
info2.diff_image_object -> getproperty, data = image_data

;help, image_data
;
;image_data[0:511,489:511] = 0

info2.diff_image_object -> setproperty, data = image_data

info2.telescope_string_object->SetProperty, hide = 1
info2.marker_Sun_position->SetProperty, hide = 0
info2.clock_angle_model->SetProperty, hide = 0
info2.ut_string_object->SetProperty, hide = 1
info2.clock_angle_deg_object->SetProperty, hide = 0
info2.leading_edge_txt_object->SetProperty, hide = 0
info2.clock_angle_marker_line->SetProperty,hide = 0
info2.clock_angle_string->SetProperty,hide = 1
info2.leading_edge_marker_line->SetProperty,hide = 0
info2.leading_edge_string->SetProperty,hide = 1
info2.goto_next_image_icon->SetProperty,hide = 1
info2.goto_previous_image_icon->SetProperty,hide = 1
info2.top_title_background->SetProperty,hide = 1

info2.main_window->Draw, info2.main_View

info2.main_window -> getproperty, image_data = image_data2

;info2.telescope_string_object->SetProperty, hide = 0
;info2.marker_Sun_position->SetProperty, hide = 0
;info2.clock_angle_model->SetProperty, hide = 0
;info2.ut_string_object->SetProperty, hide = 0
;info2.clock_angle_deg_object->SetProperty, hide = 0
;info2.leading_edge_txt_object->SetProperty, hide = 0
;
;info2.main_window->Draw, info2.main_View

snapshot2 = bytarr(3,512,512)
snapshot2[0:2,0:511,0:511] = image_data2[0:2,0:511,0:511]

info2.info_window -> getproperty, image_data = info_data2

snapshot3 = bytarr(3,512,100)
snapshot3[0:2,0:511,0:99] = info_data2[0:2,0:511,0:99]

snapshot4 = bytarr(3,512,612)
snapshot4[0:2,0:511,0:99] = snapshot3[0:2,0:511,0:99]
snapshot4[0:2,0:511,100:611] = snapshot2[0:2,0:511,0:511]

snapshot5 = bytarr(512,612,3)
for i = 0 , 2 do begin
  snapshot5[0:511,0:611,i] = snapshot4[i,0:511,0:611]
endfor
snapshot4 = snapshot5

im_info = size(snapshot4,/struc)
;help,im_info
ndimen = im_info.n_dimensions

if ndimen GT 0 then dimen = im_info.dimensions[0:ndimen-1] $
else begin
  msg1 = 'ERROR: Bad Image! Dimensions for the image are ' + string(ndimen)
  msg2 = 'which is incorrect.'
  msg3 = ' '
  msg4 = 'Please contact the SWPC Enlil software support team...'
  ok = DIALOG_MESSAGE([msg1,msg2,msg3,msg4],/ERROR,/center)
  EXIT 
endelse

idltype = im_info.type
case idltype of
     1:  bitpix = 8
     2:  bitpix = 16
     4:  bitpix = -32     
     3:  bitpix = 32
     5:  bitpix = -64
     12: bitpix = 16
     13: bitpix = 32
 endcase

;Get the name of the cme catalog file's submitter
; i.e., the name of the person logged in and running this app.
submitter_name = GETENV('USERNAME')

;Create and fill the output (CME Naming/Cataloging) FITS header
;do not separately define the header; fxaddpar will create the fits header
fxaddpar, header, 'SIMPLE','T','File conforms to FITS standard.'
fxaddpar, header, 'BITPIX',bitpix,'Number of bits per data pixel'
fxaddpar, header, 'NAXIS',ndimen,'Number of data axes'
fxaddpar, header, 'NAXIS1',dimen[0],'Length of data axis 1'
fxaddpar, header, 'NAXIS2',dimen[1],'Length of data axis 2'
fxaddpar, header, 'CME_DATE',info2.cme_catalog_datetime_for_fits_header,'Date & Time of CME image'
fxaddpar, header, 'TELESCOP',info2.which_telescope,'Name of the instrument'
fxaddpar, header, 'TELECODE',info2.telescope_code,'Unique 3 letter designation for an instrument'
fxaddpar, header, 'CME_PA',info2.clock_angle_degrees,'CME Central position angle (000 to 359 deg)'
fxaddpar, header, 'LEADEDGE',info2.leading_edge_in_rs,'Front boundary of a CME, in Solar radii'
fxaddpar, header, 'SUB_NAME',submitter_name,'Name of the person who submitted this CME catalog file'

;NOTE: FITS header variable names can be no longer than 8 CHARS !
;      Use the underscore character rather than dashes in the var names !

;print,'header=',header

;define the full-path-name of the output FITS file; full-path: where it is to be saved
;(image out location/path is defined in the input file for this program, cment.in)
image_filename = info2.image_out_location + info2.sep + info2.cme_name + '.fts'

;MAY WANT TO CHECK EXISTENCE OF LOCATION TO WRITE TO; AVAILABLE SIZE ???

;Save the output FITS file, with its "Representative Image" data, and its header
writefits,image_filename,snapshot4,header,/CHECKSUM

; png version - (not used).....
; 
;png_image_filename = info2.image_out_location + info2.sep + info2.cme_name + '.png'
;
;write_png,png_image_filename,transpose(snapshot4, [2,0,1])
;print,' '
;print,header

;print,'image_filename = ',image_filename

;msg1 = 'The CME catalog file ' + info2.cme_name + ' has been saved'
;msg2 = 'at ' + info2.image_out_location + '.'
;ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)

;msg1 = 'Done!'
;ok = DIALOG_MESSAGE(msg1,/information,/center)

; Now the cme has been saved we can name another one - so we sensitize the widget....

widget_control,info2.widget_name_another_cme, sensitive=1

widget_control, info2.widget_save_cme, sensitive=0

endif else begin
  
; you get here if you don't want to save the CME to SPI
; we just redraw and sensitize the name another CME widget....
; if the user gets here - they can have another go at saving
; the CME, or load in a new one.... 
  
  info2.Main_Window->Draw, info2.Main_View
  info2.info_Window->Draw, info2.info_View
  widget_control,info2.widget_name_another_cme, sensitive=1
  
endelse

END



PRO save_CME, event

  ;does some error checking, then calls create_cme_catalog_file
  ;to write the FITS output image file to its final location.
  
  compile_opt idl2
  
  Widget_Control, event.top, Get_UValue=info2, /No_Copy
  
  
  ;check to make sure the CME naming file does not follow the CME representative image file
  ordered_correctly = compare_image_datetimes(info2.cme_name, info2.cme_image_name)
  
  if ~ordered_correctly then begin
    msg1 = 'The CME-naming image must not be later than the CME Representative Image.'
    msg2 = 'Please select another image.'
    ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
  endif else begin
    ;create the CME catalog file, in FITS format:
    create_cme_catalog_file,info2
  endelse
  
  ;Note:  The procedure pick_and_save_file is not being used at
  ;this time; keep it in case it is needed in the future...
  ;pick_and_save_file
  
  Widget_Control, event.top, Set_UValue=info2, /No_Copy
  
END


PRO rotate_image_to_solar_north,info2,image_data

   ;rotate the image to Solar North, and place Sun center in the center of the image using IDL's ROT function
   ;First, keep image 512x512 (as set in a rebin cmd in calling routine, or sometime previously)
   ; need to adjust sun center (which specifies the sun center in the original size image, e.g.
   ;256x256, or 1024x1024) correctly in rotated 512x512 image, by applying the scale factor:
   
   compile_opt idl2
 
   image_data = rot(image_data, -1.*info2.rotation, 1.0, $
                    info2.center_of_sunX*info2.scaling_factor,  $
                    info2.center_of_sunY*info2.scaling_factor, /interp)
  

END





PRO change_current_image, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

which_image = event.value - 1

info2.current_image_number = which_image

actually_change_the_image, info2

info2.Main_Window->Draw, info2.Main_View

Widget_Control, event.top, Set_UValue=info2, /No_Copy

END



pro actually_change_the_image, info2

compile_opt idl2

all_time_strings = *info2.pointer_to_full_time_string
info2.ut_string_object -> setproperty, strings = all_time_strings[info2.current_image_number]

image_data = (info2.all_image_data_list)[info2.current_image_number]
if info2.which_diff gt 0 then begin
  if info2.which_diff eq 1 then begin
  back_image = info2.current_image_number - 1
  if back_image lt 0 then back_image = 0
  endif else begin
    back_image = info2.background_image_number
  endelse
image_data2 = (info2.all_image_data_list)[back_image]
image_data = image_data - image_data2

image_data = image_data < info2.image_saturation_value
image_data = image_data > (0. - info2.image_saturation_value)

endif

info2.diff_image_object -> setproperty, data = bytscl(image_data)

end






PRO top_level_base_events,event

;*****************************************************************
;
;  This routine is called if the main window is moved on the screen
;  The position of the window on the screen is monitored so that
;  - on exiting the application - the 'last position' is noted.
;  .... so the next time the application is launched - it comes
;  right back up in the same place on the monitor.... 
;
;               it just makes life a little nicer ;o)
;
;*****************************************************************

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

info2.tlb_position[0] = event.x
info2.tlb_position[1] = event.y

Widget_Control, event.top, Set_UValue=info2, /No_Copy

end





PRO main_window_button_click_events, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

if info2.i_representative_image_has_been_chosen eq 0 then begin
if info2.i_images_are_loaded eq 1 then begin

drawTypes = ['PRESS', 'RELEASE', 'MOTION', 'SCROLL', 'EXPOSE']
thisEvent = drawTypes[event.type]


CASE thisEvent OF

   'EXPOSE': BEGIN
       END
   'PRESS': BEGIN

       info2.latest_click_X = event.x
       info2.latest_click_Y = event.y
       info2.previous_click_X = event.x
       info2.previous_click_Y = event.y
       info2.click_and_drag  = 1

       near_clock_angle = 0
       near_leading_edge = 0   
       myLoc = [event.x, event.y]    
       oObjArr = info2.Main_Window->Select(info2.Main_View, myLoc, dimensions=[20,20])
       nSel = N_ELEMENTS(oObjArr)
       if nsel gt 1 then begin
       iobject = intarr(nSel)
       
       FOR i=0, nSel-1 DO BEGIN
          oObjArr[i]->GetProperty, NAME=name
          if name eq 'leading_edge' then begin
            near_leading_edge = 1
          endif
          if name eq 'clock_angle' then begin
            near_clock_angle = 1
          endif
       ENDFOR
       endif
       
       if near_leading_edge eq 1 then begin
       info2.the_action = 2
       endif else begin
       if near_clock_angle eq 1 then begin
       info2.the_action = 1
       endif else begin
       info2.the_action = 0
       endelse
       endelse

       END

   'RELEASE': BEGIN
    
    if info2.the_action eq 0 then begin
      
; clicking the previous and next image icons....

    if event.y le 277 and event.y ge 252 then begin
      
      if event.x le 30 then begin
      info2.current_image_number --
      if info2.current_image_number lt 0 then info2.current_image_number = 0
      endif 
      
      if event.x ge 482 then begin
      info2.current_image_number ++
      if info2.current_image_number gt info2.number_of_images - 1 then info2.current_image_number = info2.number_of_images - 1
      endif
      
    actually_change_the_image, info2
    widget_control,info2.widget_image_time_slider,set_value=info2.current_image_number+1
    info2.Main_Window->Draw, info2.Main_View
       
    endif else begin
     
    endelse
    
    endif

       info2.latest_release_X = event.x
       info2.latest_release_Y = event.y
       info2.click_and_drag  = 0
       
       END
   'MOTION': BEGIN
    
if info2.i_cme_has_been_named eq 0 then begin

   if info2.click_and_drag eq 1 then begin
   
if info2.the_action eq 2 then begin

     ;Save the current leading edge value:
     old_leading_edge = info2.leading_edge_in_rs


; Translate the leading edge_model
         prev_x = info2.previous_click_X
         prev_y = info2.previous_click_Y
         prev_dist = sqrt((prev_x - 255)^2 + (prev_y - 255)^2)
         now_x = event.x
         now_y = event.y
         now_dist = sqrt((now_x - 255)^2 + (now_y - 255)^2)
         info2.leading_edge_model -> translate, 0.0 , now_dist - prev_dist , 0.0
         info2.leading_edge_model -> GetProperty, Transform = Transform
;           info2.Ellipse_center_location[0] = transform[3,0]
;           info2.Ellipse_center_location[1] = transform[3,1]
         info2.previous_click_X = event.x
         info2.previous_click_Y = event.y
         
         ;20110715 ahn add this calc here, plus add leading_edge_to_rs to the info2 structure
         ;The pixel scale is related to the original image, before we converted the image
         ;to 512x512. So we must apply a scale_factor to the pixel scale to make the
         ;leading edge calculation correct. The pixel scale is the number of arc seconds
         ;per pixel. So, for example, if the original image is twice as large as the 512x512
         ;current image then there will be twice as many arc seconds per pixel, ( i.e. pixel_scale* 2 or
         ;pixel scale / 0.5 ) in the current image.
         
         leading_edge_float = (now_dist * (info2.pixel_scale / info2.scaling_factor)) / info2.rsun 
         info2.leading_edge_in_rs = round(leading_edge_float)   
         
         lead_edge_string = string(leading_edge_float,format='(f5.1)')
;         if info2.leading_edge_in_rs lt 10 then lead_edge_string = ' ' + lead_edge_string
         ;Displays new value of the leading edge near the yellow leading edge line
         info2.leading_edge_string->SetProperty,strings=lead_edge_string
         ;Display new value of the location of the CME leading edge
         ;in lower rt corner of image:
         lead_edge_string_integer = strcompress(string(info2.leading_edge_in_rs),/remove_all)
         if info2.leading_edge_in_rs lt 10 then lead_edge_string_integer = ' ' + lead_edge_string_integer
         info2.leading_edge_txt_object->SetProperty,strings=lead_edge_string_integer

        ; if user has changed the leading edge, and if
        ; the CME file name has already been set,
        ; then reset/erase it:
        if old_leading_edge ne info2.leading_edge_in_rs and $
           strlen(info2.cme_name) ne 0 then begin
           reset_CME_name, info2
           
           ;check the Representative Image; if set, then reset/erase it
           if strlen(info2.cme_image_name) ne 0 then begin
              reset_CME_image,info2
           endif           
                     
        endif
         
endif     ; leading edge action      
if info2.the_action eq 1 then begin

if info2.full_halo eq 0 then begin

     ;Save the current clock angle:
     old_degrees = info2.clock_angle_degrees
     
     ; The user has rotated (or at least touched) the clock angle model....
     
     ; so calculate the clock angle in degrees
     calc_clock_angle,info2,event.x,event.y
         
     ; if user has changed the clock angle, and if
     ; the CME file name has already been set,
     ; then reset/erase it:
     if old_degrees ne info2.clock_angle_degrees and $
        strlen(info2.cme_name) ne 0 then begin
        
        reset_CME_name, info2
        
        ;check the Representative Image; if set, then reset/erase it
        if strlen(info2.cme_image_name) ne 0 then begin
           reset_CME_image,info2
        endif
     endif
     
     ;save the event values, in case a set full halo is released,
     ;and the program needs to calc & display the degrees at the current
     ;red clock angle's line location
     info2.rotate_x = event.x
     info2.rotate_y = event.y
     
endif  ; info2.full_halo eq 0              

endif  ; clock angle action

if info2.the_action eq 0 then begin
  
; we get here if we are dragging in the window
; but not moving the clock or leading edge (ie, somewhere else).
  
;  print, ' on my own ', event.x
  
endif


         
   endif
   
   endif  ; if info2.i_cme_has_been_named eq 0 then begin


       END

   ELSE:


ENDCASE

info2.Main_Window->Draw, info2.Main_View

endif ; if info2.i_cme_has_been_named 
endif ; i_representative_image_has_been_chosen
Widget_Control, event.top, Set_UValue=info2, /No_Copy
END








PRO info_Window_button_click_events, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info2, /No_Copy

drawTypes = ['PRESS', 'RELEASE', 'MOTION', 'SCROLL', 'EXPOSE']
thisEvent = drawTypes[event.type]

;print,' shift click ',event.modifiers
shift_click = event.modifiers

CASE thisEvent OF

   'EXPOSE': BEGIN
       END
   'PRESS': BEGIN

;       print,'Plot Window X and Y : ',event.x , event.y

       END

   'RELEASE': BEGIN


       END
   'MOTION': BEGIN ; Trackball events


       END

   ELSE:


ENDCASE

    ; Draw the view.

info2.info_Window->Draw, info2.info_View

    ;Put the info structure back.
    
Widget_Control, event.top, Set_UValue=info2, /No_Copy
END







PRO cment

;Main procedure.
;Defines and initializes the Widget (object graphics) objects.
;Defines and initializes the info2 structure.
;Realizes the Main View object.

compile_opt idl2

; a variable to control the 'Are you sure?' popup. 
i_confirm_cme_save_to_spi = 0

color_red = [138,226,70]
color_yellow = [220,220,0]
color_fill = [50,50,50]

;get directory path associated with this IDL procedure
source_path = SourcePath()
;print,'source_path=',source_path

;get input values that can change, from a configuration file
get_config,source_path,sep,telescope_array,image_in_folder_array, $
           image_in_root,image_out_location,max_interval_in_days,input_file
           
number_of_images = 0

Xsize = 512
Ysize = 512
Ysize_info = 100

;  Need to define a main_view and a main_model......

Main_View = OBJ_NEW('IDLgrView', Color=[0,0,0], Viewplane_Rect=[0.0,0.0,Xsize,Ysize])

Main_Model = OBJ_NEW('IDLgrModel') 

Main_View -> Add, Main_Model

 

;  and a info_View and info_Model......

info_View = OBJ_NEW('IDLgrView', Color=[0,0,0], Viewplane_Rect=[0.0,0.0,xsize,ysize_info])

info_Model = OBJ_NEW('IDLgrModel') 

info_View -> Add, info_Model

app_data_path = getenv('APPDATA')
cment_preferences_file = app_data_path + sep + 'swpc_cment_prefs'
;
;  define all of the widgets in our GUI ...
;

openr,lun,cment_preferences_file, error = err, /get_lun
if err eq 0 then begin
   readf,lun, x_pos , y_pos
   close,lun
   free_lun, lun
endif else begin
   x_pos = 20 & y_pos = 75
endelse
tlb_position = [x_pos,y_pos]

date_array = intarr(10)
recent_date_array = intarr(20,10)
recent_date_string = strarr(20)

tlb2 = Widget_Base(Title='CME Naming Tool', row=1, MBar=menubase, xoffset = x_pos , yoffset = y_pos, /tlb_move_events)

filer = Widget_Button(menubase, Value='File', /Menu)

widget_select_telescope = Widget_Button(filer, Value='Select Telescope',Event_Pro='select_telescope',/Menu,sensitive = 1)
for i = 0 , n_elements(telescope_array)-1 do begin
select_telescope = Widget_Button(widget_select_telescope, Value=telescope_array[i],uvalue = i)
endfor

widget_reset = Widget_Button(filer, /Separator, Value='Reset',Event_Pro='Reset')
widget_name_another_cme = Widget_Button(filer, /Separator, Value='Name Another CME',Event_Pro='name_another_cme',sensitive=0)

quitter = Widget_Button(filer, /Separator, Value='Exit',Event_Pro='exit_the_prog')

widget_image_menu = Widget_Button(menubase, Value='Image', /Menu,sensitive=0)
difference_image = Widget_Button(widget_image_menu, Value='Normal',Uvalue='Normal',Event_Pro='diff_format')
difference_image = Widget_Button(widget_image_menu, Value='Diff: Running',Uvalue='Diff: Running',Event_Pro='diff_format')
difference_image = Widget_Button(widget_image_menu, Value='Diff: Set Current Image as Background', $
                           Uvalue='Diff: Set Current Image as Background',Event_Pro='diff_format')

vertical_base=widget_base(tlb2,column=1)

drawID_main = Widget_Draw(vertical_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, Retain=0, $
   Expose_Events=1, Event_Pro='Main_window_button_click_events', Button_Events=1,/motion_events,render = 0)
   
widget_image_time_slider = widget_slider(vertical_base, Value=1, event_pro='change_current_image',  $
   minimum = 1, maximum = 2, scroll = 1, drag = 1, sensitive=0, Suppress_Value = 1)
   
define_CME_buttons_base=widget_base(vertical_base,row=1,/align_center)

fullhalo_base = Widget_Base(define_CME_buttons_base, column=1,Frame = 0,/NonExclusive) ;, $
;                /GRID_LAYOUT,/BASE_ALIGN_RIGHT, YSIZE=1)
widget_full_halo_cme = widget_button(fullhalo_base, UNAME='full_halo', $
                VALUE='Full Halo', Event_Pro='set_full_halo_cme', $
                TOOLTIP='Toggle: Set full halo cme on/off.', sensitive=0)



widget_define_CME_name = Widget_Button(define_CME_buttons_base, Value='Define CME Name',Event_Pro='define_cme_name',sensitive=0)
widget_define_CME_representative_image = Widget_Button(define_CME_buttons_base, Value='Define Representative Image',Event_Pro='define_cme_representative_image',sensitive=0)
widget_save_cme = Widget_Button(define_CME_buttons_base, Value='Save to SPI',Event_Pro='save_CME',sensitive=0)

button_and_sliders_base=widget_base(vertical_base,row=1)

  sliderbase = Widget_Base(vertical_base, row=1, frame=1) ; , XPad=0, YPad=0, Frame=0)
  
  widget_botSlider = Widget_Slider(sliderbase, Value=0, Min=0, $
    Max=120, XSize=100,Event_Pro='color_scale_bottom_slider', $
    Title='             bot', Drag=1,sensitive=0,/Suppress_Value)
  widget_topSlider = Widget_Slider(sliderbase, Value=100, Min=0, $
    Max=120, XSize=100, Event_Pro='color_scale_top_slider', $
    Title='             top', Drag=1,sensitive=0, /Suppress_Value)
  widget_gammaSlider = Widget_Slider(sliderbase, Value=10, Min=1, Max=30, $
    Drag=1, XSize=120, /Suppress_Value, Event_Pro='color_scale_gamma_slider', $
    Title='             gamma',sensitive=0)
    
  widget_saturationSlider = Widget_Slider(sliderbase, Value=50, Min=1, Max=200, $
    Drag=1, XSize=120,Suppress_Value=1, Event_Pro='change_the_image_saturation', $
    Title='             sat',sensitive=0)
    
  bot_initial_value = 0
  top_initial_value = 100
  gamma_initial_value = 10
  sat_initial_value = 50
    
  widget_reset_image_controls = Widget_Button(sliderbase, Value='Reset',Event_Pro='reset_image_controls',sensitive=0)


drawID_info = Widget_Draw(vertical_base, XSize=Xsize, YSize=100, Graphics_Level=2, Retain=0, $
Event_Pro='info_Window_button_click_events',render = 0) 


   




;fullhalo_base = Widget_Base(button_and_sliders_base, column=1,Frame = 10,/NonExclusive, $
;                /GRID_LAYOUT,/BASE_ALIGN_RIGHT,YSIZE=10)
;widget_full_halo_cme = widget_button(fullhalo_base, UNAME='full_halo', $
;                VALUE='Full Halo', Event_Pro='set_full_halo_cme', $
;                TOOLTIP='Toggle: Set full halo cme on/off.', sensitive=0)


Widget_Control, tlb2, /Realize
Widget_Control, drawID_main, Get_Value=Main_Window
Widget_Control, drawID_info, Get_Value=info_Window

;
;  End of defining all of the GUI widget stuff.........
;

image_color_palette = OBJ_NEW('IDLgrPalette')
image_color_palette -> loadct,0

;  Process the input diff data array - which contains
;  all of the diff images in a single array.
;  Define an object array (all_diff_image_objects)
;  that will hold all of the images.....

files_directory = ''
filename = ''
pointer_to_list_of_image_files = ptr_new(/allocate_heap)



diff_Image_object = OBJ_NEW('IDLgrImage')
;diff_Image_object -> SetProperty, DATA = diff_data
diff_Image_object -> SetProperty, location = [0 , 0]
diff_Image_object -> setproperty, palette = image_color_palette
diff_Image_object -> setproperty, alpha_channel = 1.0
diff_Image_object -> setproperty, blend_function = [3,4]

Main_Model -> Add, diff_image_object

;  Calculate a Julian time for each of the images (from the filenames).

pointer_to_image_time_Julian = ptr_new(/allocate_heap)
pointer_to_full_time_string = ptr_new(/allocate_heap)

Main_Window->Draw, Main_View
info_Window->Draw, info_View


;  The center of the Sun is defined as a yellow dot.
;  The dot is made by a 9 pixel square colored yellow (3 x 3 pixels)
; make the Sun position totally wrong - to make sure that the user
; realizes it is wrong and sets it manually....

center_of_SunX = 100.
center_of_SunY = 100.

;half of 3 pixels:
half_sunsize = 1.5

xvals = [center_of_sunX-half_sunsize,center_of_sunX+half_sunsize,center_of_sunX+half_sunsize,center_of_sunX-half_sunsize]
yvals = [center_of_sunY+half_sunsize,center_of_sunY+half_sunsize,center_of_sunY-half_sunsize,center_of_sunY-half_sunsize]
marker_Sun_position = Obj_New("IDLgrPolygon", xvals, yvals, Color=color_yellow)
Main_Model -> Add, marker_Sun_position

marker_Sun_position -> SetProperty, Hide = 1

next_xvals = [0,0,20]
next_yvals = [0,20,10]

goto_next_image_icon = Obj_New("IDLgrPolygon", next_xvals + 485, next_yvals + 256, $
  Color=[255,255,255],hide=1, alpha_channel=0.5)
Main_Model -> Add, goto_next_image_icon

prev_xvals = [20,0,20]
prev_yvals = [0,10,20]

goto_previous_image_icon = Obj_New("IDLgrPolygon", prev_xvals + 10, prev_yvals + 256, $
  Color=[255,255,255],hide=1, alpha_channel=0.5)
Main_Model -> Add, goto_previous_image_icon




;be sure to use reals, unless you only want integers...
scaling_factor = 1.
rotation = 0.
datetime_string = ''



clock_angle_model = OBJ_NEW('IDLgrModel') 

clock_angle_model_location = [0.,0.]

clock_angle_model -> translate,255,255,0.0
clock_angle_model_location = [255,255]

clock_angle_model->GetProperty, transform = clock_angle_model_original_transform

; Marker for the 'clock angle' line

clock_angle_marker_line_x = [0.,0.]
clock_angle_marker_line_y = [0.,ysize/2.5]
clock_angle_marker_line = Obj_New("IDLgrPolyline", clock_angle_marker_line_x , clock_angle_marker_line_y, $
name = 'clock_angle', Color=color_red, thick=2, alpha_channel = 0.7)
clock_angle_model -> Add, clock_angle_marker_line

clock_angle_model -> SetProperty, Hide = 1

Courier_small = Obj_New('IDLgrFont', 'Courier', Size=14)
Courier_large = Obj_New('IDLgrFont', 'Courier', Size=20)

rotate_x = 0.
rotate_y = 0.

clock_angle_string_value = ' 0'
;Displays the clock angle degrees near the clock angle line
clock_angle_string = obj_new("idlgrtext", strings = clock_angle_string_value,color=color_red, $
font=Courier_small, /onglass, locations = [3,(ysize/3.) + 23],hide=0,fill_background=1,fill_color=color_fill, $
alpha_channel = 0.7)

clock_angle_model -> Add, clock_angle_string

; Model and Marker for the 'leading edge' line

leading_edge_string_value = ' 0'
leading_edge_model = OBJ_NEW('IDLgrModel')

;marker for leading edge line:
leading_edge_marker_line_x = [-xsize/20.,+xsize/20.]
leading_edge_marker_line_y = [ysize/7.,ysize/7.]
leading_edge_marker_line = Obj_New("IDLgrPolyline", leading_edge_marker_line_x , leading_edge_marker_line_y, $
name = 'leading_edge', Color=color_yellow, thick=2,alpha_channel = 0.7)
leading_edge_model -> Add, leading_edge_marker_line

leading_edge_model -> SetProperty, Hide = 0

;Displays the value of the leading edge near the leading edge line
leading_edge_string = obj_new("idlgrtext", strings = leading_edge_string_value,color=color_yellow, $
font=Courier_small, /onglass, locations = [+xsize/20.,(ysize/7.) + 10],hide=0,fill_background=1,fill_color=color_fill, $
alpha_channel = 0.7)


leading_edge_model -> Add, leading_edge_string

clock_angle_model -> Add, leading_edge_model

;  Here the clock_angle_model is added to the main_model....

Main_Model -> Add, clock_angle_model


Main_Window->Draw, Main_View

latest_click_X = 0
latest_click_Y = 0
latest_release_X = 0
latest_release_Y = 0
click_and_drag = 0
previous_click_X = 0
previous_click_Y = 0

 
current_image_number = 0

which_diff = 1 ; start out with running differences....
background_image_number = -1

cme_name = ''
CME_name_string = obj_new("idlgrtext", strings = ' ',color=[255,255,255], $
locations = [xsize/2,ysize_info - 40], alignment= 0.5, font=Courier_large)
cme_image_name = ''
cme_datetime = ''
cme_catalog_datetime_for_fits_header = ''
CME_image_string = obj_new("idlgrtext", strings = ' ',color=[255,255,255], $
locations = [xsize/2,ysize_info - 70], alignment = 0.5, font=Courier_small)

info_Model -> add, CME_name_string
info_Model -> add, CME_image_string

info_Window -> draw, info_View

which_telescope = ''
telescope_field_of_view_Rs = 30.
image_size_pixels = 512.
conversion_factor_Rs_per_pixel = telescope_field_of_view_Rs * 2. / image_size_pixels

rsun = 950.   ;set initial value, solar radius from L1 in arcseconds
pixel_scale = 1.
clock_angle_degrees = 0
leading_edge_in_rs = 0
telescope_code = ''

image_type = 'fits'
image_fits_header = ptr_new(/allocate_heap)
image_saturation_value = 50.
the_action = 0

full_halo = 0

i_images_are_loaded = 0
i_cme_has_been_named = 0
i_representative_image_has_been_chosen = 0
date_array_int = intarr(10)
start_date = intarr(5)
end_date = intarr(5)

;initialize julian start & end datetimes as doubles
start_julian = -1.0D
end_julian = -1.0D

start_str=''  ;start date/time in a string
end_str=''    ;end date/time in a string

ut_string_object = OBJ_NEW('IDLgrText','')
ut_string_object -> setproperty, location = [xsize - 10, ysize - 20]
ut_string_object -> setproperty, alignment = 1.0
ut_string_object -> setproperty, color = color_red
ut_string_object -> setproperty, font = Courier_small
Main_Model -> Add, ut_string_object

telescope_string_object = OBJ_NEW('IDLgrText','')
telescope_string_object -> setproperty, location = [10, ysize - 20]
telescope_string_object -> setproperty, alignment = 0.0
telescope_string_object -> setproperty, color = color_red
telescope_string_object -> setproperty, font = Courier_small
telescope_string_object -> setproperty, strings = which_telescope
Main_Model -> Add, telescope_string_object


top_title_background_xvals = [0,0,511,511]
top_title_background_yvals = [488,511,511,488]

top_title_background = Obj_New("IDLgrPolygon", top_title_background_xvals, top_title_background_yvals, $
  Color=[50,50,50],hide=0, alpha_channel=0.5)
Main_Model -> Add, top_title_background

;Displays clock angle degrees in lower left corner on image.
clock_angle_deg_object = OBJ_NEW('IDLgrText','')
clock_angle_deg_object -> setproperty, location = [10, ysize - 490]
clock_angle_deg_object -> setproperty, alignment = 0.0
clock_angle_deg_object -> setproperty, color = color_red
clock_angle_deg_object -> setproperty, font = Courier_small
clock_angle_deg_object -> setproperty, strings = clock_angle_string_value
Main_Model -> Add, clock_angle_deg_object
clock_angle_deg_object->SetProperty, hide = 1

;Displays CME leading edge value in lower right corner on image.
leading_edge_txt_object = OBJ_NEW('IDLgrText','')
leading_edge_txt_object -> setproperty, location = [xsize - 10, ysize - 490]
leading_edge_txt_object -> setproperty, alignment = 1.0
leading_edge_txt_object -> setproperty, color = color_yellow
leading_edge_txt_object -> setproperty, font = Courier_small
leading_edge_txt_object -> setproperty, strings = leading_edge_string_value
Main_Model -> Add, leading_edge_txt_object
leading_edge_txt_object->SetProperty, hide = 1

telescope_folder = ''
all_image_data_list = list()

info2 = $
       { tlb2:tlb2, $
         Main_Window:Main_Window, $      
         Main_View:Main_View, $            
         Main_Model:Main_Model, $            
         drawID_main:drawID_main, $
         info_Window:info_Window, $      
         info_View:info_View, $            
         info_Model:info_Model, $            
         drawID_info:drawID_info, $
         tlb_position:tlb_position, $
         clock_angle_model:clock_angle_model, $            
         clock_angle_model_location:clock_angle_model_location, $
         clock_angle_marker_line:clock_angle_marker_line, $
         clock_angle_degrees:clock_angle_degrees, $
         clock_angle_string:clock_angle_string, $
         clock_angle_deg_object:clock_angle_deg_object, $
         rotate_x:rotate_x, $
         rotate_y:rotate_y, $
         rsun:rsun, $
         pixel_scale:pixel_scale, $
         leading_edge_model:leading_edge_model, $
         leading_edge_marker_line:leading_edge_marker_line, $ 
         leading_edge_in_rs: leading_edge_in_rs, $
         leading_edge_string:leading_edge_string, $
         leading_edge_txt_object:leading_edge_txt_object, $       
         files_directory:files_directory, $
         pointer_to_list_of_image_files:pointer_to_list_of_image_files, $
         number_of_images:number_of_images, $
         image_type:image_type, $
         image_saturation_value:image_saturation_value, $
         diff_image_object:diff_image_object, $
         image_fits_header:image_fits_header, $
         latest_click_X:latest_click_X, $
         latest_click_Y:latest_click_Y, $
         latest_release_X:latest_release_X, $
         latest_release_Y:latest_release_Y, $
         previous_click_X:previous_click_X, $
         previous_click_Y:previous_click_Y, $
         image_color_palette:image_color_palette, $
         marker_Sun_position:marker_Sun_position, $                       
         click_and_drag:click_and_drag, $
;         widget_select_dates:widget_select_dates, $
;         widget_select_recent:widget_select_recent, $
         ;widget_get_input_images:widget_get_input_images, $
         widget_image_time_slider:widget_image_time_slider, $
;         widget_alpha_slider:widget_alpha_slider, $
         widget_botSlider:widget_botSlider, $
         widget_topSlider:widget_topSlider, $
         widget_gammaSlider:widget_gammaSlider, $
         widget_saturationSlider:widget_saturationSlider, $
         widget_image_menu:widget_image_menu, $
         widget_full_halo_cme:widget_full_halo_cme, $
         ;widget_set_sun_position:widget_set_sun_position, $
         full_halo: full_halo, $
         pointer_to_image_time_Julian:pointer_to_image_time_Julian, $
         pointer_to_full_time_string:pointer_to_full_time_string, $
         ut_string_object:ut_string_object, $
         current_image_number:current_image_number, $
         center_of_sunX:center_of_sunX, $
         center_of_sunY:center_of_sunY, $
         scaling_factor:scaling_factor, $
         rotation:rotation, $
         background_center_of_sunX:center_of_sunX, $
         background_center_of_sunY:center_of_sunY, $
         background_rotation:rotation, $
         which_diff:which_diff, $
         background_image_number:background_image_number, $
         which_telescope:which_telescope, $
         telescope_field_of_view_Rs:telescope_field_of_view_Rs, $
         image_size_pixels:image_size_pixels, $
         conversion_factor_Rs_per_pixel:conversion_factor_Rs_per_pixel, $
         the_action:the_action, $
         telescope_string_object:telescope_string_object, $
         telescope_code:telescope_code, $
         cme_name:cme_name, $
         cme_name_string:cme_name_string, $
         cme_image_name:cme_image_name, $
         cme_image_string:cme_image_string, $
         cme_datetime:cme_datetime, $
         widget_define_CME_name:widget_define_CME_name, $
         widget_define_CME_representative_image:widget_define_CME_representative_image, $
         widget_save_cme:widget_save_cme, $
         i_images_are_loaded:i_images_are_loaded, $
         i_cme_has_been_named:i_cme_has_been_named, $
         i_representative_image_has_been_chosen:i_representative_image_has_been_chosen, $
         date_array_int:date_array_int, $
         recent_date_array:recent_date_array, $
         start_date:start_date, $
         end_date:end_date, $
         start_julian:start_julian, $
         end_julian:end_julian, $
         start_str:start_str, $
         end_str:end_str, $
         cment_preferences_file:cment_preferences_file, $
         sep:sep, $
         telescope_array:telescope_array, $
         image_in_folder_array:image_in_folder_array, $
         image_in_root:image_in_root, $
         image_out_location:image_out_location, $
         max_interval_in_days:max_interval_in_days, $
         input_file:input_file, $
         widget_select_telescope:widget_select_telescope, $
         clock_angle_model_original_transform:clock_angle_model_original_transform, $
         telescope_folder:telescope_folder, $
         cme_catalog_datetime_for_fits_header:cme_catalog_datetime_for_fits_header, $
         widget_name_another_cme:widget_name_another_cme, $
         i_confirm_cme_save_to_spi:i_confirm_cme_save_to_spi, $
         goto_next_image_icon:goto_next_image_icon, $
         goto_previous_image_icon:goto_previous_image_icon, $
         top_title_background:top_title_background, $
         all_image_data_list:all_image_data_list, $
         bot_initial_value:bot_initial_value, $
         top_initial_value:top_initial_value, $
         gamma_initial_value:gamma_initial_value, $
         sat_initial_value:sat_initial_value, $
         widget_reset_image_controls:widget_reset_image_controls}
         

Widget_Control, tlb2, Set_UValue=info2, /No_Copy

XManager, 'Main_View', tlb2, Cleanup='stp_Cleanup2', /No_Block, $
   Event_Handler='top_level_base_events', Group_Leader=groupLeader
END

