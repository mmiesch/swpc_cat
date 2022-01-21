; $Id:build_cment.pro 2012-04-20  anewman $
;+
; The build script for cment; it makes a SAVE file that can be
; executed with the IDL Virtual Machine or in IDL Runtime.
; It also makes a runtime distribution, using make_rt. 
; 
; This file should be in same directory as the cment source code.
; At the top of this script edit the 2 locations for source_code_folder and
; executable folder, ie:
; source_code_folder='C:\Users\george.millward.SWPC\IDL_projects\cat\'
; executable_folder = 'C:\Users\george.millward.SWPC\Desktop\CMENT_DEMO\'
; 
; This build script should be run from the IDL command line (see the example)
; 
; @examples
; <pre>
; IDL> @build_cment
; </pre>
; 
;
; When the build script has created the make_rt run time distribution, you may want to include
; 2 more files: the CME.jpg file (used to create a desktop icon) and the cment.in file. These
; can be copied from the source file directory.
;
; @categories cment, build script, projects
;-

; Clear out memory.
.reset_session

source_code_folder='C:\Users\george.millward.SWPC\IDL_projects\cat\'
executable_folder = 'C:\Users\george.millward.SWPC\Desktop\CMENT_DEMO2\'

path_sep = PATH_SEP(/search_path)
!path = EXPAND_PATH('+'+source_code_folder) + path_sep + !path

start_path = !path

unresolved = ''

; Switch to the project directory.
cd, current=start_dir
cd, source_code_folder

; Compile cment this is the routine called in the SAVE file.
.compile cment

;Add other *.pro files to be compiled here ( do not include build_cment):
.compile cgcolor
.compile cgsnapshot
.compile check_dates
.compile check_fits
.compile check_image
.compile checksum32
.compile datebox_ghm_no_basedate
.compile daycnv
.compile dbltostr
.compile decomposedcolor
.compile detabify
.compile error_message
.compile fits_add_checksum
.compile fits_ascii_encode
.compile fits_test_checksum
.compile fsc_inputfield
.compile fxaddpar
.compile fxparpos
.compile get_date
.compile gettok
.compile host_to_ieee
.compile ieee_to_host
.compile is_ieee_big
.compile leapyear
.compile mkhdr
.compile month_cnv
.compile monthlen
.compile mrd_hread
.compile n_bytes
.compile pickcolorname
.compile progressbar__define
.compile readfits
.compile sxaddpar
.compile sxdelpar
.compile sxpar
.compile systimex
.compile valid_num
.compile validate_date
.compile writefits
.compile strsplit



; Resolve dependencies -- compile what the compiled routines call. Note that
; cment routines referenced by a string ( e.g. event handlers & the cleanup 
; routine are called by XMANAGER) need to be resolved, as well as any classes
; that are written in IDL source code.
resolve_all, $
   class=['progressbar', 'idlgrmodel', 'FSC_INPUTFIELD'], $
   resolve_either=['FSC_InputField_Event_Handler', $ ; these are DFanning's event handlers
                    'Example_Event','PickColorName_Select_color', $
                    'PickColorName_Buttons', $
                    'top_level_base_events', 'select_telescope', $ ;SWPC code event handlers
                    'Reset', 'exit_the_prog','diff_format', $
                    'Main_window_button_click_events','set_full_halo_cme', $
                    'define_cme_name','define_cme_representative_image','save_CME', $
                    'info_Window_button_click_events','change_current_image', $
                    'change_transparency','color_scale_bottom_slider', $
                    'color_scale_top_slider','color_scale_gamma_slider', $
                    'change_the_image_saturation', $
                    'stp_Cleanup2'], $   ;cleanup routine called by XMANAGER
   /CONTINUE_ON_ERROR, UNRESOLVED=unresolved
   
; Write the SAVE file, naming it after the routine to be called, cment.
file = 'cment.sav'
save, /routines, filename=file, /verbose
if FILE_TEST(file) eq 1 then $
   print, 'File "' + file + '" created.' $
else $
   print, 'Error: Save file not created.'

print,''
print,'unresolved in next line will be empty, if all dependencies are resolved.'
print,'unresolved=',unresolved
print,''
print,'NOTE: check for any syntax errors in this build file, especially if'
print,'it has been modified since the last run. Just scroll up in the'
print,'IDL console and look for errors in red, usually above the SAVE: sentences.'

print,''
print,'Making the Real-time distribution now. It will take a minute or so; watch for message that it finished...'

print,'executable_folder=',executable_folder
print,'source_code_folder=',source_code_folder

if FILE_TEST(executable_folder,/DIRECTORY) eq 0 then $ ;create the output dir, if not exist
   FILE_MKDIR,executable_folder

MAKE_RT,'cment',executable_folder, $
  savefile=source_code_folder + 'cment.sav', $
  /OVERWRITE,/VM,/WIN32,logfile=executable_folder + '\cment\make_rt_log.txt'
  
FILE_COPY, source_code_folder + 'cment.in', executable_folder + '\cment\cment.in'
FILE_COPY, source_code_folder + 'CME.jpg', executable_folder + '\cment\CME.jpg'

; Revert to the start directory and reset the original path.
cd, start_dir
!path = start_path


