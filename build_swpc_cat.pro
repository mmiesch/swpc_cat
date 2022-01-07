; $Id:build_swpc_cat.pro 2012-04-20  anewman $
;+
; The build script for swpc_cat; it makes a SAVE file that can be
; executed with the IDL Virtual Machine or in IDL Runtime.
; It also makes a runtime distribution, using make_rt. 
; 
; This file should be in same directory as the swpc_cat source code.
; At the top of this script edit the 2 locations for source_code_folder and 
; executable folder, ie:
; source_code_folder='C:\Users\george.millward.SWPC\IDL_projects\swpc_cat_production\'
; executable_folder = 'C:\Users\george.millward.SWPC\Desktop\CAT_DEMO'
; 
; This build script should be run from the IDL command line (see the example)
; 
; @examples
; <pre>
; IDL> @build_swpc_cat
; </pre>
; 
;
; When the build script has created the make_rt run time distribution, you may want to include
; 2 more files: the CME.jpg file (used to create a desktop icon) and the swpc_cat.in file. These
; can be copied from the source file directory.
;
; @categories swpc_cat, build script, projects
;-

; Clear out memory.
.reset_session

;source_code_folder='C:\Users\george.millward.SWPC\IDL_projects\swpc_cat_production\'
;executable_folder = 'C:\Users\george.millward.SWPC\Desktop\CAT_DEMO2\'
source_code_folder='C:\Users\george.millward\IDLWorkspace\analysis_tools\swpc_cat\'
executable_folder = 'C:\Users\george.millward\Desktop\swpc_cat_new_4\'

path_sep = PATH_SEP(/search_path)
!path = EXPAND_PATH('+'+source_code_folder) + path_sep + !path

start_path = !path

unresolved = ''

; Switch to the project directory.
cd, current=start_dir
cd, source_code_folder

.compile swpc_cat_fsc_inputfield

; Compile swpc_cat: this is the routine called in the SAVE file.
.compile swpc_cat


; Resolve dependencies -- compile what the compiled routines call. Note that
; cat routines referenced by a string ( e.g. event handlers & the cleanup 
; routine are called by XMANAGER) need to be resolved, as well as any classes
; that are written in IDL source code.
resolve_all  , $
   class=['swpc_cat_progressbar', 'idlgrmodel', 'swpc_cat_FSC_INPUTFIELD', 'trackball'], $
   resolve_either=['swpc_cat_FSC_InputField_Event_Handler'], $ ;, $ ; these are DFanning's event handlers
   /CONTINUE_ON_ERROR, UNRESOLVED=unresolved
   
; Write the SAVE file, naming it after the routine to be called, swpc_cat.
file = 'swpc_cat.sav'
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

MAKE_RT,'swpc_cat',executable_folder, $
  savefile=source_code_folder + 'swpc_cat.sav', $
  /OVERWRITE,/VM,/WIN32,logfile=executable_folder + '\swpc_cat\make_rt_log.txt'

FILE_COPY, source_code_folder + 'swpc_cat.in', executable_folder + '\swpc_cat\swpc_cat.in'
FILE_MKDIR, executable_folder + '\swpc_cat\swpc_cat_data'
FILE_COPY, source_code_folder + '\swpc_cat_data\Earth_position_2011_to_2019.txt', $
  executable_folder + '\swpc_cat\swpc_cat_data\Earth_position_2011_to_2019.txt'
FILE_COPY, source_code_folder + '\swpc_cat_data\Earth_position_2011_to_2030.txt', $
  executable_folder + '\swpc_cat\swpc_cat_data\Earth_position_2011_to_2030.txt'

; Revert to the start directory and reset the original path.
cd, start_dir
!path = start_path


