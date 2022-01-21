
;+
; NAME:
;  TEXTBOX
;
; PURPOSE:
;
;  This function allows the user to type some text in a
;  pop-up dialog widget and have it returned to the program.
;  This is an example of a Pop-Up Dialog Widget.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;  Utility, Widgets
;
; CALLING SEQUENCE:
;
;  thetext = TextBox()
;
; INPUTS:
;
;  None.
;
; KEYWORD PARAMETERS:
;
;  CANCEL: An output parameter. If the user kills the widget or clicks the Cancel
;       button this keyword is set to 1. It is set to 0 otherwise. It
;       allows you to determine if the user canceled the dialog without
;       having to check the validity of the answer.
;
;       theText = TextBox(Title='Provide Phone Number...', Label='Number:', Cancel=cancelled)
;       IF cancelled THEN Return
;
;  GROUP_LEADER: The widget ID of the group leader of this pop-up
;       dialog. This should be provided if you are calling
;       the program from within a widget program:
;
;          thetext = TextBox(Group_Leader=event.top)
;
;       If a group leader is not provided, an unmapped top-level base widget
;       will be created as a group leader.
;
;  LABEL: A string the appears to the left of the text box.
;
;  TITLE:  The title of the top-level base. If not specified, the
;       string 'Provide Input:' is used by default.
;
;  VALUE: A string variable that is the intial value of the textbox. By default, a null string.
;
;  XSIZE: The size of the text widget in pixel units. By default, 200.
;
; OUTPUTS:
;
;  theText: The string of characters the user typed in the
;       text widget. No error checking is done.
;
; RESTRICTIONS:
;
;  The widget is destroyed if the user clicks on either button or
;  if they hit a carriage return (CR) in the text widget. The
;  text is recorded if the user hits the ACCEPT button or hits
;  a CR in the text widget.
;
; MODIFICATION HISTORY:
;
;  Written by: David W. Fanning, December 20, 2001.
;  Added VALUE keyword to set the initial value of the text box. 4 Nov 2002. DWF.
;-
;
;******************************************************************************************;
;  Copyright (c) 2008, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;


PRO datebox_ghm_CenterTLB, tlb

   ; This utility routine centers the TLB.

Device, Get_Screen_Size=screenSize
IF screenSize[0] GT 2000 THEN screenSize[0] = screenSize[0]/2 ; Dual monitors.
xCenter = screenSize(0) / 2
yCenter = screenSize(1) / 2

geom = Widget_Info(tlb, /Geometry)
xHalfSize = geom.Scr_XSize / 2
yHalfSize = geom.Scr_YSize / 2

Widget_Control, tlb, XOffset = xCenter-xHalfSize, $
   YOffset = yCenter-yHalfSize

END ;-----------------------------------------------------

PRO set_directory, event

; This event handler responds to the droplist widget where the USER is specifying the data location

Widget_Control, event.top, Get_UValue=info

CASE event.index OF
  0: BEGIN
      setenv,'STEREO_IMG_ROOT=\\Nas-bes\wsa-enlil\Images\'
      setenv,'STEREO_ROOT_STR=ops'
     END
  1: BEGIN
      setenv,'STEREO_IMG_ROOT=\\Nas-bes\v_v\STEREO_data\beacon'
      setenv,'STEREO_ROOT_STR=v_v'
     END
  2: BEGIN
      setenv,'STEREO_IMG_ROOT=\\Nas-bes\wsa-enlil\Images\beacon'
      setenv,'STEREO_ROOT_STR=sci'
     END
ENDCASE

END

PRO datebox_ghm_no_basedate_Event, event

   ; This event handler responds to all events. Widget
   ; is always destoyed. The text is recorded if ACCEPT
   ; button is selected or user hits CR in text widget.

Widget_Control, event.top, Get_UValue=info

CASE event.ID OF
   info.cancelID: Widget_Control, event.top, /Destroy
   ELSE: BEGIN

         ; Get the text and store it in the pointer location.

;      Widget_Control, info.textID, Get_Value=theText
      start_year = info.start_yearID -> Get_Value()
      start_month = info.start_monthID -> Get_Value()
      start_day = info.start_dayID -> Get_Value()
      start_hour = info.start_hourID -> Get_Value()
      start_minute = info.start_minuteID -> Get_Value()
      end_year = info.end_yearID -> Get_Value()
      end_month = info.end_monthID -> Get_Value()
      end_day = info.end_dayID -> Get_Value()
      end_hour = info.end_hourID -> Get_Value()
      end_minute = info.end_minuteID -> Get_Value()
      
      IF N_Elements(start_hour) eq 0 then start_hour = 0
      IF N_Elements(end_hour) eq 0 then end_hour = 0
      IF N_Elements(start_minute) eq 0 then start_minute = 0
      IF N_Elements(end_minute) eq 0 then end_minute = 0
      
      if N_Elements(start_year) eq 0 or N_Elements(start_month) eq 0 or $
      N_Elements(start_day) eq 0 or $
      N_Elements(end_year) eq 0 or N_Elements(end_month) eq 0 or $
      N_Elements(end_day) eq 0 then begin
      
      msg1 = 'Start And End Dates Cannot Contain Blanks'
      
      ok = DIALOG_MESSAGE(msg1,/ERROR,/center)
      
      endif else begin

      (*info.ptr).start_year = start_year[0]
      (*info.ptr).start_month = start_month[0]
      (*info.ptr).start_day = start_day[0]
      (*info.ptr).start_hour = start_hour[0]
      (*info.ptr).start_minute = start_minute[0]
      (*info.ptr).end_year = end_year[0]
      (*info.ptr).end_month = end_month[0]
      (*info.ptr).end_day = end_day[0]
      (*info.ptr).end_hour = end_hour[0]
      (*info.ptr).end_minute = end_minute[0]
      (*info.ptr).cancel = 0
      Widget_Control, event.top, /Destroy
      endelse
      ENDCASE
ENDCASE
END ;-----------------------------------------------------



FUNCTION datebox_ghm_no_basedate, Title=title, Label=label, Cancel=cancel, $
   Group_Leader=groupleader, XSize=xsize, Value=value, $
   date_array = date_array


   ; Return to caller if there is an error. Set the cancel
   ; flag and destroy the group leader if it was created.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = Dialog_Message(!Error_State.Msg)
   IF destroy_groupleader THEN Widget_Control, groupleader, /Destroy
   cancel = 1
   RETURN, ""
ENDIF

   ; Check parameters and keywords.

IF N_Elements(title) EQ 0 THEN title = 'Provide Input:'
IF N_Elements(label) EQ 0 THEN label = ""
IF N_Elements(value) EQ 0 THEN value = ""
IF N_Elements(xsize) EQ 0 THEN xsize = 200

   ; Provide a group leader if not supplied with one. This
   ; is required for modal operation of widgets. Set a flag
   ; for destroying the group leader widget before returning.

IF N_Elements(groupleader) EQ 0 THEN BEGIN
   groupleader = Widget_Base(Map=0)
   Widget_Control, groupleader, /Realize
   destroy_groupleader = 1
ENDIF ELSE destroy_groupleader = 0

IF N_Elements(date_array) eq 0 then begin

; figure out the current time so we can provide some decent defaults for the coming
; date/time text boxes.............

current_julian_time = systime(/julian,/utc)
CALDAT, current_julian_time, current_Month, current_Day, current_Year, current_Hour, current_Minute
;print,'right now ', current_Month, current_Day, current_Year, current_Hour, current_Minute

julian_this_hour = JULDAY(current_Month, current_Day, current_Year, current_Hour) 
julian_this_hour_ending = julian_this_hour
julian_this_hour_starting = julian_this_hour - 1.d

CALDAT, julian_this_hour_ending, current_Month_ending, current_Day_ending, current_Year_ending, current_Hour_ending
CALDAT, julian_this_hour_starting, current_Month_starting, current_Day_starting, current_Year_starting, current_Hour_starting

current_year_starting = strcompress(string(current_year_starting),/remove_all)
current_month_starting = strcompress(string(current_month_starting),/remove_all)
current_day_starting = strcompress(string(current_day_starting),/remove_all)
current_hour_starting = strcompress(string(current_hour_starting),/remove_all)
current_minute_starting = value

current_year_ending = strcompress(string(current_year_ending),/remove_all)
current_month_ending = strcompress(string(current_month_ending),/remove_all)
current_day_ending = strcompress(string(current_day_ending),/remove_all)
current_hour_ending = strcompress(string(current_hour_ending),/remove_all)
current_minute_ending = value

endif else begin

;date_array = [start_year,start_month,start_day,start_hour,start_minute, $
;              end_year,end_month,end_day,end_hour,end_minute, $
;              base_year,base_month,base_day,base_hour,base_minute]

current_year_starting = date_array[0]
current_month_starting = date_array[1]
current_day_starting = date_array[2]
current_hour_starting = date_array[3]
current_minute_starting = date_array[4]
current_year_ending = date_array[5]
current_month_ending = date_array[6]
current_day_ending = date_array[7]
current_hour_ending = date_array[8]
current_minute_ending = date_array[9]

endelse

;

   ; Create modal base widget.

tlb = Widget_Base(Title=title, Column=1, /Modal, $
   /Base_Align_Center, Group_Leader=groupleader)

   ; Create the rest of the widgets.

labelbase = Widget_Base(tlb, Row=6)
row1base = Widget_Base(labelbase, Row=1)
row2base = Widget_Base(labelbase, Row=1)
row3base = Widget_Base(labelbase, Row=1)
row4base = Widget_Base(labelbase, Row=1)


start_label = Widget_Label(row1base, Value='Start Date [Y M D  H M]')
;start_yearID = Widget_Text(row2base, /Editable, XSize=4, Value=current_Year)
;start_monthID = Widget_Text(row2base, /Editable, XSize=2, Value=current_Month)
;start_dayID = Widget_Text(row2base, /Editable, XSize=2, Value=current_Day)
;start_hourID = Widget_Text(row2base, /Editable, XSize=2, Value=current_Hour)
;start_minuteID = Widget_Text(row2base, /Editable, XSize=2, Value=current_minute)
start_yearID = FSC_INPUTFIELD(row2base, Value=current_Year_starting, /IntegerValue, digits = 4 , xsize = 4 , title='')
start_monthID = FSC_INPUTFIELD(row2base, Value=current_Month_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')
start_dayID = FSC_INPUTFIELD(row2base, Value=current_Day_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')
start_hourID = FSC_INPUTFIELD(row2base, Value=current_Hour_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')
start_minuteID = FSC_INPUTFIELD(row2base, Value=current_Minute_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')

end_label = Widget_Label(row3base, Value='End Date [Y M D  H M]')
;end_yearID = Widget_Text(row4base, /Editable, XSize=4, Value=current_Year_ending)
;end_monthID = Widget_Text(row4base, /Editable, XSize=2, Value=current_Month_ending)
;end_dayID = Widget_Text(row4base, /Editable, XSize=2, Value=current_Day_ending)
;end_hourID = Widget_Text(row4base, /Editable, XSize=2, Value=current_Hour_ending)
;end_minuteID = Widget_Text(row4base, /Editable, XSize=2, Value=current_minute_ending)
end_yearID = FSC_INPUTFIELD(row4base, Value=current_Year_ending, /IntegerValue, digits = 4 , xsize = 4 , title='')
end_monthID = FSC_INPUTFIELD(row4base, Value=current_Month_ending, /IntegerValue, digits = 2 , xsize = 2 , title='')
end_dayID = FSC_INPUTFIELD(row4base, Value=current_Day_ending, /IntegerValue, digits = 2 , xsize = 2 , title='')
end_hourID = FSC_INPUTFIELD(row4base, Value=current_Hour_ending, /IntegerValue, digits = 2 , xsize = 2 , title='')
end_minuteID = FSC_INPUTFIELD(row4base, Value=current_Minute_ending, /IntegerValue, digits = 2 , xsize = 2 , title='')

buttonBase = Widget_Base(tlb, Row=1)
cancelID = Widget_Button(buttonBase, Value='Cancel')
acceptID = Widget_Button(buttonBase, Value='Accept')

  ; Center the widgets on display.

datebox_ghm_CenterTLB, tlb
Widget_Control, tlb, /Realize

   ; Create a pointer for the text the user will type into the program.
   ; The cancel field is set to 1 to indicate that the user canceled
   ; the operation. Only if a successful conclusion is reached (i.e.,
   ; a Carriage Return or Accept button selection) is the cancel field
   ; set to 0.

ptr = Ptr_New({start_year:"", $
               start_month:"", $
               start_day:"", $
               start_hour:"", $
               start_minute:"", $
               end_year:"", $
               end_month:"", $
               end_day:"", $
               end_hour:"", $
               end_minute:"", $
               cancel:1})

   ; Store the program information:

info = {ptr:ptr, $
        start_yearID:start_yearID, $
        start_monthID:start_monthID, $
        start_dayID:start_dayID, $
        start_hourID:start_hourID, $
        start_minuteID:start_minuteID, $
        end_yearID:end_yearID, $
        end_monthID:end_monthID, $
        end_dayID:end_dayID, $
        end_hourID:end_hourID, $
        end_minuteID:end_minuteID, $
        cancelID:cancelID}

Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Blocking or modal widget, depending upon group leader.

XManager, 'datebox_ghm_no_basedate', tlb

   ; Return from block. Return the text to the caller of the program,
   ; taking care to clean up pointer and group leader, if needed.
   ; Set the cancel keyword.

start_year = (*ptr).start_year
start_month = (*ptr).start_month
start_day = (*ptr).start_day
start_hour = (*ptr).start_hour
start_minute = (*ptr).start_minute
end_year = (*ptr).end_year
end_month = (*ptr).end_month
end_day = (*ptr).end_day
end_hour = (*ptr).end_hour
end_minute = (*ptr).end_minute

if start_minute eq 'NULLVALUE' then start_minute = '0'
if end_minute eq 'NULLVALUE' then end_minute = '0'


date_array = [start_year,start_month,start_day,start_hour,start_minute, $
              end_year,end_month,end_day,end_hour,end_minute]

theText = date_array
cancel = (*ptr).cancel
Ptr_Free, ptr
IF destroy_groupleader THEN Widget_Control, groupleader, /Destroy

RETURN, theText
END ;-----------------------------------------------------
