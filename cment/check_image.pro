FUNCTION check_image, fits_file

   ;Identify bad input data (e.g., input FITS files that are corrupted)
   ;If the file is corrupted, flag it, by returning 0
   ;
   ;This is "Trust-but-verify code; bad images should have been deleted
   ;upstream of this program, so that this program would never see them.
   ;But just in case that did not happen...we do some checking here.
   
   compile_opt idl2


   im_data = READFITS(fits_file,header,/silent)
   ;if im_data eq -1  i.e, if error encountered with reading the FITS file
   if size(im_data,/n_dimensions) eq 0 then begin
      If im_data eq -1 then begin
         msg1 = 'There was a problem reading FITS file ' + fits_file + '.'
         msg2 = 'Excluding this image from the list of images to view...'
         ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
         return,0
      endif
   endif

   ;check for correct size of the image; the 2 dimensions 
   ;must be the same size, or the file is corrupted
   ;Example bad images of this type, used in testing:
   ;   X:\Images\20110819\LASCO_C2\20110819_143322_C2.fts
   ;   X:\Images\20110819\LASCO_C2\20110819_144128_C2.fts
   
   naxis1 = fix(sxpar(header, 'NAXIS1'))
   naxis2 = fix(sxpar(header, 'NAXIS2'))
   
   if naxis1 ne naxis2 then begin
      dim2 = '(' + strtrim(string(naxis1),2) + 'x' + strtrim(string(naxis2),2) + ')'
      msg1 = 'There is a dimension problem ' + dim2 + ' with image ' + fits_file
      msg2 = 'Excluding this image from the list of images to view...'
      ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
      return,0
   endif
   
   return,1  ;good image
 
END