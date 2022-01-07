;----------------------------------------------------------------
;
;       NOAA/SWPC CME ANALYSIS TOOL (CAT)
;       Developed as part of the WSA-Enlil transition project
;       Released to the Solar Physics community via Solarsoft
;       February 2013 
;       
;----------------------------------------------------------------
;
; This software was developed at the NOAA Space Weather Prediction 
; Center by employees of the Federal Government in the course 
; of their official duties. Pursuant to title 17 Section 105 of the 
; United States Code this software is not subject to copyright 
; protection and is in the public domain. It is an experimental system. 
; NOAA assumes no responsibility whatsoever for its use by other parties,
; and makes no guarantees, expressed or implied, about its quality, 
; reliability, or any other characteristic.
; 
; We would appreciate acknowledgement if the software is used. 
; This software can be redistributed and/or modified freely provided 
; that any derivative works bear some notice that they are derived from it, 
; and any modified versions bear some notice that they have been modified.
; 
; We would like to acknowledge use of several routines taken from the wider
; Physics and idl programming communities - most notably software developed
; by Craig Markwardt and David Fanning:
; 
; http://cow.physics.wisc.edu/~craigm/idl/idl.html
; http://www.idlcoyote.com/
; 
; The development of CAT and it's use within the context of the WSA-Enlil
; project has been published in Space Weather Journal:
; http://onlinelibrary.wiley.com/doi/10.1002/swe.20024/abstract
;
; George Millward
; Curt de Koning
; Vic Pizzo
; Doug Biesecker
;
; NOAA Space Weather Prediction Center
; 325 Broadway, Boulder, CO 80305
;
; george.millward@noaa.gov
; 
;----------------------------------------------------------------
;
;
;+
; NAME:
;   QTCOMPOSE
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Convert a rotation angle and axis into quaternion
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   Q = QTCOMPOSE(VAXIS, PHI)
;
; DESCRIPTION:
;
;  The function QTCOMPOSE accepts a unit vector rotation axis VAXIS
;  and a rotation angle PHI, and returns the corresponding quaternion.
;
;  The user must take care to pass the same number of axes as rotation
;  angles.
;
;  Use QTAXIS and QTANG to extract the properties of an existing
;  quaternion.  Use QTCOMPOSE to combine a rotation axis and angle
;  into a new quaternion.
;
;  Conventions for storing quaternions vary in the literature and from
;  library to library.  This library uses the convention that the
;  first three components of each quaternion are the 3-vector axis of
;  rotation, and the 4th component is the rotation angle.  Expressed
;  in formulae, a single quaternion is given by:
;
;     Q(0:2) = [VX, VY, VZ]*SIN(PHI/2)
;     Q(3)   =              COS(PHI/2)
;
;  where PHI is the rotation angle, and VAXIS = [VX, VY, VZ] is the
;  rotation eigen axis expressed as a unit vector.  This library
;  accepts quaternions of both signs, but by preference returns
;  quaternions with a positive 4th component.
;
; INPUTS:
;
;  VAXIS - array of one or more unit vectors specifying the rotation
;          axes.  For a single rotation, VAXIS should be a 3-vector.
;          For N vectors, VAXIS should be a 3xN array.
;
;  PHI - one or more rotation angles, in radians.  For a single
;        rotation, PHI should be a scalar.  For N rotations, PHI
;        should be an N-vector.
;
; RETURNS:
;
;  For a single rotation, returns a quaternion as a 4-vector.  For N
;  rotations, returns a 4xN vector of quaternions.
;
;
; KEYWORD PARAMETERS:
;
;  NONE
;
; EXAMPLE:
;
;   IDL> print, qtcompose([0d,1,0], !dpi/4)
;          0.0000000      0.38268343       0.0000000      0.92387953
;
;   Prints the quaternion composed of a rotation of !dpi/4 radians
;   around the axis [0,1,0]
;
;
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, July 2001, CM
;   Documented, Dec 2001, CM
;   Allow output to be DOUBLE, 27 Jan 2002, CM
;   Allow vector vs scalar arguments, 28 Jan 2002, CM
;   Usage message, error checking, 15 Mar 2002, CM
;
;  $Id: qtcompose.pro,v 1.11 2008/12/14 20:00:31 craigm Exp $
;
;-
; Copyright (C) 2001, 2002, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
function swpc_cat_qtcompose, axis, phi

  if n_params() EQ 0 then begin
      info = 1
      USAGE:
      message, 'USAGE:', /info
      message, 'Q = QTCOMPOSE(AXIS, PHI)', info=1
      return, 0
  endif

  nph = n_elements(phi)
  nv  = n_elements(axis)/3
  if nph LT 1 OR nv LT 1 then goto, USAGE

  nq  = nv > nph
  q = make_array(value=axis(0)*phi(0)*0., 4,nq)

  sph = sin(phi/2) & cph = cos(phi/2)
  if nph EQ 1 AND nv EQ 1 then return, [ axis(0:2) * sph(0), cph(0) ]
  if nph GT 1 AND nv EQ 1 then begin
      ;; Single axis, multiple rotation angles
      q(0,*) = axis(0)*sph
      q(1,*) = axis(1)*sph
      q(2,*) = axis(2)*sph
      q(3,*) = cph
  endif else if nph EQ 1 AND nv GT 1 then begin
      ;; Multiple axis, single rotation
      q(0:2,*) = axis*sph(0)
      q(3,*)   = cph(0)
  endif else if nph EQ nv then begin
      ;; Multiple axes, multiple rotations
      q(0:2,*) = axis*rebin(reform(temporary(sph),1,nq),3,nq)
      q(3,*)   = temporary(cph)
  endif else begin
      message, 'ERROR: number of axes and angles do not match'
  endelse

  return, q
end


function swpc_cat_qtmult, aqt, bqt, inv1=inverse1, inv2=inverse2

; THIS ROUTINE MULTIPLIES QUATERNIONS
; CQT CORRESPONDS TO THE ROTATION AQT FOLLOWED BY BQT
; ASSUMING S/C COORDINATES ARE INITIALLY ALIGN WITH INERTIAL COORD.
; THEN ROTATION AQT DESCRIBES ROTATION SUCH THAT THE SUBROUTINE
;   QTXRA GIVES THE INERTIAL COORDINATES OF THE S/C X-AXIS
;   THE FIRST 3 COMPONENTS OF AQT GIVE THE EIGENAXIS EXPRESSED
;   IN S/C COORDINATES BEFORE THE ROTATION (=INTERTIAL COORD.).
; THE BQT ROTATION FOLLOWS THE AQT ROTATION. CQT THEN DESCRIBES
;   THIS COMBINATION SUCH THAT QTXRA GIVES THE INERTIAL COORDINATES
;   OF THE S/C X-AXIS AFTER BOTH ROTATIONS. 
;   THE FIRST 3 COMPONENTS OF BQT GIVE THE EIGENAXIS EXPRESSED
;   IN S/C COORDINATES AFTER THE AQT ROTATION.

  if n_params() EQ 0 then begin
      info = 1
      USAGE:
      message, 'USAGE:', /info
      message, 'QNEW = QTMULT(Q1, Q2)', info=info
      return, 0
  endif

  sz1 = size(aqt)
  sz2 = size(bqt)
  if sz1(0) LT 1 OR sz2(0) LT 1 then $
    message, 'ERROR: Q1 and Q2 must be quaternions'
  if sz1(1) NE 4 OR sz2(1) NE 4 then $
    message, 'ERROR: Q1 and Q2 must be quaternions'
  n1 = n_elements(aqt)/4
  n2 = n_elements(bqt)/4
  if n1 NE n2 AND n1 NE 1 AND n2 NE 1 then $
    message, 'ERROR: Q1 and Q2 must both have the same number of quaternions'

  nq = n1>n2
  cqt = make_array(value=aqt(0)*bqt(0)*0, dimension=[4,nq])

  if n1 GT 1 then begin
      aqt0 = aqt(0,*) & aqt1 = aqt(1,*) & aqt2 = aqt(2,*) & aqt3 = aqt(3,*)
  endif else begin
      aqt0 = aqt(0) & aqt1 = aqt(1) & aqt2 = aqt(2) & aqt3 = aqt(3)
  endelse
  if n2 GT 1 then begin
      bqt0 = bqt(0,*) & bqt1 = bqt(1,*) & bqt2 = bqt(2,*) & bqt3 = bqt(3,*)
  endif else begin
      bqt0 = bqt(0) & bqt1 = bqt(1) & bqt2 = bqt(2) & bqt3 = bqt(3)
  endelse
  if keyword_set(inverse1) then begin
      aqt0 = -aqt0 & aqt1 = -aqt1 & aqt2 = -aqt2
  endif
  if keyword_set(inverse2) then begin
      bqt0 = -bqt0 & bqt1 = -bqt1 & bqt2 = -bqt2
  endif

  CQT(0,0) = AQT0*BQT3 + AQT1*BQT2 - AQT2*BQT1 + AQT3*BQT0
  CQT(1,0) =-AQT0*BQT2 + AQT1*BQT3 + AQT2*BQT0 + AQT3*BQT1
  CQT(2,0) = AQT0*BQT1 - AQT1*BQT0 + AQT2*BQT3 + AQT3*BQT2
  CQT(3,0) =-AQT0*BQT0 - AQT1*BQT1 - AQT2*BQT2 + AQT3*BQT3
  
  return, cqt
end



;+
; NAME:
;   QTEULER
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Compose a series of euler-type rotations into a single quaternion
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   Q = QTEULER(AXES, ANG0, ANG1, ... )
;
; DESCRIPTION:
;
;  The function QTEULER composes a series of Euler-type rotations into
;  a single set of quaternion representations.
;
;  The user specifies a set of axes, and the angles to rotation about
;  those axes, and QTEULER constructs the corresponding quaternion.
;
;  There must be a one-to-one correspondence between the elements of
;  AXES and the number of rotations.  AXES specifies the rotation axes
;  as an string, which must be one of 'X', 'Y', or 'Z'.  Other axes
;  are invalid.  For example, the following call:
;
;    QTEULER(['X','Z'], THETA, PHI)
;
;  will rotate first about the *Z* axis by the angle PHI, and then
;  around the *resulting X* axis by angle THETA.
;
;  Several things are worth noting here.  First, rotations are applied
;  first from the right, not the left.  This conforms to the usual
;  matrix notation for applying rotations to a vector on the right
;  hand side.  For example, in matrix notation,
;
;      XNEW = A3 A2 A1 XOLD
;
;  applies first A1, then A2 and finally A3 to the XOLD vector,
;  resulting in the new vector XNEW.  The same semantics apply here.
;
;  A second thing to bear in mind is that the axes themselves change
;  during the rotations.  Thus, the coordinates specified in AXES
;  should be considered attached to the "body" and not the inertial
;  frame.
;
;
; INPUTS:
;
;  AXES - a string array, specifies the rotation axes.  Rotations are
;         applied last element first.  Each element of AXES must be
;         one of 'X', 'Y' or 'Z'.
;
;  ANG0, ..., ANGi - the successive rotation angles.  Angle ANGi
;         corresponds to axis AXES(i).
;
;         If ANGi is a scalar, then it will be promoted to a vector
;         the same size as the other rotation angles being performed.
;         Otherwise, if the angles ANGi are vectors, then they must
;         all be of the same size.
;
; RETURNS:
;
;  The resulting quaternion (or, if ANGi are vectors, array of
;  quaternions), which represent the requested rotations.
;
;
; KEYWORD PARAMETERS:
;
;  NONE
;
; EXAMPLE:
;
;    ;;              Precession        Nutation    
;    qtot = qteuler(['z','y','z',      'x','z','x'        ], $
;                    -zeta, +theta, -z, +eps0, -dpsi, -eps)
;
;   Applies a series of rotations to correct for earth nutation and
;   precession.  The order of rotations on a vector would be
;   X-Z-X-Z-Y-Z (i.e., the reverse order printed).
;
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, 27 Jan 2002, CM
;   More error checking, 03 Mar 2002, CM
;
;  $Id: qteuler.pro,v 1.4 2002/05/09 23:03:27 craigm Exp $
;
;-
; Copyright (C) 2002, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

;; Extract axis ei and angle angi
pro swpc_cat_qteuler_extract, ax, i, ei, angi, $
                     ang0, ang1, ang2, ang3, ang4, $
                     ang5, ang6, ang7, ang8, ang9, $
                     status=status, errmsg=errmsg
                     
  status = 0
  zero = ang0(0)*0
  ex = [1,zero,zero] & ey = [zero,1,zero] & ez = [zero,zero,1]

  ei = [0D, 0D, 0D]


; the following line will not work under the idl virtual machine.
; commenting it out for now because I don't think it is ever used.
; ghgm 7/22/2013

if i eq 0 then begin
ei = ez
angi = ang0 
endif
if i eq 1 then begin
ei = ey
angi = ang1
endif
if i eq 2 then begin
ei = ex
angi = ang2
endif

;print, ' here 1'
;
;print,'i = *',i,'*'
;print, ax(i)
;print, strtrim(i,2)
;
;void = execute('ei = e'+ax(i)+' & angi = ang'+strtrim(i,2))
;
;print, void
;
;print, ' here 2'


;  if execute('ei = e'+ax(i)+' & angi = ang'+strtrim(i,2)) NE 1 then begin
;      stop
;      errmsg = 'Invalid axis specification'
;      return
;  endif

  status = 1
  return
end

function swpc_cat_qteuler, axes, block=block, $
            ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9, $
            ang10, ang11, ang12, ang13, ang14, ang15

  if n_params() EQ 0 then begin
      info = 1
      USAGE_ERR:
      message, 'USAGE: Q = QTEULER(AXES, ANG0, ...)', /info
      message, '  AXES = ["X",...]    ("X" or "Y" or "Z")', /info
      message, '  ANGn = rotation angle (radians)', info=info
      return, 0
  endif
  if n_elements(axes) LT 1 OR n_elements(ang0) LT 1 then $
    goto, USAGE_ERR
  nang = n_params()-1

  ;; Check to be sure each axis label is 'X' 'Y' or 'Z'
  ax = strupcase(strmid(strtrim(axes,2),0,1))
  wh = where(ax NE 'X' AND ax NE 'Y' AND ax NE 'Z', ct)
  if ct GT 0 then begin
      errmsg = 'AXES must be one of "X", "Y" or "Z"'
      goto, BAD_AXIS
  endif
  if n_elements(ax) NE nang then begin
      errmsg = 'Number of AXES and rotations ANGi must agree'
      goto, BAD_AXIS
  endif

  swpc_cat_qteuler_extract, ax, 0, ev, angv, status=status, errmsg=errmsg, $
    ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9

  if status EQ 0 then begin
      BAD_AXIS:
      message, 'ERROR: '+errmsg, /info
      goto, USAGE_ERR
  endif

  qq = swpc_cat_qtcompose(ev, angv)
  for i = 1, nang-1 do begin
      swpc_cat_qteuler_extract, ax, i, ev, angv, status=status, errmsg=errmsg, $
        ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9
      if status EQ 0 then goto, BAD_AXIS
      
      qq = swpc_cat_qtmult(qq, swpc_cat_qtcompose(ev, angv))
  endfor

  return, qq
end


;+
; NAME:
;   QTMULT
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Multiply quaternions
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   QRESULT = QTMULT(Q1, [/INV1,]   Q2, [/INV2])
;
; DESCRIPTION:
;
;   The function QTMULT performs multiplication of quaternions.
;   Quaternion multiplication is not component-by-component, but
;   rather represents the composition of two rotations, namely Q2
;   followed by Q1.
;
;   More than one multiplication can be performed at one time if Q1
;   and Q2 are 4xN arrays.  In that case both input arrays must be of
;   the same dimension.
;
;   If INV1 is set, then the inverse of Q1 is used.  This is a
;   convenience, to avoid the call QTINV(Q1).  Of course, INV2 can
;   be set to use the inverse of Q2.
;
;   Note that quaternion multiplication is not commutative.
;
;  Conventions for storing quaternions vary in the literature and from
;  library to library.  This library uses the convention that the
;  first three components of each quaternion are the 3-vector axis of
;  rotation, and the 4th component is the rotation angle.  Expressed
;  in formulae, a single quaternion is given by:
;
;     Q(0:2) = [VX, VY, VZ]*SIN(PHI/2)
;     Q(3)   =              COS(PHI/2)
;
;  where PHI is the rotation angle, and VAXIS = [VX, VY, VZ] is the
;  rotation eigen axis expressed as a unit vector.  This library
;  accepts quaternions of both signs, but by preference returns
;  quaternions with a positive 4th component.
;
;
; INPUTS:
;
;  Q1 - array of one or more unit quaternions, the first operand in
;       the multiplication.  For a single quaternion, Q1 should be a
;       4-vector.  For N quaternions, Q1 should be a 4xN array.
;       If INV1 is set, then the inverse of Q1 is used.
;
;  Q2 - same as Q1, for the second operand.
;       If INV2 is set, then the inverse of Q2 is used.
;
; RETURNS:
;
;   The resulting multiplied unit quaternions.  For a single inputs,
;   returns a 4-vector.  For N input quaternions, returns N
;   quaternions as a 4xN array.
;
;
; KEYWORD PARAMETERS:
;
;  INV1 - if set, use QTINV(Q1) in place of Q1.
;
;  INV2 - if set, use QTINV(Q2) in place of Q2.
;
; EXAMPLE:
;
;   Q1 = qtcompose([0,0,1],  32d*!dpi/180d)
;   Q2 = qtcompose([1,0,0], 116d*!dpi/180d)
;
;   IDL> print, qtmult(q1, q2)
;        0.81519615      0.23375373      0.14606554      0.50939109
;
;   Form a rotation quaternion of 32 degrees around the Z axis, and 
;   116 degrees around the X axis, then multiply the two quaternions.
;   
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTMULTN, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, July 2001, CM
;   Documented, Dec 2001, CM
;   Documentation, allow 1xN or Nx1 multiplies, 27 Jan 2002, CM
;   Usage message, error checking, 15 Mar 2002, CM
;   Add the INV1 and INV2 keywords, 30 Aug 2007, CM
;
;  $Id: qtmult.pro,v 1.8 2007/09/03 07:18:25 craigm Exp $
;
;-
; Copyright (C) 2001, 2002, 2007, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-




;+
; NAME:
;   QTVROT
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Apply quaternion rotation to a 3-vector
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   VNEW = QTVROT(V, Q, [/INVERT])
;
; DESCRIPTION:
;
;   The function QTVROT applies a quaternion rotation (or its inverse)
;   to a 3-vector V to produce a new vector VNEW.
;
;   If both V and VNEW are vector components measured in the same
;   inertial coordinate system, then VNEW returns the components of
;   the vector V rotated by quaternion Q.  I.e., the AXES stay fixed
;   and the VECTOR rotates.  Replace Q by QTINV(Q) in the case of
;   /INVERT.
;
;   If V are components of a vector measured in the "body" coordinate
;   frame, and Q represents the orientation of the body frame
;   w.r.t. the inertial frame, then VNEW are the components of the
;   same vector in the inertial frame.  I.e., the VECTOR stays fixed
;   and the AXES rotate.  For /INVERT, the coordinate transformation
;   is from inertial frame to body frame.
;
;   If either Q is a single quaternion, or V is a single 3-vector,
;   then QTVROT will expand the single to the number of elements of
;   the other operand.  Otherwise, the number of quaternions and
;   vectors must be equal.
;
;  Conventions for storing quaternions vary in the literature and from
;  library to library.  This library uses the convention that the
;  first three components of each quaternion are the 3-vector axis of
;  rotation, and the 4th component is the rotation angle.  Expressed
;  in formulae, a single quaternion is given by:
;
;     Q(0:2) = [VX, VY, VZ]*SIN(PHI/2)
;     Q(3)   =              COS(PHI/2)
;
;  where PHI is the rotation angle, and VAXIS = [VX, VY, VZ] is the
;  rotation eigen axis expressed as a unit vector.  This library
;  accepts quaternions of both signs, but by preference returns
;  quaternions with a positive 4th component.
;
;
; INPUTS:
;
;  V - array of one or more 3-vectors.  For a single vector, V should
;      be a 3-vector.  For N vectors, V should be a 3xN array.
;
;  Q - array of one or more unit quaternions.  For a single
;      quaternion, Q should be a 4-vector.  For N quaternions, Q
;      should be a 4xN array.
;
;
; RETURNS:
;
;   The resulting rotated vectors.  For single inputs, returns a
;   3-vector.  For N inputs, returns N vectors as a 3xN array.
;
;
; KEYWORD PARAMETERS:
;
;   INVERT - if set, then the antirotation represented by QTINV(Q) is
;            performed.
;
;
; EXAMPLE:
;
;   Q1 = qtcompose([0,0,1],  32d*!dpi/180d)
;   Q2 = qtcompose([1,0,0], 116d*!dpi/180d)
;   Q = qtmult(Q1, Q2)
;
;   V = [[1d,0,0],[0,1,0],[0,0,1]]
;
;   IDL> print, qtvrot(v, q)
;         0.84804810      0.52991926       0.0000000
;         0.23230132     -0.37175982      0.89879405
;         0.47628828     -0.76222058     -0.43837115
;
;
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, July 2001, CM
;   Documented, Dec 2001, CM
;   Small changes, 28 Jan 2002, CM
;   Usage message, error checking, 15 Mar 2002, CM
;
;  $Id: qtvrot.pro,v 1.7 2002/05/09 23:03:27 craigm Exp $
;
;-
; Copyright (C) 2001, 2002, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

;; QVROT
;;
;; The FORWARD (default) transform: 
;;
;; * takes a vector vin (components given in inertial coordinates) and
;;  returns the components of the rotated vector vout (components
;;  given in inertial coordinates) -- ie, the AXES stay fixed and the
;;  VECTOR rotates; OR, equivalently,
;;
;; * takes a fixed vector vin (components given in body coordinates)
;;   and returns the components of the vector in inertial coordinates,
;;   where the body system is described by quaternion q -- ie, the
;;   VECTOR stays fixed and the AXES rotate.
;;
;;
;; The INVERSE transform (gotten by setting /INVERT):
;;
;; * takes a vector vin (components given in inertial coordinates) and
;;   returns the components of the anti-rotated vector vout
;;   (components given in inertial coordinates) -- ie, the AXES stay
;;   fixed and the VECTOR rotates.  Anti-rotated here means rotated in
;;   the opposite direction of q; OR, equivalently,
;;
;; * takes a fixed vector vin (components given in inertial
;;   coordinates) and returns the components of the vector in body
;;   coordinates, where the body system is described by quaternion q
;;   -- ie, the VECTOR stays fixed and the AXES rotate.
;;

function swpc_cat_qtvrot, vin, q, invert=invert

  if n_params() EQ 0 then begin
      info = 1
      USAGE:
      message, 'USAGE:', /info
      message, 'VNEW = QTVROT(V, Q)', info=info
      return, 0
  endif
  nq = n_elements(q)/4
  nv = n_elements(vin)/3
  if nq LT 1 OR nv LT 1 then goto, USAGE

  if n_elements(q) GT 4 AND n_elements(vin) GT 3 then begin
      if n_elements(q)/4 NE n_elements(vin)/3 then begin
          message, 'ERROR: incompatible number of quaternions & vectors'
          return, -1L
      end
      vout = vin*q(0)*0.
      nq = n_elements(q)/4
      nv = nq
  endif else if n_elements(q) GT 4 then begin
      nq = n_elements(q)/4
      nv = 1L
      vout = vin(*) # (fltarr(nq)+1) * q(0)*0.
  endif else begin
      nq = 1L
      nv = n_elements(vin)/3
      vout = vin*q(0)*0.
  endelse
  vout = reform(vout, 3, max([nv,nq]), /overwrite)
  q1 = q(0,*) & q2 = q(1,*) & q3 = q(2,*) & q4 = q(3,*)
  if n_elements(q1) EQ 1 then begin
      q1 = q1(0)  & q2 = q2(0)  & q3 = q3(0)  & q4 = q4(0)
  endif else begin
      q1 = q1(*)  & q2 = q2(*)  & q3 = q3(*)  & q4 = q4(*)
  endelse
  v0 = vin(0,*) & v1 = vin(1,*) & v2 = vin(2,*)
  if n_elements(v0) EQ 1 then begin
      v0 = v0(0)  & v1 = v1(0)  & v2 = v2(0)
  endif else begin
      v0 = v0(*)  & v1 = v1(*)  & v2 = v2(*)
  endelse

  if NOT keyword_set(INVERT) then begin

      ;; FORWARD TRANSFORMATION
      VOUT(0,*)=((Q1*Q1-Q2*Q2-Q3*Q3+Q4*Q4)*V0 $
                 + 2.D0*(Q1*Q2-Q3*Q4)*V1 $
                 + 2.D0*(Q1*Q3+Q2*Q4)*V2)
      VOUT(1,*)=(2.D0*(Q1*Q2+Q3*Q4)*V0 $
                 + (-Q1*Q1+Q2*Q2-Q3*Q3+Q4*Q4)*V1 $
                 + 2.D0*(Q2*Q3-Q1*Q4)*V2)
      VOUT(2,*)=(2.D0*(Q1*Q3-Q2*Q4)*V0 $
                 + 2.D0*(Q2*Q3+Q1*Q4)*V1 $
                 + (-Q1*Q1-Q2*Q2+Q3*Q3+Q4*Q4)*V2)
  endif else begin
      ;; INVERSE TRANSFORMATION
      VOUT(0,*)=((Q1*Q1-Q2*Q2-Q3*Q3+Q4*Q4)*V0 $
                 + 2.D0*(Q1*Q2+Q3*Q4)*V1 $
                 + 2.D0*(Q1*Q3-Q2*Q4)*V2)
      VOUT(1,*)=(2.D0*(Q1*Q2-Q3*Q4)*V0 $
               + (-Q1*Q1+Q2*Q2-Q3*Q3+Q4*Q4)*V1 $
               + 2.D0*(Q2*Q3+Q1*Q4)*V2)
      VOUT(2,*)=(2.D0*(Q1*Q3+Q2*Q4)*V0 $
               + 2.D0*(Q2*Q3-Q1*Q4)*V1 $
               + (-Q1*Q1-Q2*Q2+Q3*Q3+Q4*Q4)*V2)
  endelse
  vout = vout

  return, vout
end





pro swpc_cat_GENLEM_CALC_VERTICES, params, vert_data,calc_info, lemniscate_style $
                        , xverts=xverts, yverts=yverts

; BUILD THE LEMNISCATE IN THE STANDARD FRAME
; in the standard frame the central axis of the lemniscate is just the
; x-axis

if not(keyword_set(xverts))    then xverts = 30
if not(keyword_set(yverts))    then yverts = 90

nrows = long(xverts) > 2
ncols = long(yverts) > 2

dblp = double(params)

; build lemniscate used parametric equation
; tttt is a parameter ranging from [0,pi/2] (excluding endpoints)
; ssss is a parameter of revolution, ranging from [0,2pi]
tttt = !dpi/2.0d0 * (1.0d0+dindgen(nrows))/double(nrows+1L)
ssss = 2.0d0*!dpi * (dindgen(ncols))/double(ncols-1L)

xmax = dblp[0]
tmpx = xmax * cos(tttt)
if lemniscate_style then tmpx /= (1.0d0+sin(tttt)^2)

ymax = dblp[1] * cos(tttt)*sin(tttt)
tmpy = cos(ssss)
if lemniscate_style then ymax /= (1.0d0+sin(tttt)^2)

tmpz = -sin(ssss)

nverts = nrows*ncols + 2L
tmpdat = dblarr(nverts,3)

for rrr=0L,nrows-1L do begin
   x1d = replicate(tmpx[rrr],ncols)
   y1d = ymax[rrr]*tmpy
   z1d = dblp[2]*ymax[rrr]*tmpz

; roundoff error can occur when y=ymax --- this should result in z=0,
;                                          but sometimes it results in
;                                          z^2=-1.0e-17, for example;
;                                          in that case sqrt(z^2) will
;                                          return a NaN
   z1d[0] = 0.0d0
   z1d[ncols-1] = 0.0d0

   tmpdat[rrr*ncols:(rrr+1L)*ncols-1L,*] = [[x1d],[y1d],[z1d]]
endfor
tmpdat[nverts-2L,*] = [0.0d0,0.0d0,0.0d0]
tmpdat[nverts-1L,*] = [xmax,0.0d0,0.0d0]

vert_data = transpose(tmpdat)

calc_info = { info,            $
              n_xverts: nrows, $
              n_yverts: ncols, $
              numverts: nverts $
            }

end



pro swpc_cat_GENLEM_ROTATE, angs,vert_orig, vert_rot

; ROTATE LEMNISCATE FROM STANDARD FRAME TO DATA FRAME

npts = n_elements(vert_orig) / 3L

; set up quaternion to rotate points from model-space to data-space
;   1st rotation: adjust for tilt
;   2nd rotation: rotate from inertial xy-plane to central plane of
;                 point cloud 
;   3rd rotation: rotate from inertial x-axis to central axis of the
;                 point cloud 
; see, also, wikipedia article on euler angles: rotation XYZ of points
; note that the complete 3-step rotation from Arfken [p 200] is not
; consistent with the physical concept of tilt
qrot = swpc_cat_qteuler(['Z','Y','X'],angs[0],-angs[1],-angs[2])

; perform rotation
vert_rot = swpc_cat_qtvrot(vert_orig,rebin(qrot,4L,npts))

END



pro swpc_cat_GENLEM_CALC_CONNECTIONS, calc_info, connex

; FILL IN THE CONNECTIVITY ARRAY

nrows = calc_info.n_xverts
ncols = calc_info.n_yverts
nconn = (ncols*(nrows-1L)*5L) + (2L*ncols*4L)
connex = lonarr(nconn)

i = 0L
for k=0L,nrows-2L do begin
   for j=0L,ncols-1L do begin
      connex[i] = 4L

      connex[i+1L] = (k+1L)*ncols + j

      w = ((j + 1L) mod ncols) + (k+1L)*ncols
      connex[i+2L] = w

      w = ((j + 1L) mod ncols) + k*ncols
      connex[i+3L] = w

      connex[i+4L] = k*ncols + j

      i = i + 5L
   endfor
endfor

; connect base point to first strip
for j=0L,ncols-1L do begin
   connex[i] = 3L
   connex[i+1L] = calc_info.numverts-1L
   connex[i+2L] = j
   connex[i+3L] = (j+1L) mod ncols
   i = i + 4L
endfor

; connect xmax to last strip
for j=0L,ncols-1L do begin
   connex[i] = 3L
   connex[i+1L] = calc_info.numverts-2L
   connex[i+2L] = j+(nrows-1L)*ncols
   connex[i+3L] = ((j+1L) mod ncols) + (nrows-1L)*ncols
   i = i + 4L
endfor

END



;-------------------------------------------------------------------------
; GENERAL_LEMNISCATE_DEFINE
;
; Purpose:
;   Defines and populates a data structure for a generalized lemniscate.
;   In the so-called standard frame, a generalized lemniscate is based
;   on a quartic surface.
;   For the lemniscate of Gerono, the surface is parameterized as
;      x = c1*cos(t)
;      y = c2*cos(t) sin(t) cos(s)
;      z = c3*sin(s)
;   For the lemniscate of Bernoulli, the surface is parameterized as
;      x = c1*cos(t)/(1+sin(t)^2)
;      y = c2*cos(t) sin(t) cos(s) / (1+sin(t)^2)
;      z = c3*sin(s)
;   In both cases, t=[0,pi/2], s=[0,2pi]
;   To fully generalized the surface, it will be rotated such that its
;   central axis (the x-axis in the standard frame) has an arbitrary
;   orientaion; in addition, the generalized lemniscate will have an
;   arbitrary tilt
;-------------------------------------------------------------------------

function swpc_cat_GENERAL_LEMNISCATE_DEFINE, params, in_angles, lemniscate_style $
                                  , deg_or_rad=deg_or_rad               $
                                  , xverts=xverts, yverts=yverts

if not(keyword_set(xverts))    then xverts = 0
if not(keyword_set(yverts))    then yverts = 0

if not(keyword_set(deg_or_rad)) then deg_or_rad = 'rad'
if strlowcase(deg_or_rad) eq 'deg' then begin
   pro_angs = !dpi/180.0d0 * in_angles
endif else pro_angs = in_angles

swpc_cat_genlem_calc_vertices,params,verts_stdfr,calc_info,lemniscate_style, $
                     xverts=xverts,yverts=yverts
swpc_cat_genlem_rotate,pro_angs,verts_stdfr,verts_datfr
swpc_cat_genlem_calc_connections,calc_info,connex

genlem = { lemniscate,                  $
           parameters: params,          $
           orientation: in_angles[0:1], $
           tilt: in_angles[2],          $
           vertices: verts_datfr,       $
           connections: connex          $
         }

return, genlem

END


;
;+
; NAME:
;       FSC_NORMALIZE
;
; PURPOSE:
;
;       This is a utility routine to calculate the scaling vector
;       required to position a graphics primitive of specified range
;       at a specific position in an arbitray coordinate system. The
;       scaling vector is given as a two-element array like this:
;
;          scalingVector = [translationFactor, scalingFactor]
;
;       The scaling vector should be used with the [XYZ]COORD_CONV
;       keywords of a graphics object or model. For example, if you
;       wanted to scale an X axis into the coordinate range of -0.5 to 0.5,
;       you might type something like this:
;
;          xAxis->GetProperty, Range=xRange
;          xScale = FSC_Normalize(xRange, Position=[-0.5, 0.5])
;          xAxis, XCoord_Conv=xScale
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

;       Object Graphics
;
; CALLING SEQUENCE:
; 
;       xscaling = FSC_NORMALIZE(xrange, POSITION=position)
;
; INPUTS:
; 
;       XRANGE: A two-element vector specifying the data range.
;
; KEYWORD PARAMETERS:
; 
;       POSITION: A two-element vector specifying the location
;       in the coordinate system you are scaling into. The vector [0,1]
;       is used by default if POSITION is not specified.
;
; COMMON BLOCKS:
; 
;       None.
;
; EXAMPLE:
; 
;       See above.
;
; MODIFICATION HISTORY:
;       Written by:  David W. Fanning, OCT 1997.
;       Fixed a problem with illegal divide by zero. 21 April 2005. DWF.
;       Fixed a problem when range[0] is greater than range[1]. 11 July 2006. DWF.
;       Renamed to FSC_Normalize to avoid conflicts with 10,000 other programs named NORMALIZE. 17 October 2008. DWF.
;-
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

function swpc_cat_FSC_Normalize, range, Position=position

On_Error, 1
IF N_Params() EQ 0 THEN Message, 'Please pass range vector as argument.'

IF (N_Elements(position) EQ 0) THEN position = [0.0D, 1.0D] ELSE $
    position=Double(position)
range = Double(range)

IF range[1] GE range[0] THEN BEGIN
   scale = [((position[0]*range[1])-(position[1]*range[0])) / $
       ((range[1]-range[0]) > 1e-12), (position[1]-position[0])/((range[1]-range[0]) > 1e-12)]
ENDIF ELSE BEGIN
   scale = [((position[1]*range[0])-(position[0]*range[1])) / $
       ((range[0]-range[1]) > 1e-12), (position[1]-position[0])/((range[0]-range[1]) > 1e-12)]
   scale[1] = -scale[1]
ENDELSE
RETURN, scale
END
;-------------------------------------------------------------------------





; ECDP
;pro swpc_cat_get_config,source_path,sep,telescope_array,image_in_folder_ary, $
;              image_in_root,export_location_root, $
;              max_interval_in_days,OPS_or_VnV,input_file

;compile_opt idl2

; CASE STRUPCASE(!Version.OS_FAMILY)OF
;   'WINDOWS' : sep = '\'
;   'UNIX'    : sep = '/'   ;if running on linux, IDL returns this...
;   ELSE      : sep = '/'
; ENDCASE

;telescopes = ''           ; Coronagraph telescopes having images to be processed.
;image_in_folders = ''     ; Names of the folders containing the input FITS images.
;image_in_root = ''        ; Root of the location of the input FITS image folders on the server.
;export_location_root = '' ; Root location for export of results.
;max_interval_in_days = '' ; Upper limit of the datetime interval when selecting input images.
;OPS_or_VnV = ''
;line = ''

;get configuration info from input file
;input_file = source_path + sep + 'swpc_cat.in'
;openr,lun,input_file,/GET_LUN
;WHILE NOT EOF(lun) DO BEGIN
;  READF, lun, line
;  pos = STRPOS(line,'=')
;  if STRCMP( line,'telescopes', pos, /FOLD_CASE ) then telescopes = STRMID(line,pos+1)
;  if STRCMP( line,'image_in_folders', pos, /FOLD_CASE ) then image_in_folders = STRMID(line,pos+1)
;  if STRCMP( line,'image_in_root', pos, /FOLD_CASE ) then image_in_root = STRMID(line,pos+1)
;  if STRCMP( line,'export_location_root', pos, /FOLD_CASE ) then export_location_root = STRMID(line,pos+1)
;  if STRCMP( line,'max_interval_in_days', pos, /FOLD_CASE ) then max_interval_in_days = STRMID(line,pos+1)
;  if STRCMP( line,'run_as_OPS_or_VnV', pos, /FOLD_CASE ) then OPS_or_VnV = STRMID(line,pos+1)
     
;ENDWHILE
;close,lun
;free_lun, lun

;convert some strings to arrays:
;telescope_array = strsplit(telescopes,',',/EXTRACT)
;image_in_folder_ary = strsplit(image_in_folders,',',/EXTRACT)

;convert max interval to integer
;max_interval_in_days = FIX(max_interval_in_days)

;END

pro swpc_cat_get_config,source_path,sep,telescope_array,image_in_folder_ary, $
              image_in_root,export_location_root, preferences_folder, $
              max_interval_in_days,OPS_or_VnV,earth_position_filename, $
              input_file

compile_opt idl2

 CASE STRUPCASE(!Version.OS_FAMILY)OF
   'WINDOWS' : sep = '\'
   'UNIX'    : sep = '/'   ;if running on linux, IDL returns this...
   ELSE      : sep = '/'
 ENDCASE

telescopes = ''           ; Coronagraph telescopes having images to be processed.
image_in_folders = ''     ; Names of the folders containing the input FITS images.
image_in_root = ''        ; Root of the location of the input FITS image folders on the server.
export_location_root = '' ; Root location for export of results.
preferences_folder = ''   ; local folder to hold preferences
max_interval_in_days = '' ; Upper limit of the datetime interval when selecting input images.
OPS_or_VnV = ''
earth_position_filename = '' ; Name of the Earth position ephemeris file (Earth_position_1972_to_2030.txt)
line = ''

;get configuration info from input file
input_file = source_path + sep + 'swpc_cat.in'
openr,lun,input_file,/GET_LUN
WHILE NOT EOF(lun) DO BEGIN
  READF, lun, line
  pos = STRPOS(line,'=')
  if STRCMP( line,'telescopes', pos, /FOLD_CASE ) then telescopes = STRMID(line,pos+1)
  if STRCMP( line,'image_in_folders', pos, /FOLD_CASE ) then image_in_folders = STRMID(line,pos+1)
  if STRCMP( line,'image_in_root', pos, /FOLD_CASE ) then image_in_root = STRMID(line,pos+1)
  if STRCMP( line,'export_location_root', pos, /FOLD_CASE ) then export_location_root = STRMID(line,pos+1)
;  if STRCMP( line,'preferences_folder', pos, /FOLD_CASE ) then preferences_folder = STRMID(line,pos+1)
  if STRCMP( line,'max_interval_in_days', pos, /FOLD_CASE ) then max_interval_in_days = STRMID(line,pos+1)
  if STRCMP( line,'run_as_OPS_or_VnV', pos, /FOLD_CASE ) then OPS_or_VnV = STRMID(line,pos+1)
  if STRCMP( line,'earth_position_filename', pos, /FOLD_CASE ) then earth_position_filename = STRMID(line,pos+1)
  
     
ENDWHILE
close,lun

;convert some strings to arrays:
telescope_array = strsplit(telescopes,',',/EXTRACT)
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

function swpc_cat_SourcePath, Base_Name = BaseName, _Extra = Extra
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



function swpc_cat_which_style_lemniscate, event
WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

if info.lemniscate_style eq 0 then begin
info.lemniscate_style = 1
endif else begin
info.lemniscate_style = 0
endelse

if info.debug_mode eq 1 then print, info.lemniscate_style

swpc_cat_define_cme_lemniscate, info.radial_distance_lemniscate, info.angular_width_lemniscate, info.lemniscate_style, fitted_cme_info

info.L_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.L_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.L_Window_copy->Draw, info.L_second_view_copy
swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline
info.L_Window->Draw, info.L_both_views

info.C_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.C_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.C_Window_copy->Draw, info.C_second_view_copy
swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
info.C_Window->Draw, info.C_both_views

info.R_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.R_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.R_Window_copy->Draw, info.R_second_view_copy
swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline
info.R_Window->Draw, info.R_both_views


WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
return,1
END




pro swpc_cat_show_line_plot, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

CASE event.id OF

   info.L_widget_show_line_plot : BEGIN
   
   info.L_image_lineplot_model -> SetProperty, hide = 0
   info.L_clock_angle_model -> SetProperty, hide = 0
   info.show_image_line_plot = 1
      swpc_cat_replot_image_line_plot, $
      info.L_clock_angle_degrees, info.L_coronagraph_image_object, info.L_image_lineplot, $
      info.position_image_lineplot, info.L_cme_outline
   info.L_Window->Draw, info.L_both_views

   ENDCASE
      
   info.C_widget_show_line_plot : BEGIN
   
   info.C_image_lineplot_model -> SetProperty, hide = 0
   info.C_clock_angle_model -> SetProperty, hide = 0
   info.show_image_line_plot = 1
      swpc_cat_replot_image_line_plot, $
      info.C_clock_angle_degrees, info.C_coronagraph_image_object, info.C_image_lineplot, $
      info.position_image_lineplot, info.C_cme_outline
   info.C_Window->Draw, info.C_both_views

   ENDCASE
      
   info.R_widget_show_line_plot : BEGIN
   
   info.R_image_lineplot_model -> SetProperty, hide = 0
   info.R_clock_angle_model -> SetProperty, hide = 0
   info.show_image_line_plot = 1
      swpc_cat_replot_image_line_plot, $
      info.R_clock_angle_degrees, info.R_coronagraph_image_object, info.R_image_lineplot, $
      info.position_image_lineplot, info.R_cme_outline
   info.R_Window->Draw, info.R_both_views

   ENDCASE
     
ENDCASE


WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

end











pro swpc_cat_show_3d_view, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

;         info.Earth_pos_AU
;         info.Earth_pos_HG_LAT_deg

;info.L_HEEQ_coords = [0.,0.]
;info.C_HEEQ_coords = [0.,0.]
;info.R_HEEQ_coords = [0.,0.]

;print, info.L_HEEQ_coords
;print, info.C_HEEQ_coords
;print, info.R_HEEQ_coords

;stereo_b_phi = info.L_HEEQ_coords[0] * !dtor + !pi
;stereo_b_theta = !pi/2 - info.L_HEEQ_coords[1] * !dtor
;
;earth_phi = info.C_HEEQ_coords[0] * !dtor + !pi
;earth_theta = !pi/2 - info.C_HEEQ_coords[1] * !dtor
;
;stereo_a_phi = info.R_HEEQ_coords[0] * !dtor + !pi
;stereo_a_theta = !pi/2 - info.R_HEEQ_coords[1] * !dtor

stereo_b_phi = info.L_HEEQ_coords[0]
stereo_b_theta = info.L_HEEQ_coords[1]

earth_phi = info.C_HEEQ_coords[0]
earth_theta = info.C_HEEQ_coords[1]

stereo_a_phi = info.R_HEEQ_coords[0]
stereo_a_theta = info.R_HEEQ_coords[1]

;swpc_cat_3d_view, earth_phi=earth_phi, earth_theta=earth_theta, stereo_a_phi=stereo_a_phi, stereo_a_theta=stereo_a_theta, $
;             stereo_b_phi=stereo_b_phi, stereo_b_theta=stereo_b_theta, tlb_3D, Three_D_view_event_widget
             
;swpc_cat_3d_view2, tlb_3D, Three_D_view_event_widget, $
;                     earth_theta, earth_phi, $
;                     stereo_a_theta, stereo_a_phi, $
;                     stereo_b_theta, stereo_b_phi, $
;                     info.latitude_degrees, $
;                     info.longitude_degrees, $
;                     info.radial_distance_lemniscate, $
;                     info.angular_width_lemniscate                  
             
info.Three_D_view_tlb_ID = tlb_3D
info.Three_D_view_event_widget = Three_D_view_event_widget

;print, ' in CAT, tlb_3D = ', info.Three_D_view_tlb_ID
             

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

end





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


pro swpc_cat_datebox_add_manual_point_CenterTLB, tlb

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

pro swpc_cat_datebox_add_manual_point_Event, event

   ; This event handler responds to all events. Widget
   ; is always destoyed. The text is recorded if ACCEPT
   ; button is selected or user hits CR in text widget.

Widget_Control, event.top, Get_UValue=info

CASE event.ID OF
   info.cancelID: Widget_Control, event.top, /Destroy
   ELSE: BEGIN

         ; Get the text and store it in the pointer location.

;      Widget_Control, info.textID, Get_Value=theText
      year = info.yearID -> Get_Value()
      month = info.monthID -> Get_Value()
      day = info.dayID -> Get_Value()
      hour = info.hourID -> Get_Value()
      minute = info.minuteID -> Get_Value()
      distance = info.distanceID -> Get_Value()
      
      IF N_Elements(hour) eq 0 then hour = 0
      IF N_Elements(minute) eq 0 then minute = 0
      
      if N_Elements(year) eq 0 or N_Elements(month) eq 0 or $
      N_Elements(day) eq 0 then begin
      
      msg1 = 'Date Cannot Contain Blanks'
      
      ok = DIALOG_MESSAGE(msg1,/ERROR,/center)
      
      endif else begin

      (*info.ptr).year = year[0]
      (*info.ptr).month = month[0]
      (*info.ptr).day = day[0]
      (*info.ptr).hour = hour[0]
      (*info.ptr).minute = minute[0]
      (*info.ptr).distance = distance[0]
      (*info.ptr).cancel = 0
      Widget_Control, event.top, /Destroy
      endelse
      ENDCASE
ENDCASE
END ;-----------------------------------------------------





function swpc_cat_datebox_add_manual_point, Title=title, Label=label, Cancel=cancel, $
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


label = Widget_Label(row1base, Value='Date/Time [Y M D  H M]')
yearID = swpc_cat_FSC_InputField(row2base, Value=current_Year_starting, /IntegerValue, digits = 4 , xsize = 4 , title='')
monthID = swpc_cat_FSC_InputField(row2base, Value=current_Month_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')
dayID = swpc_cat_FSC_InputField(row2base, Value=current_Day_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')
hourID = swpc_cat_FSC_InputField(row2base, Value=current_Hour_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')
minuteID = swpc_cat_FSC_InputField(row2base, Value=current_Minute_starting, /IntegerValue, digits = 2 , xsize = 2 , title='')

distance_label = Widget_Label(row3base, Value='Distance')
distanceID = swpc_cat_FSC_InputField(row4base, Value='1.0', /FloatValue, digits = 4 , xsize = 4 , title='')

buttonBase = Widget_Base(tlb, Row=1)
cancelID = Widget_Button(buttonBase, Value='Cancel')
acceptID = Widget_Button(buttonBase, Value='Accept')

  ; Center the widgets on display.

swpc_cat_datebox_add_manual_point_CenterTLB, tlb
Widget_Control, tlb, /Realize

   ; Create a pointer for the text the user will type into the program.
   ; The cancel field is set to 1 to indicate that the user canceled
   ; the operation. Only if a successful conclusion is reached (i.e.,
   ; a Carriage Return or Accept button selection) is the cancel field
   ; set to 0.

ptr = Ptr_New({year:"", $
               month:"", $
               day:"", $
               hour:"", $
               minute:"", $
               distance:"", $
               cancel:1})

   ; Store the program information:

info = {ptr:ptr, $
        yearID:yearID, $
        monthID:monthID, $
        dayID:dayID, $
        hourID:hourID, $
        minuteID:minuteID, $
        distanceID:distanceID, $
        cancelID:cancelID}

Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Blocking or modal widget, depending upon group leader.

XManager, 'swpc_cat_datebox_add_manual_point', tlb

   ; Return from block. Return the text to the caller of the program,
   ; taking care to clean up pointer and group leader, if needed.
   ; Set the cancel keyword.

year = (*ptr).year
month = (*ptr).month
day = (*ptr).day
hour = (*ptr).hour
minute = (*ptr).minute
distance = (*ptr).distance

if minute eq 'NULLVALUE' then minute = '0'


date_array = list(year,month,day,hour,minute,distance)

theText = date_array
cancel = (*ptr).cancel
Ptr_Free, ptr
IF destroy_groupleader THEN Widget_Control, groupleader, /Destroy

RETURN, theText
END ;-----------------------------------------------------









pro swpc_cat_manually_add_or_remove_point, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY
   
widget_control,info.plot_window_manually_add_point,get_value=the_text
   
if the_text eq 'Manually Add Point' then begin

output = swpc_cat_datebox_add_manual_point(Title='Add Point', Cancel=cancelled, date_array=info.start_date)
info.manual_point = output
use_manual_point = 1
year = output[0]
month = output[1]
day = output[2]
hour = output[3]
minute = output[4]
distance = float(output[5])

datetime_Julian = JULDAY(Month, Day, Year, Hour, Minute)

(info.CME_matches_image_Julian).add, datetime_Julian
(info.CME_matches_image_DateTime_string).add, '-1'
(info.CME_matches_image_telescope).add, 'Manual'
(info.CME_matches_image_Rs_leading_edge).add, distance
(info.CME_matches_image_Image_number).add, -1
info.CME_matches_image_CME_outline.add,-1

widget_control,info.plot_window_manually_add_point,set_value='Remove Manual Point'

endif else begin

the_index_to_be_removed = (info.CME_matches_image_telescope).FindValue('Manual')   

(info.CME_matches_image_Julian).remove, the_index_to_be_removed
(info.CME_matches_image_DateTime_string).remove, the_index_to_be_removed
(info.CME_matches_image_telescope).remove, the_index_to_be_removed
(info.CME_matches_image_Rs_leading_edge).remove, the_index_to_be_removed
(info.CME_matches_image_Image_number).remove, the_index_to_be_removed
info.CME_matches_image_CME_outline.remove, the_index_to_be_removed

widget_control,info.plot_window_manually_add_point,set_value='Manually Add Point'

endelse

swpc_cat_sort_out_the_line_plot_and_widgets, info

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

end







pro swpc_cat_remove_this_image, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

CASE event.id OF

   info.L_widget_remove_this_image : BEGIN

the_index_to_be_removed = info.L_current_image_number
   
(info.L_list_of_image_names).remove, the_index_to_be_removed
(info.L_list_of_image_data).remove, the_index_to_be_removed
(info.L_list_of_datetime_strings).remove, the_index_to_be_removed
(info.L_list_of_datetime_Julian).remove, the_index_to_be_removed
(info.L_list_of_full_time_strings).remove, the_index_to_be_removed
(info.L_list_of_image_exposure_times).remove, the_index_to_be_removed
(info.L_list_of_image_offsets).remove, the_index_to_be_removed
(info.L_list_of_image_scaling_factors).remove, the_index_to_be_removed
(info.L_list_of_HEEQ_coords).remove, the_index_to_be_removed
(info.L_list_of_pixel_scales).remove, the_index_to_be_removed
(info.L_list_of_rsuns).remove, the_index_to_be_removed
(info.L_list_of_Sun_satellite_distances).remove, the_index_to_be_removed
info.L_number_of_images = info.L_number_of_images - 1
info.L_current_image_number = info.L_current_image_number - 1
widget_control, info.L_widget_image_sequence_slider,set_slider_max = info.L_number_of_images
widget_control,info.L_widget_image_sequence_slider,set_value = info.L_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.L_current_image_number,info.L_background_image_number,info.L_difference_imaging, $
    info.L_list_of_image_data,info.L_image_saturation_value,info.L_coronagraph_image_object,info.L_border_image_object, $
    info.CME_matches_image_L_Image_number,info.L_current_background_color, $
    info.background_color,info.L_current_text_color,info.color_stereo_B,info.L_cme_outline,info.L_cme_MATCH_outline, $
    info.L_widget_outline_matches_image,info.CME_matches_image_L_CME_outline, $
    info.L_ut_string_object,info.L_list_of_full_time_strings,info.L_title_object,info.L_Window,info.L_both_views,0,0, info.i_log_scale
    
L_yvals = fltarr(n_elements(info.L_list_of_datetime_Julian)) + 0.2
info.L_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.L_plot->SetProperty, dataX = (info.L_list_of_datetime_Julian).toarray() - info.start_julian, dataY = L_yvals
info.L_plot->SetProperty, color=[100,100,255]
info.L_plot->GetProperty, symbol=L_thisSymbol
    
swpc_cat_sort_out_the_timeline_symbols, info.L_number_of_images, info.CME_matches_image_L_Image_number, $
                                  info.L_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline

swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color

info.images_timeline_window -> draw, info.images_timeline_view
   
   ENDCASE
   
   info.C_widget_remove_this_image : BEGIN
   
CASE info.currently_showing_lasco_C2_or_C3 OF

   'C3' : BEGIN

the_index_to_be_removed = info.C_current_image_number
   
(info.C_list_of_image_names).remove, the_index_to_be_removed
(info.C_list_of_image_data).remove, the_index_to_be_removed
(info.C_list_of_datetime_strings).remove, the_index_to_be_removed
(info.C_list_of_datetime_Julian).remove, the_index_to_be_removed
(info.C_list_of_full_time_strings).remove, the_index_to_be_removed
(info.C_list_of_image_exposure_times).remove, the_index_to_be_removed
(info.C_list_of_image_offsets).remove, the_index_to_be_removed
(info.C_list_of_image_scaling_factors).remove, the_index_to_be_removed
(info.C_list_of_HEEQ_coords).remove, the_index_to_be_removed
(info.C_list_of_pixel_scales).remove, the_index_to_be_removed
(info.C_list_of_rsuns).remove, the_index_to_be_removed
(info.C_list_of_Sun_satellite_distances).remove, the_index_to_be_removed
info.C_number_of_images = info.C_number_of_images - 1
info.C_current_image_number = info.C_current_image_number - 1
widget_control, info.C_widget_image_sequence_slider,set_slider_max = info.C_number_of_images
widget_control,info.C_widget_image_sequence_slider,set_value = info.C_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.C_current_image_number,info.C_background_image_number,info.C_difference_imaging, $
    info.C_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c3,info.C_cme_outline,info.C_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C_CME_outline, $
    info.C_ut_string_object,info.C_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale
    
C_yvals = fltarr(n_elements(info.C_list_of_datetime_Julian)) + 0.65
info.C_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.C_plot->SetProperty, dataX = (info.C_list_of_datetime_Julian).toarray() - info.start_julian, dataY = C_yvals
info.C_plot->SetProperty, color=[100,255,100]
info.C_plot->GetProperty, symbol=C_thisSymbol
    
swpc_cat_sort_out_the_timeline_symbols, info.C_number_of_images, info.CME_matches_image_C_Image_number, $
                                  info.C_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color

info.images_timeline_window -> draw, info.images_timeline_view

ENDCASE

   'C2' : BEGIN
   
the_index_to_be_removed = info.C2_current_image_number
   
(info.C2_list_of_image_names).remove, the_index_to_be_removed
(info.C2_list_of_image_data).remove, the_index_to_be_removed
(info.C2_list_of_datetime_strings).remove, the_index_to_be_removed
(info.C2_list_of_datetime_Julian).remove, the_index_to_be_removed
(info.C2_list_of_full_time_strings).remove, the_index_to_be_removed
(info.C2_list_of_image_exposure_times).remove, the_index_to_be_removed
(info.C2_list_of_image_offsets).remove, the_index_to_be_removed
(info.C2_list_of_image_scaling_factors).remove, the_index_to_be_removed
(info.C2_list_of_HEEQ_coords).remove, the_index_to_be_removed
(info.C2_list_of_pixel_scales).remove, the_index_to_be_removed
(info.C2_list_of_rsuns).remove, the_index_to_be_removed
(info.C2_list_of_Sun_satellite_distances).remove, the_index_to_be_removed
info.C2_number_of_images = info.C2_number_of_images - 1
info.C2_current_image_number = info.C2_current_image_number - 1
widget_control, info.C_widget_image_sequence_slider,set_slider_max = info.C2_number_of_images
widget_control,info.C_widget_image_sequence_slider,set_value = info.C2_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.C2_current_image_number,info.C2_background_image_number,info.C2_difference_imaging, $
    info.C2_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C2_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c2,info.C_cme_outline,info.C2_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C2_CME_outline, $
    info.C_ut_string_object,info.C2_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale
    
C2_yvals = fltarr(n_elements(info.C2_list_of_datetime_Julian)) + 0.4
info.C2_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.C2_plot->SetProperty, dataX = (info.C2_list_of_datetime_Julian).toarray() - info.start_julian, dataY = C2_yvals
info.C2_plot->SetProperty, color=info.color_c2
info.C2_plot->GetProperty, symbol=C2_thisSymbol
    
swpc_cat_sort_out_the_timeline_symbols, info.C2_number_of_images, info.CME_matches_image_C2_Image_number, $
                                  info.C2_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline

swpc_cat_set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, info.C2_current_image_number, info.color_C2, info.cme_outline_color

info.images_timeline_window -> draw, info.images_timeline_view
   
ENDCASE

ENDCASE  

   
   ENDCASE
   
   info.R_widget_remove_this_image : BEGIN
   
the_index_to_be_removed = info.R_current_image_number
   
(info.R_list_of_image_names).remove, the_index_to_be_removed
(info.R_list_of_image_data).remove, the_index_to_be_removed
(info.R_list_of_datetime_strings).remove, the_index_to_be_removed
(info.R_list_of_datetime_Julian).remove, the_index_to_be_removed
(info.R_list_of_full_time_strings).remove, the_index_to_be_removed
(info.R_list_of_image_exposure_times).remove, the_index_to_be_removed
(info.R_list_of_image_offsets).remove, the_index_to_be_removed
(info.R_list_of_image_scaling_factors).remove, the_index_to_be_removed
(info.R_list_of_HEEQ_coords).remove, the_index_to_be_removed
(info.R_list_of_pixel_scales).remove, the_index_to_be_removed
(info.R_list_of_rsuns).remove, the_index_to_be_removed
(info.R_list_of_Sun_satellite_distances).remove, the_index_to_be_removed
info.R_number_of_images = info.R_number_of_images - 1
info.R_current_image_number = info.R_current_image_number - 1
widget_control, info.R_widget_image_sequence_slider,set_slider_max = info.R_number_of_images
widget_control,info.R_widget_image_sequence_slider,set_value = info.R_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.R_current_image_number,info.R_background_image_number,info.R_difference_imaging, $
    info.R_list_of_image_data,info.R_image_saturation_value,info.R_coronagraph_image_object,info.R_border_image_object, $
    info.CME_matches_image_R_Image_number,info.R_current_background_color, $
    info.background_color,info.R_current_text_color,info.color_stereo_B,info.R_cme_outline,info.R_cme_MATCH_outline, $
    info.R_widget_outline_matches_image,info.CME_matches_image_R_CME_outline, $
    info.R_ut_string_object,info.R_list_of_full_time_strings,info.R_title_object,info.R_Window,info.R_both_views,0,0, info.i_log_scale
    
R_yvals = fltarr(n_elements(info.R_list_of_datetime_Julian)) + 0.9
info.R_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.R_plot->SetProperty, dataX = (info.R_list_of_datetime_Julian).toarray() - info.start_julian, dataY = R_yvals
info.R_plot->SetProperty, color=[255,100,100]
info.R_plot->GetProperty, symbol=R_thisSymbol
    
swpc_cat_sort_out_the_timeline_symbols, info.R_number_of_images, info.CME_matches_image_R_Image_number, $
                                  info.R_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_A, info.cme_outline_color

info.images_timeline_window -> draw, info.images_timeline_view

    
   
   ENDCASE


ENDCASE


WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

end









pro swpc_cat_define_representative_image, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

CASE event.id OF

   info.L_widget_representative_image : BEGIN
   
   info.L_Window -> getproperty, image_data = image_data
   
   ENDCASE
   
   info.C_widget_representative_image : BEGIN
   
   info.C_Window -> getproperty, image_data = image_data
   
   ENDCASE
   
   info.R_widget_representative_image : BEGIN
   
   info.R_Window -> getproperty, image_data = image_data
   
   ENDCASE


ENDCASE

;for ix = 0, 511 do begin
;for iy = 0, 511 do begin
;for iz = 0, 2 do begin
;info.representative_image_data[ix,iy,iz] = image_data[iz,ix,iy]
;endfor
;endfor
;endfor

; image_data is ordered z,x,y - but we need the representative image to
; be x,y,z - this transposition is done efficiently with the transpose
; function.....

info.representative_image_data = transpose(image_data,[1,2,0])
info.representative_image_has_been_defined = 1

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

end





pro swpc_cat_REDRAW_THE_IMAGE, $
    current_image_number,background_image_number,difference_imaging, $
    list_of_image_data,image_saturation_value,coronagraph_image_object, border_image_object, $
    CME_matches_image_TELESCOPE_Image_number,current_background_color, $
    background_color,current_text_color,telescope_text_color,cme_outline,cme_MATCH_outline, $
    widget_outline_matches_image,CME_matches_image_TELESCOPE_CME_outline, $
    ut_string_object,list_of_full_time_strings,title_object,Window,view,hide_the_cme,animating, i_log_scale

; Redraw the image taking into account whether this image is matched or not

if animating eq 0 then begin

if float(!version.release) gt 8.05 then begin
   Result = (CME_matches_image_TELESCOPE_Image_number).Count(current_image_number)
   Result2 = (CME_matches_image_TELESCOPE_Image_number).FindValue(current_image_number)
endif else begin
   Result2 = where(CME_matches_image_TELESCOPE_Image_number eq current_image_number, Result)
endelse
if result eq 0 then begin
current_background_color = background_color
current_text_color = telescope_text_color
cme_outline -> setProperty, hide = hide_the_cme
cme_MATCH_outline->SetProperty, hide = 1
widget_control,widget_outline_matches_image,set_value='CME Matches Image'
endif else begin
current_background_color = telescope_text_color
current_text_color = background_color
cme_outline -> setProperty, hide = 1
if float(!version.release) gt 8.15 then result2 = result2[0] ; fix for idl8.2 which creates an array(1) for this variable
cme_MATCH_outline->SetProperty, data = (CME_matches_image_TELESCOPE_CME_outline)[result2]
cme_MATCH_outline->SetProperty, hide = 0
widget_control,widget_outline_matches_image,set_value='Unmatch'
endelse

endif

swpc_cat_image_difference_and_scaling, current_background_color,  current_image_number, background_image_number, difference_imaging, $
                 list_of_image_data, image_saturation_value, coronagraph_image_object, border_image_object, i_log_scale                               
ut_string_object->SetProperty, strings = list_of_full_time_strings[current_image_number]
ut_string_object->setProperty, color = current_text_color
title_object->setProperty, color = current_text_color
                          
if animating eq 0 then Window->Draw, view

end











pro swpc_cat_reset_cme_analysis, event
WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

;       Widget_Control, event.id, Draw_Motion_Events=1
	Widget_Control, info.draw_C, Draw_Motion_Events=1

       reset_for_new_images_as_well_as_the_cme_analysis = 0
       swpc_cat_full_reset, info, reset_for_new_images_as_well_as_the_cme_analysis
       
; Resetting the analysis also means greying out the velocity and export buttons....
       
       widget_control, info.widget_calculate_velocity, sensitive = 0
       widget_control, info.widget_export_cone, sensitive = 0
; ECDP, but resetting the right click buttons, regardless of the mode we were in prior to the reset      
       widget_control, info.widget_calculate_pos_velocity, sensitive = 0
       widget_control, info.widget_use_emp_cone, sensitive = 0       
       widget_control, info.C_widget_pos_outline_matches_image, sensitive = 1
       widget_control, info.C_widget_outline_matches_image, sensitive = 1
       widget_control, info.C_widget_representative_image, sensitive = 1
       widget_control, info.L_widget_outline_matches_image, sensitive = 1
       widget_control, info.L_widget_representative_image, sensitive = 1
       widget_control, info.R_widget_outline_matches_image, sensitive = 1
       widget_control, info.R_widget_representative_image, sensitive = 1
;       Widget_Control, event.id, Draw_Motion_Events=1
; ECDP       
WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
END




pro swpc_cat_sort_out_the_timeline_symbols, number_of_images, CME_matches_list, $
                                  the_plot, xSymbolSize_timeline, ySymbolSize_timeline  
   
   symarray = objarr(number_of_images)
   for i = 0 , number_of_images - 1 do begin
   
   if float(!version.release) gt 8.05 then begin
      Result = (CME_matches_list).Count(i)
   endif else begin
      !null = where(CME_matches_list eq i, Result)
   endelse
   
   if result eq 0 then begin
   num = 6
   filled = 0
   endif else begin
   num = 6
   filled = 1
   endelse
   
   thisSymbol = obj_new("IDLgrsymbol",data=num, Size=[xSymbolSize_timeline, ySymbolSize_timeline],filled = filled)
   symarray[i] = thisSymbol
   
   endfor
   the_plot->SetProperty, symbol=symarray

end



pro swpc_cat_set_image_border, color, border_image_object

border_image_object -> getproperty, data = image_data2
       
image_data2[0,*,0:18] = color[0]
image_data2[1,*,0:18] = color[1]
image_data2[2,*,0:18] = color[2]
image_data2[0,*,493:511] = color[0]
image_data2[1,*,493:511] = color[1]
image_data2[2,*,493:511] = color[2]
image_data2[0,0:18,*] = color[0]
image_data2[1,0:18,*] = color[1]
image_data2[2,0:18,*] = color[2]
image_data2[0,493:511,*] = color[0]
image_data2[1,493:511,*] = color[1]
image_data2[2,493:511,*] = color[2]

border_image_object -> setproperty, data = image_data2

end





pro swpc_cat_define_cme_lemniscate, radial_distance, angular_width, lemniscate_style, fitted_cme_info
                            
xmax = radial_distance
x_pt = xmax/5.0
tatx = acos(0.2)
ymax = x_pt*tan((angular_width/2.0)*!dtor)
zmax = x_pt*tan((angular_width/2.0)*!dtor)

c1 = xmax
c2 = ymax/(cos(tatx)*sin(tatx))
c3 = 1.0
; c3=1 will produce a symmetric, or circular, front
; c3 ne 1 will produce an elliptical front
;-----------------------------

lon = 0.*!dtor
lat = 90.*!dtor
tilt = 0.

polysurf_params = [lon,lat,tilt,c1,c2,c3]
polysurf_const = polysurf_params[3:5] ; [c1,c2,c3]
polysurf_dirs  = polysurf_params[0:2] ; [lon,lat,tilt]

fitted_cme_info = swpc_cat_general_lemniscate_define(polysurf_const,polysurf_dirs, $
                                            lemniscate_style,             $
                                            xverts=30,yverts=90)

end





pro swpc_cat_full_reset, info, reset_for_new_images_as_well_as_the_cme_analysis

; ECDP
if (info.compact eq 0) then begin
widget_control, info.c_mode_name, set_value = 'SELECT MODE'
endif
info.which_mode = 0
; hide the clock
info.C_clock_angle_degrees = 0.
info.C_leading_edge_in_rs = 0.
swpc_cat_set_clock_angle,info,info.C_clock_angle_degrees

info.C_clock_angle_model->SetProperty, transform = info.C_clock_angle_model_original_transform
;print, info.C_clock_angle_model_original_transform 
info.C_clock_angle_model -> SetProperty, Hide = 1
info.C_leading_edge_model -> SetProperty, Hide = 1

info.representative_image_has_been_defined = 0

info.cone_radius = 0.

info.rotation = 0.
info.pixel_scale = 0.
info.scaling_factor = 0.


info.rotation = 0.
info.rsun = 950.   ;set initial value, solar radius from L1 in arcseconds
info.pixel_scale = 1.
info.scaling_factor = 1.0
info.rotate_3D_cme =  1   ; trackball won't work.

info.initial_transform = dblarr(4,4)
info.initial_transform[0,0] = 1.
info.initial_transform[1,1] = 1.
info.initial_transform[2,2] = 1.
info.initial_transform[3,3] = 1.

info.latitude_degrees = 0.
info.longitude_degrees = 0.

info.CME_matches_image_Julian = list()
info.CME_matches_image_DateTime_string = list()
info.CME_matches_image_telescope = list()
info.CME_matches_image_Rs_leading_edge = list()
info.CME_matches_image_Image_number = list()
info.CME_matches_image_CME_outline = list()
info.CME_matches_image_L_Image_number = list()
info.CME_matches_image_C_Image_number = list()
info.CME_matches_image_C2_Image_number = list()
info.CME_matches_image_R_Image_number = list()
info.CME_matches_image_L_CME_outline = list()
info.CME_matches_image_C_CME_outline = list()
info.CME_matches_image_C2_CME_outline = list()
info.CME_matches_image_R_CME_outline = list()

info.POS_Julian = list()
info.POS_leading_edge = list()
info.POS_matches_image_DateTime_string = list()
info.POS_matches_image_telescope = list()

info.L_current_background_color = info.background_color
info.L_current_text_color = info.color_stereo_B
;info.L_cme_outline -> setProperty, hide = 0
info.L_cme_outline -> setProperty, hide = 1
info.L_cme_MATCH_outline->SetProperty, hide = 1

swpc_cat_set_image_border, info.background_color, info.L_border_image_object

info.C_current_background_color = info.background_color

info.C_cme_outline -> setProperty, hide = 1
info.C_cme_MATCH_outline->SetProperty, hide = 1
info.C2_cme_MATCH_outline->SetProperty, hide = 1

; ECDP - need to set this to zero so that 'motion' without click is switched on again.
info.C_the_action = 0
info.pressed_the_right_button = 0
info.show_image_line_plot = 1
info.C_latest_click_X = 0
info.C_latest_click_Y = 0
info.C_previous_click_X = 0
info.C_previous_click_Y = 0
info.C_click_and_drag  = 0
info.C_rotate_x = 0
info.C_rotate_y = 0
info.timeline_left_mouse_button_being_pressed = 0
info.clicked_L = 0
info.clicked_C = 0
info.clicked_C2 = 0
info.clicked_R = 0
;Widget_Control, Draw_Motion_Events=1
; ECDP

info.R_current_background_color = info.background_color
info.R_current_text_color = info.color_stereo_A
;info.R_cme_outline -> setProperty, hide = 0
info.R_cme_outline -> setProperty, hide = 1
info.R_cme_MATCH_outline->SetProperty, hide = 1

swpc_cat_set_image_border, info.background_color, info.R_border_image_object

if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
info.C_current_text_color = info.color_c3
endif else begin
info.C_current_text_color = info.color_c2
endelse

swpc_cat_set_image_border, info.background_color, info.C_border_image_object

info.L_ut_string_object->setProperty, color = info.L_current_text_color
info.L_title_object->setProperty, color = info.L_current_text_color
info.C_ut_string_object->setProperty, color = info.C_current_text_color
info.C_title_object->setProperty, color = info.C_current_text_color
info.R_ut_string_object->setProperty, color = info.R_current_text_color
info.R_title_object->setProperty, color = info.R_current_text_color


info.anim_start_frame = 1
info.anim_end_frame = 1
info.i_anim_back_and_forth = 0
info.i_move = 1

info.which_window_to_animate = 1

; ECDP - in the new version, we want the sliders to be insensitive until a mode is selected
; (we also want to hide the clock and lemniscate until a mode has been selected as this forces the user
; to actively select the mode they require) 
widget_control,info.widget_latitude_slider, sensitive = 0
widget_control,info.widget_longitude_slider, sensitive = 0
widget_control,info.widget_angular_width_slider, sensitive = 0
widget_control,info.widget_radial_distance_slider, sensitive = 0
; ECDP

info.thisPlot_B2->SetProperty, hide = 1
info.thisPlot2_c3->SetProperty, hide = 1
info.thisPlot2_c2->SetProperty, hide = 1
info.thisPlot2_stereo_a->SetProperty, hide = 1
info.thisPlot2_stereo_b->SetProperty, hide = 1
info.latitude_text->SetProperty, hide = 1
info.longitude_text->SetProperty, hide = 1
info.cone_angle_text->SetProperty, hide = 1
info.velocity_text->SetProperty, hide = 1
; ECDP
info.pos_velocity_text->SetProperty, hide = 1
info.emp_cone_angle_text->SetProperty, hide = 1
; ECDP
info.T21_5_text->SetProperty, hide = 1

info.radial_distance_lemniscate = 8.0
info.angular_width_lemniscate = 90.

swpc_cat_define_cme_lemniscate, info.radial_distance_lemniscate, info.angular_width_lemniscate, info.lemniscate_style, fitted_cme_info

; ECDP we have re-drawn the lemnisacte, but now we want to hide it until the mode is selected
info.L_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.L_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.L_Window_copy->Draw, info.L_second_view_copy
swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline
info.L_Window->Draw, info.L_both_views

info.C_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.C_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.C_Window_copy->Draw, info.C_second_view_copy
swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
info.C_Window->Draw, info.C_both_views

info.R_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.R_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.R_Window_copy->Draw, info.R_second_view_copy
swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline
info.R_Window->Draw, info.R_both_views

info.latitude_degrees = 0.
info.longitude_degrees = 0.

info.L_cme_model->SetProperty, transform = info.initial_transform
info.L_cme_model_copy->SetProperty, transform = info.initial_transform

info.C_cme_model->SetProperty, transform = info.initial_transform
info.C_cme_model_copy->SetProperty, transform = info.initial_transform

info.R_cme_model->SetProperty, transform = info.initial_transform
info.R_cme_model_copy->SetProperty, transform = info.initial_transform

swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline
info.L_Window->Draw, info.L_both_views

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
info.C_Window->Draw, info.C_both_views

swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline
info.R_Window->Draw, info.R_both_views

;info.C_cme_outline->SetProperty, Hide = 1
;info.L_cme_outline->SetProperty, Hide = 1
;info.R_cme_outline->SetProperty, Hide = 1

info.L_cme_MATCH_outline->SetProperty, hide = 1
info.C_cme_MATCH_outline->SetProperty, hide = 1
info.C2_cme_MATCH_outline->SetProperty, hide = 1
info.R_cme_MATCH_outline->SetProperty, hide = 1

widget_control,info.L_widget_outline_matches_image,set_value='CME Matches Image'
widget_control,info.C_widget_outline_matches_image,set_value='CME Matches Image'
widget_control,info.R_widget_outline_matches_image,set_value='CME Matches Image'

widget_control, info.widget_latitude_slider, set_value = 0

widget_control, info.widget_longitude_slider, set_value = 0

widget_control, info.widget_angular_width_slider, set_value = 90

widget_control, info.widget_radial_distance_slider, set_value = 80

info.lat_string = string(info.latitude_degrees,format='(f6.1)')
info.lat_string_object -> setproperty, strings = 'lat  :' + info.lat_string

info.lon_string = string(info.longitude_degrees,format='(f6.1)')
info.lon_string_object -> setproperty, strings = 'lon  :' + info.lon_string

radial_distance_string = 'dist :' + string(info.radial_distance_lemniscate, format='(f6.1)')
info.radial_distance_string_object -> setproperty, strings = radial_distance_string

angular_width_string = 'cone :' + string(info.angular_width_lemniscate, format='(f6.1)')
info.angular_width_string_object -> setproperty, strings = angular_width_string

info.plot_window->erase, color=info.background_color
info.cme_info_Window->erase, color=info.background_color
info.enlil_info_window->erase, color=info.background_color

;info.plot_window->Draw, info.plot_View
info.cme_info_Window->Draw, info.cme_info_view

if info.L_number_of_images gt 0 then begin
swpc_cat_sort_out_the_timeline_symbols, info.L_number_of_images, info.CME_matches_image_L_Image_number, $
                                  info.L_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline
endif

if info.C_number_of_images gt 0 then begin
swpc_cat_sort_out_the_timeline_symbols, info.C_number_of_images, info.CME_matches_image_C_Image_number, $
                                  info.C_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline
endif
                                  
if info.C2_number_of_images gt 0 then begin
swpc_cat_sort_out_the_timeline_symbols, info.C2_number_of_images, info.CME_matches_image_C2_Image_number, $
                                  info.C2_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline                 
endif

if info.R_number_of_images gt 0 then begin
swpc_cat_sort_out_the_timeline_symbols, info.R_number_of_images, info.CME_matches_image_R_Image_number, $
                                  info.R_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline
endif
                  
info.images_timeline_window->Draw, info.images_timeline_view

; ECDP reset the title to the standard one (even though it doesn't make much difference)
if (info.compact eq 0) then begin
widget_control, info.time_plot_title_widget, set_value = 'CME LEADING EDGE (Rs) vs TIME (UT)'
endif

if reset_for_new_images_as_well_as_the_cme_analysis eq 1 then begin

info.L_list_of_image_names = list()
info.L_list_of_image_data = list()
info.L_list_of_datetime_strings = list()
info.L_list_of_datetime_Julian = list()
info.L_list_of_full_time_strings = list()
info.L_list_of_image_exposure_times = list()
info.L_list_of_image_offsets = list()
info.L_list_of_image_scaling_factors = list()
info.C_list_of_image_names = list()
info.C_list_of_image_data = list()
info.C_list_of_datetime_strings = list()
info.C_list_of_datetime_Julian = list()
info.C_list_of_full_time_strings = list()
info.C_list_of_image_exposure_times = list()
info.C_list_of_image_offsets = list()
info.C_list_of_image_scaling_factors = list()
info.C2_list_of_image_names = list()
info.C2_list_of_image_data = list()
info.C2_list_of_datetime_strings = list()
info.C2_list_of_datetime_Julian = list()
info.C2_list_of_full_time_strings = list()
info.C2_list_of_image_exposure_times = list()
info.C2_list_of_image_offsets = list()
info.C2_list_of_image_scaling_factors = list()
info.R_list_of_image_names = list()
info.R_list_of_image_data = list()
info.R_list_of_datetime_strings = list()
info.R_list_of_datetime_Julian = list()
info.R_list_of_full_time_strings = list()
info.R_list_of_image_exposure_times = list()
info.R_list_of_image_offsets = list()
info.R_list_of_image_scaling_factors = list()
info.L_list_of_HEEQ_coords = list()
info.C_list_of_HEEQ_coords = list()
info.C2_list_of_HEEQ_coords = list()
info.R_list_of_HEEQ_coords = list()
info.L_HEEQ_coords = [0.,0.]
info.C_HEEQ_coords = [0.,0.]
info.C2_HEEQ_coords = [0.,0.]
info.R_HEEQ_coords = [0.,0.]
info.master_list = list()
info.L_indexes = list()
info.C_indexes = list()
info.C2_indexes = list()
info.R_indexes = list()
info.master_list_size = 0
info.L_current_image_number = 1  ; start with the second image in the sequence
info.C_current_image_number = 1  ; (ie, the 1st one that shows anything with
info.C2_current_image_number = 1 ; difference imaging)....
info.R_current_image_number = 1
info.L_list_of_pixel_scales = list()
info.C_list_of_pixel_scales = list()
info.C2_list_of_pixel_scales = list()
info.R_list_of_pixel_scales = list()
info.L_list_of_rsuns = list()
info.C_list_of_rsuns = list()
info.C2_list_of_rsuns = list()
info.R_list_of_rsuns = list()
info.L_list_of_Sun_satellite_distances = list()
info.C_list_of_Sun_satellite_distances = list()
info.C2_list_of_Sun_satellite_distances = list()
info.R_list_of_Sun_satellite_distances = list()
info.L_telescope_code = ''
info.C_telescope_code = ''
info.C2_telescope_code = ''
info.R_telescope_code = ''
info.date_array_int = intarr(10)
info.start_date = intarr(5)
info.end_date = intarr(5)
info.start_julian = -1.0D
info.end_julian = -1.0D

info.L_window->erase, color=info.background_color
info.C_window->erase, color=info.background_color
info.R_window->erase, color=info.background_color
info.images_timeline_window->erase, color=info.background_color

info.L_difference_imaging = 'running'
info.C_difference_imaging = 'running'
info.C2_difference_imaging = 'running'
info.R_difference_imaging = 'running'

info.stopflag=1
widget_control,info.play_pause_button,set_value='Play '
widget_control,info.widget_anim_back_and_forth,set_value=0
info.i_anim_back_and_forth=0
widget_control,info.widget_animation_delay,set_value=50
info.delay = 0.25
widget_control,info.widget_which_window_to_animate,set_value=1

endif


end



pro swpc_cat_set_timeline_highlight_block, plot, nvals, index, regular_color, highlight_color
vert_colors = intarr(3,nvals)
vert_colors[0,0:nvals-1] = regular_color[0]
vert_colors[1,0:nvals-1] = regular_color[1]
vert_colors[2,0:nvals-1] = regular_color[2]
vert_colors[0,index] = highlight_color[0]
vert_colors[1,index] = highlight_color[1]
vert_colors[2,index] = highlight_color[2]
plot->SetProperty, vert_colors=vert_colors
end



pro swpc_cat_update_cme_outline,Window_copy,second_view_copy,cme_outline
Window_copy->Draw, second_view_copy
Window_copy->GetProperty, Image_Data=snapshot
CONTOUR, snapshot[0,0:511,0:511], /PATH_DATA_COORDS, PATH_XY=data, nlevels=1
cme_outline->setProperty, data = data
end







pro swpc_cat_show_C2_or_C3, event
WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

widget_control,info.widget_show_C2_or_C3,get_value=the_text

if the_text eq 'Show LASCO C2' then begin

info.currently_showing_lasco_C2_or_C3 = 'C2'
widget_control,info.widget_show_C2_or_C3,set_value='Show LASCO C3'
swpc_cat_calculate_c2_c3_lead_edge, info
;print, 'Here! =', info.C_leading_edge_in_rs
; make sure, as we switch from C3 to C2, that any C3 match (green line)
; is hidden.... 
info.C_cme_MATCH_outline->SetProperty, hide = 1

widget_control, info.C_widget_image_sequence_slider,set_slider_max = n_elements(info.C2_list_of_datetime_Julian)

info.C_title_object -> setproperty, strings = 'SOHO LASCO C2'

swpc_cat_REDRAW_THE_IMAGE, $
    info.C2_current_image_number,info.C2_background_image_number,info.C2_difference_imaging, $
    info.C2_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C2_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c2,info.C_cme_outline,info.C2_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C2_CME_outline, $
    info.C_ut_string_object,info.C2_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_Calculate_Earth_B_Angle,(info.C2_list_of_datetime_Julian)[0],B_angle_degrees
info.C2_HEEQ_coords[1] = B_angle_degrees

info.C_coords_offset_model -> SetProperty,transform = info.initial_transform
info.C_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.C_coords_offset_model->rotate,[1,0,0], info.C2_HEEQ_coords[1]
info.C_coords_offset_model_copy->rotate,[1,0,0], info.C2_HEEQ_coords[1]

info.C2_telescope_field_of_view_Rs = (256. * ((info.C2_list_of_pixel_scales)[0] / (info.C2_list_of_image_scaling_factors)[0])) / (info.C2_list_of_rsuns)[0]

if info.debug_mode eq 1 then print, 'C2 ', info.C2_telescope_field_of_view_Rs, (info.C2_list_of_pixel_scales)[0], (info.C2_list_of_image_scaling_factors)[0], (info.C2_list_of_rsuns)[0]

info.C_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.C2_telescope_field_of_view_Rs,2.0*info.C2_telescope_field_of_view_Rs]
info.C_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.C2_telescope_field_of_view_Rs,2.0*info.C2_telescope_field_of_view_Rs]

the_day = long((info.C2_list_of_datetime_Julian)[0])
i_day = where(the_day lt info.Julian_day_for_Earth_pos)
i_day = i_day[0]

print, i_day, info.Earth_pos_AU[i_day]

if i_day eq -1 then begin
  print, 'Date not found in the Earth_position file, stopping!'
  print, 'comment out the stop statement in the line below for an emergency fix (SWPC_CAT will work with a slightly incorrect Earth_pos_AU)'
  stop ; << comment out this stop statement for an emegency fix (as in, you don't have a new Earth_pos file, but need SWPC_CAT to work ok)
endif

info.C_second_view -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99  ; 0.99 factor is for L1 as opposed to Earth.
info.C_second_view_copy -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99

info.C_Fixed_Sun_Satellite_model -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model -> translate, info.C2_telescope_field_of_view_Rs , info.C2_telescope_field_of_view_Rs, 0.0
info.C_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model_copy -> translate, info.C2_telescope_field_of_view_Rs , info.C2_telescope_field_of_view_Rs, 0.0
;info.C_cme_outline -> setProperty, hide = 0
; ECDP
if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
  info.C_cme_outline -> setProperty, hide = 1
endif else begin
  info.C_cme_outline -> setProperty, hide = 0
endelse
; ECDP  

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
info.C_Window->Draw, info.C_both_views

endif else begin

info.currently_showing_lasco_C2_or_C3 = 'C3'
widget_control,info.widget_show_C2_or_C3,set_value='Show LASCO C2'
swpc_cat_calculate_c2_c3_lead_edge, info
;print, 'Here!! =', info.C_leading_edge_in_rs
; make sure, as we switch from C2 to C3, that any C2 match (green line)
; is hidden.... 
info.C2_cme_MATCH_outline->SetProperty, hide = 1

widget_control, info.C_widget_image_sequence_slider,set_slider_max = n_elements(info.C_list_of_datetime_Julian)

info.C_title_object -> setproperty, strings = 'SOHO LASCO C3'

swpc_cat_REDRAW_THE_IMAGE, $
    info.C_current_image_number,info.C_background_image_number,info.C_difference_imaging, $
    info.C_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c3,info.C_cme_outline,info.C_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C_CME_outline, $
    info.C_ut_string_object,info.C_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_Calculate_Earth_B_Angle,(info.C_list_of_datetime_Julian)[0],B_angle_degrees
info.C_HEEQ_coords[1] = B_angle_degrees

info.C_coords_offset_model -> SetProperty,transform = info.initial_transform
info.C_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.C_coords_offset_model->rotate,[1,0,0], info.C_HEEQ_coords[1]
info.C_coords_offset_model_copy->rotate,[1,0,0], info.C_HEEQ_coords[1]

info.C_telescope_field_of_view_Rs = (256. * ((info.C_list_of_pixel_scales)[0] / (info.C_list_of_image_scaling_factors)[0])) / (info.C_list_of_rsuns)[0]

if info.debug_mode eq 1 then print, 'C3 ', info.C_telescope_field_of_view_Rs, (info.C_list_of_pixel_scales)[0], (info.C_list_of_image_scaling_factors)[0], (info.C2_list_of_rsuns)[0]

info.C_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]
info.C_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]

the_day = long((info.C_list_of_datetime_Julian)[0])
i_day = where(the_day lt info.Julian_day_for_Earth_pos)
i_day = i_day[0]

info.C_second_view -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99  ; 0.99 factor is for L1 as opposed to Earth.
info.C_second_view_copy -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99

info.C_Fixed_Sun_Satellite_model -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
info.C_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model_copy -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
;info.C_cme_outline -> setProperty, hide = 0
; ECDP
if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
  info.C_cme_outline -> setProperty, hide = 1
endif else begin
  info.C_cme_outline -> setProperty, hide = 0
endelse
; ECDP  
swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline

;swpc_cat_Window_button_click_events, event
info.C_Window->Draw, info.C_both_views

endelse

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
END






pro swpc_cat_make_end_date_Xhours_after_start, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

CASE event.id OF

   info.widget_plus_12h_button : delta_t = 0.5d
   info.widget_plus_24h_button : delta_t = 1.0d
   
ENDCASE


start_year = info.date_array[0]
start_month = info.date_array[1]
start_day = info.date_array[2]
start_hour = info.date_array[3]
start_minute = info.date_array[4]

start_julday = JULDAY(start_month, start_day, start_year, start_hour, start_minute)

end_julday = start_julday + delta_t

CALDAT, end_julday, end_month, end_day, end_year, end_hour, end_minute

info.date_array[5] = end_year
info.date_array[6] = end_month
info.date_array[7] = end_day
info.date_array[8] = end_hour
info.date_array[9] = end_minute

info.end_yearID -> Set_Value,strcompress(string(end_year),/remove_all)
info.end_monthID -> Set_Value,strcompress(string(end_month),/remove_all)
info.end_dayID -> Set_Value,strcompress(string(end_day),/remove_all)
info.end_hourID -> Set_Value,strcompress(string(end_hour),/remove_all)
info.end_minuteID -> Set_Value,strcompress(string(end_minute),/remove_all)

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

END







pro swpc_cat_datebox_Event, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

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
     
      if start_minute eq 'NULLVALUE' then start_minute = '0'
      if end_minute eq 'NULLVALUE' then end_minute = '0'


info.date_array = [start_year,start_month,start_day,start_hour,start_minute, $
                    end_year,end_month,end_day,end_hour,end_minute]

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

END




function swpc_cat_which_window_to_animate, event

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

CASE event.value OF

   0 : BEGIN

   info.which_window_to_animate = 0
   array_of_julians = info.L_list_of_datetime_Julian->toarray()

   ENDCASE
   1 : BEGIN

   info.which_window_to_animate = 1
   array_of_julians = info.C_list_of_datetime_Julian->toarray()

   ENDCASE
   2 : BEGIN

   info.which_window_to_animate = 2
   array_of_julians = info.R_list_of_datetime_Julian->toarray()

   ENDCASE
   
ENDCASE

l = where(array_of_julians gt info.anim_start_julian)
info.anim_start_frame = l[0]
info.thisFrame = info.anim_start_frame
l = where(array_of_julians lt info.anim_end_julian)
info.anim_end_frame = l[-1]

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY

Return, 1
End





pro swpc_cat_set_animation_start, event
WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

this_julian = (((info.timeline_normalized_x - info.position_timeline[0]) / (info.position_timeline[2] - info.position_timeline[0])) * (info.end_julian - info.start_julian)) + info.start_julian

animation_start_position = ((this_julian - info.start_julian)/(info.end_julian - info.start_julian) * 0.9) + 0.05
line_data=fltarr(2,2)
line_data[0,0] = animation_start_position
line_data[1,0] = 0.
line_data[0,1] = animation_start_position
line_data[1,1] = 1.
info.animation_start_time_marker->setproperty,data=line_data

xpos = animation_start_position
info.animation_start_time_marker_handle_data=[[xpos,0.],[xpos+0.005,0.],[xpos+0.005,0.08],[xpos,0.08]]
info.animation_start_time_marker_handle->setproperty,data=info.animation_start_time_marker_handle_data

info.images_timeline_window->Draw, info.images_timeline_view

CASE info.which_window_to_animate OF
   0 : BEGIN
array_of_julians = info.L_list_of_datetime_Julian->toarray()
   ENDCASE
   1 : BEGIN
array_of_julians = info.C_list_of_datetime_Julian->toarray()
   ENDCASE
   2 : BEGIN
array_of_julians = info.R_list_of_datetime_Julian->toarray()
   ENDCASE
   
ENDCASE



l = where(array_of_julians gt this_julian)
info.anim_start_frame = l[0]
info.anim_start_julian = this_julian
info.thisFrame = info.anim_start_frame

if info.debug_mode eq 1 then print, this_julian, format='(f20.10)'
CALDAT, this_julian, Month, Day, Year, Hour, Minute
if info.debug_mode eq 1 then print, year, month, day, hour, minute
if info.debug_mode eq 1 then print, array_of_julians, format='(3f20.10)'
if info.debug_mode eq 1 then print, l[0]
if info.debug_mode eq 1 then print, array_of_julians[l[0]],format='(f20.10)'

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
END


pro swpc_cat_set_animation_end, event
WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

this_julian = (((info.timeline_normalized_x - info.position_timeline[0]) / (info.position_timeline[2] - info.position_timeline[0])) * (info.end_julian - info.start_julian)) + info.start_julian

animation_end_position = ((this_julian - info.start_julian)/(info.end_julian - info.start_julian) * 0.9) + 0.05
line_data=fltarr(2,2)
line_data[0,0] = animation_end_position
line_data[1,0] = 0.
line_data[0,1] = animation_end_position
line_data[1,1] = 1.
info.animation_end_time_marker->setproperty,data=line_data

xpos = animation_end_position
info.animation_end_time_marker_handle_data=[[xpos,0.],[xpos-0.005,0.],[xpos-0.005,0.08],[xpos,0.08]]
info.animation_end_time_marker_handle->setproperty,data=info.animation_end_time_marker_handle_data

info.images_timeline_window->Draw, info.images_timeline_view

CASE info.which_window_to_animate OF

   0 : BEGIN
array_of_julians = info.L_list_of_datetime_Julian->toarray()
   ENDCASE
   1 : BEGIN
array_of_julians = info.C_list_of_datetime_Julian->toarray()
   ENDCASE
   2 : BEGIN
array_of_julians = info.R_list_of_datetime_Julian->toarray()
   ENDCASE
   
ENDCASE

l = where(array_of_julians lt this_julian)
info.anim_end_frame = l[-1]
info.anim_end_julian = this_julian

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
END



;------------------------------------------------------------------------------
function swpc_cat_xmovie_back_and_forth, event
WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY
if info.i_anim_back_and_forth eq 0 then begin
info.i_anim_back_and_forth = 1
endif else begin
info.i_anim_back_and_forth = 0
endelse
WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
return,1
END
;------------------------------------------------------------------------------
pro swpc_cat_XMOVIE_DELAY_EVENTS, event

   ; Get the INFO structure.

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

   ; Set the delay.

info.delay = float((100 - event.value)) / 200.0 ; Max delay 0.5 seconds.

      ; Put INFO structure back.

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY


END ; ************************ of XMOVIE_DELAY_EVENTS ***********************
;------------------------------------------------------------------------------


pro swpc_cat_XMOVIE_EVENT, event

   ; Get the INFO structure.

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

   ; What kind of event is this?

eventName = TAG_NAMES(event, /STRUCTURE_NAME)

   ; Code for BUTTON events.

IF eventName EQ 'WIDGET_BUTTON' THEN BEGIN

widget_control,info.play_pause_button,get_value = play_or_pause_text

   CASE play_or_pause_text OF

      'Play ' : BEGIN
widget_control,info.play_pause_button,set_value = 'Pause'

            ; If stop flag is 0, then display next frame.

         IF info.stopflag EQ 0 THEN BEGIN

p = info.thisFrame


CASE info.which_window_to_animate OF

   0 : BEGIN
   
info.L_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
                                
info.L_ut_string_object->SetProperty, strings = (info.L_list_of_full_time_strings)[info.L_current_image_number]                           
info.L_Window->Draw, info.L_both_views

swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color

widget_control,info.L_widget_image_sequence_slider,set_value = info.L_current_image_number + 1

   ENDCASE
   1 : BEGIN
   
if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
   
info.C_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
                                
info.C_ut_string_object->SetProperty, strings = (info.C_list_of_full_time_strings)[info.C_current_image_number]                           
info.C_Window->Draw, info.C_both_views

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color

widget_control,info.C_widget_image_sequence_slider,set_value = info.C_current_image_number + 1

endif else begin

info.C2_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
                                
info.C_ut_string_object->SetProperty, strings = (info.C2_list_of_full_time_strings)[info.C2_current_image_number]                           
info.C_Window->Draw, info.C_both_views

swpc_cat_set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, info.C2_current_image_number, info.color_C2, info.cme_outline_color

widget_control,info.C_widget_image_sequence_slider,set_value = info.C2_current_image_number + 1

endelse

   ENDCASE
   2 : BEGIN
   
info.R_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
                                
info.R_ut_string_object->SetProperty, strings = (info.R_list_of_full_time_strings)[info.R_current_image_number]                           
info.R_Window->Draw, info.R_both_views

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_A, info.cme_outline_color

widget_control,info.R_widget_image_sequence_slider,set_value = info.R_current_image_number + 1

   ENDCASE
   
ENDCASE


info.images_timeline_window -> draw, info.images_timeline_view

         ENDIF

            ; Update frame number.
         if info.i_anim_back_and_forth eq 0 then begin
         
         info.thisFrame = info.thisFrame + 1

         IF info.thisFrame GT (info.anim_end_frame) THEN info.thisFrame = info.anim_start_frame
         
         endif else begin
         
         info.thisFrame = info.thisFrame + info.i_move
         
         IF info.thisFrame GT (info.anim_end_frame) THEN begin
         info.thisFrame = info.thisFrame -1
         info.i_move = -1        
         endif
         IF info.thisFrame LT (info.anim_start_frame) THEN begin
         info.thisFrame = info.thisFrame +1
         info.i_move = 1
         endif
         
         endelse

            ; Set a timer event to get back into this event handler.

         WIDGET_CONTROL, event.ID, TIMER=0

            ; Update stop flag.

         info.stopflag = 0
      END ; of Play button CASE.

      'Pause' : BEGIN
      
widget_control,info.play_pause_button,set_value = 'Play '

         info.stopflag = 1
;         info.thisFrame = info.thisFrame - 1  ; Subtract a frame.
; make sure that when we stop the animation the current image number
; matches the image shown....
p = info.thisFrame

CASE info.which_window_to_animate OF

   0 : BEGIN
   
info.L_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
                                
info.L_ut_string_object->SetProperty, strings = (info.L_list_of_full_time_strings)[info.L_current_image_number]                           
info.L_Window->Draw, info.L_both_views

swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color

   ENDCASE
   1 : BEGIN
   
if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
   
info.C_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
                                
info.C_ut_string_object->SetProperty, strings = (info.C_list_of_full_time_strings)[info.C_current_image_number]                           
info.C_Window->Draw, info.C_both_views

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color

widget_control,info.C_widget_image_sequence_slider,set_value = info.C_current_image_number + 1

endif else begin

info.C2_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
                                
info.C_ut_string_object->SetProperty, strings = (info.C2_list_of_full_time_strings)[info.C2_current_image_number]                           
info.C_Window->Draw, info.C_both_views

swpc_cat_set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, info.C2_current_image_number, info.color_C2, info.cme_outline_color

widget_control,info.C_widget_image_sequence_slider,set_value = info.C2_current_image_number + 1

endelse

   ENDCASE
   2 : BEGIN
   
info.R_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
                                
info.R_ut_string_object->SetProperty, strings = (info.R_list_of_full_time_strings)[info.R_current_image_number]                           
info.R_Window->Draw, info.R_both_views

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_A, info.cme_outline_color

   ENDCASE
   
ENDCASE
info.images_timeline_window -> draw, info.images_timeline_view

         END

   ENDCASE

      ; Put INFO structure back.

   IF WIDGET_INFO(event.top, /VALID_ID) THEN $
      WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
   RETURN
ENDIF ; of BUTTON code

   ; Code for TIMER events.

IF eventName EQ 'WIDGET_TIMER' THEN BEGIN

   IF info.stopflag EQ 0 THEN BEGIN

         ; Delay for a moment.

      Wait, info.delay

p = info.thisFrame

CASE info.which_window_to_animate OF

   0 : BEGIN
   
info.L_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
                                
info.L_ut_string_object->SetProperty, strings = (info.L_list_of_full_time_strings)[info.L_current_image_number]                           
info.L_Window->Draw, info.L_both_views

swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color

widget_control,info.L_widget_image_sequence_slider,set_value = info.L_current_image_number + 1

   ENDCASE
   1 : BEGIN
   
if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
   
info.C_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
                                
info.C_ut_string_object->SetProperty, strings = (info.C_list_of_full_time_strings)[info.C_current_image_number]                           
info.C_Window->Draw, info.C_both_views

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color

widget_control,info.C_widget_image_sequence_slider,set_value = info.C_current_image_number + 1

endif else begin

info.C2_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
                                
info.C_ut_string_object->SetProperty, strings = (info.C2_list_of_full_time_strings)[info.C2_current_image_number]                           
info.C_Window->Draw, info.C_both_views

swpc_cat_set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, info.C2_current_image_number, info.color_C2, info.cme_outline_color

widget_control,info.C_widget_image_sequence_slider,set_value = info.C2_current_image_number + 1

endelse


   ENDCASE
   2 : BEGIN
   
info.R_current_image_number = p

swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
                                
info.R_ut_string_object->SetProperty, strings = (info.R_list_of_full_time_strings)[info.R_current_image_number]                           
info.R_Window->Draw, info.R_both_views

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_A, info.cme_outline_color

widget_control,info.R_widget_image_sequence_slider,set_value = info.R_current_image_number + 1

   ENDCASE
   
ENDCASE
info.images_timeline_window -> draw, info.images_timeline_view

            ; Update frame number.
         if info.i_anim_back_and_forth eq 0 then begin
         
         info.thisFrame = info.thisFrame + 1
         IF info.thisFrame GT (info.anim_end_frame) THEN info.thisFrame = info.anim_start_frame
         
         endif else begin
         
         info.thisFrame = info.thisFrame + info.i_move
         
         IF info.thisFrame GT (info.anim_end_frame) THEN begin
         info.thisFrame = info.thisFrame -2
         info.i_move = -1        
         endif
         IF info.thisFrame LT (info.anim_start_frame) THEN begin
         info.thisFrame = info.thisFrame +2
         info.i_move = 1
         endif
         
         endelse

         ; Set a timer event to get back into this event handler.

      IF info.stopflag EQ 0 THEN WIDGET_CONTROL, event.ID, TIMER=0
   ENDIF

      ; Put INFO structure back.

   WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
ENDIF ; of TIMER code

END ; ******************************** of XMOVIE_EVENT ***********************












pro swpc_cat_change_image_tab, event

; process comes here when changing the image adjust TAB -
; but nothing is defined at present.

END





pro swpc_cat_Reset_L, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.L_image_color_palette -> SetProperty, top_stretch = 100
info.L_image_color_palette -> SetProperty, bottom_stretch = 0
info.L_image_color_palette -> SetProperty, gamma = 1.
info.L_image_saturation_value = 200
swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
info.L_Window->Draw, info.L_both_views

widget_control, info.L_widget_topSlider, set_value = 100
widget_control, info.L_widget_botSlider, set_value = 0
widget_control, info.L_widget_gammaSlider, set_value = 10
widget_control, info.L_widget_saturationSlider, set_value = 200

Widget_Control, event.top, Set_UValue=info, /No_Copy

END

pro swpc_cat_Reset_C, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.C_image_color_palette -> SetProperty, top_stretch = 100
info.C_image_color_palette -> SetProperty, bottom_stretch = 0
info.C_image_color_palette -> SetProperty, gamma = 1.
info.C_image_saturation_value = 50

if info.currently_showing_lasco_C2_or_C3 eq 'C2' then begin
swpc_cat_image_difference_and_scaling, info.color_c2,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
endif else begin
swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
endelse


info.C_Window->Draw, info.C_both_views

widget_control, info.C_widget_topSlider, set_value = 100
widget_control, info.C_widget_botSlider, set_value = 0
widget_control, info.C_widget_gammaSlider, set_value = 10
widget_control, info.C_widget_saturationSlider, set_value = 50

Widget_Control, event.top, Set_UValue=info, /No_Copy

END

pro swpc_cat_Reset_R, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.R_image_color_palette -> SetProperty, top_stretch = 100
info.R_image_color_palette -> SetProperty, bottom_stretch = 0
info.R_image_color_palette -> SetProperty, gamma = 1.
info.R_image_saturation_value = 200
swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
info.R_Window->Draw, info.R_both_views

widget_control, info.R_widget_topSlider, set_value = 100
widget_control, info.R_widget_botSlider, set_value = 0
widget_control, info.R_widget_gammaSlider, set_value = 10
widget_control, info.R_widget_saturationSlider, set_value = 200

Widget_Control, event.top, Set_UValue=info, /No_Copy

END







pro swpc_cat_copy_L_to_C, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.L_image_color_palette -> GetProperty, top_stretch = top_stretch
info.L_image_color_palette -> GetProperty, bottom_stretch = bottom_stretch
info.L_image_color_palette -> GetProperty, gamma = gamma
info.C_image_color_palette -> SetProperty, top_stretch = top_stretch
info.C_image_color_palette -> SetProperty, bottom_stretch = bottom_stretch
info.C_image_color_palette -> SetProperty, gamma = gamma
info.C_image_saturation_value = info.L_image_saturation_value
swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
info.C_Window->Draw, info.C_both_views
widget_control, info.L_widget_topSlider, get_value = topslider
widget_control, info.L_widget_botSlider, get_value = botslider
widget_control, info.L_widget_gammaSlider, get_value = gammaslider
widget_control, info.L_widget_saturationSlider, get_value = saturationSlider
widget_control, info.C_widget_topSlider, set_value = topslider
widget_control, info.C_widget_botSlider, set_value = botslider
widget_control, info.C_widget_gammaSlider, set_value = gammaslider
widget_control, info.C_widget_saturationSlider, set_value = saturationSlider

Widget_Control, event.top, Set_UValue=info, /No_Copy

END



pro swpc_cat_copy_L_to_R, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.L_image_color_palette -> GetProperty, top_stretch = top_stretch
info.L_image_color_palette -> GetProperty, bottom_stretch = bottom_stretch
info.L_image_color_palette -> GetProperty, gamma = gamma
info.R_image_color_palette -> SetProperty, top_stretch = top_stretch
info.R_image_color_palette -> SetProperty, bottom_stretch = bottom_stretch
info.R_image_color_palette -> SetProperty, gamma = gamma
info.R_image_saturation_value = info.L_image_saturation_value
swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
info.R_Window->Draw, info.R_both_views
widget_control, info.L_widget_topSlider, get_value = topslider
widget_control, info.L_widget_botSlider, get_value = botslider
widget_control, info.L_widget_gammaSlider, get_value = gammaslider
widget_control, info.L_widget_saturationSlider, get_value = saturationSlider
widget_control, info.R_widget_topSlider, set_value = topslider
widget_control, info.R_widget_botSlider, set_value = botslider
widget_control, info.R_widget_gammaSlider, set_value = gammaslider
widget_control, info.R_widget_saturationSlider, set_value = saturationSlider

Widget_Control, event.top, Set_UValue=info, /No_Copy

END

pro swpc_cat_copy_C_to_L, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.C_image_color_palette -> GetProperty, top_stretch = top_stretch
info.C_image_color_palette -> GetProperty, bottom_stretch = bottom_stretch
info.C_image_color_palette -> GetProperty, gamma = gamma
info.L_image_color_palette -> SetProperty, top_stretch = top_stretch
info.L_image_color_palette -> SetProperty, bottom_stretch = bottom_stretch
info.L_image_color_palette -> SetProperty, gamma = gamma
info.L_image_saturation_value = info.C_image_saturation_value
swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
info.L_Window->Draw, info.L_both_views
widget_control, info.C_widget_topSlider, get_value = topslider
widget_control, info.C_widget_botSlider, get_value = botslider
widget_control, info.C_widget_gammaSlider, get_value = gammaslider
widget_control, info.C_widget_saturationSlider, get_value = saturationSlider
widget_control, info.L_widget_topSlider, set_value = topslider
widget_control, info.L_widget_botSlider, set_value = botslider
widget_control, info.L_widget_gammaSlider, set_value = gammaslider
widget_control, info.L_widget_saturationSlider, set_value = saturationSlider

Widget_Control, event.top, Set_UValue=info, /No_Copy

END

pro swpc_cat_copy_C_to_R, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.C_image_color_palette -> GetProperty, top_stretch = top_stretch
info.C_image_color_palette -> GetProperty, bottom_stretch = bottom_stretch
info.C_image_color_palette -> GetProperty, gamma = gamma
info.R_image_color_palette -> SetProperty, top_stretch = top_stretch
info.R_image_color_palette -> SetProperty, bottom_stretch = bottom_stretch
info.R_image_color_palette -> SetProperty, gamma = gamma
info.R_image_saturation_value = info.C_image_saturation_value
swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
info.R_Window->Draw, info.R_both_views
widget_control, info.C_widget_topSlider, get_value = topslider
widget_control, info.C_widget_botSlider, get_value = botslider
widget_control, info.C_widget_gammaSlider, get_value = gammaslider
widget_control, info.C_widget_saturationSlider, get_value = saturationSlider
widget_control, info.R_widget_topSlider, set_value = topslider
widget_control, info.R_widget_botSlider, set_value = botslider
widget_control, info.R_widget_gammaSlider, set_value = gammaslider
widget_control, info.R_widget_saturationSlider, set_value = saturationSlider

Widget_Control, event.top, Set_UValue=info, /No_Copy

END

pro swpc_cat_copy_R_to_L, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.R_image_color_palette -> GetProperty, top_stretch = top_stretch
info.R_image_color_palette -> GetProperty, bottom_stretch = bottom_stretch
info.R_image_color_palette -> GetProperty, gamma = gamma
info.L_image_color_palette -> SetProperty, top_stretch = top_stretch
info.L_image_color_palette -> SetProperty, bottom_stretch = bottom_stretch
info.L_image_color_palette -> SetProperty, gamma = gamma
info.L_image_saturation_value = info.R_image_saturation_value
swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
info.L_Window->Draw, info.L_both_views
widget_control, info.R_widget_topSlider, get_value = topslider
widget_control, info.R_widget_botSlider, get_value = botslider
widget_control, info.R_widget_gammaSlider, get_value = gammaslider
widget_control, info.R_widget_saturationSlider, get_value = saturationSlider
widget_control, info.L_widget_topSlider, set_value = topslider
widget_control, info.L_widget_botSlider, set_value = botslider
widget_control, info.L_widget_gammaSlider, set_value = gammaslider
widget_control, info.L_widget_saturationSlider, set_value = saturationSlider

Widget_Control, event.top, Set_UValue=info, /No_Copy

END

pro swpc_cat_copy_R_to_C, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.R_image_color_palette -> GetProperty, top_stretch = top_stretch
info.R_image_color_palette -> GetProperty, bottom_stretch = bottom_stretch
info.R_image_color_palette -> GetProperty, gamma = gamma
info.C_image_color_palette -> SetProperty, top_stretch = top_stretch
info.C_image_color_palette -> SetProperty, bottom_stretch = bottom_stretch
info.C_image_color_palette -> SetProperty, gamma = gamma
info.C_image_saturation_value = info.R_image_saturation_value
swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
info.C_Window->Draw, info.C_both_views
widget_control, info.R_widget_topSlider, get_value = topslider
widget_control, info.R_widget_botSlider, get_value = botslider
widget_control, info.R_widget_gammaSlider, get_value = gammaslider
widget_control, info.R_widget_saturationSlider, get_value = saturationSlider
widget_control, info.C_widget_topSlider, set_value = topslider
widget_control, info.C_widget_botSlider, set_value = botslider
widget_control, info.C_widget_gammaSlider, set_value = gammaslider
widget_control, info.C_widget_saturationSlider, set_value = saturationSlider

Widget_Control, event.top, Set_UValue=info, /No_Copy

END










pro swpc_cat_change_latitude, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.latitude_degrees = float(event.value) / 10.

;print,'latitude ', info.latitude_degrees

info.L_cme_model->SetProperty, transform = info.initial_transform
info.L_cme_model_copy->SetProperty, transform = info.initial_transform
info.L_cme_model->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.L_cme_model_copy->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.L_cme_model->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply
info.L_cme_model_copy->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply

info.C_cme_model->SetProperty, transform = info.initial_transform
info.C_cme_model_copy->SetProperty, transform = info.initial_transform
info.C_cme_model->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.C_cme_model_copy->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.C_cme_model->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply
info.C_cme_model_copy->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply

info.R_cme_model->SetProperty, transform = info.initial_transform
info.R_cme_model_copy->SetProperty, transform = info.initial_transform
info.R_cme_model->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.R_cme_model_copy->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.R_cme_model->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply
info.R_cme_model_copy->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply



swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline

if ((info.which_mode eq 0) and (info.L_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.L_clock_angle_degrees, info.L_coronagraph_image_object, info.L_image_lineplot, $
                            info.position_image_lineplot, info.L_cme_outline                           
endif
endif
info.L_Window->Draw, info.L_both_views

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline

if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.C_clock_angle_degrees, info.C_coronagraph_image_object, info.C_image_lineplot, $
                            info.position_image_lineplot, info.C_cme_outline                           
endif

info.C_Window->Draw, info.C_both_views

swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline

if ((info.which_mode eq 0) and (info.R_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.R_clock_angle_degrees, info.R_coronagraph_image_object, info.R_image_lineplot, $
                            info.position_image_lineplot, info.R_cme_outline                           
endif
endif
info.R_Window->Draw, info.R_both_views

info.lat_string = string(info.latitude_degrees,format='(f6.1)')
info.lat_string_object -> setproperty, strings = 'lat  :' + info.lat_string
info.cme_info_Window->Draw, info.cme_info_view

;send_event_to_3D_view, info.allow_show_3D_view, info.Three_D_view_tlb_ID, info.Three_D_view_event_widget, $
;                       info.latitude_degrees, info.longitude_degrees, $
;                       info.radial_distance_lemniscate, info.angular_width_lemniscate 

;if info.allow_show_3D_view eq 1 and info.Three_D_view_tlb_ID ne 0L then begin
;send_event = {id:0L, top:0L, handler:0L, $
;              latitude:info.latitude_degrees, $
;              longitude:info.longitude_degrees, $
;              radial_distance:info.radial_distance_lemniscate, $
;              angular_width:info.angular_width_lemniscate}
;widget_control,info.Three_D_view_event_widget, send_event = send_event
;endif

if info.allow_show_3D_view eq 1 and info.Three_D_view_tlb_ID ne 0L then begin
send_event = {id:0L, top:0L, handler:0L, which_parameter:1, parameter:info.latitude_degrees}
widget_control,info.Three_D_view_event_widget, send_event = send_event
endif




Widget_Control, event.top, Set_UValue=info, /No_Copy

END





pro swpc_cat_change_longitude, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.longitude_degrees = float(event.value) / 10.

;print,'longitude_degrees ', info.longitude_degrees

info.L_cme_model->SetProperty, transform = info.initial_transform
info.L_cme_model_copy->SetProperty, transform = info.initial_transform
info.L_cme_model->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.L_cme_model_copy->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.L_cme_model->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply
info.L_cme_model_copy->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply

info.C_cme_model->SetProperty, transform = info.initial_transform
info.C_cme_model_copy->SetProperty, transform = info.initial_transform
info.C_cme_model->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.C_cme_model_copy->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.C_cme_model->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply
info.C_cme_model_copy->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply

info.R_cme_model->SetProperty, transform = info.initial_transform
info.R_cme_model_copy->SetProperty, transform = info.initial_transform
info.R_cme_model->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.R_cme_model_copy->rotate,[0,1,0], info.longitude_degrees, /premultiply
info.R_cme_model->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply
info.R_cme_model_copy->rotate,[1,0,0], 0.0 - info.latitude_degrees, /premultiply


swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline

if ((info.which_mode eq 0) and (info.L_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.L_clock_angle_degrees, info.L_coronagraph_image_object, info.L_image_lineplot, $
                            info.position_image_lineplot, info.L_cme_outline                           
endif
endif

info.L_Window->Draw, info.L_both_views

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline

if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.C_clock_angle_degrees, info.C_coronagraph_image_object, info.C_image_lineplot, $
                            info.position_image_lineplot, info.C_cme_outline                           
endif
info.C_Window->Draw, info.C_both_views

swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline

if ((info.which_mode eq 0) and (info.R_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.R_clock_angle_degrees, info.R_coronagraph_image_object, info.R_image_lineplot, $
                            info.position_image_lineplot, info.R_cme_outline                           
endif
endif
info.R_Window->Draw, info.R_both_views


info.lon_string = string(info.longitude_degrees,format='(f6.1)')
info.lon_string_object -> setproperty, strings = 'lon  :' + info.lon_string
info.cme_info_Window->Draw, info.cme_info_view

if info.allow_show_3D_view eq 1 and info.Three_D_view_tlb_ID ne 0L then begin
send_event = {id:0L, top:0L, handler:0L, which_parameter:2, parameter:info.longitude_degrees}
widget_control,info.Three_D_view_event_widget, send_event = send_event
endif

Widget_Control, event.top, Set_UValue=info, /No_Copy

END











pro swpc_cat_change_lemniscate_radial_distance, event

Widget_Control, event.top, Get_UValue=info, /No_Copy


radial_distance_lemniscate = float(event.value) / 10.


cone_axis_data = fltarr(3,2)
cone_axis_data[2,1] = radial_distance_lemniscate *3.
info.L_cone_Z_axis-> SetProperty, data = cone_axis_data
info.C_cone_Z_axis-> SetProperty, data = cone_axis_data
info.R_cone_Z_axis-> SetProperty, data = cone_axis_data
if info.debug_mode eq 1 then print,' cone_axis_data ',cone_axis_data

info.C_cone_Z_axis-> GetProperty, XCOORD_CONV = xconv
info.C_cone_Z_axis-> GetProperty, YCOORD_CONV = yconv
info.C_cone_Z_axis-> GetProperty, ZCOORD_CONV = zconv 


info.radial_distance_lemniscate = radial_distance_lemniscate

swpc_cat_define_cme_lemniscate, info.radial_distance_lemniscate, info.angular_width_lemniscate, info.lemniscate_style, fitted_cme_info

info.L_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.L_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.L_Window_copy->Draw, info.L_second_view_copy
swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline

if ((info.which_mode eq 0) and (info.L_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.L_clock_angle_degrees, info.L_coronagraph_image_object, info.L_image_lineplot, $
                            info.position_image_lineplot, info.L_cme_outline                           
endif
endif

info.L_Window->Draw, info.L_both_views


info.C_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.C_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.C_Window_copy->Draw, info.C_second_view_copy

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline

if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.C_clock_angle_degrees, info.C_coronagraph_image_object, info.C_image_lineplot, $
                            info.position_image_lineplot, info.C_cme_outline                           
endif

;if (info.which_mode eq 0) then begin
info.C_Window->Draw, info.C_both_views
;endif

info.R_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.R_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.R_Window_copy->Draw, info.R_second_view_copy
swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline

if ((info.which_mode eq 0) and (info.R_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.R_clock_angle_degrees, info.R_coronagraph_image_object, info.R_image_lineplot, $
                            info.position_image_lineplot, info.R_cme_outline                           
endif
endif
info.R_Window->Draw, info.R_both_views


radial_distance_string = 'dist :' + string(radial_distance_lemniscate, format='(f6.1)')
info.radial_distance_string_object -> setproperty, strings = radial_distance_string
info.cme_info_Window->Draw, info.cme_info_view

if info.allow_show_3D_view eq 1 and info.Three_D_view_tlb_ID ne 0L then begin
send_event = {id:0L, top:0L, handler:0L, which_parameter:3, parameter:info.radial_distance_lemniscate}
widget_control,info.Three_D_view_event_widget, send_event = send_event
endif


Widget_Control, event.top, Set_UValue=info, /No_Copy

END





pro swpc_cat_change_lemniscate_angular_width, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

angular_width = float(event.value)

info.angular_width_lemniscate = angular_width

swpc_cat_define_cme_lemniscate, info.radial_distance_lemniscate, info.angular_width_lemniscate, info.lemniscate_style, fitted_cme_info

info.L_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.L_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.L_Window_copy->Draw, info.L_second_view_copy
swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline

if ((info.which_mode eq 0) and (info.L_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.L_clock_angle_degrees, info.L_coronagraph_image_object, info.L_image_lineplot, $
                            info.position_image_lineplot, info.L_cme_outline                           
endif
endif
info.L_Window->Draw, info.L_both_views

info.C_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.C_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.C_Window_copy->Draw, info.C_second_view_copy

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline


if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.C_clock_angle_degrees, info.C_coronagraph_image_object, info.C_image_lineplot, $
                            info.position_image_lineplot, info.C_cme_outline                           
endif
info.C_Window->Draw, info.C_both_views

info.R_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.R_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.R_Window_copy->Draw, info.R_second_view_copy
swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline

if ((info.which_mode eq 0) and (info.R_number_of_images GT 0)) then begin
if info.show_image_line_plot eq 1 then begin
swpc_cat_replot_image_line_plot, info.R_clock_angle_degrees, info.R_coronagraph_image_object, info.R_image_lineplot, $
                            info.position_image_lineplot, info.R_cme_outline                           
endif
endif
info.R_Window->Draw, info.R_both_views

angular_width_string = 'cone :' + string(angular_width, format='(f6.1)')
info.angular_width_string_object -> setproperty, strings = angular_width_string
info.cme_info_Window->Draw, info.cme_info_view

if info.allow_show_3D_view eq 1 and info.Three_D_view_tlb_ID ne 0L then begin
send_event = {id:0L, top:0L, handler:0L, which_parameter:4, parameter:info.angular_width_lemniscate / 2.0 }
widget_control,info.Three_D_view_event_widget, send_event = send_event
endif

Widget_Control, event.top, Set_UValue=info, /No_Copy

END






pro swpc_cat_send_event_to_3D_view, allow_show_3D_view, Three_D_view_tlb_ID, Three_D_view_event_widget, $
                       latitude_degrees, longitude_degrees, $
                       radial_distance_lemniscate, angular_width_lemniscate
if allow_show_3D_view eq 1 and Three_D_view_tlb_ID ne 0L then begin
send_event = {id:0L, top:0L, handler:0L, $
              latitude:latitude_degrees, $
              longitude:longitude_degrees, $
              radial_distance:radial_distance_lemniscate, $
              angular_width:angular_width_lemniscate}
widget_control,Three_D_view_event_widget, send_event = send_event
endif
end


pro swpc_cat_Calculate_Earth_B_Angle,Julian_date,B_angle_degrees

MJD = Julian_date - 2400000.5

ii = 7.25

nn = MJD - 51544.5
gg = 357.528 + (0.9856003 * nn)
ll = 280.460 + (0.9856474 * nn)

lambda = ll + (1.915*sin(gg*!pi/180.)) + 0.020*sin(2.*gg*!pi/180.)

omega = 73.667 + (0.013958 * (MJD + 3242) / 365.25)   ; Hapgood,  http://sspg1.bnsc.rl.ac.uk/Share/Coordinates/angles.htm

sin_b0 = sin((lambda - omega)*!pi/180.)*sin(ii*!pi/180.)

b0 = asin(sin_b0)

B_angle_degrees = b0 * 180. / !pi

end


; ECDP
pro swpc_cat_pos_cme_match_or_unmatch, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

;CASE event.id OF
case info.currently_showing_lasco_C2_or_C3 of   

   'C3' : begin
   
   widget_control,info.C_widget_pos_outline_matches_image,get_value=the_text
   
   (info.POS_Julian).add, (info.C_list_of_datetime_Julian)[info.C_current_image_number]
   (info.POS_matches_image_DateTime_string).add, (info.C_list_of_full_time_strings)[info.C_current_image_number]
   (info.POS_leading_edge).add, info.C_leading_edge_in_rs
   (info.POS_matches_image_telescope).add, info.C_telescope_code
   info.C_Window->Draw, info.C_both_views
;    print, info.POS_Julian, info.POS_leading_edge

  endcase
  
 'C2' : begin

   widget_control,info.C_widget_pos_outline_matches_image,get_value=the_text
   
  (info.POS_Julian).add, (info.C2_list_of_datetime_Julian)[info.C2_current_image_number]
   (info.POS_matches_image_DateTime_string).add, (info.C2_list_of_full_time_strings)[info.C2_current_image_number]
   (info.POS_leading_edge).add, info.C_leading_edge_in_rs
   (info.POS_matches_image_telescope).add, info.C2_telescope_code
   info.C_Window->Draw, info.C_both_views
;    print, info.POS_Julian, info.POS_leading_edge
  
    endcase
    
endcase    


swpc_cat_set_pos_speed_calc_button, info

swpc_cat_sort_out_the_pos_line_plot_and_widgets, info   
   
Widget_Control, event.top, Set_UValue=info, /No_Copy
end
; ECDP





pro swpc_cat_CME_match_or_unmatch, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

CASE event.id OF

   info.L_widget_outline_matches_image : BEGIN
   
   widget_control,info.L_widget_outline_matches_image,get_value=the_text
  
   
   if the_text eq 'CME Matches Image' then begin 
   
   (info.CME_matches_image_Julian).add, (info.L_list_of_datetime_Julian)[info.L_current_image_number]
   (info.CME_matches_image_DateTime_string).add, (info.L_list_of_full_time_strings)[info.L_current_image_number]
   (info.CME_matches_image_telescope).add, info.L_telescope_code
   (info.CME_matches_image_Rs_leading_edge).add, info.radial_distance_lemniscate
   (info.CME_matches_image_Image_number).add, info.L_current_image_number
   (info.CME_matches_image_L_Image_number).add, info.L_current_image_number

    swpc_cat_image_difference_and_scaling, info.color_stereo_B,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
    info.L_title_object->setProperty, color = info.background_color
    info.L_ut_string_object->setProperty, color = info.background_color

   info.L_cme_outline->GetProperty, data = data
   (info.CME_matches_image_CME_outline).add,data
   (info.CME_matches_image_L_CME_outline).add,data
   info.L_cme_outline -> setProperty, hide = 1
   info.L_cme_MATCH_outline->SetProperty, data = data
   info.L_cme_MATCH_outline->SetProperty, hide = 0
   info.L_Window->Draw, info.L_both_views
   
   widget_control,info.L_widget_outline_matches_image,set_value='Unmatch' 
   
   endif else begin
   
   if float(!version.release) gt 8.05 then begin
      the_index_to_be_removed = (info.CME_matches_image_Image_number).FindValue(info.L_current_image_number)
;   print, ' the_index_to_be_removed ', the_index_to_be_removed
      the_index_to_be_removed_L = (info.CME_matches_image_L_Image_number).FindValue(info.L_current_image_number)
;   print, ' the_index_to_be_removed_L ', the_index_to_be_removed_L
   endif else begin
      the_index_to_be_removed = where(info.CME_matches_image_Image_number eq info.L_current_image_number)
      the_index_to_be_removed_L = where(info.CME_matches_image_L_Image_number eq info.L_current_image_number)
   endelse
   
   
   (info.CME_matches_image_Julian).remove, the_index_to_be_removed
   (info.CME_matches_image_DateTime_string).remove, the_index_to_be_removed
   (info.CME_matches_image_telescope).remove, the_index_to_be_removed
   (info.CME_matches_image_Rs_leading_edge).remove, the_index_to_be_removed
   (info.CME_matches_image_Image_number).remove, the_index_to_be_removed
   (info.CME_matches_image_L_Image_number).remove, the_index_to_be_removed_L

    swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
    info.L_title_object->setProperty, color = info.color_stereo_B
    info.L_ut_string_object->setProperty, color = info.color_stereo_B

   (info.CME_matches_image_CME_outline).remove, the_index_to_be_removed
   (info.CME_matches_image_L_CME_outline).remove, the_index_to_be_removed_L
   info.L_cme_outline -> setProperty, hide = 0
   info.L_cme_MATCH_outline->SetProperty, hide = 1
   info.L_Window->Draw, info.L_both_views
   
   widget_control,info.L_widget_outline_matches_image,set_value='CME Matches Image'
   
   endelse
   
swpc_cat_sort_out_the_timeline_symbols, info.L_number_of_images, info.CME_matches_image_L_Image_number, $
                                  info.L_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline

    
   ENDCASE
      
   info.C_widget_outline_matches_image : BEGIN
   
   case info.currently_showing_lasco_C2_or_C3 of
   
   'C3' : begin
   
   widget_control,info.C_widget_outline_matches_image,get_value=the_text
   
   if the_text eq 'CME Matches Image' then begin 
   
   (info.CME_matches_image_Julian).add, (info.C_list_of_datetime_Julian)[info.C_current_image_number]
   (info.CME_matches_image_DateTime_string).add, (info.C_list_of_full_time_strings)[info.C_current_image_number]
   (info.CME_matches_image_telescope).add, info.C_telescope_code
   (info.CME_matches_image_Rs_leading_edge).add, info.radial_distance_lemniscate
   (info.CME_matches_image_Image_number).add, info.C_current_image_number
   (info.CME_matches_image_C_Image_number).add, info.C_current_image_number

    swpc_cat_image_difference_and_scaling, info.color_c3,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
    info.C_title_object->setProperty, color = info.background_color
    info.C_ut_string_object->setProperty, color = info.background_color

   info.C_cme_outline->GetProperty, data = data
   (info.CME_matches_image_CME_outline).add,data
   (info.CME_matches_image_C_CME_outline).add,data
   info.C_cme_outline -> setProperty, hide = 1
;     if (info.which_mode eq 1) then begin
;   info.C_cme_outline -> setProperty, hide = 1
;   endif else begin
;   info.C_cme_outline -> setProperty, hide = 0 
;   endelse 
; ECDP
   info.C_cme_MATCH_outline->SetProperty, data = data
   info.C_cme_MATCH_outline->SetProperty, hide = 0
   info.C_Window->Draw, info.C_both_views
   
   widget_control,info.C_widget_outline_matches_image,set_value='Unmatch'
   
   endif else begin
   
   if float(!version.release) gt 8.05 then begin
      the_index_to_be_removed = (info.CME_matches_image_Image_number).FindValue(info.C_current_image_number)
;   print, ' the_index_to_be_removed ', the_index_to_be_removed
      the_index_to_be_removed_C = (info.CME_matches_image_C_Image_number).FindValue(info.C_current_image_number)
;   print, ' the_index_to_be_removed_C ', the_index_to_be_removed_C
   endif else begin
      the_index_to_be_removed = where(info.CME_matches_image_Image_number eq info.C_current_image_number)
      the_index_to_be_removed_C = where(info.CME_matches_image_C_Image_number eq info.C_current_image_number)
   endelse
   
   
   (info.CME_matches_image_Julian).remove, the_index_to_be_removed
   (info.CME_matches_image_DateTime_string).remove, the_index_to_be_removed
   (info.CME_matches_image_telescope).remove, the_index_to_be_removed
   (info.CME_matches_image_Rs_leading_edge).remove, the_index_to_be_removed
   (info.CME_matches_image_Image_number).remove, the_index_to_be_removed
   (info.CME_matches_image_C_Image_number).remove, the_index_to_be_removed_C

    swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
    info.C_title_object->setProperty, color = info.color_c3
    info.C_ut_string_object->setProperty, color = info.color_c3

   (info.CME_matches_image_CME_outline).remove, the_index_to_be_removed
   (info.CME_matches_image_C_CME_outline).remove, the_index_to_be_removed_C
   info.C_cme_outline -> setProperty, hide = 0
;   if (info.which_mode eq 1) then begin
;   info.C_cme_outline -> setProperty, hide = 1
;   endif else begin
;   info.C_cme_outline -> setProperty, hide = 0 
;   endelse 
; ECDP 
   info.C_cme_MATCH_outline->SetProperty, hide = 1
   info.C_Window->Draw, info.C_both_views
   
   widget_control,info.C_widget_outline_matches_image,set_value='CME Matches Image'
   
   endelse
   
   
; whether we are matching or unmatching we need to sort out the timeline symbols appropriately                                   
swpc_cat_sort_out_the_timeline_symbols, info.C_number_of_images, info.CME_matches_image_C_Image_number, $
                                  info.C_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline
    
   endcase
   
   'C2' : begin

   widget_control,info.C_widget_outline_matches_image,get_value=the_text
   
   if the_text eq 'CME Matches Image' then begin 
   
   (info.CME_matches_image_Julian).add, (info.C2_list_of_datetime_Julian)[info.C2_current_image_number]
   (info.CME_matches_image_DateTime_string).add, (info.C2_list_of_full_time_strings)[info.C2_current_image_number]
   (info.CME_matches_image_telescope).add, info.C2_telescope_code
   (info.CME_matches_image_Rs_leading_edge).add, info.radial_distance_lemniscate
   (info.CME_matches_image_Image_number).add, info.C2_current_image_number
   (info.CME_matches_image_C2_Image_number).add, info.C2_current_image_number

    swpc_cat_image_difference_and_scaling, info.color_c2,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
    info.C_title_object->setProperty, color = info.background_color
    info.C_ut_string_object->setProperty, color = info.background_color

   info.C_cme_outline->GetProperty, data = data
   (info.CME_matches_image_CME_outline).add,data
   (info.CME_matches_image_C2_CME_outline).add,data
   info.C_cme_outline -> setProperty, hide = 1
    ; ECDP only draw if we are in mode 0 (i.e. normal mode)   
;   if (info.which_mode eq 1) then begin
;   info.C_cme_outline -> setProperty, hide = 1
;   endif else begin
;   info.C_cme_outline -> setProperty, hide = 0 
;   endelse 
; ECDP
   info.C2_cme_MATCH_outline->SetProperty, data = data
   info.C2_cme_MATCH_outline->SetProperty, hide = 0
   info.C_Window->Draw, info.C_both_views
   
   widget_control,info.C_widget_outline_matches_image,set_value='Unmatch'
   
   endif else begin
   
   if float(!version.release) gt 8.05 then begin
      the_index_to_be_removed = (info.CME_matches_image_Image_number).FindValue(info.C2_current_image_number)
;   print, ' the_index_to_be_removed ', the_index_to_be_removed
      the_index_to_be_removed_C2 = (info.CME_matches_image_C2_Image_number).FindValue(info.C2_current_image_number)
;   print, ' the_index_to_be_removed_C ', the_index_to_be_removed_C2
   endif else begin
      the_index_to_be_removed = where(info.CME_matches_image_Image_number eq info.C2_current_image_number)
      the_index_to_be_removed_C2 = where(info.CME_matches_image_C2_Image_number eq info.C2_current_image_number)
   endelse
   
   
   (info.CME_matches_image_Julian).remove, the_index_to_be_removed
   (info.CME_matches_image_DateTime_string).remove, the_index_to_be_removed
   (info.CME_matches_image_telescope).remove, the_index_to_be_removed
   (info.CME_matches_image_Rs_leading_edge).remove, the_index_to_be_removed
   (info.CME_matches_image_Image_number).remove, the_index_to_be_removed
   (info.CME_matches_image_C2_Image_number).remove, the_index_to_be_removed_C2

    swpc_cat_image_difference_and_scaling, info.background_color,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
    info.C_title_object->setProperty, color = info.color_c2
    info.C_ut_string_object->setProperty, color = info.color_c2

   (info.CME_matches_image_CME_outline).remove, the_index_to_be_removed
   (info.CME_matches_image_C2_CME_outline).remove, the_index_to_be_removed_C2
   
; ECDP only draw if we are in mode 0 (i.e. normal mode)   
;   if (info.which_mode eq 1) then begin
;   info.C_cme_outline -> setProperty, hide = 1
;   endif else begin
   info.C_cme_outline -> setProperty, hide = 0 
;   endelse 
; ECDP
   info.C2_cme_MATCH_outline->SetProperty, hide = 1
   info.C_Window->Draw, info.C_both_views
   
   widget_control,info.C_widget_outline_matches_image,set_value='CME Matches Image'
   
   endelse
                                  
; whether we are matching or unmatching we need to sort out the timeline symbols appropriately                                   
swpc_cat_sort_out_the_timeline_symbols, info.C2_number_of_images, info.CME_matches_image_C2_Image_number, $
                                  info.C2_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline


   
   endcase
   endcase
   
   info.C_Window->Draw, info.C_both_views

   ENDCASE
   
   
      
   info.R_widget_outline_matches_image : BEGIN
   
   widget_control,info.R_widget_outline_matches_image,get_value=the_text
   
   if the_text eq 'CME Matches Image' then begin 
   
   (info.CME_matches_image_Julian).add, (info.R_list_of_datetime_Julian)[info.R_current_image_number]
   (info.CME_matches_image_DateTime_string).add, (info.R_list_of_full_time_strings)[info.R_current_image_number]
   (info.CME_matches_image_telescope).add, info.R_telescope_code
   (info.CME_matches_image_Rs_leading_edge).add, info.radial_distance_lemniscate
   (info.CME_matches_image_Image_number).add, info.R_current_image_number
   (info.CME_matches_image_R_Image_number).add, info.R_current_image_number
   
    swpc_cat_image_difference_and_scaling, info.color_stereo_A,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
    info.R_title_object->setProperty, color = info.background_color
    info.R_ut_string_object->setProperty, color = info.background_color

   info.R_cme_outline->GetProperty, data = data
   (info.CME_matches_image_CME_outline).add,data
   (info.CME_matches_image_R_CME_outline).add,data
   info.R_cme_outline -> setProperty, hide = 1
   info.R_cme_MATCH_outline->SetProperty, data = data
   info.R_cme_MATCH_outline->SetProperty, hide = 0
    info.R_Window->Draw, info.R_both_views
    
   widget_control,info.R_widget_outline_matches_image,set_value='Unmatch'
    
   endif else begin
   
   if float(!version.release) gt 8.05 then begin
      the_index_to_be_removed = (info.CME_matches_image_Image_number).FindValue(info.R_current_image_number)
;   print, ' the_index_to_be_removed ', the_index_to_be_removed
      the_index_to_be_removed_R = (info.CME_matches_image_R_Image_number).FindValue(info.R_current_image_number)
;   print, ' the_index_to_be_removed_R ', the_index_to_be_removed_R
   endif else begin
      the_index_to_be_removed = where(info.CME_matches_image_Image_number eq info.R_current_image_number)
      the_index_to_be_removed_R = where(info.CME_matches_image_R_Image_number eq info.R_current_image_number)
   endelse
   
   
   (info.CME_matches_image_Julian).remove, the_index_to_be_removed
   (info.CME_matches_image_DateTime_string).remove, the_index_to_be_removed
   (info.CME_matches_image_telescope).remove, the_index_to_be_removed
   (info.CME_matches_image_Rs_leading_edge).remove, the_index_to_be_removed
   (info.CME_matches_image_Image_number).remove, the_index_to_be_removed
   (info.CME_matches_image_R_Image_number).remove, the_index_to_be_removed_R

    swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
    info.R_title_object->setProperty, color = info.color_stereo_A
    info.R_ut_string_object->setProperty, color = info.color_stereo_A

;   info.L_cme_outline->GetProperty, data = data
   (info.CME_matches_image_CME_outline).remove, the_index_to_be_removed
   (info.CME_matches_image_R_CME_outline).remove, the_index_to_be_removed_R
   info.R_cme_outline -> setProperty, hide = 0
;   info.L_cme_MATCH_outline->SetProperty, data = data
   info.R_cme_MATCH_outline->SetProperty, hide = 1
   info.R_Window->Draw, info.R_both_views
   
   widget_control,info.R_widget_outline_matches_image,set_value='CME Matches Image'
   
   endelse
   
; whether we are matching or unmatching we need to sort out the timeline symbols appropriately 
swpc_cat_sort_out_the_timeline_symbols, info.R_number_of_images, info.CME_matches_image_R_Image_number, $
                                  info.R_plot, info.xSymbolSize_timeline, info.ySymbolSize_timeline

   ENDCASE
     
ENDCASE  

                                  
info.images_timeline_window->Draw, info.images_timeline_view




swpc_cat_sort_out_the_line_plot_and_widgets, info





Widget_Control, event.top, Set_UValue=info, /No_Copy
END


; ECDP
pro swpc_cat_sort_out_the_pos_line_plot_and_widgets, info


dataX = (info.POS_Julian).toarray()
dataY = (info.POS_leading_edge).toarray()

elements = n_elements(dataY)

if elements gt 0 then begin

vert_colors = intarr(3,elements)

for i = 0 , elements - 1 do begin

;  case (info.POS_matches_image_telescope)[i] of

;  'SC3' : BEGIN
;  vert_colors[0,i] = 200
;  vert_colors[1,i] = 200
;  vert_colors[2,i] = 200
;  endcase

;  'SC2' : BEGIN
  vert_colors[0,i] = 200
  vert_colors[1,i] = 200
  vert_colors[2,i] = 200
;  endcase

;endcase

endfor

endif   ; if elements gt 0

swpc_cat_update_the_line_plot,info.position_B,info.xmin,info.xmax,info.ymin,info.ymax, $
                                 info.xaxis1_B,info.xaxis2_B,info.yaxis1_B,info.yaxis2_B, $
                                 info.xtickinterval, $                                 
                                 info.LE_plot_matched_CMEs,info.plot_view,info.plot_window,dataX,dataY, $
                                 vert_colors, $
                                 xs, ys

;endif

; check the total number of matches.  If this is less than 2 then the calculate velocity
; button needs to be greyed out......

if float(!version.release) gt 8.05 then begin
   number_of_matches = (info.POS_Julian).count()
endif else begin
   number_of_matches = n_elements(info.POS_Julian)
endelse

any_manuals = (info.CME_matches_image_telescope).FindValue('Manual')

number_of_points_for_velocity = number_of_matches
number_of_image_matches = number_of_matches
if any_manuals ne !null then number_of_image_matches = number_of_image_matches - 1

if number_of_points_for_velocity ge 2 then begin
widget_control,info.widget_calculate_pos_velocity,sensitive=1
;if (info.half_cone_angle_emp ne 0) then BEGIN
;widget_control,info.widget_use_emp_cone,sensitive=1
;endif
endif else begin
widget_control,info.widget_calculate_pos_velocity,sensitive=0
;widget_control,info.widget_use_emp_cone,sensitive=0
endelse

end

; ECDP




pro swpc_cat_sort_out_the_line_plot_and_widgets, info


dataX = (info.CME_matches_image_Julian).toarray()
dataY = (info.CME_matches_image_Rs_leading_edge).toarray()

elements = n_elements(dataY)

if elements gt 0 then begin

vert_colors = intarr(3,elements)

for i = 0 , elements - 1 do begin

case (info.CME_matches_image_telescope)[i] of

'BC2' : BEGIN
vert_colors[0,i] = 100
vert_colors[1,i] = 100
vert_colors[2,i] = 255
endcase

'SC3' : BEGIN
vert_colors[0,i] = 100
vert_colors[1,i] = 255
vert_colors[2,i] = 100
endcase

'SC2' : BEGIN
vert_colors[0,i] = 0
vert_colors[1,i] = 100
vert_colors[2,i] = 0
endcase

'AC2' : BEGIN
vert_colors[0,i] = 255
vert_colors[1,i] = 100
vert_colors[2,i] = 100
endcase

'Manual' : BEGIN
vert_colors[0,i] = 255
vert_colors[1,i] = 255
vert_colors[2,i] = 255
endcase

endcase

endfor

endif else begin  ; if elements gt 0

; we get here if we've done some cme matches and then 
; unmatched them again until no matches are left.

dataX=dblarr(1)
dataY=fltarr(1)
dataX[0] = (info.xmin + info.xmax) / 2.
dataY[0] = 15.
vert_colors = info.background_color

endelse  ; if elements gt 0
; hide any lines from previous plotting as this can be confusing
info.thisPlot_B2->SetProperty, hide = 1
; plot new points
swpc_cat_update_the_line_plot,info.position_B,info.xmin,info.xmax,info.ymin,info.ymax, $
                                 info.xaxis1_B,info.xaxis2_B,info.yaxis1_B,info.yaxis2_B, $
                                 info.xtickinterval, $                                 
                                 info.LE_plot_matched_CMEs,info.plot_view,info.plot_window,dataX,dataY, $
                                 vert_colors, $
                                 xs, ys
; thisPlot_B2 = info.thisPlot_B2,dataX2 = POS_Jul_times2 , dataY2 = POS_yvals

if n_elements(info.CME_matches_image_Rs_leading_edge) eq 1 then begin

;print, ' FIRST TIME'

lat_string = strcompress(string(round(info.latitude_degrees)),/remove_all)
lon_string = strcompress(string(round(info.longitude_degrees)),/remove_all)
cone_angle_string = strcompress(string(round(info.angular_width_lemniscate / 2.)),/remove_all)

info.latitude_text -> SetProperty, strings = 'lat ' + lat_string
info.latitude_text -> SetProperty, hide = 0

info.longitude_text -> SetProperty, strings = 'lon ' + lon_string
info.longitude_text -> SetProperty, hide = 0

info.cone_angle_text -> SetProperty, strings = 'cone ' + cone_angle_string
info.cone_angle_text -> SetProperty, hide = 0

info.enlil_info_window->Draw, info.enlil_info_View

endif


; check the total number of matches.  If this is less than 2 then the calculate velocity
; button needs to be greyed out......

if float(!version.release) gt 8.05 then begin
   number_of_matches = (info.CME_matches_image_Julian).count()
endif else begin
   number_of_matches = n_elements(info.CME_matches_image_Julian)
endelse

any_manuals = (info.CME_matches_image_telescope).FindValue('Manual')

number_of_points_for_velocity = number_of_matches
number_of_image_matches = number_of_matches
if any_manuals ne !null then number_of_image_matches = number_of_image_matches - 1

if number_of_points_for_velocity ge 2 then begin
widget_control,info.widget_calculate_velocity,sensitive=1
endif else begin
widget_control,info.widget_calculate_velocity,sensitive=0
endelse

if number_of_image_matches ge 1 then begin
widget_control,info.widget_angular_width_slider, sensitive=0
widget_control,info.widget_latitude_slider, sensitive=0
widget_control,info.widget_longitude_slider, sensitive=0
widget_control,info.L_widget_remove_this_image,sensitive=0
widget_control,info.C_widget_remove_this_image,sensitive=0
widget_control,info.R_widget_remove_this_image,sensitive=0
endif else begin
widget_control,info.widget_angular_width_slider, sensitive=1
widget_control,info.widget_latitude_slider, sensitive=1
widget_control,info.widget_longitude_slider, sensitive=1
endelse

end



; ECDP
pro swpc_cat_set_pos_speed_calc_button, info

if float(!version.release) gt 8.05 then begin
   number_of_pos_matches = (info.POS_Julian).count()
endif else begin
   number_of_pos_matches = n_elements(info.POS_Julian)
endelse

;any_manuals = (info.CME_matches_image_telescope).FindValue('Manual')

number_of_points_for_pos_velocity = number_of_pos_matches
number_of_pos_image_matches = number_of_pos_matches
;if any_manuals ne !null then number_of_image_matches = number_of_image_matches - 1

if number_of_points_for_pos_velocity ge 2 then begin
widget_control,info.widget_calculate_pos_velocity,sensitive=1
endif else begin
widget_control,info.widget_calculate_pos_velocity,sensitive=0
endelse

end
; ECDP






pro swpc_cat_xcolors_top_slider, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

CASE event.id OF

   info.L_widget_topSlider : BEGIN
   info.L_image_color_palette -> setproperty,top_stretch = event.value
   info.L_Window->Draw, info.L_both_views
   ENDCASE
      
   info.C_widget_topSlider : BEGIN
   info.C_image_color_palette -> setproperty,top_stretch = event.value
   info.C_Window->Draw, info.C_both_views
   ENDCASE
      
   info.R_widget_topSlider : BEGIN
   info.R_image_color_palette -> setproperty,top_stretch = event.value
   info.R_Window->Draw, info.R_both_views
   ENDCASE
     
ENDCASE

Widget_Control, event.top, Set_UValue=info, /No_Copy
END




pro swpc_cat_xcolors_bottom_slider, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

CASE event.id OF

   info.L_widget_botSlider : BEGIN
   info.L_image_color_palette -> setproperty,bottom_stretch = event.value
   info.L_Window->Draw, info.L_both_views
   ENDCASE
      
   info.C_widget_botSlider : BEGIN
   info.C_image_color_palette -> setproperty,bottom_stretch = event.value
   info.C_Window->Draw, info.C_both_views
   ENDCASE
      
   info.R_widget_botSlider : BEGIN
   info.R_image_color_palette -> setproperty,bottom_stretch = event.value
   info.R_Window->Draw, info.R_both_views
   ENDCASE
   
ENDCASE

Widget_Control, event.top, Set_UValue=info, /No_Copy
END 




pro swpc_cat_xcolors_gamma_slider, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

CASE event.id OF

   info.L_widget_gammaSlider : BEGIN
   info.L_image_color_palette -> setproperty,gamma = float(event.value)/ 10.
   info.L_Window->Draw, info.L_both_views
   ENDCASE
      
   info.C_widget_gammaSlider : BEGIN
   info.C_image_color_palette -> setproperty,gamma = float(event.value)/ 10.
   info.C_Window->Draw, info.C_both_views
   ENDCASE
      
   info.R_widget_gammaSlider : BEGIN
   info.R_image_color_palette -> setproperty,gamma = float(event.value)/ 10.
   info.R_Window->Draw, info.R_both_views
   ENDCASE
   
ENDCASE

Widget_Control, event.top, Set_UValue=info, /No_Copy
END






pro swpc_cat_change_the_image_saturation, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

CASE event.id OF

   info.L_widget_saturationSlider : BEGIN
   info.L_image_saturation_value = event.value
swpc_cat_image_difference_and_scaling, info.L_current_background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
info.L_ut_string_object->setProperty, color = info.L_current_text_color
info.L_title_object->setProperty, color = info.L_current_text_color
   info.L_Window->Draw, info.L_both_views
   ENDCASE
      
   info.C_widget_saturationSlider : BEGIN
   
   info.C_image_saturation_value = event.value
   
case info.currently_showing_lasco_C2_or_C3 of 

'C2' : begin 
swpc_cat_image_difference_and_scaling, info.C_current_background_color,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
info.C_ut_string_object->setProperty, color = info.C_current_text_color
info.C_title_object->setProperty, color = info.C_current_text_color
endcase
'C3' : begin 
swpc_cat_image_difference_and_scaling, info.C_current_background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
info.C_ut_string_object->setProperty, color = info.C_current_text_color
info.C_title_object->setProperty, color = info.C_current_text_color
endcase
endcase
                                
                                
   info.C_Window->Draw, info.C_both_views
   
   ENDCASE
      
   info.R_widget_saturationSlider : BEGIN
   info.R_image_saturation_value = event.value
swpc_cat_image_difference_and_scaling, info.R_current_background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
info.R_ut_string_object->setProperty, color = info.R_current_text_color
info.R_title_object->setProperty, color = info.R_current_text_color 
   info.R_Window->Draw, info.R_both_views
   ENDCASE
   
ENDCASE

Widget_Control, event.top, Set_UValue=info, /No_Copy

END


pro swpc_cat_change_lemniscate_transparency, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

ellipse_alpha = float(event.value) / 100. 

info.cme_fitted_surf -> setProperty, alpha_channel= ellipse_alpha

info.C_Window->Draw, info.C_both_views                                

Widget_Control, event.top, Set_UValue=info, /No_Copy

END






pro swpc_cat_set_xmin_and_xmax,Julian_start,Julian_end,xmin,xmax,xtickinterval

;help , X_values

;xmax=max(X_values,min=xmin)

xmin = Julian_start
xmax = Julian_end

xmin_fix = long(xmin)
xmax_fix = long(xmax)

xmin_rem = xmin - xmin_fix
xmax_rem = xmax - xmax_fix

scale_to_nearest = 6

if scale_to_nearest eq 6 then begin

if xmin_rem lt 0.25 then begin
xmin_rem = 0.0
endif
if xmin_rem ge 0.25 and xmin_rem lt 0.5 then begin
xmin_rem = 0.25
endif
if xmin_rem ge 0.5 and xmin_rem lt 0.75 then begin
xmin_rem = 0.5
endif
if xmin_rem ge 0.75 then begin
xmin_rem = 0.75
endif

if xmax_rem ge 0.75 then begin
xmax_rem = 1.0
endif
if xmax_rem ge 0.5 and xmax_rem lt 0.75 then begin
xmax_rem = 0.75
endif
if xmax_rem ge 0.25 and xmax_rem lt 0.5 then begin
xmax_rem = 0.5
endif
if xmax_rem lt 0.25 then begin
xmax_rem = 0.25
endif

xtickinterval = 0.25

endif

if scale_to_nearest eq 1 then begin

xmin_rem = (float(fix(10. * xmin_rem)))/10.d

xmax_rem = (float(fix(10. * xmax_rem) + 1))/10.d

xtickinterval = 1./24.



endif




xmin = double(xmin_fix) + double(xmin_rem) - 0.0001

xmax = double(xmax_fix) + double(xmax_rem) + 0.0001

end




pro swpc_cat_calculate_ymin_ymax_to_nearest_5, ymin, ymax

ymin2 = fix(ymin) / 5
ymin = ymin2 * 5.

ymax2 = (fix(ymax) / 5) + 1
ymax = ymax2 * 5.

end





pro swpc_cat_update_the_line_plot,position_B,xmin,xmax,ymin,ymax,xaxis1_B,xaxis2_B,yaxis1_B,yaxis2_B, $
                                       xtickinterval , $
                                       LE_plot_matched_CMEs,plot_view,plot_window,dataX,dataY,vert_colors, $
                                       thisplot_B2=thisplot_B2,dataX2=dataX2,dataY2=dataY2, $
                                       xs, ys
                                       
;help, dataX
;print, dataX , dataX - xmin , n_elements(dataY)

if n_elements(dataY) ge 2 then begin
xmax = max(dataX,min = xmin)
CALDAT, xmin, min_Month, min_Day, min_Year, min_Hour, min_Minute
xmin = JULDAY(min_Month, min_Day, min_Year, min_Hour)
CALDAT, xmax, max_Month, max_Day, max_Year, max_Hour, max_Minute
xmax = JULDAY(max_Month, max_Day, max_Year, max_Hour) + (1.d/24.d)
xtickinterval = 1.d/24.d
ymax=max(dataY, min = ymin)
swpc_cat_calculate_ymin_ymax_to_nearest_5,ymin,ymax
ytickinterval = 5.
endif else begin
xtickinterval = 0.25d
ytickinterval = 10.
endelse


; update the line plot with this info.......

result = label_date(date_format='%H',offset=xmin)
xrange = [0.d,xmax-xmin]
an_hour = 1./24.
six_hours = 0.25

xs = swpc_cat_FSC_Normalize(xrange, Position=[position_B[0], position_B[2]])

xaxis1_B->SetProperty, exact = 1
xaxis1_B->SetProperty, range=xrange
xAxis1_B->SetProperty, XCoord_Conv=xs
xaxis1_B->SetProperty, tickinterval = xtickinterval
xaxis1_B->SetProperty, tickformat='label_date'

xaxis2_B->SetProperty, exact = 1
xaxis2_B->SetProperty, range=xrange
xAxis2_B->SetProperty, XCoord_Conv=xs
xaxis2_B->SetProperty, tickinterval = xtickinterval
xaxis2_B->SetProperty, tickformat='label_date'

yrange = [ymin,ymax]
ys = swpc_cat_FSC_Normalize(yrange, Position=[position_B[1], position_B[3]])

yaxis1_B->SetProperty, exact = 1
yaxis1_B->SetProperty, range=yrange
yaxis1_B->SetProperty, YCoord_Conv=ys
yaxis1_B->SetProperty, tickinterval = ytickinterval

yaxis2_B->SetProperty, exact = 1
yaxis2_B->SetProperty, range=yrange
yaxis2_B->SetProperty, YCoord_Conv=ys
yaxis1_B->SetProperty, tickinterval = ytickinterval

LE_plot_matched_CMEs->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
LE_plot_matched_CMEs->SetProperty, dataX = dataX - xmin, dataY = dataY
;help, dataY
LE_plot_matched_CMEs->GetProperty, symbol=thisSymbol

;help,thisSymbol

;help, vert_colors
;print, vert_colors

if n_elements(vert_colors) eq 3 then begin
LE_plot_matched_CMEs->SetProperty, vert_colors=0  ; make sure if the plot has previously had vert_colors set
                                                  ; then this is now removed (ie, set to zero).
LE_plot_matched_CMEs->SetProperty, color=vert_colors
endif else begin
LE_plot_matched_CMEs->SetProperty, vert_colors=vert_colors
endelse

    ; Size the symbols appropriately for the plot.
    
symSize = 1.0
plot_position_scaling_factor = (position_B[2] - position_B[0]) / (position_B[3] - position_B[1]) 

xSymSize = (xmax - xmin) * 0.015 * symSize
ySymSize = (ymax - ymin) * 0.015 * plot_position_scaling_factor * symSize
thisSymbol->SetProperty, Size=[xSymSize, ySymSize]

;thisSymbol ->GetProperty, size = symbolsize
;print, 'symbol size ',symbolsize
;print, 'xrange ',xmax - xmin
;print, 'yrange ',ymax - ymin
;symbolsize[1] = symbolsize[0] * ((ymax - ymin) / (xmax - xmin))
;thisSymbol ->SetProperty, size = symbolsize
;thisplot_B->SetProperty, symbol=thisSymbol


if keyword_set(thisplot_B2) then begin
thisPlot_B2->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
thisPlot_B2->SetProperty, dataX = dataX2 - xmin, dataY = dataY2
thisPlot_B2->SetProperty, color = [255,0,0]
thisPlot_B2->SetProperty, hide = 0
endif

plot_window->Draw, plot_View




end





pro swpc_cat_create_plot_objects,position,xTitle,yTitle,thisPlot,thisPlot2,xAxis1,xAxis2,yAxis1,yAxis2, $
                        thisSymbol,symsize,exact,helvetica10pt,helvetica14pt,x,y, $
                        thisPlot2_c3, thisPlot2_c2, thisPlot2_stereo_a, thisPlot2_stereo_b


xTitle = Obj_New('IDLgrText', xtitle, Color=[255,255,0])
yTitle = Obj_New('IDLgrText', ytitle, Color=[255,255,0])


    ; Create a plot object. The plot will be in the coordinate
    ; space 0->1. The view will be in the range -0.35->1.25 so
    ; that the plot axis annotation will be visable. Make the plot
    ; a green color.
thisSymbol = obj_new("IDLgrsymbol",data=6)

thisPlot = Obj_New("IDLgrPLOT", x, y, _Extra=extra, $
   Color=[255,255,0], Symbol=thisSymbol, Thick=1 , linestyle = 6)

thisPlot2 = Obj_New("IDLgrPLOT", x, y, _Extra=extra, $
   Color=[255,255,0], Thick=1, linestyle = 0)

thisPlot2_c3 = Obj_New("IDLgrPLOT", x, y, _Extra=extra, $
   Color=[255,255,0], Thick=1, linestyle = 0) 
thisPlot2_c2 = Obj_New("IDLgrPLOT", x, y, _Extra=extra, $
   Color=[255,255,0], Thick=1, linestyle = 0)  
thisPlot2_stereo_a = Obj_New("IDLgrPLOT", x, y, _Extra=extra, $
   Color=[255,255,0], Thick=1, linestyle = 0)
thisPlot2_stereo_b = Obj_New("IDLgrPLOT", x, y, _Extra=extra, $
   Color=[255,255,0], Thick=1, linestyle = 0)

    ; Get the data ranges from the Plot Object.

thisPlot->GetProperty, XRange=xrange, YRange=yrange

    ; Create plot box style axes. Make the axes yellow.
    ; The large values in the LOCATION keyword indicates which
    ; values are NOT used. The axes text is set to Helvetica
    ; 10 point font.

xAxis1 = Obj_New("IDLgrAxis", 0, Color=[255,255,255], Ticklen=0.01, $
    Minor=0, Range=xrange, Title=xtitle, $
    Location=[1000, position[1] ,0], Exact=exact[0])
xAxis1->GetProperty, Ticktext=xAxisText
xAxisText->SetProperty, Recompute_Dimensions=2    ; ensure that when axes are changed, the text size is recomputed properly.
xAxisText->SetProperty, Font=helvetica10pt

xAxis2 = Obj_New("IDLgrAxis", 0, Color=[255,255,255], Ticklen=0.01, $
    Minor=0, /NoText, Range=xrange, TickDir=1, $
    Location=[1000, position[3], 0], Exact=exact[0])

yAxis1 = Obj_New("IDLgrAxis", 1, Color=[255,255,255], Ticklen=0.01, $
    Minor=0, Title=ytitle, Range=yrange, $
    Location=[position[0], 1000, 0], Exact=exact[1])
yAxis1->GetProperty, Ticktext=yAxisText
yAxisText->SetProperty, Recompute_Dimensions=2    ; ensure that when axes are changed, the text size is recomputed properly.
yAxisText->SetProperty, Font=helvetica10pt

yAxis2 = Obj_New("IDLgrAxis", 1, Color=[255,255,255], Ticklen=0.01, $
    Minor=0, /NoText, Range=yrange, TickDir=1, $
    Location=[position[2], 1000, 0], Exact=exact[1])

    ; Because we may not be using exact axis ranging, the axes
    ; may extend further than the xrange and yrange. Get the
    ; actual axis range so that the plot, etc. can be scaled
    ; appropriately.
xAxis1->GetProperty, CRange=xrange
yAxis1->GetProperty, CRange=yrange

    ; Set up the scaling so that the axes for the plot and the
    ; plot data extends from 0->1 in the X and Y directions.

xs = swpc_cat_FSC_Normalize(xrange, Position=[position[0], position[2]])
ys = swpc_cat_FSC_Normalize(yrange, Position=[position[1], position[3]])

    ; Scale the plot data and axes into 0->1.

thisPlot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
xAxis1->SetProperty, XCoord_Conv=xs
xAxis2->SetProperty, XCoord_Conv=xs
yAxis1->SetProperty, YCoord_Conv=ys
yAxis2->SetProperty, YCoord_Conv=ys

    ; Size the symbols appropriately for the plot.

xSymSize = (xrange[1] - xrange[0]) * 0.015 * symSize
ySymSize = (yrange[1] - yrange[0]) * 0.015 * symSize
IF Obj_Valid(thisSymbol) THEN thisSymbol->SetProperty, Size=[xSymSize, ySymSize]
   
xAxis1 -> SetProperty, hide = 1
xAxis2 -> SetProperty, hide = 1
yAxis1 -> SetProperty, hide = 1
yAxis2 -> SetProperty, hide = 1
thisPlot -> SetProperty, hide = 1
thisPlot2 -> SetProperty, hide = 1

end






pro swpc_cat_diff_format, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

Widget_Control, event.id, Get_UValue=difference_imaging_format

L_C_or_R = strmid(difference_imaging_format,0,1)

CASE difference_imaging_format OF

   'L Normal': BEGIN
   
      info.L_difference_imaging = 'none'
      widget_control,info.L_widget_remove_this_image,sensitive=1
      palette_top = info.palette_top_normal
      palette_bottom = info.palette_bottom_normal
      palette_gamma = info.palette_gamma_normal

      ENDCASE
      
   'L Diff: Running': BEGIN
   
      info.L_difference_imaging = 'running'
      widget_control,info.L_widget_remove_this_image,sensitive=0
      palette_top = info.palette_top_diff
      palette_bottom = info.palette_bottom_diff
      palette_gamma = info.palette_gamma_diff

      ENDCASE
      
   'L Diff: Set Current Image as Background': BEGIN
   
      info.L_difference_imaging = 'base'
      widget_control,info.L_widget_remove_this_image,sensitive=0
      info.L_background_image_number = info.L_current_image_number
      palette_top = info.palette_top_diff
      palette_bottom = info.palette_bottom_diff
      palette_gamma = info.palette_gamma_diff

      ENDCASE
      
   'C Normal': BEGIN
   
      if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
            info.C_difference_imaging = 'none'
            widget_control,info.C_widget_remove_this_image,sensitive=1
      endif else begin
            info.C2_difference_imaging = 'none'
            widget_control,info.C_widget_remove_this_image,sensitive=1
      endelse
      
      palette_top = info.palette_top_normal
      palette_bottom = info.palette_bottom_normal
      palette_gamma = info.palette_gamma_normal

      ENDCASE
      
   'C Diff: Running': BEGIN
   
      if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
            info.C_difference_imaging = 'running'
            widget_control,info.C_widget_remove_this_image,sensitive=0
      endif else begin
            info.C2_difference_imaging = 'running'
            widget_control,info.C_widget_remove_this_image,sensitive=0
      endelse
      
      palette_top = info.palette_top_diff
      palette_bottom = info.palette_bottom_diff
      palette_gamma = info.palette_gamma_diff

      ENDCASE
      
   'C Diff: Set Current Image as Background': BEGIN

      if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
            info.C_difference_imaging = 'base'
            widget_control,info.C_widget_remove_this_image,sensitive=0
            info.C_background_image_number = info.C_current_image_number
      endif else begin
            info.C2_difference_imaging = 'base'
            widget_control,info.C_widget_remove_this_image,sensitive=0
            info.C2_background_image_number = info.C2_current_image_number
      endelse
      
      palette_top = info.palette_top_diff
      palette_bottom = info.palette_bottom_diff
      palette_gamma = info.palette_gamma_diff

      ENDCASE
      
   'R Normal': BEGIN
   
      info.R_difference_imaging = 'none'
      widget_control,info.R_widget_remove_this_image,sensitive=1
      palette_top = info.palette_top_normal
      palette_bottom = info.palette_bottom_normal
      palette_gamma = info.palette_gamma_normal

      ENDCASE
      
   'R Diff: Running': BEGIN
   
      info.R_difference_imaging = 'running'
      widget_control,info.R_widget_remove_this_image,sensitive=0
      palette_top = info.palette_top_diff
      palette_bottom = info.palette_bottom_diff
      palette_gamma = info.palette_gamma_diff

      ENDCASE
      
   'R Diff: Set Current Image as Background': BEGIN
   
      info.R_difference_imaging = 'base'
      widget_control,info.R_widget_remove_this_image,sensitive=0
      info.R_background_image_number = info.R_current_image_number
      palette_top = info.palette_top_diff
      palette_bottom = info.palette_bottom_diff
      palette_gamma = info.palette_gamma_diff

      ENDCASE

ENDCASE

CASE L_C_or_R OF

   'L': BEGIN

swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
                                
info.L_Window->Draw, info.L_image_view

        ENDCASE


   'C': BEGIN
      
      if info.currently_showing_lasco_C2_or_C3 eq 'C3' then begin
swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
      endif else begin
swpc_cat_image_difference_and_scaling, info.background_color,  info.C2_current_image_number, info.C2_background_image_number, info.C2_difference_imaging, $
                 info.C2_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale
      endelse
                                
info.C_Window->Draw, info.C_image_view

        ENDCASE
        
   'R': BEGIN
                
swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
                                
info.R_Window->Draw, info.R_image_view

        ENDCASE
        
ENDCASE                          

info.L_Window->Draw, info.L_image_view


Widget_Control, event.top, Set_UValue=info, /No_Copy

END





pro swpc_cat_export_cme_analysis_data, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

if info.representative_image_has_been_defined eq 1 then begin

output_latitude_solar_equatorial = info.latitude_degrees
output_longitude = info.longitude_degrees
output_cone_angle = info.angular_width_lemniscate / 2.
output_radial_velocity = info.radial_velocity
; ECDP 
output_pos_velocity = info.pos_velocity
output_emp_cone_angle = (-13.0 + 18.6*alog10(output_pos_velocity))
output_emp_cone_angle_unc = 6.5

t_enlil_Year = info.time_at_Enlil_boundary[0]
t_enlil_Month = info.time_at_Enlil_boundary[1]
t_enlil_Day = info.time_at_Enlil_boundary[2]
t_enlil_Hour = info.time_at_Enlil_boundary[3]
t_enlil_Minute = info.time_at_Enlil_boundary[4]
t_enlil_Second = info.time_at_Enlil_boundary[5]

t_enlil_year_string = strcompress(string(t_enlil_Year),/remove_all)
t_enlil_month_string = strcompress(string(t_enlil_Month),/remove_all)
if t_enlil_month_string lt 10 then t_enlil_month_string = '0' + t_enlil_month_string
t_enlil_day_string = strcompress(string(t_enlil_day),/remove_all)
if t_enlil_day_string lt 10 then t_enlil_day_string = '0' + t_enlil_day_string
t_enlil_hour_string = strcompress(string(t_enlil_hour),/remove_all)
if t_enlil_hour_string lt 10 then t_enlil_hour_string = '0' + t_enlil_hour_string
t_enlil_minute_string = strcompress(string(t_enlil_minute),/remove_all)
if t_enlil_minute_string lt 10 then t_enlil_minute_string = '0' + t_enlil_minute_string

Enlil_boundary_string = t_enlil_year_string + t_enlil_month_string + t_enlil_day_string + 'T' + t_enlil_hour_string + t_enlil_minute_string

output_latitude_nint = round(output_latitude_solar_equatorial)
output_longitude_nint = round(output_longitude)
output_cone_angle_nint = round(output_cone_angle)
output_radial_velocity_nint = round(output_radial_velocity)

; ECDP
output_pos_velocity_nint = round(output_pos_velocity)
output_emp_cone_angle_nint = round(output_emp_cone_angle)


export_filename = 'swpc_cat_' + $
                Enlil_boundary_string + '_' + $
                strcompress(string(output_latitude_nint),/remove_all) + '_' + $
                strcompress(string(output_longitude_nint),/remove_all) + '_' + $
                strcompress(string(output_cone_angle_nint),/remove_all) + '_' + $
                strcompress(string(output_radial_velocity_nint),/remove_all)

if info.OPS_or_VnV ne 'OPS' then begin

suffix = ['','_a','_b','_c','_d','_e','_f','_g','_h','_i','_j','_k','_l','_m' $
          ,'_n','_o','_p','_q','_r','_s','_t','_u','_v','_w','_x','_y','_z']
iii = 0
directory_test_result = 1
while directory_test_result eq 1 do begin
;  print, iii
export_directory = info.export_location_root + info.sep + strcompress(t_enlil_year_string $
                   + t_enlil_month_string + t_enlil_day_string + suffix[iii],/remove_all)
directory_test_result = FILE_TEST(export_directory , /DIRECTORY)
iii ++
endwhile

file_mkdir, export_directory
             
if export_directory ne '' then begin
               
filename_minus_extension = export_directory

openw,export_lun, export_directory + info.sep + export_filename + '.txt',/get_lun, error = err

IF (err NE 0) then begin

PRINT,'error writing '

endif else begin



printf,export_lun,info.cone_tool_version
printf,export_lun,strcompress(string(Enlil_boundary_string),/remove_all)
printf,export_lun,strcompress(string(output_latitude_nint),/remove_all)
printf,export_lun,strcompress(string(output_longitude_nint),/remove_all)
printf,export_lun,strcompress(string(output_cone_angle_nint),/remove_all)
printf,export_lun,strcompress(string(output_radial_velocity_nint),/remove_all)
; ECDP
printf,export_lun,strcompress(string(output_pos_velocity_nint),/remove_all)
printf,export_lun,strcompress(string(output_emp_cone_angle_nint),/remove_all)
;
printf,export_lun,''

printf,export_lun,'Time at 21.5Rs (yyyy-mm-dd hh:mm) : ',info.T21_5_string
printf,export_lun,'CME axis latitude       (degrees) : ',output_latitude_solar_equatorial
printf,export_lun,'CME axis longitude      (degrees) : ',output_longitude
printf,export_lun,'half angle              (degrees) : ',output_cone_angle
printf,export_lun,'Radial velocity            (km/s) : ',output_radial_velocity
; ECDP
printf,export_lun,'POS velocity               (km/s) : ',output_pos_velocity
printf,export_lun,'EMP half angle          (degrees) : ',output_emp_cone_angle
;
printf,export_lun,''
printf,export_lun,'************************************'
printf,export_lun,''

n_matches = n_elements(info.CME_matches_image_telescope)
which_manual = (info.CME_matches_image_telescope).FindValue('Manual')
if which_manual eq !null then begin
any_manuals = 0
endif else begin
any_manuals = 1
endelse

if any_manuals eq 1 then begin
printf,export_lun,'Number of CME matches : ', strcompress(string(n_matches),/remove_all), ', Including 1 manual data point'
endif else begin
printf,export_lun,'Number of CME matches : ', n_matches
endelse

printf,export_lun,''
for i = 0 , n_matches - 1 do begin

if any_manuals eq 1 then begin
if i eq which_manual then begin
CALDAT, (info.CME_matches_image_Julian)[which_manual], Month, Day, Year, Hour, Minute
printf,export_lun,'CME match    : ', i + 1
printf,export_lun,'*********** Manual data point *************'
printf,export_lun,'Julian       : ', (info.CME_matches_image_Julian)[which_manual], format = '(A15,F20.10)'
printf,export_lun,'y m d h m    : ', year, month, day, hour, minute
printf,export_lun,'Leading Edge : ', (info.CME_matches_image_Rs_leading_edge)[which_manual]
printf,export_lun,''
endif else begin

printf,export_lun,'CME match    : ', i + 1
printf,export_lun,'Telescope    : ', (info.CME_matches_image_telescope)[i]
printf,export_lun,'Date-Time    : ', (info.CME_matches_image_DateTime_string)[i]
printf,export_lun,'Julian       : ', (info.CME_matches_image_Julian)[i], format = '(A15,F20.10)'
printf,export_lun,'Leading Edge : ', (info.CME_matches_image_Rs_leading_edge)[i]
;printf,export_lun,''
;printf,export_lun,'************************************'



if info.output_matched_line_data_in_txt_file eq 1 then begin
  
n_points = (n_elements((info.CME_matches_image_CME_outline)[i])) / 3
printf,export_lun,'N points     : ', n_points
for j = 0 , n_points - 1 do begin
printf,export_lun,j+1,((info.CME_matches_image_CME_outline)[i])[0,j] , ((info.CME_matches_image_CME_outline)[i])[1,j]
endfor

endif

printf,export_lun,''

endelse ; any_manuals
endif else begin

; regular CME matches (with no manual match)

printf,export_lun,'CME match    : ', i + 1
printf,export_lun,'Telescope    : ', (info.CME_matches_image_telescope)[i]
printf,export_lun,'Date-Time    : ', (info.CME_matches_image_DateTime_string)[i]
printf,export_lun,'Julian       : ', (info.CME_matches_image_Julian)[i], format = '(A15,F20.10)'
printf,export_lun,'Leading Edge : ', (info.CME_matches_image_Rs_leading_edge)[i]
printf,export_lun,''

if info.output_matched_line_data_in_txt_file eq 1 then begin

n_points = (n_elements((info.CME_matches_image_CME_outline)[i])) / 3
printf,export_lun,'N points     : ', n_points
for j = 0 , n_points - 1 do begin
printf,export_lun,j+1,((info.CME_matches_image_CME_outline)[i])[0,j] , ((info.CME_matches_image_CME_outline)[i])[1,j]
endfor

endif

;printf,export_lun,''

endelse

endfor

; ECDP - output the leading edge information, just in case it's useful 
  n_pos_matches = n_elements(info.POS_Julian)
  
if (n_pos_matches gt 0) then begin
;  printf,export_lun,''
  printf,export_lun,'************************************'
  printf,export_lun,''
  printf,export_lun,'Number of CME POS matches : ', n_pos_matches
  printf,export_lun,''
  
  for k = 0, n_pos_matches - 1 do begin
    printf,export_lun,'CME POS match    : ', k + 1
    printf,export_lun,'Telescope    : ', (info.POS_matches_image_telescope)[k]
    printf,export_lun,'Date-Time    : ', (info.POS_matches_image_DateTime_string)[k]
    printf,export_lun,'Julian       : ', (info.POS_Julian)[k], format = '(A15,F20.10)'
    printf,export_lun,'POS Leading Edge : ', (info.POS_leading_edge)[k]
    printf,export_lun,''   
  endfor   

endif

; ECDP
printf,export_lun,'************************************'

close,export_lun
free_lun,export_lun

;

julian_matches_array = (info.CME_matches_image_Julian).toarray()
jul_max = max(julian_matches_array,min = jul_min)

jul_min = jul_min - 1.d/24.d
jul_max = jul_max + 1.d/24.d

; loop over all images and save as pngs................

L_png_directory = export_directory + info.sep + 'images_STEREO_B'
file_mkdir, L_png_directory
for i = 0, info.L_number_of_images - 1 do begin

if (info.L_list_of_datetime_Julian)[i] ge jul_min and (info.L_list_of_datetime_Julian)[i] le jul_max then begin

info.L_current_image_number = i

swpc_cat_REDRAW_THE_IMAGE, $
    info.L_current_image_number,info.L_background_image_number,info.L_difference_imaging, $
    info.L_list_of_image_data,info.L_image_saturation_value,info.L_coronagraph_image_object,info.L_border_image_object, $
    info.CME_matches_image_L_Image_number,info.L_current_background_color, $
    info.background_color,info.L_current_text_color,info.color_stereo_B,info.L_cme_outline,info.L_cme_MATCH_outline, $
    info.L_widget_outline_matches_image,info.CME_matches_image_L_CME_outline, $
    info.L_ut_string_object,info.L_list_of_full_time_strings,info.L_title_object,info.L_Window,info.L_both_views,1,0, info.i_log_scale
   
t_str = (info.L_list_of_full_time_strings)[info.L_current_image_number]
image_name = strmid(t_str,0,4) + strmid(t_str,5,2) + strmid(t_str,8,2) + '_' + strmid(t_str,11,2) + strmid(t_str,14,2)
info.L_Window -> getproperty, image_data = image_data
filename = L_png_directory + info.sep + image_name + '.png'
WRITE_PNG, Filename, Image_data

endif
    
endfor

info.C_cme_outline -> setProperty, hide = 1
info.C_cme_MATCH_outline-> setProperty, hide = 1
info.C2_cme_MATCH_outline-> setProperty, hide = 1
info.currently_showing_lasco_C2_or_C3 = 'C3'
info.C_title_object -> setproperty, strings = 'SOHO LASCO C3'

C_png_directory = export_directory + info.sep + 'images_LASCO_C3'
file_mkdir, C_png_directory
for i = 0 , info.C_number_of_images - 1 do begin

if (info.C_list_of_datetime_Julian)[i] ge jul_min and (info.C_list_of_datetime_Julian)[i] le jul_max then begin

info.C_current_image_number = i

swpc_cat_REDRAW_THE_IMAGE, $
    info.C_current_image_number,info.C_background_image_number,info.C_difference_imaging, $
    info.C_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c3,info.C_cme_outline,info.C_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C_CME_outline, $
    info.C_ut_string_object,info.C_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,1,0, info.i_log_scale
    
t_str = (info.C_list_of_full_time_strings)[info.C_current_image_number]
image_name = strmid(t_str,0,4) + strmid(t_str,5,2) + strmid(t_str,8,2) + '_' + strmid(t_str,11,2) + strmid(t_str,14,2)
info.C_Window -> getproperty, image_data = image_data
filename = C_png_directory + info.sep + image_name + '.png'
WRITE_PNG, Filename, Image_data

endif
    
endfor

info.C_cme_outline -> setProperty, hide = 1
info.C_cme_MATCH_outline-> setProperty, hide = 1
info.C2_cme_MATCH_outline-> setProperty, hide = 1
info.currently_showing_lasco_C2_or_C3 = 'C2'
info.C_title_object -> setproperty, strings = 'SOHO LASCO C2'
widget_control,info.widget_show_C2_or_C3,set_value='Show LASCO C3'

C2_png_directory = export_directory + info.sep + 'images_LASCO_C2'
file_mkdir, C2_png_directory
for i = 0 , info.C2_number_of_images - 1 do begin

if (info.C2_list_of_datetime_Julian)[i] ge jul_min and (info.C2_list_of_datetime_Julian)[i] le jul_max then begin
  
info.C2_current_image_number = i

swpc_cat_REDRAW_THE_IMAGE, $
    info.C2_current_image_number,info.C2_background_image_number,info.C2_difference_imaging, $
    info.C2_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C2_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c2,info.C_cme_outline,info.C2_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C2_CME_outline, $
    info.C_ut_string_object,info.C2_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,1,0, info.i_log_scale
    
t_str = (info.C2_list_of_full_time_strings)[info.C2_current_image_number]
image_name = strmid(t_str,0,4) + strmid(t_str,5,2) + strmid(t_str,8,2) + '_' + strmid(t_str,11,2) + strmid(t_str,14,2)
info.C_Window -> getproperty, image_data = image_data
filename = C2_png_directory + info.sep + image_name + '.png'
WRITE_PNG, Filename, Image_data

endif
    
endfor


R_png_directory = export_directory + info.sep + 'images_STEREO_A'
file_mkdir, R_png_directory
for i = 0 , info.R_number_of_images - 1 do begin

if (info.R_list_of_datetime_Julian)[i] ge jul_min and (info.R_list_of_datetime_Julian)[i] le jul_max then begin

info.R_current_image_number = i
  
swpc_cat_REDRAW_THE_IMAGE, $
    info.R_current_image_number,info.R_background_image_number,info.R_difference_imaging, $
    info.R_list_of_image_data,info.R_image_saturation_value,info.R_coronagraph_image_object,info.R_border_image_object, $
    info.CME_matches_image_R_Image_number,info.R_current_background_color, $
    info.background_color,info.R_current_text_color,info.color_stereo_A,info.R_cme_outline,info.R_cme_MATCH_outline, $
    info.R_widget_outline_matches_image,info.CME_matches_image_R_CME_outline, $
    info.R_ut_string_object,info.R_list_of_full_time_strings,info.R_title_object,info.R_Window,info.R_both_views,1,0, info.i_log_scale
    
t_str = (info.R_list_of_full_time_strings)[info.R_current_image_number]
image_name = strmid(t_str,0,4) + strmid(t_str,5,2) + strmid(t_str,8,2) + '_' + strmid(t_str,11,2) + strmid(t_str,14,2)
info.R_Window -> getproperty, image_data = image_data
filename = R_png_directory + info.sep + image_name + '.png'
WRITE_PNG, Filename, Image_data

endif
    
endfor

; end of image loop.


; line plots
if (info.output_line_plots eq 1) then begin
; re-calculate velocity
swpc_cat_do_velocity, info
                                 
info.plot_window -> getproperty, image_data = image_data
filename = export_directory + info.sep + export_filename + 'CME_LinePlot.png'
write_png, filename, image_data

; re-calculate POS velocity
swpc_cat_do_pos_velocity, info

info.plot_window -> getproperty, image_data = image_data
filename = export_directory + info.sep + export_filename + 'CME_POS_LinePlot.png'
write_png, filename, image_data
endif


endelse ; error writing

endif  ; filename = ''

endif ; if info.OPS_or_VnV ne 'OPS'


;   fits file output....

Enlil_boundary_string = strcompress(string(Enlil_boundary_string),/remove_all)
output_latitude_nint = strcompress(string(output_latitude_nint),/remove_all)
output_longitude_nint = strcompress(string(output_longitude_nint),/remove_all)
output_cone_angle_nint = strcompress(string(output_cone_angle_nint),/remove_all)
output_radial_velocity_nint = strcompress(string(output_radial_velocity_nint),/remove_all)
; ECDP
output_pos_velocity_nint = strcompress(string(output_pos_velocity_nint),/remove_all)
output_emp_cone_angle_nint = strcompress(string(output_emp_cone_angle_nint), /remove_all)
;
snapshot = info.representative_image_data

im_info = size(snapshot,/struc)
;help,im_info
ndimen = im_info.n_dimensions
dimen = im_info.dimensions[0:ndimen-1]

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

swpc_cat_fxaddpar, header, 'SIMPLE','T','File conforms to FITS standard.'
swpc_cat_fxaddpar, header, 'BITPIX',bitpix,'Number of bits per data pixel'
swpc_cat_fxaddpar, header, 'NAXIS',ndimen,'Number of data axes'
swpc_cat_fxaddpar, header, 'NAXIS1',dimen[0],'Length of data axis 1'
swpc_cat_fxaddpar, header, 'NAXIS2',dimen[1],'Length of data axis 2'
swpc_cat_fxaddpar, header, 'NAXIS3',dimen[2],'Length of data axis 3'
swpc_cat_fxaddpar, header, 'SUB_NAME',submitter_name,'Name of the submitter'
swpc_cat_fxaddpar, header, 'CME-DATE-TIME',Enlil_boundary_string,'Date & Time at Enlil boundary'
swpc_cat_fxaddpar, header, 'CME-LAT',output_latitude_nint,'CME Latitude'
swpc_cat_fxaddpar, header, 'CME-LON',output_longitude_nint,'CME Longitude'
swpc_cat_fxaddpar, header, 'CME-CONE-ANGLE',output_cone_angle_nint,'CME Cone Angle'
swpc_cat_fxaddpar, header, 'CME-VELOCITY',output_radial_velocity_nint,'CME Radial Velocity'
; ECDP
swpc_cat_fxaddpar, header, 'CME-POSVELOCITY',output_pos_velocity_nint,'CME POS Velocity'
swpc_cat_fxaddpar, header, 'CME-EMP-CONE-ANGLE', output_emp_cone_angle_nint, 'CME EMP Cone Angle'
;

if info.OPS_or_VnV eq 'OPS' then begin
fits_directory = info.export_location_root
endif else begin
fits_directory = export_directory
endelse

the_fits_file = fits_directory + info.sep + export_filename + '.fts'
                
swpc_cat_writefits,the_fits_file,snapshot,header,/CHECKSUM

; also write a .png version...

the_png_file = fits_directory + info.sep + export_filename + '.png'

write_png,the_png_file,transpose(snapshot, [2,0,1])


if info.OPS_or_VnV eq 'OPS' then begin
   msg1 = 'You have exported your CME analysis as file:'
   msg2 = ''
   msg3 = export_filename
   msg4 = 'This is now available for upload into SPI.'
   ok = DIALOG_MESSAGE([msg1,msg2,msg3,msg2,msg4],/information,/center)
endif else begin
   msg1 = 'Your results have been exported to folder:'
   msg2 = ''
   msg3 = export_directory
   ok = DIALOG_MESSAGE([msg1,msg2,msg3],/information,/center)  
endelse

endif else begin  ;  representative image has not been defined....

msg1 = 'You do not have a Representative Image defined.'
msg2 = 'This is needed before you can export your result'
ok = DIALOG_MESSAGE([msg1,msg2],/error,/center)

endelse


Widget_Control, event.top, Set_UValue=info, /No_Copy

END

;
pro swpc_cat_set_opening_angle, event

compile_opt idl2
Widget_Control, event.top, Get_UValue=info, /No_Copy
; switch mode back to NORMAL
info.which_mode = 0
; since the output of the fitting procedure is the half opening angle, we need to multiply this by 2.
value = 2.0*info.half_cone_angle_emp
info.angular_width_lemniscate = value
widget_control, info.widget_angular_width_slider, set_value = info.angular_width_lemniscate
swpc_cat_define_cme_lemniscate, info.radial_distance_lemniscate, info.angular_width_lemniscate, info.lemniscate_style, fitted_cme_info
info.C_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.C_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.C_Window_copy->Draw, info.C_second_view_copy
info.C_cme_outline->SetProperty, Hide = 0
info.L_cme_outline->SetProperty, Hide = 0
info.R_cme_outline->SetProperty, Hide = 0
swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
info.C_Window->Draw, info.C_both_views
info.R_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.R_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.R_Window_copy->Draw, info.R_second_view_copy
swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline
info.R_Window->Draw, info.R_both_views
info.L_cme_fitted_surf->SetProperty, data = fitted_cme_info.vertices
info.L_cme_fitted_surf_copy->SetProperty, data = fitted_cme_info.vertices
info.L_Window_copy->Draw, info.L_second_view_copy
swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline
info.L_Window->Draw, info.L_both_views


angular_width_string = 'cone :' + string(value, format='(f6.1)')
info.angular_width_string_object -> setproperty, strings = angular_width_string
info.cme_info_Window->Draw, info.cme_info_view

info.LE_plot_matched_CMEs ->SetProperty, hide = 1
info.thisPlot_B2->SetProperty, hide = 1
info.plot_window->Draw, info.plot_View
info.LE_plot_matched_CMEs ->SetProperty, hide = 0
info.thisPlot_B2->SetProperty, hide = 0
;swpc_cat_change_lemniscate_angular_width
;swpc_cat_define_cme_lemniscate, info.radial_distance_lemniscate, info.angular_width_lemniscate, info.lemniscate_style, fitted_cme_info

;C_widget_pos_outline_matches_image
; switch mode back to NORMAL
;info.which_mode = 0
widget_control,info.widget_latitude_slider, /sensitive 
widget_control,info.widget_longitude_slider, /sensitive 
widget_control,info.widget_angular_width_slider, /sensitive 
widget_control,info.widget_radial_distance_slider, /sensitive 
widget_control, info.C_widget_outline_matches_image, sensitive = 1
widget_control, info.C_widget_pos_outline_matches_image, sensitive = 0
widget_control, info.C_widget_representative_image, sensitive = 1
widget_control, info.L_widget_representative_image, sensitive = 1
widget_control, info.R_widget_representative_image, sensitive = 1
if (info.compact eq 0) then begin
widget_control, info.c_mode_name, set_value = 'NORMAL MODE'
widget_control, info.time_plot_title_widget, set_value = 'CME LEADING EDGE (Rs) vs TIME (UT)'
endif
;widget_control, info.mode_name, set_value = 'NORMAL (SOHO + STEREO)'
; remove the clock
info.C_clock_angle_model -> SetProperty, Hide = 1
info.C_leading_edge_model -> SetProperty, Hide = 1
; re-draw everything
info.C_Window->Draw, info.C_both_views
;info.L_Window->Draw, info.L_both_views
;info.R_Window->Draw, info.R_both_views

Widget_Control, event.top, Set_UValue=info, /No_Copy
end
;

pro swpc_cat_do_pos_velocity, info

compile_opt idl2
;Widget_Control, event.top, Get_UValue=info, /No_Copy

number_of_matches = n_elements(info.POS_Julian)

POS_Jul_times = (info.POS_Julian).toarray()
POS_le = (info.POS_leading_edge).toarray()

POS_sorted_list = sort(POS_Jul_times)

; use Jul_times2 and R_le2 for the sorted list.
; Don't want to mess with the original Jul_times
; and R_le above.

POS_Jul_times2 = POS_Jul_times[POS_sorted_list]
POS_le2 = POS_le[POS_sorted_list]

result_POS = LINFIT(POS_Jul_times2,POS_le2)

POS_velocity = result_POS[1] * info.Solar_radius_km / 86400.d
POS_t_zero = -(result_POS[0] / result_POS[1])

;if abs(radial_velocity) gt 10.0 then begin

POS_a = result_POS[0]
POS_b = result_POS[1]
POS_yvals = POS_a + (POS_b * POS_Jul_times2)

info.LE_plot_matched_CMEs->GetProperty, vert_colors=vert_colors
;print, POS_a, POS_b, POS_yvals
;print, POS_velocity
info.POS_velocity = POS_velocity
vert_colors = [200,200,200]
if info.debug_mode eq 1 then print, 'info.xmin 1 ',info.xmin
this_xmin = info.xmin
swpc_cat_update_the_line_plot,info.position_B,this_xmin,info.xmax,info.ymin,info.ymax, $
                                 info.xaxis1_B,info.xaxis2_B,info.yaxis1_B,info.yaxis2_B, $
                                 info.xtickinterval, $                                 
                                 info.LE_plot_matched_CMEs,info.plot_view,info.plot_window, $
                                 POS_Jul_times, POS_le, vert_colors, $
                                 thisPlot_B2 = info.thisPlot_B2,dataX2 = POS_Jul_times2 , dataY2 = POS_yvals, $
                                 xs, ys
if info.debug_mode eq 1 then print, 'info.xmin 2 ',this_xmin

if (info.compact eq 0) then begin
widget_control, info.time_plot_title_widget, set_value = 'POS LEADING EDGE (Rs) vs TIME (UT)'
endif
info.LE_plot_matched_CMEs ->SetProperty, hide = 0
      info.thisPlot_B2->SetProperty, hide = 0
      info.plot_window->Draw, info.plot_View


; Now calculate the empirical opening angle
half_cone_angle_emp = (-13.0 + 18.6*alog10(POS_velocity))
info.half_cone_angle_emp = half_cone_angle_emp
half_cone_angle_emp_int = round(half_cone_angle_emp)
str_half_cone_angle_emp = strcompress(string(half_cone_angle_emp_int),/remove_all)
; 1-sigma uncertainty on half angle = 6.5 degrees\
half_cone_angle_emp_unc = 6.5
half_cone_angle_emp_unc_int = round(half_cone_angle_emp_unc)
str_half_cone_angle_emp_unc = strcompress(string(half_cone_angle_emp_unc_int),/remove_all)

; ECDP; pos velocity and emp cone angle
pvel_string = strcompress(string(fix(pos_velocity)),/remove_all)
info.pos_velocity_text -> SetProperty, strings = 'POS Vel ' + pvel_string
info.pos_velocity_text -> SetProperty, hide = 0

; ECDP define the plusminus symbol
plusminus = string(177b)
emp_cone_angle_string = strcompress(string(fix(half_cone_angle_emp)),/remove_all)
emp_cone_angle_unc_string = strcompress(string(fix(half_cone_angle_emp_unc)),/remove_all)
info.emp_cone_angle_text -> SetProperty, strings = 'Emp cone angle ' + emp_cone_angle_string + plusminus + emp_cone_angle_unc_string
info.emp_cone_angle_text -> SetProperty, hide = 0

; enable the user to share the information
widget_control,info.widget_use_emp_cone,sensitive=1

; write the values to the screen
info.enlil_info_window->Draw, info.enlil_info_View
;print, 'Here'
;print, info.background_color
;print, emp_cone_angle_unc_string
end

; ECDP



; ECDP
pro swpc_cat_calculate_pos_velocity, event

;*****************************************************************
;
;  This is executed when the calculate_pos_velocity button is pressed
;  With at least 2 defined positions for the CME the POS velocity
;  can be calculated...
;
;*****************************************************************

compile_opt idl2
Widget_Control, event.top, Get_UValue=info, /No_Copy
; ECDP - go to the code which actually calculates the POS velocity (above)
swpc_cat_do_pos_velocity, info

Widget_Control, event.top, Set_UValue=info, /No_Copy

end


; Code which calculates the velocity
pro swpc_cat_do_velocity, info

compile_opt idl2

number_of_matches = n_elements(info.CME_matches_image_Julian)

Jul_times = (info.CME_matches_image_Julian).toarray()
R_le = (info.CME_matches_image_Rs_leading_edge).toarray()

sorted_list = sort(Jul_times)

; use Jul_times2 and R_le2 for the sorted list.
; Don't want to mess with the original Jul_times
; and R_le above.

Jul_times2 = Jul_times[sorted_list]
R_le2 = R_le[sorted_list]

result_R = LINFIT(Jul_times2,R_le2)

radial_velocity = result_R[1] * info.Solar_radius_km / 86400.d
t_zero = - (result_R[0] / result_R[1])

if abs(radial_velocity) gt 10.0 then begin

a = result_R[0]
b = result_R[1]
yvals = a + (b * Jul_times2)

info.LE_plot_matched_CMEs->GetProperty, vert_colors=vert_colors

vert_colors = [100,255,100]

; ECDP make sure that the title is always correct
if (info.compact eq 0) then begin
widget_control, info.time_plot_title_widget, set_value = 'CME LEADING EDGE (Rs) vs TIME (UT)'
endif

if info.debug_mode eq 1 then print, 'info.xmin 1 ',info.xmin
this_xmin = info.xmin
swpc_cat_update_the_line_plot,info.position_B,this_xmin,info.xmax,info.ymin,info.ymax, $
                                 info.xaxis1_B,info.xaxis2_B,info.yaxis1_B,info.yaxis2_B, $
                                 info.xtickinterval, $                                 
                                 info.LE_plot_matched_CMEs,info.plot_view,info.plot_window, $
                                 Jul_times, R_le, vert_colors, $
                                 thisPlot_B2 = info.thisPlot_B2, dataX2 = Jul_times2 , dataY2 = yvals, $
                                 xs, ys
                                
if info.debug_mode eq 1 then print, 'info.xmin 2 ',this_xmin

;print, ' Radial velocity (km/s) ',  radial_velocity
;print, t_zero, format = '(f20.10)'
CALDAT, t_zero, t_zero_Month , t_zero_Day , t_zero_Year , t_zero_Hour , t_zero_Minute , t_zero_Second
;print, 'Time at Sun',t_zero_Year,t_zero_Month,t_zero_Day,t_zero_Hour,t_zero_Minute,t_zero_Second

time_from_Sun_to_Enlil_seconds = (21.5 * info.Solar_radius_km) / radial_velocity
time_from_Sun_to_Enlil_days = time_from_Sun_to_Enlil_seconds / 86400.d

Time_at_Enlil_boundary = t_zero + time_from_Sun_to_Enlil_days
CALDAT, Time_at_Enlil_boundary, t_enlil_Month , t_enlil_Day , t_enlil_Year , t_enlil_Hour , t_enlil_Minute , t_enlil_Second
;print, 'Time at Enlil (21.5Rs) ',t_enlil_Year,t_enlil_Month,t_enlil_Day,t_enlil_Hour,t_enlil_Minute,t_enlil_Second

info.radial_velocity = radial_velocity
info.time_at_Enlil_boundary_Julian = Time_at_Enlil_boundary
info.time_at_Enlil_boundary = [t_enlil_Year,t_enlil_Month,t_enlil_Day,t_enlil_Hour,t_enlil_Minute,t_enlil_Second]



; ********** write the results in the format needed for input to Enlil *************

half_radius = info.angular_width_lemniscate / 2.0
half_radius_int = round(half_radius)
str_half_radius = strcompress(string(half_radius_int),/remove_all)
 
latitude_degrees = info.latitude_degrees
lat_deg_int = round(latitude_degrees)
str_latitude = strcompress(string(lat_deg_int),/remove_all)

longitude_degrees = info.longitude_degrees
lon_deg_int = round(longitude_degrees)
str_longitude = strcompress(string(lon_deg_int),/remove_all)

V_radial_int = round(info.radial_velocity)
str_V_radial = strcompress(string(V_radial_int),/remove_all)
 
str_year = strcompress(string(t_enlil_year),/remove_all)
str_month = strcompress(string(t_enlil_month),/remove_all)
if t_enlil_month lt 10 then str_month = '0' + str_month
str_date = strcompress(string(t_enlil_day),/remove_all)
if t_enlil_day lt 10 then str_date = '0' + str_date
 
str_hour = strcompress(string(t_enlil_hour),/remove_all)
if t_enlil_hour lt 10 then str_hour = '0' + str_hour
str_minute = strcompress(string(t_enlil_minute),/remove_all)
if t_enlil_minute lt 10 then str_minute = '0' + str_minute


; print the results out to the plot_window......

vel_string = strcompress(string(fix(radial_velocity)),/remove_all)
                                 
info.velocity_text -> SetProperty, strings = 'Vel ' + vel_string
info.velocity_text -> SetProperty, hide = 0

info.T21_5_string = str_year + "-" + str_month + "-" + str_date + ' ' + str_hour + ":" + str_minute
                                 
info.T21_5_text -> SetProperty, strings =    'T ' + info.T21_5_string
info.T21_5_text -> SetProperty, hide = 0

info.enlil_info_window->Draw, info.enlil_info_View


; **************write the Enlil cone input file stuff *******************

; once we've done the calculation we can un-grey the export Enlil data button....
Widget_Control,info.widget_export_cone,sensitive = 1

if info.calculate_individual_velocities_for_each_telescope eq 1 then begin
; 
; calculate velocity for each coronagraph independently......
;

sc3_indexes = (info.CME_matches_image_telescope).FindValue('SC3')
sc2_indexes = (info.CME_matches_image_telescope).FindValue('SC2')
ac2_indexes = (info.CME_matches_image_telescope).FindValue('AC2')
bc2_indexes = (info.CME_matches_image_telescope).FindValue('BC2')

if info.debug_mode eq 1 then print, 'SC3 ',n_elements(sc3_indexes) ,': ', sc3_indexes
if info.debug_mode eq 1 then print, 'SC2 ',n_elements(sc2_indexes) ,': ', sc2_indexes
if info.debug_mode eq 1 then print, 'AC2 ',n_elements(ac2_indexes) ,': ', ac2_indexes
if info.debug_mode eq 1 then print, 'BC2 ',n_elements(bc2_indexes) ,': ', bc2_indexes

for i_tele = 1 , 4 do begin
  
if i_tele eq 1 then begin
  tele_indexes = sc3_indexes
endif
if i_tele eq 2 then begin
  tele_indexes = sc2_indexes
endif
if i_tele eq 3 then begin
  tele_indexes = ac2_indexes
endif
if i_tele eq 4 then begin
  tele_indexes = bc2_indexes
endif

if n_elements(tele_indexes) ge 2 then begin
  
  Jul_times_tele = ((info.CME_matches_image_Julian)[tele_indexes]).ToArray()
  R_le_tele = ((info.CME_matches_image_Rs_leading_edge)[tele_indexes]).ToArray()
  sorted_list_tele = sort(Jul_times_tele)
  Jul_times2_tele = Jul_times_tele[sorted_list_tele]
  R_le2_tele = R_le_tele[sorted_list_tele]
  result_R_tele = LINFIT(Jul_times2_tele,R_le2_tele)
  radial_velocity_tele = result_R_tele[1] * info.Solar_radius_km / 86400.d
  t_zero_tele = - (result_R_tele[0] / result_R_tele[1])
  
  if abs(radial_velocity_tele) gt 10.0 then begin
    
  time_from_Sun_to_Enlil_days_tele = (21.5 * info.Solar_radius_km) / (radial_velocity_tele * 86400.d)
  Time_at_Enlil_boundary_tele = t_zero_tele + time_from_Sun_to_Enlil_days_tele
  CALDAT, Time_at_Enlil_boundary_tele, t_enlil_Month_tele , t_enlil_Day_tele , t_enlil_Year_tele , $
        t_enlil_Hour_tele , t_enlil_Minute_tele , t_enlil_Second_tele

    if i_tele eq 1 then begin
      if info.debug_mode eq 1 then print, '*********** LASCO C3 *************'
  thisPlot_tele = info.thisplot2_c3
    endif
    if i_tele eq 2 then begin
      if info.debug_mode eq 1 then print, '*********** LASCO C2 *************'
  thisPlot_tele = info.thisplot2_c2
    endif
    if i_tele eq 3 then begin
      if info.debug_mode eq 1 then print, '*********** STEREO A *************'
  thisPlot_tele = info.thisplot2_stereo_a
    endif
    if i_tele eq 4 then begin
      if info.debug_mode eq 1 then print, '*********** STEREO B *************'
  thisPlot_tele = info.thisplot2_stereo_b
    endif

    if info.debug_mode eq 1 then print, 'velocity: ', radial_velocity_tele
    if info.debug_mode eq 1 then print, 'Enlil time: ',t_enlil_Year_tele, t_enlil_Month_tele, t_enlil_Day_tele, 'T ', $
                              t_enlil_Hour_tele , t_enlil_Minute_tele

    a_tele = result_R_tele[0]
    b_tele = result_R_tele[1]
    yvals_tele = a_tele + (b_tele * Jul_times2_tele)
    thisPlot_tele->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
    thisPlot_tele->SetProperty, dataX = Jul_times2_tele - this_xmin, dataY = yvals_tele
    thisPlot_tele->SetProperty, hide = 0
    
    vel_string_tele = strcompress(string(fix(radial_velocity_tele)),/remove_all)
    
    if i_tele eq 1 then begin
       info.velocity_text_c3 -> SetProperty, strings = vel_string_tele
    endif
    if i_tele eq 2 then begin
       info.velocity_text_c2 -> SetProperty, strings = vel_string_tele
    endif
    if i_tele eq 3 then begin
       info.velocity_text_stereo_a -> SetProperty, strings = vel_string_tele
    endif
    if i_tele eq 4 then begin
       info.velocity_text_stereo_b -> SetProperty, strings = vel_string_tele
    endif
    
    info.enlil_info_window->Draw, info.enlil_info_View
    
    info.plot_window->Draw, info.plot_View

  endif else begin  ; velocity less than 10
    
  if info.debug_mode eq 1 then print, 'tele velocity less than 10 kms'

  endelse          ; velocity less than 10
  
endif else begin   ; less than 2 points for this telescope
  
    if i_tele eq 1 then begin
      if info.debug_mode eq 1 then print, '*********** LASCO C3 *************'
  thisPlot_tele = info.thisplot2_c3
    endif
    if i_tele eq 2 then begin
      if info.debug_mode eq 1 then print, '*********** LASCO C2 *************'
  thisPlot_tele = info.thisplot2_c2
    endif
    if i_tele eq 3 then begin
      if info.debug_mode eq 1 then print, '*********** STEREO A *************'
  thisPlot_tele = info.thisplot2_stereo_a
    endif
    if i_tele eq 4 then begin
      if info.debug_mode eq 1 then print, '*********** STEREO B *************'
  thisPlot_tele = info.thisplot2_stereo_b
    endif  
  
  
    thisPlot_tele->SetProperty, hide = 1
  
endelse            ; less than 2 points for this telescope

endfor

endif ; calculate_individual_velocities_for_each_telescope



endif else begin

Result = DIALOG_MESSAGE("Velocity is less than 10kms-1. Something's amiss",/center)

endelse

end


; ***************************************************************

pro swpc_cat_calculate_velocity, event

;*****************************************************************
;
;  This is executed when the calculate_velocity button is pressed
;  A primary lemniscate has already been defined (locked) and at least 
;  one other lemniscate has been matched to an image.
;  With at least 2 defined positions for the CME the radial velocity
;  can be calculated...
;
;*****************************************************************

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy
; go to the code where we actually calculate the velocity
swpc_cat_do_velocity, info
Widget_Control, event.top, Set_UValue=info, /No_Copy

END



pro swpc_cat_change_L_C_or_R_image_with_slider, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

which_image = event.value - 1



       CASE event.id  OF
          info.L_widget_image_sequence_slider : BEGIN


info.L_current_image_number = which_image

swpc_cat_REDRAW_THE_IMAGE, $
    info.L_current_image_number,info.L_background_image_number,info.L_difference_imaging, $
    info.L_list_of_image_data,info.L_image_saturation_value,info.L_coronagraph_image_object,info.L_border_image_object, $
    info.CME_matches_image_L_Image_number,info.L_current_background_color, $
    info.background_color,info.L_current_text_color,info.color_stereo_B,info.L_cme_outline,info.L_cme_MATCH_outline, $
    info.L_widget_outline_matches_image,info.CME_matches_image_L_CME_outline, $
    info.L_ut_string_object,info.L_list_of_full_time_strings,info.L_title_object,info.L_Window,info.L_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color

      ; ECDP only draw if we are in mode 0 (i.e. normal mode)   
   if (info.which_mode eq 1) then begin
; if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin 
   info.L_cme_outline -> setProperty, hide = 1
   endif else begin
   info.L_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP

    info.L_Window->Draw, info.L_both_views



          END
          
          info.C_widget_image_sequence_slider : BEGIN
          
case info.currently_showing_lasco_C2_or_C3 of

'C3' : begin

info.C_current_image_number = which_image

swpc_cat_REDRAW_THE_IMAGE, $
    info.C_current_image_number,info.C_background_image_number,info.C_difference_imaging, $
    info.C_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c3,info.C_cme_outline,info.C_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C_CME_outline, $
    info.C_ut_string_object,info.C_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color

if info.show_image_line_plot eq 1 then begin

swpc_cat_replot_image_line_plot, info.C_clock_angle_degrees, info.C_coronagraph_image_object, info.C_image_lineplot, $
                            info.position_image_lineplot, info.C_cme_outline
                            
endif

; ECDP only draw if we are in mode 0 (i.e. normal mode)   
   if (info.which_mode eq 1) then begin
;  if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.C_cme_outline -> setProperty, hide = 1
   endif else begin
   info.C_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP

  info.C_Window->Draw, info.C_both_views


endcase


'C2' : begin

info.C2_current_image_number = which_image

swpc_cat_REDRAW_THE_IMAGE, $
    info.C2_current_image_number,info.C2_background_image_number,info.C2_difference_imaging, $
    info.C2_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C2_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c2,info.C_cme_outline,info.C2_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C2_CME_outline, $
    info.C_ut_string_object,info.C2_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, info.C2_current_image_number, info.color_C2, info.cme_outline_color

     C_xp = 255 - (255*sin(info.C_clock_angle_degrees * !dtor))  
     C_yp = 255 + (255*cos(info.C_clock_angle_degrees * !dtor))
     info.C_coronagraph_image_object -> getproperty, data = image_data
     xvals = round(((findgen(255)/254.) * (C_xp - 255.)) + 255.)
     yvals = round(((findgen(255)/254.) * (C_yp - 255.)) + 255.)
     line_Y = image_data[xvals,yvals]
     info.C_image_lineplot->SetProperty, datax = findgen(255) , datay = line_Y
     maxy = max(line_Y,min=minY)
     yrange = [minY,maxY]
     C_image_xs = swpc_cat_FSC_Normalize([0.,254.], Position=[info.position_image_lineplot[0], info.position_image_lineplot[2]])
     C_image_ys = swpc_cat_FSC_Normalize(yrange, Position=[info.position_image_lineplot[1], info.position_image_lineplot[3]])
     info.C_image_lineplot->SetProperty, XCoord_Conv=C_image_xs, YCoord_Conv=C_image_ys

      ; ECDP only draw if we are in mode 0 (i.e. normal mode)   
   if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.C_cme_outline -> setProperty, hide = 1
   endif else begin
   info.C_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP

info.C_Window->Draw, info.C_both_views

endcase

endcase

         END
         
          info.R_widget_image_sequence_slider : BEGIN


info.R_current_image_number = which_image

swpc_cat_REDRAW_THE_IMAGE, $
    info.R_current_image_number,info.R_background_image_number,info.R_difference_imaging, $
    info.R_list_of_image_data,info.R_image_saturation_value,info.R_coronagraph_image_object,info.R_border_image_object, $
    info.CME_matches_image_R_Image_number,info.R_current_background_color, $
    info.background_color,info.R_current_text_color,info.color_stereo_A,info.R_cme_outline,info.R_cme_MATCH_outline, $
    info.R_widget_outline_matches_image,info.CME_matches_image_R_CME_outline, $
    info.R_ut_string_object,info.R_list_of_full_time_strings,info.R_title_object,info.R_Window,info.R_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_a, info.cme_outline_color

      ; ECDP only draw if we are in mode 0 (i.e. normal mode)   
   if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.R_cme_outline -> setProperty, hide = 1
   endif else begin
   info.R_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP

info.R_Window->Draw, info.R_both_views


         END
         
ENDCASE

info.images_timeline_window->Draw, info.images_timeline_view


Widget_Control, event.top, Set_UValue=info, /No_Copy


END







pro swpc_cat_image_difference_and_scaling, background_color, current_image_number, background_image_number, difference_imaging, $
                 list_of_image_data, image_saturation_value, coronagraph_image_object, border_image_object, i_log_scale

;print, ' current_image_number ',current_image_number

       CASE difference_imaging  OF
          'none' : BEGIN

image_data = bytscl((list_of_image_data)[current_image_number])

          END
          'running' : BEGIN
          
image_data = (list_of_image_data)[current_image_number]
image_data_background = (list_of_image_data)[current_image_number - 1 > 0]
if i_log_scale eq 0 then begin
image_data = (image_data - image_data_background) < image_saturation_value
image_data = bytscl((0. - image_saturation_value) > image_data)
endif else begin
image_data = (alog10(image_data) - alog10(image_data_background)) < alog10(image_saturation_value)
image_data = bytscl((0. - alog10(image_saturation_value)) > image_data)
endelse

          END
          'base' : BEGIN

image_data = (list_of_image_data)[current_image_number]
image_data_background = (list_of_image_data)[background_image_number]
image_data = (image_data - image_data_background) < image_saturation_value
image_data = bytscl((0. - image_saturation_value) > image_data)

          END
       ENDCASE
       
coronagraph_image_object -> setproperty, data = image_data

swpc_cat_set_image_border, background_color, border_image_object  

end










pro swpc_cat_set_start_and_end_time, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

       reset_for_new_images_as_well_as_the_cme_analysis = 1
       swpc_cat_full_reset, info, reset_for_new_images_as_well_as_the_cme_analysis 

       date_array = info.date_array

       ;check the validity of the dates entered
       checking = 1
       while checking eq 1 do begin

          swpc_cat_check_dates,date_array,start_julian,end_julian,info.max_interval_in_days, $
                      dates_are_fine,problem_string

          if dates_are_fine eq 0 then begin  ;there's a problem with one or both dates...
          Result = DIALOG_MESSAGE(problem_string,/center)
          Widget_Control, event.top, Set_UValue=info, /No_Copy
          return
          endif  else begin
             checking = 0  ;have valid dates; so, done checking them
          endelse
       endwhile

; At this point we should have valid dates.......

       ;print, 'VALID DATES :'
;       print,' date_array ', date_array
       info.date_array_int = fix(date_array)
       ;print, 'date_array_int ',info.date_array_int

;       Result = DIALOG_MESSAGE(date_array,/center)

       info.start_date = info.date_array_int[0:4]
       info.end_date = info.date_array_int[5:9]
       info.start_julian = start_julian
       info.end_julian = end_julian

          info.which_telescope = 'LASCO C3'
          info.C_telescope_code = 'SC3'
          info.C_telescope_field_of_view_Rs = 30.
          swpc_cat_get_images,info, 1
          
          info.which_telescope = 'LASCO C2'
          info.C2_telescope_code = 'SC2'
          info.C2_telescope_field_of_view_Rs = 6.
          swpc_cat_get_images,info, 3

          info.which_telescope = 'STEREO A COR2'
          info.R_telescope_code = 'AC2'
          info.R_telescope_field_of_view_Rs = 15.
          swpc_cat_get_images,info, 2

          info.which_telescope = 'STEREO B COR2'
          info.L_telescope_code = 'BC2'
          info.L_telescope_field_of_view_Rs = 15.
          swpc_cat_get_images,info, 0

        
               
master_list = info.L_list_of_datetime_strings + info.C_list_of_datetime_strings + info.R_list_of_datetime_strings

if n_elements(master_list) eq 0 then begin
  
  msg1 = 'There are no images for those dates. Try again!!'
  ok = DIALOG_MESSAGE(msg1,/ERROR,/center)
  
endif else begin

;help, master_list
;print, master_list

master_array = master_list.ToArray()

;print, master_array

master_array2 = master_array[sort(master_array)]
master_array = master_array2[uniq(master_array2)]

(info.master_list).add, master_array,/extract

info.master_list_size = n_elements(master_array)
count_L = 0
count_C = 0
count_R = 0
for i = 0 , info.master_list_size - 1 do begin
Result_L = (info.L_list_of_datetime_strings).Count(master_array[i])
count_L = count_L + Result_L
Result_C = (info.C_list_of_datetime_strings).Count(master_array[i])
count_C = count_C + Result_C
Result_R = (info.R_list_of_datetime_strings).Count(master_array[i])
count_R = count_R + Result_R
if result_L gt 0 then (info.L_indexes).add, count_L else (info.L_indexes).add, result_L
if result_C gt 0 then (info.C_indexes).add, count_C else (info.C_indexes).add, result_C
if result_R gt 0 then (info.R_indexes).add, count_R else (info.R_indexes).add, result_R
endfor

;print, ' *********************'
;print, info.master_list
;print, ' *********************'
;print, info.L_indexes
;print, ' *********************'
;print, info.C_indexes
;print, ' *********************'
;print, info.R_indexes
;print, ' *********************'

;widget_control,info.widget_image_sequence_slider, set_slider_min = 1
;widget_control,info.widget_image_sequence_slider, set_slider_max = info.master_list_size
;widget_control,info.widget_image_sequence_slider, sensitive=1


; line plot stuff..............   

;set_xmin_and_xmax, info.start_julian, info.end_julian, info.xmin, info.xmax, info.xtickinterval

result = label_date(date_format='%H',offset=info.start_julian)
xrange = [0.d,info.end_julian - info.start_julian]
yrange = [0.,1.0]
an_hour = 1./24.
six_hours = 0.25

xs = swpc_cat_FSC_Normalize(xrange, Position=[info.position_timeline[0], info.position_timeline[2]])
ys = swpc_cat_FSC_Normalize(yrange, Position=[info.position_timeline[1], info.position_timeline[3]])

info.xaxis_images_timeline->SetProperty, exact = 1
info.xaxis_images_timeline->SetProperty, range=xrange
info.xaxis_images_timeline->SetProperty, XCoord_Conv=xs
info.xaxis_images_timeline->SetProperty, tickinterval = an_hour
info.xaxis_images_timeline->SetProperty, tickformat='label_date'

;print, ' info.L_number_of_images ', info.L_number_of_images


symSize = 1.0
plot_position_scaling_factor = (info.position_timeline[2] - info.position_timeline[0]) / (info.position_timeline[3] - info.position_timeline[1]) 
info.xSymbolSize_timeline = (info.end_julian - info.start_julian) * 0.015 * symSize / 4.
info.ySymbolSize_timeline = (1.0 - 0.0) * 0.015 * plot_position_scaling_factor * symSize * 4.

; ECDP
if info.L_number_of_images gt 0 then begin
;if ((info.L_number_of_images gt 0) and (info.which_mode eq 0 )) then begin
; if we are in NORMAL mode

L_yvals = fltarr(n_elements(info.L_list_of_datetime_Julian)) + 0.2
info.L_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.L_plot->SetProperty, dataX = (info.L_list_of_datetime_Julian).toarray() - info.start_julian, dataY = L_yvals
info.L_plot->SetProperty, color=[100,100,255]
info.L_plot->GetProperty, symbol=L_thisSymbol

L_points = n_elements(L_yvals)
symarray = objarr(L_points)
for i = 0 , L_points - 1 do begin
num = 6
filled = 0
thisSymbol_L = obj_new("IDLgrsymbol",data=num, Size=[info.xSymbolSize_timeline, info.ySymbolSize_timeline],filled = filled)
symarray[i] = thisSymbol_L
endfor
info.L_plot->SetProperty, symbol=symarray
; ECDP - in the new set up we don't want to see the lemniscate until the mode has been selected
info.L_cme_outline -> setProperty, hide = 1
; ECDP
swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color
widget_control, info.L_widget_image_sequence_slider,set_slider_max = info.L_number_of_images
widget_control,info.L_widget_image_sequence_slider,set_value = info.L_current_image_number + 1
info.L_ut_string_object->SetProperty, strings = (info.L_list_of_full_time_strings)[info.L_current_image_number]
endif else begin
info.L_ut_string_object->SetProperty, strings = 'NO IMAGES'
info.L_cme_outline -> setProperty, hide = 1
widget_control, info.L_widget_image_sequence_slider,sensitive=0
endelse
info.L_Window->Draw, info.L_both_views

if info.C_number_of_images gt 0 then begin

C_yvals = fltarr(n_elements(info.C_list_of_datetime_Julian)) + 0.65
info.C_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.C_plot->SetProperty, dataX = (info.C_list_of_datetime_Julian).toarray() - info.start_julian, dataY = C_yvals
info.C_plot->SetProperty, color=[100,255,100]
info.C_plot->GetProperty, symbol=C_thisSymbol

C_points = n_elements(C_yvals)
symarray = objarr(C_points)
for i = 0 , C_points - 1 do begin
num = 6
filled = 0
thisSymbol_C = obj_new("IDLgrsymbol",data=num, Size=[info.xSymbolSize_timeline, info.ySymbolSize_timeline],filled = filled)
symarray[i] = thisSymbol_C
endfor
info.C_plot->SetProperty, symbol=symarray

; ECDP - in the new set up we don't want to see the lemniscate until the mode has been selected
info.C_cme_outline -> setProperty, hide = 1
; ECDP

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color
widget_control, info.C_widget_image_sequence_slider,set_slider_max = info.C_number_of_images
widget_control,info.C_widget_image_sequence_slider,set_value = info.C_current_image_number + 1
info.C_ut_string_object->SetProperty, strings = (info.C_list_of_full_time_strings)[info.C_current_image_number]
endif else begin
info.C_ut_string_object->SetProperty, strings = 'NO IMAGES'
info.C_cme_outline -> setProperty, hide = 1
widget_control, info.C_widget_image_sequence_slider,sensitive=0
endelse
info.C_Window->Draw, info.C_both_views

if info.C2_number_of_images gt 0 then begin
; ECDP  
info.C_cme_outline -> setProperty, hide = 1
; ECDP
C2_yvals = fltarr(n_elements(info.C2_list_of_datetime_Julian)) + 0.4
info.C2_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.C2_plot->SetProperty, dataX = (info.C2_list_of_datetime_Julian).toarray() - info.start_julian, dataY = C2_yvals
info.C2_plot->SetProperty, color=info.color_c2
info.C2_plot->GetProperty, symbol=C2_thisSymbol

C2_points = n_elements(C2_yvals)
symarray = objarr(C2_points)
for i = 0 , C2_points - 1 do begin
num = 6
filled = 0
thisSymbol_C2 = obj_new("IDLgrsymbol",data=num, Size=[info.xSymbolSize_timeline, info.ySymbolSize_timeline],filled = filled)
symarray[i] = thisSymbol_C2
endfor
info.C2_plot->SetProperty, symbol=symarray

endif
info.C_Window->Draw, info.C_both_views


; ECDP
if info.R_number_of_images gt 0 then begin
;if ((info.R_number_of_images gt 0) and (info.which_mode eq 0 )) then begin
; ECDP

R_yvals = fltarr(n_elements(info.R_list_of_datetime_Julian)) + 0.9
info.R_plot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
info.R_plot->SetProperty, dataX = (info.R_list_of_datetime_Julian).toarray() - info.start_julian, dataY = R_yvals
info.R_plot->SetProperty, color=[255,100,100]
info.R_plot->GetProperty, symbol=R_thisSymbol

R_points = n_elements(R_yvals)
symarray = objarr(R_points)
for i = 0 , R_points - 1 do begin
num = 6
filled = 0
thisSymbol_R = obj_new("IDLgrsymbol",data=num, Size=[info.xSymbolSize_timeline, info.ySymbolSize_timeline],filled = filled)
symarray[i] = thisSymbol_R
endfor
info.R_plot->SetProperty, symbol=symarray

; ECDP - in the new set up we don't want to see the lemniscate until the mode has been selected
info.R_cme_outline -> setProperty, hide = 1
; ECDP

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_A, info.cme_outline_color
widget_control, info.R_widget_image_sequence_slider,set_slider_max = info.R_number_of_images
widget_control,info.R_widget_image_sequence_slider,set_value = info.R_current_image_number + 1
info.R_ut_string_object->SetProperty, strings = (info.R_list_of_full_time_strings)[info.R_current_image_number]
endif else begin
info.R_ut_string_object->SetProperty, strings = 'NO IMAGES'
info.R_cme_outline -> setProperty, hide = 1
widget_control, info.R_widget_image_sequence_slider,sensitive=0
endelse
info.R_Window->Draw, info.R_both_views

line_data=fltarr(2,2)
line_data[0,0] = 0. + 1./24.
line_data[1,0] = 0.
line_data[0,1] = 0. + 1./24.
line_data[1,1] = 1.
info.animation_start_time_marker->setproperty,data=line_data

xpos = 1./24.
info.animation_start_time_marker_handle_data=[[xpos,0.],[xpos+0.005,0.],[xpos+0.005,0.08],[xpos,0.08]]
info.animation_start_time_marker_handle->setproperty,data=info.animation_start_time_marker_handle_data

line_data=fltarr(2,2)
line_data[0,0] = info.end_julian - info.start_julian - 1./24.
line_data[1,0] = 0.
line_data[0,1] = info.end_julian - info.start_julian - 1./24.
line_data[1,1] = 1.
info.animation_end_time_marker->setproperty,data=line_data

xpos = info.end_julian - info.start_julian - 1./24.
info.animation_end_time_marker_handle_data=[[xpos,0.],[xpos-0.005,0.],[xpos-0.005,0.08],[xpos,0.08]]
info.animation_end_time_marker_handle->setproperty,data=info.animation_end_time_marker_handle_data


line_data=fltarr(2,2)
line_data[0,0] = (info.end_julian - info.start_julian) / 2.
line_data[1,0] = 0.
line_data[0,1] = (info.end_julian - info.start_julian) / 2.
line_data[1,1] = 1.
info.animation_current_time_marker->setproperty,data=line_data


info.images_timeline_window->Draw, info.images_timeline_view

info.xmin = info.start_julian
info.xmax = info.end_julian
info.ymin = 0.
info.ymax = 60.

dataX=dblarr(1)
dataY=fltarr(1)
dataX[0] = (info.xmin + info.xmax) / 2.
dataY[0] = 15.

info.xAxis1_B -> SetProperty, hide = 0
info.yAxis1_B -> SetProperty, hide = 0
info.xAxis2_B -> SetProperty, hide = 0
info.yAxis2_B -> SetProperty, hide = 0
info.LE_plot_matched_CMEs -> SetProperty, hide = 0
vert_colors = info.background_color

swpc_cat_update_the_line_plot,info.position_B,info.xmin,info.xmax,info.ymin,info.ymax, $
                                 info.xaxis1_B,info.xaxis2_B,info.yaxis1_B,info.yaxis2_B, $
                                 info.xtickinterval, $                                 
                                 info.LE_plot_matched_CMEs,info.plot_view,info.plot_window,dataX,dataY,vert_colors, $
                                 xs, ys
                                 
                                 
info.images_are_loaded = 1

widget_control, info.animation_controls_outer_base, sensitive=1

widget_control, info.image_controls_base, sensitive=1
widget_control, info.cme_controls_base, sensitive=1
widget_control, info.plot_window_base, sensitive=1
widget_control, info.enlil_info_base, sensitive=1

widget_control, info.L_widget_image_sequence_slider, sensitive=1
widget_control, info.C_widget_image_sequence_slider, sensitive=1
widget_control, info.R_widget_image_sequence_slider, sensitive=1
  
endelse
      
Widget_Control, event.top, Set_UValue=info, /No_Copy

END






pro swpc_cat_get_images, info , i 

compile_opt idl2 

 CASE i OF
   0 : Begin
list_of_image_names = info.L_list_of_image_names
list_of_image_data = info.L_list_of_image_data
list_of_datetime_strings = info.L_list_of_datetime_strings
list_of_datetime_Julian = info.L_list_of_datetime_Julian
list_of_full_time_strings = info.L_list_of_full_time_strings
list_of_image_exposure_times = info.L_list_of_image_exposure_times
list_of_image_offsets = info.L_list_of_image_offsets
list_of_image_scaling_factors = info.L_list_of_image_scaling_factors
list_of_HEEQ_coords = info.L_list_of_HEEQ_coords
list_of_pixel_scales = info.L_list_of_pixel_scales
list_of_rsuns = info.L_list_of_rsuns
list_of_Sun_satellite_distances = info.L_list_of_Sun_satellite_distances
   End
   1    : Begin
list_of_image_names = info.C_list_of_image_names
list_of_image_data = info.C_list_of_image_data
list_of_datetime_strings = info.C_list_of_datetime_strings
list_of_datetime_Julian = info.C_list_of_datetime_Julian
list_of_full_time_strings = info.C_list_of_full_time_strings
list_of_image_exposure_times = info.C_list_of_image_exposure_times
list_of_image_offsets = info.C_list_of_image_offsets
list_of_image_scaling_factors = info.C_list_of_image_scaling_factors
list_of_HEEQ_coords = info.C_list_of_HEEQ_coords
list_of_pixel_scales = info.C_list_of_pixel_scales
list_of_rsuns = info.C_list_of_rsuns
list_of_Sun_satellite_distances = info.C_list_of_Sun_satellite_distances
   End  
   2      : Begin
list_of_image_names = info.R_list_of_image_names
list_of_image_data = info.R_list_of_image_data
list_of_datetime_strings = info.R_list_of_datetime_strings
list_of_datetime_Julian = info.R_list_of_datetime_Julian
list_of_full_time_strings = info.R_list_of_full_time_strings
list_of_image_exposure_times = info.R_list_of_image_exposure_times
list_of_image_offsets = info.R_list_of_image_offsets
list_of_image_scaling_factors = info.R_list_of_image_scaling_factors
list_of_HEEQ_coords = info.R_list_of_HEEQ_coords
list_of_pixel_scales = info.R_list_of_pixel_scales
list_of_rsuns = info.R_list_of_rsuns
list_of_Sun_satellite_distances = info.R_list_of_Sun_satellite_distances
   End
   3    : Begin
list_of_image_names = info.C2_list_of_image_names
list_of_image_data = info.C2_list_of_image_data
list_of_datetime_strings = info.C2_list_of_datetime_strings
list_of_datetime_Julian = info.C2_list_of_datetime_Julian
list_of_full_time_strings = info.C2_list_of_full_time_strings
list_of_image_exposure_times = info.C2_list_of_image_exposure_times
list_of_image_offsets = info.C2_list_of_image_offsets
list_of_image_scaling_factors = info.C2_list_of_image_scaling_factors
list_of_HEEQ_coords = info.C2_list_of_HEEQ_coords
list_of_pixel_scales = info.C2_list_of_pixel_scales
list_of_rsuns = info.C2_list_of_rsuns
list_of_Sun_satellite_distances = info.C2_list_of_Sun_satellite_distances
   End 
 ENDCASE

;get 1st image from the start date folder

   swpc_cat_progressbar = Obj_New('swpc_cat_progressbar', Color='Orange', Text='Loading Image files',/nocancel, background = 'White', $
                    xsize = 500 , ysize = 20 , title='Getting Images') ;, group_leader = info.tlb)
   swpc_cat_progressbar -> Start
   
   progress_base_percentage = 10.

swpc_cat_find_images_by_date,info.start_date,info,list_of_image_names,list_of_image_data, $
                    list_of_datetime_strings, list_of_image_exposure_times, $
                    list_of_image_offsets,list_of_image_scaling_factors,swpc_cat_progressbar, $
                    progress_base_percentage, list_of_HEEQ_coords, list_of_pixel_scales, $
                    list_of_rsuns,list_of_Sun_satellite_distances

;get beginning of start_date,end_date, in julian date format
start_0000 = JULDAY(info.start_date[1],info.start_date[2],info.start_date[0])
end_0000 = JULDAY(info.end_date[1],info.end_date[2],info.end_date[0])

number_of_day_folders = end_0000 - start_0000 + 1

;print,' number_of_day_folders ', number_of_day_folders

;get remaining images, if the date/time interval spans multiple days
next_date = start_0000
if (start_0000 NE end_0000) then begin  ;start date is not the same as end date
    while next_date LT end_0000 do begin
      ;get next date
      next_date = next_date + 1.0D
      CALDAT,next_date,next_month,next_day,next_year,next_hour,next_minute,next_sec
      next_date_array = [next_year,next_month,next_day]
      swpc_cat_find_images_by_date,next_date_array,info,list_of_image_names,list_of_image_data, $
                          list_of_datetime_strings, list_of_image_exposure_times, $
                          list_of_image_offsets,list_of_image_scaling_factors,swpc_cat_progressbar, $
                          progress_base_percentage, list_of_HEEQ_coords, list_of_pixel_scales, $
                          list_of_rsuns,list_of_Sun_satellite_distances
    endwhile    

endif

if progress_base_percentage lt 80. then swpc_cat_progressbar -> Update, 90., Text=''

number_of_images = n_elements(list_of_image_names)
if number_of_images eq 0 then begin
;   help,info.start_date
;   help,info.end_date
;   print, info.start_date , info.end_date
   st = string(info.start_date)
   en = string(info.end_date)
   start_str = strcompress(st[0]+st[1]+st[2]+st[3]+st[4])
   end_str = strcompress(en[0]+en[1]+en[2]+en[3]+en[4])
   msg1='No images found for ' + info.which_telescope + ' for interval'
   msg2=start_str + ' through ' + end_str + '.'
;   msg3='Please use "Set Dates" to select another date/time interval.'
;   ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
endif else begin

image_time_Julian = dblarr(number_of_images)
full_time_string = strarr(number_of_images)

for it = 0 , number_of_images - 1 do begin
  datetime_string = list_of_datetime_strings[it]
  year_string = strmid(datetime_string,0,4)
  month_string = strmid(datetime_string,5,2)
  day_string = strmid(datetime_string,8,2)
  hour_string = strmid(datetime_string,11,2)
  min_string = strmid(datetime_string,14,2)
  sec_string = strmid(datetime_string,17,2)
  year = fix(year_string)
  month = fix(month_string)
  day = fix(day_string)
  hours = fix(hour_string)
  mins = fix(min_string)
  secs = fix(sec_string)
  list_of_datetime_Julian.add, JULDAY(month, day, year, hours, mins, secs)
  date_string = year_string + '-' + month_string + '-' + day_string
  list_of_full_time_strings.add, date_string + ' ' + hour_string + ':' + min_string + ' UT'
endfor

CASE i OF
   0 : Begin

;info.anim_end_frame = number_of_images
info.L_number_of_images = number_of_images
info.L_list_of_image_names = list_of_image_names
info.L_list_of_image_data = list_of_image_data
info.L_list_of_datetime_strings = list_of_datetime_strings
info.L_list_of_datetime_Julian = list_of_datetime_Julian
info.L_list_of_full_time_strings = list_of_full_time_strings
info.L_list_of_image_exposure_times = list_of_image_exposure_times
info.L_list_of_image_offsets = list_of_image_offsets
info.L_list_of_image_scaling_factors = list_of_image_scaling_factors
info.L_list_of_HEEQ_coords = list_of_HEEQ_coords
info.L_list_of_pixel_scales = list_of_pixel_scales
info.L_list_of_rsuns = list_of_rsuns
info.L_list_of_Sun_satellite_distances = list_of_Sun_satellite_distances

; special case for where we just have 1 image:
if info.L_number_of_images eq 1 then begin
info.L_current_image_number = 0
info.L_difference_imaging = 'none'
endif

swpc_cat_image_difference_and_scaling, info.background_color,  info.L_current_image_number, info.L_background_image_number, info.L_difference_imaging, $
                 info.L_list_of_image_data, info.L_image_saturation_value, info.L_coronagraph_image_object, info.L_border_image_object, info.i_log_scale
                                                                
info.L_HEEQ_coords = list_of_HEEQ_coords[0]
;print, 'L HEEQ longitude = ', info.L_HEEQ_coords

info.L_coords_offset_model -> SetProperty,transform = info.initial_transform
info.L_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.L_coords_offset_model->rotate,[0,1,0], 0.0 - info.L_HEEQ_coords[0]
info.L_coords_offset_model_copy->rotate,[0,1,0], 0.0 - info.L_HEEQ_coords[0]
info.L_coords_offset_model->rotate,[1,0,0], info.L_HEEQ_coords[1]
info.L_coords_offset_model_copy->rotate,[1,0,0], info.L_HEEQ_coords[1]

;print, 'pixel scale L ', (info.L_list_of_pixel_scales)[0]
;print, ' L rsun ',(info.L_list_of_rsuns)[0]
;print, 'Field of view L', (256. * ((info.L_list_of_pixel_scales)[0] / (info.L_list_of_image_scaling_factors)[0])) / (info.L_list_of_rsuns)[0]

info.L_telescope_field_of_view_Rs = (256. * ((info.L_list_of_pixel_scales)[0] / (info.L_list_of_image_scaling_factors)[0])) / (info.L_list_of_rsuns)[0]

info.L_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.L_telescope_field_of_view_Rs,2.0*info.L_telescope_field_of_view_Rs]
info.L_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.L_telescope_field_of_view_Rs,2.0*info.L_telescope_field_of_view_Rs]
info.L_second_view -> setproperty, eye = (info.L_list_of_Sun_satellite_distances)[0]
info.L_second_view_copy -> setproperty, eye = (info.L_list_of_Sun_satellite_distances)[0]


info.L_Fixed_Sun_Satellite_model -> SetProperty,transform = info.L_Fixed_Sun_Satellite_model_original_transform
info.L_Fixed_Sun_Satellite_model -> translate,  info.L_telescope_field_of_view_Rs , info.L_telescope_field_of_view_Rs , 0.0
info.L_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.L_Fixed_Sun_Satellite_model_original_transform
info.L_Fixed_Sun_Satellite_model_copy -> translate,  info.L_telescope_field_of_view_Rs , info.L_telescope_field_of_view_Rs, 0.0
;info.L_cme_outline -> setProperty, hide = 0
info.L_cme_outline -> setProperty, hide = 1
swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline
;if (info.which_mode eq 1) then begin
;   info.L_cme_outline -> setProperty, hide = 1
;   endif else begin
;   info.L_cme_outline -> setProperty, hide = 0 
;   endelse 
; ECDP
info.L_Window->Draw, info.L_both_views


   End
   1    : Begin

info.anim_end_frame = number_of_images - 1  
info.C_number_of_images = number_of_images
info.C_list_of_image_names = list_of_image_names
info.C_list_of_image_data = list_of_image_data
info.C_list_of_datetime_strings = list_of_datetime_strings
info.C_list_of_datetime_Julian = list_of_datetime_Julian
info.C_list_of_full_time_strings = list_of_full_time_strings
info.C_list_of_image_exposure_times = list_of_image_exposure_times
info.C_list_of_image_offsets = list_of_image_offsets
info.C_list_of_image_scaling_factors = list_of_image_scaling_factors
info.C_list_of_HEEQ_coords = list_of_HEEQ_coords
info.C_list_of_pixel_scales = list_of_pixel_scales
info.C_list_of_rsuns = list_of_rsuns
info.C_list_of_Sun_satellite_distances = list_of_Sun_satellite_distances

; special case for where we just have 1 image:
if info.C_number_of_images eq 1 then begin
info.C_current_image_number = 0
info.C_difference_imaging = 'none'
endif

swpc_cat_image_difference_and_scaling, info.background_color,  info.C_current_image_number, info.C_background_image_number, info.C_difference_imaging, $
                 info.C_list_of_image_data, info.C_image_saturation_value, info.C_coronagraph_image_object, info.C_border_image_object, info.i_log_scale

heeq_coords = list_of_HEEQ_coords[0]
info.C_HEEQ_coords[0] = heeq_coords[0]
swpc_cat_Calculate_Earth_B_Angle,list_of_datetime_Julian[0],B_angle_degrees
info.C_HEEQ_coords[1] = B_angle_degrees
;print, 'C HEEQ coords = ', info.C_HEEQ_coords

info.C_coords_offset_model -> SetProperty,transform = info.initial_transform
info.C_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.C_coords_offset_model->rotate,[1,0,0], info.C_HEEQ_coords[1]
info.C_coords_offset_model_copy->rotate,[1,0,0], info.C_HEEQ_coords[1]

;print, 'pixel scale C ', (info.C_list_of_pixel_scales)[0]
;print, ' C rsun ',(info.C_list_of_rsuns)[0]
;print, 'C image scaling factors: ', (info.C_list_of_image_scaling_factors)[0]
;print, 'Field of view C', (256. * ((info.C_list_of_pixel_scales)[0] / (info.C_list_of_image_scaling_factors)[0])) / (info.C_list_of_rsuns)[0]

info.C_telescope_field_of_view_Rs = (256. * ((info.C_list_of_pixel_scales)[0] / (info.C_list_of_image_scaling_factors)[0])) / (info.C_list_of_rsuns)[0]

info.C_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]
info.C_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]

the_day = long((info.C_list_of_datetime_Julian)[0])
i_day = where(the_day lt info.Julian_day_for_Earth_pos)
i_day = i_day[0]
;print,  info.Earth_pos_AU[i_day] , info.Earth_pos_HG_LAT_deg[i_day]

info.C_second_view -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99  ; 0.99 factor is for L1 as opposed to Earth.
info.C_second_view_copy -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99

info.C_Fixed_Sun_Satellite_model -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
info.C_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model_copy -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
;info.C_cme_outline -> setProperty, hide = 0
info.C_cme_outline -> setProperty, hide = 1
swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
info.C_Window->Draw, info.C_both_views

   End  
   2      : Begin
   
info.R_number_of_images = number_of_images
info.R_list_of_image_names = list_of_image_names
info.R_list_of_image_data = list_of_image_data
info.R_list_of_datetime_strings = list_of_datetime_strings
info.R_list_of_datetime_Julian = list_of_datetime_Julian
info.R_list_of_full_time_strings = list_of_full_time_strings
info.R_list_of_image_exposure_times = list_of_image_exposure_times
info.R_list_of_image_offsets = list_of_image_offsets
info.R_list_of_image_scaling_factors = list_of_image_scaling_factors
info.R_list_of_HEEQ_coords = list_of_HEEQ_coords
info.R_list_of_pixel_scales = list_of_pixel_scales
info.R_list_of_rsuns = list_of_rsuns
info.R_list_of_Sun_satellite_distances = list_of_Sun_satellite_distances

; special case for where we just have 1 image:
if info.R_number_of_images eq 1 then begin
info.R_current_image_number = 0
info.R_difference_imaging = 'none'
endif

swpc_cat_image_difference_and_scaling, info.background_color,  info.R_current_image_number, info.R_background_image_number, info.R_difference_imaging, $
                 info.R_list_of_image_data, info.R_image_saturation_value, info.R_coronagraph_image_object, info.R_border_image_object, info.i_log_scale
                                
info.R_HEEQ_coords = list_of_HEEQ_coords[0]
;print, 'R HEEQ longitude = ', info.R_HEEQ_coords

info.R_coords_offset_model -> SetProperty,transform = info.initial_transform
info.R_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.R_coords_offset_model->rotate,[0,1,0], 0.0 - info.R_HEEQ_coords[0]
info.R_coords_offset_model_copy->rotate,[0,1,0], 0.0 - info.R_HEEQ_coords[0]
info.R_coords_offset_model->rotate,[1,0,0], info.R_HEEQ_coords[1]
info.R_coords_offset_model_copy->rotate,[1,0,0], info.R_HEEQ_coords[1]

;print, 'pixel scale R ', (info.R_list_of_pixel_scales)[0]
;print, ' R rsun ',(info.R_list_of_rsuns)[0]
;print, 'Field of view R', (256. * ((info.R_list_of_pixel_scales)[0] / (info.R_list_of_image_scaling_factors)[0])) / (info.R_list_of_rsuns)[0]

info.R_telescope_field_of_view_Rs = (256. * ((info.R_list_of_pixel_scales)[0] / (info.R_list_of_image_scaling_factors)[0])) / (info.R_list_of_rsuns)[0]

info.R_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.R_telescope_field_of_view_Rs,2.0*info.R_telescope_field_of_view_Rs]
info.R_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.R_telescope_field_of_view_Rs,2.0*info.R_telescope_field_of_view_Rs]
info.R_second_view -> setproperty, eye = (info.R_list_of_Sun_satellite_distances)[0]
info.R_second_view_copy -> setproperty, eye = (info.R_list_of_Sun_satellite_distances)[0]

info.R_Fixed_Sun_Satellite_model -> SetProperty,transform = info.R_Fixed_Sun_Satellite_model_original_transform
info.R_Fixed_Sun_Satellite_model -> translate, info.R_telescope_field_of_view_Rs , info.R_telescope_field_of_view_Rs, 0.0
info.R_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.R_Fixed_Sun_Satellite_model_original_transform
info.R_Fixed_Sun_Satellite_model_copy -> translate, info.R_telescope_field_of_view_Rs , info.R_telescope_field_of_view_Rs, 0.0
;info.R_cme_outline -> setProperty, hide = 0
info.R_cme_outline -> setProperty, hide = 1
swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline

info.R_Window->Draw, info.R_both_views

   End
   3    : Begin
   
info.C2_number_of_images = number_of_images
info.C2_list_of_image_names = list_of_image_names
info.C2_list_of_image_data = list_of_image_data
info.C2_list_of_datetime_strings = list_of_datetime_strings
info.C2_list_of_datetime_Julian = list_of_datetime_Julian
info.C2_list_of_full_time_strings = list_of_full_time_strings
info.C2_list_of_image_exposure_times = list_of_image_exposure_times
info.C2_list_of_image_offsets = list_of_image_offsets
info.C2_list_of_image_scaling_factors = list_of_image_scaling_factors
info.C2_list_of_HEEQ_coords = list_of_HEEQ_coords
info.C2_list_of_pixel_scales = list_of_pixel_scales
info.C2_list_of_rsuns = list_of_rsuns
info.C2_list_of_Sun_satellite_distances = list_of_Sun_satellite_distances

if info.debug_mode eq 1 then print, '******************'
if info.debug_mode eq 1 then print, 'C2 lists rsuns ',list_of_rsuns
if info.debug_mode eq 1 then print, '******************'
if info.debug_mode eq 1 then print, 'C2 lists pixel scales ',list_of_pixel_scales
if info.debug_mode eq 1 then print, '******************'
if info.debug_mode eq 1 then print, 'C2 lists scaling ',list_of_image_scaling_factors
if info.debug_mode eq 1 then print, '******************'

; special case for where we just have 1 image:
if info.C2_number_of_images eq 1 then begin
info.C2_current_image_number = 0
info.C2_difference_imaging = 'none'
endif

   End  
ENDCASE

; now we have moved the initial ellipse into view - also show the default cone parameters
; in the info window.....

radial_distance_string = 'dist :' + string(info.radial_distance_lemniscate, format='(f6.1)')
info.radial_distance_string_object -> setproperty, strings = radial_distance_string

angular_width_string = 'cone :' + string(info.angular_width_lemniscate, format='(f6.1)')
info.angular_width_string_object -> setproperty, strings = angular_width_string

info.lat_string = string(info.latitude_degrees,format='(f6.1)')
info.lon_string = string(info.longitude_degrees,format='(f6.1)')
info.lat_string_object -> setproperty, strings = 'lat  :' + info.lat_string
info.lon_string_object -> setproperty, strings = 'lon  :' + info.lon_string

info.cme_info_model->SetProperty, hide = 0
info.cme_info_Window->Draw, info.cme_info_view


endelse

widget_control, info.L_widget_botSlider, sensitive=1
widget_control, info.L_widget_topSlider, sensitive=1
widget_control, info.L_widget_gammaSlider, sensitive=1
widget_control, info.L_widget_saturationSlider, sensitive=1
widget_control, info.C_widget_botSlider, sensitive=1
widget_control, info.C_widget_topSlider, sensitive=1
widget_control, info.C_widget_gammaSlider, sensitive=1
widget_control, info.C_widget_saturationSlider, sensitive=1
widget_control, info.R_widget_botSlider, sensitive=1
widget_control, info.R_widget_topSlider, sensitive=1
widget_control, info.R_widget_gammaSlider, sensitive=1
widget_control, info.R_widget_saturationSlider, sensitive=1

;widget_control,info.widget_image_menu, sensitive=1
widget_control, info.widget_radial_distance_slider, sensitive = 0
widget_control, info.widget_angular_width_slider, sensitive = 0
; 

; 


   swpc_cat_progressbar -> Destroy

   
;print, '*****************************************************'
;print, ' STEREO B '
;print, info.L_list_of_datetime_strings
;print, info.L_list_of_full_time_strings
;;print, ' STEREO B ',info.L_list_of_image_names
;print, '*****************************************************'
;print, ' LASCO    '
;print, info.C_list_of_datetime_strings
;;print, info.C_list_of_full_time_strings
;;print, ' LASCO ',info.C_list_of_image_names
;print, '*****************************************************'
;print, ' STEREO A '
;print, info.R_list_of_datetime_strings
;;print, info.R_list_of_full_time_strings
;;print, ' STEREO A ',info.R_list_of_image_names
;print, '*****************************************************'

   
;   jumphere :




END







pro swpc_cat_find_images_by_date, ddate,info,list_of_image_names,list_of_image_data, $
                         list_of_datetime_strings, list_of_image_exposure_times, $
                         list_of_image_offsets,list_of_image_scaling_factors,swpc_cat_progressbar, $
                         progress_base_percentage, list_of_HEEQ_coords, list_of_pixel_scales, $
                         list_of_rsuns, list_of_Sun_satellite_distances
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
;info.telescope_array as in info.image_in_folder_array
res = STRCMP(info.which_telescope,info.telescope_array)
ele = where(res,count)

if count eq 1 then begin
  telescope_folder = info.image_in_folder_array[ele] 
endif else begin
  telescope_folder = ''
  msg1 = 'ERROR: folder for telescope ' + info.which_telescope +' is not defined.'
  msg2 = 'There should be matching entries for each telescope '
  msg3 = 'in lines "telescopes" and "image_in_folders" in file' 
  msg4 = info.input_file + '.'
  ok = DIALOG_MESSAGE([msg1,msg2,msg3,msg4],/ERROR,/center)
  EXIT   
endelse  

;assemble the location to look for files:
image_location = info.image_in_root + info.sep + date_folder + info.sep + telescope_folder

result = file_test(image_location,/directory)
if result eq 0 then begin
  msg1 = 'Warning: The subdir ' + image_location + ' does not exist;'
  msg2 = 'this means there are no images available for telescope ' + info.which_telescope
  msg3 = 'on date ' + date_folder  + ' [yet].'
;  ok = DIALOG_MESSAGE([msg1,msg2,msg3],/information,/center)
  RETURN
  
endif

   swpc_cat_progressbar -> Update, progress_base_percentage, Text='Searching for image files'

;get the files from this folder:
;  ( do not get the STEREO polarization images, for now...)
if (STRCMP( info.which_telescope, 'STEREO', 6 )) then begin
   names_of_all_files_in_folder = file_search(image_location,'*d7*.fts',count = number_of_files_in_folder)
;   print, 'names_of_all_files_in_folder ', names_of_all_files_in_folder
endif else begin
   names_of_all_files_in_folder = file_search(image_location,'*.fts',count = number_of_files_in_folder)
endelse
if number_of_files_in_folder lt 1 then begin
   msg1 = 'Warning...No FITS files in this folder: ' + image_location
;   ok = DIALOG_MESSAGE(msg1,/information,/center)
   RETURN
   
endif else begin
   progress_base_percentage = progress_base_percentage + 10.
   
   swpc_cat_progressbar -> Update, progress_base_percentage, Text='Searching for image files'
   ;Keep only the files between the user-selected 
   ;start and end date interval, inclusive
   
;   help,names_of_all_files_in_folder
   datetime_strings = STRMID(FILE_BASENAME(names_of_all_files_in_folder,'.fts'),0,15)
   
   
;   print, 'number_of_files_in_folder ',number_of_files_in_folder
;   help, datetime_strings
;   print, 'datetime strings ', datetime_strings
   
   i=0
   while (i lt number_of_files_in_folder) do begin
   
      data_year = FIX(STRMID(datetime_strings[i],0,4))
      data_month = FIX(STRMID(datetime_strings[i],4,2))
      data_day = FIX(STRMID(datetime_strings[i],6,2))
      data_hour = FIX(STRMID(datetime_strings[i],9,2))
      data_minute = FIX(STRMID(datetime_strings[i],11,2))
      data_second = FIX(STRMID(datetime_strings[i],13,2))
      data_dt = JULDAY(data_month,data_day,data_year,data_hour,data_minute,data_second)
      
      progress_percent = progress_base_percentage + (30.*(float(i+1)/float(number_of_files_in_folder)))
      
      swpc_cat_progressbar -> Update, progress_percent, Text='Loading Image file : ' + file_basename(names_of_all_files_in_folder[i])    ; + STRMID(datetime_strings[i],0,4)+ '-' + STRMID(datetime_strings[i],4,2)+'-' + STRMID(datetime_strings[i],6,2)

      ;Add to the list of images that are between start and end dates, inclusive
      if ((data_dt GE info.start_julian) && (data_dt LE info.end_julian)) then begin
;         print,'loading ',i
         good_image = swpc_cat_check_image_and_load_data(info.which_telescope, telescope_folder, names_of_all_files_in_folder[i], $
                                                scaling_factor, rotation, center_of_sunX, center_of_sunY, $
                                                rsun, pixel_scale, exposure_time, offset, datetime_string, $ 
                                                HEEQ_coords, Sun_satellite_distance, image_data)         
;         print,'good_image=',good_image
         if good_image then begin
            list_of_image_names.Add,names_of_all_files_in_folder[i],/NO_COPY
            list_of_image_data.Add,image_data,/NO_COPY
            list_of_datetime_strings.add,datetime_string,/no_copy
            list_of_image_exposure_times.Add,exposure_time,/NO_COPY
            list_of_image_offsets.Add,offset,/NO_COPY
            list_of_image_scaling_factors.Add,scaling_factor,/NO_COPY
            list_of_HEEQ_coords.Add,HEEQ_coords,/NO_COPY
            list_of_pixel_scales.Add,pixel_scale,/NO_COPY
            list_of_rsuns.Add,rsun,/NO_COPY
            list_of_Sun_satellite_distances.Add,Sun_satellite_distance,/NO_COPY
      endif
      endif
      i=i+1
   endwhile
   
;   print, 'number of files ', i
   
   progress_base_percentage = progress_percent
     
endelse

END




pro swpc_cat_check_dates,date_array,start_julian,end_julian,max_interval_in_days, $
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
valid = swpc_cat_validate_date(start_year,start_month,start_day,start_hour,start_minute)
if (valid eq 0) then begin
 dates_are_fine = 0
 problem_string = 'Start date is not a valid date'
endif else begin

   ;validate end date:
   valid = swpc_cat_validate_date(end_year,end_month,end_day,end_hour,end_minute)
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



function swpc_cat_validate_date,year,month,day,hour,minute

;20110628 anewman - started

compile_opt idl2

if year lt 1000 or year gt 2999 then return,0
if month lt 1 or month gt 12 then return,0
mlen = swpc_cat_MonthLen(month,year4=year)
if day lt 1 or day gt mlen then return,0
if hour lt 0 or hour gt 23 then return,0
if minute lt 0 or minute gt 59 then return,0

return,1
end



function swpc_cat_MonthLen,month,year4=year4
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

;  if not keyword_set(year4) then begin
;    year4 = SystimeX(/year,/str)
;  endif
  year4 = fix(year4)
  if size(month,/type) eq 7 then begin
    month = strmid(month,0,3)
    mstr = ["xx","Jan","Feb","Mar","Apr","May","Jun", $
        "Jul","Aug","Sep","Oct","Nov","Dec"]
    mnum = where(month eq mstr,cnt)
    month = mnum
  endif
  mlen = [0,31,28,31,30,31,30,31,31,30,31,30,31]
  if swpc_cat_LeapYear(year4) then mlen[2] = 29
  return,mlen[month]
end



function swpc_cat_LeapYear,Year
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




pro swpc_cat_get_fits_header_data, fits_header, which_telescope, telescope_folder, scaling_factor, rotation, $
                          sun_X, sun_Y, rsun, pixel_scale, exposure_time, offset, $ 
                          datetime_string, HEEQ_coords, Sun_satellite_distance, $
                          OBS_ID_divide_by_2_factor,corrupt_file

;background:  1 => this is a background image
;             0 => this is a normal image

compile_opt idl2

   HEEQ_coords = fltarr(2)
   OBS_ID_divide_by_2_factor = 1.


;The FITS header for LASCO and STEREO differ; STEREO does not have a TIME-OBS item
if (STRMID(which_telescope,0,5) eq 'LASCO') then begin

   date_string = swpc_cat_FXPAR(fits_header,'DATE-OBS')
   ;use dashes, not slashes, between YYYY MM DD
   date_string = strjoin(strsplit(date_string,'/',/EXTRACT),'-')
   time_string = swpc_cat_FXPAR(fits_header,'TIME-OBS')
   time_string = strmid(time_string,0,8) ;strip off msec
   datetime_string = date_string + 'T' + time_string ;assemble datetime string
   naxis1 = fix(swpc_cat_sxpar(fits_header, 'NAXIS1'))
   ;thanks to Doug Biesecker for the following:
   rotation = float(swpc_cat_sxpar(fits_header, 'CROTA1')); specified in degrees CW
   pixel_scale = float(swpc_cat_sxpar(fits_header, 'CDELT1')) ; assumes CDELT1 and CDELT2 are identical
   binning1 = fix(swpc_cat_sxpar(fits_header, 'SUMROW')); 0 signifies no binning
   binning2 = fix(swpc_cat_sxpar(fits_header, 'LEBXSUM')); 1 signifies no binning
   if binning2 mod 2 eq 1 then binning2 = binning2 - 1   
   binning = max([binning1, binning2]); binning isn't needed so long as CDELT1 is used
   sun_x = float(swpc_cat_sxpar(fits_header, 'CRPIX1')); confirm this is x-value
   sun_y = float(swpc_cat_sxpar(fits_header, 'CRPIX2')); confirm this is y-value
   rsun = 941.; solar radius from L1 in arcseconds, a default value 
   exposure_time = float(swpc_cat_FXPAR(fits_header,'EXPTIME'))
   offset = float(swpc_cat_FXPAR(fits_header,'OFFSET'))
   HEEQ_coords[0] = 0.0
   HEEQ_coords[1] = 0.0
   Sun_satellite_distance = 215.0
   corrupt_file = 0
;   
; At SWPC, when we pull down LASCO C2 and C3 files we rebin them 
; to be 512 * 512, from the original 1024 * 1024.  We do this because
; we only need them to be 512 and the rebin makes them 1/4 the size and
; much faster to load.  Our rebinned files have a new fits header keyword
; 'REBINVER' which is set to 1.  The following code tests for REBINVER and
; sets sun position and pixel_scale appropriately.  Original 1024 * 1024
; files do not have the REBINVER keyword and thus the resulting swpc_rebin_param 
; is zero and the code skips over the following 'if block'.
; Original 1024 files work fine - they just load much more slowly.
;     
   swpc_rebin_param = fix(swpc_cat_sxpar(fits_header, 'REBINVER'))
   
   if swpc_rebin_param eq 1 then begin
     SWPC_rebin_version1_factor = 0.5
     sun_x = sun_x * SWPC_rebin_version1_factor
     sun_y = sun_y * SWPC_rebin_version1_factor
     pixel_scale = pixel_scale / SWPC_rebin_version1_factor
   endif else begin
     SWPC_rebin_version1_factor = 1.0
   endelse


endif else if (STRMID(which_telescope,0,6) eq 'STEREO') then begin

   date_string = swpc_cat_FXPAR(fits_header,'DATE-CMD')
   ;time is stored in DATE-CMD, for STEREO FITS files
   ;want this format: YYYY-MM-DDTHH:MM:SS
   time_string = strmid(date_string,11,8) ;strip off msec
   date_string = strmid(date_string,0,10) ;strip off time 
   datetime_string = date_string + 'T' + time_string ;reassemble datetime string
   naxis1 = fix(swpc_cat_sxpar(fits_header, 'NAXIS1'))
   ;thanks to Doug Biesecker for the following:
   rotation = float(swpc_cat_sxpar(fits_header, 'CROTA'))
   pixel_scale = float(swpc_cat_sxpar(fits_header, 'CDELT1'))
   binning = fix(swpc_cat_sxpar(fits_header, 'SUMMED')); dimension = original/(2^(SUMMED-1)); 1 signifies no binning
   ; NOTE: binning isn't needed since CDELT1 is updated properly
   sun_x = float(swpc_cat_sxpar(fits_header, 'CRPIX1'))
   sun_y = float(swpc_cat_sxpar(fits_header, 'CRPIX2'))
   rsun = float(swpc_cat_sxpar(fits_header, 'RSUN')); The fits files contain the solar radius in arcsec
;   print, which_telescope, ' ',pixel_scale, sun_x , sun_y , rsun
   exposure_time = float(swpc_cat_fxpar(fits_header,'EXPTIME'))
   offset = float(swpc_cat_fxpar(fits_header,'OFFSET'))
   HEEQ_X = float(swpc_cat_fxpar(fits_header,'HEQX_OBS'))
   HEEQ_Y = float(swpc_cat_fxpar(fits_header,'HEQY_OBS'))
   HEEQ_Z = float(swpc_cat_fxpar(fits_header,'HEQZ_OBS'))
   HEEQ_coords[0] = atan(HEEQ_Y,HEEQ_X) * 180. / !pi
   HEEQ_coords[1] = atan(HEEQ_Z,sqrt(HEEQ_X^2 + HEEQ_Y^2)) * 180. / !pi
   OBS_ID = fix(swpc_cat_fxpar(fits_header,'OBS_ID'))
   
   biassdev = float(swpc_cat_fxpar(fits_header,'BIASSDEV'))
   ip_time = fix(swpc_cat_fxpar(fits_header,'IP_TIME'))
   NMISSING = fix(swpc_cat_fxpar(fits_header,'NMISSING'))
   misslist = swpc_cat_fxpar(fits_header,'MISSLIST')
   datazer = fix(swpc_cat_fxpar(fits_header,'DATAZER'))
   evcount = swpc_cat_fxpar(fits_header,'EVCOUNT')
   datamin = float(swpc_cat_fxpar(fits_header,'DATAMIN'))
   
   corrupt_file = nmissing
   
;   print, which_telescope, ' ',date_string,' ',time_string, ip_time,biassdev,nmissing,misslist, datazer,evcount,datamin
   
   if OBS_ID eq 1694 then OBS_ID_divide_by_2_factor = 2.
   Sun_satellite_distance = sqrt(HEEQ_X^2 + HEEQ_Y^2 + HEEQ_Z^2) / 695500000.d
;   print, ' Sun_satellite_distance ', Sun_satellite_distance
;   print, fits_header

        
endif

;The input FITS files we get from NOC do not have 512x512 images, but
; we convert each image to 512x512 (see rebin command elsewhere in this
; program). We need to keep track of the
; scaling factor for the image conversion,
; for use in rotation and leading edge calculations


scaling_factor = 1.
if (naxis1 eq 256) then scaling_factor = 2. $         ;256x256 to 512x512: multiply by 2
else if (naxis1 eq 1024) then scaling_factor = 0.5   ;1024x1024 to 512x512: divide by 2 (i.e. multiply by 0.5)


END







pro swpc_cat_rotate_image_to_solar_north,center_of_sunX,center_of_sunY,rotation,scaling_factor,image_data

   ;rotate the image to Solar North, and place Sun center in the center of the image using IDL's ROT function
   ;First, keep image 512x512 (as set in a rebin cmd in calling routine, or sometime previously)
   ; need to adjust sun center (which specifies the sun center in the original size image, e.g.
   ;256x256, or 1024x1024) correctly in rotated 512x512 image, by applying the scale factor:
   
   compile_opt idl2
 
   image_data = rot(image_data, -1.*rotation, 1.0, $
                    center_of_sunX*scaling_factor,  $
                    center_of_sunY*scaling_factor, /interp)
  

END







function swpc_cat_check_image_and_load_data, which_telescope, telescope_folder, fits_file, scaling_factor, rotation, $
                                    center_of_sunX, center_of_sunY, rsun, pixel_scale, exposure_time, $
                                    offset, datetime_string, HEEQ_coords, Sun_satellite_distance, image_data 

   ;Identify bad input data (e.g., input FITS files that are corrupted)
   ;If the file is corrupted, flag it, by returning 0
   ;
   ;This is "Trust-but-verify code; bad images should have been deleted
   ;upstream of this program, so that this program would never see them.
   ;But just in case that did not happen...we do some checking here.
   
   compile_opt idl2

   im_data = swpc_cat_READFITS(fits_file,fits_header,/silent)
   ;if im_data eq -1  i.e, if error encountered with reading the FITS file
   if size(im_data,/n_dimensions) eq 0 then begin
      If im_data eq -1 then begin
         msg1 = 'There was a problem reading FITS file ' + fits_file + '.'
         msg2 = 'Excluding this image from the list of images to view...'
;         ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
         return,0
      endif
   endif

   ;check for correct size of the image; the 2 dimensions 
   ;must be the same size, or the file is corrupted
   ;Example bad images of this type, used in testing:
   ;   X:\Images\20110819\LASCO_C2\20110819_143322_C2.fts
   ;   X:\Images\20110819\LASCO_C2\20110819_144128_C2.fts
   
   naxis1 = fix(swpc_cat_sxpar(fits_header, 'NAXIS1'))
   naxis2 = fix(swpc_cat_sxpar(fits_header, 'NAXIS2'))
   
   if naxis1 ne naxis2 then begin
      dim2 = '(' + strtrim(string(naxis1),2) + 'x' + strtrim(string(naxis2),2) + ')'
      msg1 = 'There is a dimension problem ' + dim2 + ' with image ' + fits_file
      msg2 = 'Excluding this image from the list of images to view...'
;      ok = DIALOG_MESSAGE([msg1,msg2],/information,/center)
      return,0
   endif else begin
      swpc_cat_get_fits_header_data, fits_header, which_telescope, telescope_folder, scaling_factor, rotation, $
      center_of_sunX, center_of_sunY, rsun, pixel_scale, exposure_time, offset, $
      datetime_string, HEEQ_coords, Sun_satellite_distance, OBS_ID_divide_by_2_factor,corrupt_file
      
      image_data = rebin(im_data,512,512)
      image_data = float(image_data)
      swpc_cat_rotate_image_to_solar_north,center_of_sunX,center_of_sunY,rotation,scaling_factor,image_data
      image_data = image_data/OBS_ID_divide_by_2_factor - offset
      image_data = image_data * 20. / exposure_time
   endelse
   
   ; sometimes we have a corrupt 'noisy' STEREO file but the nmissing header variable does not reflect this.
   ; so we also need to check to see if the data is noisy and discard the image if it is.
   ; code provided by Bill Thompson.... william.t.thompson@nasa.gov
   
   if which_telescope eq 'STEREO A COR2' or which_telescope eq 'STEREO B COR2' then begin
   a = im_data
   stat = fltarr(4)
   sz = size(a)
   stat[0] = stddev(median(a[1:10, 1:10],3))
   stat[1] = stddev(median(a[1:10, sz[2]-11:sz[2]-2],3))
   stat[2] = stddev(median(a[sz[1]-11:sz[1]-2, 1:10],3))
   stat[3] = stddev(median(a[sz[1]-11:sz[1]-2,sz[2]-11:sz[2]-2],3))
;   if keyword_set(beacon) then mstat=10 else mstat=5
   mstat = 10
   if median(stat) gt mstat then corrupt_file = 1
   endif
   
;   print, 'return code ', 1 - corrupt_file
   
   return,1 - corrupt_file  ;  return a 1 if a good image, otherwise a zero
 
END



PRO swpc_cat_Cleanup, tlb

;*****************************************************************
;
;  cleanup - cleans up - Frees all created objects.
;
;****************************************************************
compile_opt idl2

Widget_Control, tlb, Get_UValue=info2
IF N_Elements(info) NE 0 THEN begin
help,/mem
heap_free,info
help,/mem
heap_gc
help,/mem
ENDIF
END





pro swpc_cat_top_level_base_events,event

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

Widget_Control, event.top, Get_UValue=info, /No_Copy

info.tlb_position[0] = event.x
info.tlb_position[1] = event.y

;openw,lun,info.swpc_cat_preferences_file,/get_lun
;printf,lun,info.tlb_position
;free_lun,lun

;directory_result = file_test(info.preferences_folder,/directory)
;if directory_result ne 1 then begin
;file_mkdir, info.preferences_folder
;endif
;
;openw,lun,info.preferences_folder+info.sep+'last_window_position.txt',/get_lun
;printf,lun,info.tlb_position
;free_lun,lun


Widget_Control, event.top, Set_UValue=info, /No_Copy

end












pro swpc_cat_plot_Window_button_click_events, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

drawTypes = ['PRESS', 'RELEASE', 'MOTION', 'SCROLL', 'EXPOSE']
thisEvent = drawTypes[event.type]

CASE thisEvent OF

   'EXPOSE':  ; Nothing required except to draw the view.
   'PRESS': BEGIN
       Widget_Control, event.id, Draw_Motion_Events=1 ; Motion events ON.

       END
   'RELEASE': BEGIN
       
       Widget_Control, event.id, Draw_Motion_Events=0 ; Motion events OFF.
       
       
       if event.release eq 4 then begin
       
    if info.allow_manual_data_point eq 1 then begin
      contextBase = WIDGET_INFO(event.TOP, FIND_BY_UNAME = 'plot_window_contextBase')
      WIDGET_DISPLAYCONTEXTMENU, event.ID, event.X, event.Y, contextBase
    endif
       
       endif
           
       END
   'MOTION': BEGIN ; Trackball events
  

       END
   ELSE:


ENDCASE

Widget_Control, event.top, Set_UValue=info, /No_Copy
END











pro swpc_cat_Window_button_click_events, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

;case event.id of 
;    info.draw_L : print, 'L'
;    info.draw_C : print, 'C'
;    info.draw_R : print, 'R'
;endcase

drawTypes = ['PRESS', 'RELEASE', 'MOTION', 'SCROLL', 'EXPOSE']
thisEvent = drawTypes[event.type]


CASE thisEvent OF

   'EXPOSE':  ; Nothing required except to draw the view.
   'PRESS': BEGIN
       Widget_Control, event.id, Draw_Motion_Events=1 ; Motion events ON.


if info.debug_mode eq 1 then print, 'event.press ', event.press
if event.press eq 4 then begin

info.pressed_the_right_button = 1 

endif else begin

info.pressed_the_right_button = 0

if info.show_image_line_plot eq 1 then begin
       
case event.id of 
    info.draw_L : Begin
    
       info.L_latest_click_X = event.x
       info.L_latest_click_Y = event.y
       info.L_previous_click_X = event.x
       info.L_previous_click_Y = event.y
       info.L_click_and_drag  = 1

       near_clock_angle = 0
       near_leading_edge = 0   
       myLoc = [event.x, event.y]    
       oObjArr = info.L_Window->Select(info.L_image_view, myLoc, dimensions=[20,20])
       nSel = N_ELEMENTS(oObjArr)
       if nsel gt 1 then begin
       iobject = intarr(nSel)
       
       FOR i=0, nSel-1 DO BEGIN
          oObjArr[i]->GetProperty, NAME=name
          if name eq 'leading_edge' then begin
            near_leading_edge = 1
          endif
          if name eq 'L_clock_angle' then begin
            near_clock_angle = 1
          endif
       ENDFOR
       endif
       
       if near_leading_edge eq 1 then begin
       info.L_the_action = 2
       endif else begin
       if near_clock_angle eq 1 then begin
       info.L_the_action = 1
       endif else begin
       info.L_the_action = 0
       endelse
       endelse

       END
       

    info.draw_C : Begin
    
       info.C_latest_click_X = event.x
       info.C_latest_click_Y = event.y
       info.C_previous_click_X = event.x
       info.C_previous_click_Y = event.y
       info.C_click_and_drag  = 1


       near_clock_angle = 0
       near_leading_edge = 0   
       myLoc = [event.x, event.y]    
;       oObjArr = info.C_Window->Select(info.C_image_view, myLoc, dimensions=[20,20])
; ECDP - we need a restricted window around the leading edge marker to ensure that we don't get errors due to different mouse click points.
       oObjArr = info.C_Window->Select(info.C_image_view, myLoc, dimensions=[4,4])
       nSel = N_ELEMENTS(oObjArr)
       if nsel gt 1 then begin
       iobject = intarr(nSel)
       
       FOR i=0, nSel-1 DO BEGIN
          oObjArr[i]->GetProperty, NAME=name
          if name eq 'leading_edge' then begin
            near_leading_edge = 1
;            print, 'Near leading edge'
          endif
          if name eq 'C_clock_angle' then begin
            near_clock_angle = 1
;            print, 'Near clock angle'  
          endif
       ENDFOR
       endif
       
       if near_leading_edge eq 1 then begin
       info.C_the_action = 2
       endif else begin
; only allowed to move the clock angle if it's the first attempt. Otherwise we want to keep it fixed, for consistency        
       if ((near_clock_angle eq 1) and (n_elements(info.POS_Julian) eq 0))then begin
       info.C_the_action = 1
       endif else begin
       info.C_the_action = 0
       endelse
       endelse

       END
    
    info.draw_R : Begin
    
       info.R_latest_click_X = event.x
       info.R_latest_click_Y = event.y
       info.R_previous_click_X = event.x
       info.R_previous_click_Y = event.y
       info.R_click_and_drag  = 1

       near_clock_angle = 0
       near_leading_edge = 0   
       myLoc = [event.x, event.y]    
       oObjArr = info.R_Window->Select(info.R_image_view, myLoc, dimensions=[20,20])
       nSel = N_ELEMENTS(oObjArr)
       if nsel gt 1 then begin
       iobject = intarr(nSel)
       
       FOR i=0, nSel-1 DO BEGIN
          oObjArr[i]->GetProperty, NAME=name
          if name eq 'leading_edge' then begin
            near_leading_edge = 1
          endif
          if name eq 'R_clock_angle' then begin
            near_clock_angle = 1
          endif
       ENDFOR
       endif
       
       if near_leading_edge eq 1 then begin
       info.R_the_action = 2
       endif else begin
       if(near_clock_angle eq 1) then begin
       info.R_the_action = 1
       endif else begin
       info.R_the_action = 0
       endelse
       endelse

       END

endcase

endif  ;  if info.show_image_line_plot eq 1
       
endelse

       END
       
   'RELEASE': BEGIN
       
       Widget_Control, event.id, Draw_Motion_Events=0 ; Motion events OFF.
       
if event.release eq 4 then begin
       
case event.id of 
    info.draw_L : Begin
    contextBase = WIDGET_INFO(event.TOP, FIND_BY_UNAME = 'L_drawContext')
    WIDGET_DISPLAYCONTEXTMENU, event.ID, event.X, event.Y, contextBase
    END

    info.draw_C : Begin
    contextBase = WIDGET_INFO(event.TOP, FIND_BY_UNAME = 'C_drawContext')
    WIDGET_DISPLAYCONTEXTMENU, event.ID, event.X, event.Y, contextBase
    
    END
    info.draw_R : Begin
    contextBase = WIDGET_INFO(event.TOP, FIND_BY_UNAME = 'R_drawContext')
    WIDGET_DISPLAYCONTEXTMENU, event.ID, event.X, event.Y, contextBase
    END
endcase
       
endif
       

           
       END
       
   'MOTION': BEGIN ; motion events

if (info.pressed_the_right_button eq 1) then begin

endif else begin

if info.show_image_line_plot eq 1 then begin
   
case event.id of 
    info.draw_L : Begin
    
if info.L_the_action eq 1 then begin
    
     ;Save the current clock angle:
     old_degrees = info.L_clock_angle_degrees
     
     ; The user has rotated (or at least touched) the clock angle model....
     
     ; so calculate the clock angle in degrees
     L_previous_click_X = info.L_previous_click_X
     L_previous_click_Y = info.L_previous_click_Y 
;     swpc_cat_calc_clock_angle, event.x,event.y, L_previous_click_X, L_previous_click_Y, $
;       info.L_clock_angle_model, L_clock_angle

; ECDP we don't need the clock for L window
    swpc_cat_calc_clock_angle, event.x,event.y
       
     info.L_previous_click_X = L_previous_click_X 
     info.L_previous_click_Y = L_previous_click_Y    
     info.L_rotate_x = event.x
     info.L_rotate_y = event.y
     
     info.L_clock_angle_degrees = L_clock_angle
     
     if info.images_are_loaded eq 1 then swpc_cat_replot_image_line_plot, $
      info.L_clock_angle_degrees, info.L_coronagraph_image_object, info.L_image_lineplot, $
      info.position_image_lineplot, info.L_cme_outline
      
endif

    END

    info.draw_C : Begin
    
;print, 'C_the_action', info.C_the_action
;stop

; if near the clock arm
if info.C_the_action eq 1 then begin

     ;Save the current clock angle:
     old_degrees = info.C_clock_angle_degrees
     
     ; The user has rotated (or at least touched) the clock angle model....
     
     ; so calculate the clock angle in degrees
;     C_previous_click_X = info.C_previous_click_X
;     C_previous_click_Y = info.C_previous_click_Y 
;     swpc_cat_calc_clock_angle, event.x,event.y, C_previous_click_X, C_previous_click_Y, info.C_clock_angle_model, C_clock_angle
     swpc_cat_calc_clock_angle, info, event.x,event.y
;     info.C_previous_click_X = C_previous_click_X 
;     info.C_previous_click_Y = C_previous_click_Y    
     info.C_rotate_x = event.x
     info.C_rotate_y = event.y
     
;     info.C_clock_angle_degrees = C_clock_angle

     if info.images_are_loaded eq 1 then swpc_cat_replot_image_line_plot, $
      info.C_clock_angle_degrees, info.C_coronagraph_image_object, info.C_image_lineplot, $
      info.position_image_lineplot, info.C_cme_outline
 
; print, 'motion1', event.x, event.y

endif

; ECDP
;  print, 'Here1'
;  print, 'motion ', event.x, event.y
;  print, 'C_the_action', info.C_the_action
;  stop

; if near the leading edge
if info.C_the_action eq 2 then begin

     ;Save the current leading edge value:
     C_old_leading_edge = info.C_leading_edge_in_rs
     
; Translate the leading edge_model
         C_prev_x = info.C_previous_click_X
         C_prev_y = info.C_previous_click_Y
;         C_prev_x = info.C_now_x
;         C_prev_y = info.C_now_y
         C_prev_dist = sqrt((C_prev_x - 255)^2 + (C_prev_y - 255)^2)
         C_now_x = event.x
         C_now_y = event.y
         C_now_dist = sqrt((C_now_x - 255)^2 + (C_now_y - 255)^2)
         info.c_now_dist = C_now_dist
         info.C_leading_edge_model -> translate, 0.0 , C_now_dist - C_prev_dist , 0.0
         info.C_leading_edge_model -> GetProperty, Transform = Transform
;           info2.Ellipse_center_location[0] = transform[3,0]
;           info2.Ellipse_center_location[1] = transform[3,1]
         info.C_previous_click_X = event.x
         info.C_previous_click_Y = event.y
        
;         if (C_prev_y ne info.C_previous_click_Y) then print, 'Error'
         ;20110715 ahn add this calc here, plus add leading_edge_to_rs to the info2 structure
         ;The pixel scale is related to the original image, before we converted the image
         ;to 512x512. So we must apply a scale_factor to the pixel scale to make the
         ;leading edge calculation correct. The pixel scale is the number of arc seconds
         ;per pixel. So, for example, if the original image is twice as large as the 512x512
         ;current image then there will be twice as many arc seconds per pixel, ( i.e. pixel_scale* 2 or
         ;pixel scale / 0.5 ) in the current image.
                
        swpc_cat_calculate_c2_c3_lead_edge, info
         ;Display new value of the location of the CME leading edge
         ;in lower rt corner of image:
         C_lead_edge_string_integer = strcompress(string(info.C_leading_edge_in_rs),/remove_all)
         if info.C_leading_edge_in_rs lt 10 then C_lead_edge_string_integer = ' ' + C_lead_edge_string_integer
;         info.C_leading_edge_txt_object->SetProperty,strings=C_lead_edge_string_integer

        ; if user has changed the leading edge, and if
        ; the CME file name has already been set,
        ; then reset/erase it:
;        if C_old_leading_edge ne info.C_leading_edge_in_rs and $
;           strlen(info.cme_name) ne 0 then begin
;           reset_CME_name, info2
           
           ;check the Representative Image; if set, then reset/erase it
;           if strlen(info.cme_image_name) ne 0 then begin
;              reset_CME_image,info
;           endif           
                     
;        endif
        ; ECDP
;        if (info.which_mode eq 1) then begin
;          info.C_cme_outline -> setProperty, hide = 1
;        endif else begin
;          info.C_cme_outline -> setProperty, hide = 0
;        endelse
; ECDP  
;	print, 'motion2', event.x, event.y
;	Widget_Control, info.draw_C, Draw_Motion_Events=1
        info.C_Window->Draw, info.C_both_views

endif     ; leading edge action      

;  ECDP
;	info.C_Window->Draw, info.C_both_views
    
    END
    info.draw_R : Begin
    
if info.R_the_action eq 1 then begin

     ;Save the current clock angle:
     old_degrees = info.R_clock_angle_degrees
     
     ; The user has rotated (or at least touched) the clock angle model....
     
     ; so calculate the clock angle in degrees
     R_previous_click_X = info.R_previous_click_X
     R_previous_click_Y = info.R_previous_click_Y 
     swpc_cat_calc_clock_angle, event.x,event.y, R_previous_click_X, R_previous_click_Y, $
       info.R_clock_angle_model, R_clock_angle
       
     info.R_previous_click_X = R_previous_click_X 
     info.R_previous_click_Y = R_previous_click_Y    
     info.R_rotate_x = event.x
     info.R_rotate_y = event.y
     
     info.R_clock_angle_degrees = R_clock_angle
     
     if info.images_are_loaded eq 1 then swpc_cat_replot_image_line_plot, $
      info.R_clock_angle_degrees, info.R_coronagraph_image_object, info.R_image_lineplot, $
      info.position_image_lineplot, info.R_cme_outline        

endif

    END
    
endcase

endif   ; if info.show_image_line_plot eq 1
   
endelse ; event.motion

       END
       


   ELSE:


ENDCASE


case event.id of 
    info.draw_L : info.L_Window->Draw, info.L_both_views
    info.draw_C : info.C_Window->Draw, info.C_both_views
    info.draw_R : info.R_Window->Draw, info.R_both_views
endcase

Widget_Control, event.top, Set_UValue=info, /No_Copy
END


;pro swpc_cat_replot_image_line_plot, clock_angle, coronagraph_image_object, image_lineplot, position_image_lineplot, cme_outline
pro swpc_cat_replot_image_line_plot, C_clock_angle, coronagraph_image_object, image_lineplot, position_image_lineplot, cme_outline
;     xp = 255 - (255*sin(clock_angle * !dtor))  
;     yp = 255 + (255*cos(clock_angle * !dtor))
     xp = 255 - (255*sin(C_clock_angle * !dtor))  
     yp = 255 + (255*cos(C_clock_angle * !dtor))
     coronagraph_image_object -> getproperty, data = image_data
     xvals = round(((findgen(255)/254.) * (xp - 255.)) + 255.)
     yvals = round(((findgen(255)/254.) * (yp - 255.)) + 255.)
     line_Y = image_data[xvals,yvals]
     image_lineplot->SetProperty, datax = findgen(255) , datay = line_Y
     maxy = max(line_Y,min=minY)
     yrange = [minY,maxY]
     image_xs = swpc_cat_FSC_Normalize([0.,254.], Position=[position_image_lineplot[0], position_image_lineplot[2]])
     image_ys = swpc_cat_FSC_Normalize(yrange, Position=[position_image_lineplot[1], position_image_lineplot[3]])
     image_lineplot->SetProperty, XCoord_Conv=image_xs, YCoord_Conv=image_ys
     
     
;     Result = cme_outline.ContainsPoints( xvals, yvals )
;     vert_colors = intarr(3,255)
;     for i = 0 , 254 do begin
;     if result[i] eq 1 then begin
;     vert_colors[0,i] = 255
;     endif else begin
;     vert_colors[1,i] = 255
;     endelse
;     endfor
     
     Result = cme_outline.ContainsPoints( xvals, yvals )
     exterior_locations = where(result eq 0)
     vert_colors = intarr(3,255)
     vert_colors[0,*] = 255
     vert_colors[0,exterior_locations] = 0
     vert_colors[1,exterior_locations] = 255
;     vert_colors = intarr(3,255)
;     for i = 0 , 254 do begin
;     if result[i] eq 1 then begin
;     vert_colors[0,i] = 255
;     endif else begin
;     vert_colors[1,i] = 255
;     endelse
;     endfor
     
     image_lineplot->SetProperty, vert_colors = vert_colors
     
     
end







;pro swpc_cat_calc_clock_angle, x2, y2, C_previous_click_X, C_previous_click_Y, C_clock_angle_model, the_clock_angle
PRO swpc_cat_calc_clock_angle, info, x2, y2
   compile_opt idl2

;   if info.full_halo then $
;      ;When the CME is full halo,
;      ; the degrees are always 360...
;      set_clock_angle,info,360 $
;   else  begin
      ;Get the clock angle, based on previous location
      ;and current location of the red clock angle model
      ;line
      x1 = info.C_previous_click_X
      y1 = info.C_previous_click_Y
      xc = 255
      yc = 255
      theta2 = atan((y2 - yc) , (x2 - xc))
      theta1 = atan((y1 - yc) , (x1 - xc))
      rotation_angle_degrees = (theta2 - theta1) * 180. / !pi
      info.C_previous_click_X = x2
      info.C_previous_click_Y = y2
         
      info.C_clock_angle_model -> rotate, [0,0,1], rotation_angle_degrees, /premultiply
      info.C_clock_angle_model -> getProperty, transform=transform
      
      val1=round(( atan(transform[1,0], transform[0,0]) * 180. / !pi ))
      if val1 gt 0 then begin
         the_clock_angle = 360 - val1
      endif else begin
         the_clock_angle = 0 - val1
      endelse
; set the clock angle
      swpc_cat_set_clock_angle,info,the_clock_angle
  
;  endelse

END

PRO swpc_cat_set_clock_angle,info,degrees

   compile_opt idl2
;   print, 'Degrees =', degrees
;   print, info.C_clock_angle_degrees
   ;save the clock angle, in degrees, and display it on the image
   info.C_clock_angle_degrees = degrees
   angle_string = strcompress(string(info.C_clock_angle_degrees))
   if info.C_clock_angle_degrees lt 100 then angle_string = ' ' + angle_string
   if info.C_clock_angle_degrees lt 10 then angle_string = ' ' + angle_string
   ;Displays new value of the clock angle degrees near the outer part of the rotating line
   info.C_clock_angle_string->SetProperty,strings=angle_string
   ;Display new value of the clock angle degrees in lower left corner of image:
;    info.C_clock_angle_deg_object->SetProperty,strings=angle_string

END



pro swpc_cat_images_timeline_window_button_click_events, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

drawTypes = ['PRESS', 'RELEASE', 'MOTION', 'SCROLL', 'EXPOSE']
thisEvent = drawTypes[event.type]

;print,' shift click ',event.modifiers
shift_click = event.modifiers

;print, event.x , event.y , float(event.y)/float(info.ysize_timeline)

if (info.end_julian gt 0.) then begin

normalized_x = float(event.x) / info.xsize_timeline
this_julian = (((normalized_x - info.position_timeline[0]) / (info.position_timeline[2] - info.position_timeline[0])) * (info.end_julian - info.start_julian)) + info.start_julian

CALDAT, this_julian, this_Month, this_day, this_year, this_hour, this_minute, this_Second

this_year_str = string(this_year)
this_month_str = string(this_month)
if this_month lt 10 then this_month_str = '0' + this_month_str
this_day_str = string(this_day)
if this_day lt 10 then this_day_str = '0' + this_day_str
this_hour_str = string(this_hour)
if this_hour lt 10 then this_hour_str = '0' + this_hour_str
this_minute_str = string(this_minute)
if this_minute lt 10 then this_minute_str = '0' + this_minute_str


this_time_string = strcompress(this_year_str + '-' + this_month_str + '-' + this_day_str,/remove_all) + ' ' + $
                   strcompress(this_hour_str + ':' + this_minute_str,/remove_all) + ' UT'

;print, this_time_string



CASE thisEvent OF

   'EXPOSE': BEGIN
       END
   'PRESS': BEGIN
   
   info.timeline_left_mouse_button_being_pressed = 1

;       print,'Plot Window X and Y : ',event.x , event.y

       END

   'RELEASE': BEGIN
   
       if event.release eq 4 then begin
       
    contextBase = WIDGET_INFO(event.TOP, FIND_BY_UNAME = 'timeline_contextBase')
    WIDGET_DISPLAYCONTEXTMENU, event.ID, event.X, event.Y, contextBase
    
info.timeline_normalized_x = float(event.x) / info.xsize_timeline
info.timeline_left_mouse_button_being_pressed = 0
       
       endif else begin
   
    info.timeline_left_mouse_button_being_pressed = 0
   
   info.clicked_L = 0
   info.clicked_C = 0
   info.clicked_C2 = 0
   info.clicked_R = 0
   
   if shift_click eq 1 then begin
   info.clicked_L = 1
   info.clicked_C = 1
   info.clicked_C2 = 0
   info.clicked_R = 1
   endif else begin
   fraction_y = float(event.y)/float(info.ysize_timeline)
   if fraction_y ge 0.34375 and fraction_y le 0.4375 then info.clicked_L = 1
   if fraction_y ge 0.453125 and fraction_y le 0.5625 then info.clicked_C2 = 1
   if fraction_y ge 0.559375 and fraction_y le 0.703125 then info.clicked_C = 1
   if fraction_y ge 0.734375 and fraction_y le 0.84375 then info.clicked_R = 1      
   endelse
   

if info.L_number_of_images gt 0 and info.clicked_L eq 1 then begin

; ECDP
;if (info.which_mode eq 1) then begin
;  info.L_cme_outline -> setProperty, hide = 1
;endif else begin
;  info.L_cme_outline -> setProperty, hide = 0
;endelse
; ECDP  

L_julian = (info.L_list_of_datetime_Julian).toarray()
L_index = (where(this_julian-L_julian lt 0.0))[0]

if L_index gt 0 then begin
   if abs(L_julian[L_index - 1] - this_julian) lt abs(L_julian[L_index] - this_julian) then L_index --
endif

if abs(L_julian[L_index] - this_julian) lt (1./48.) then begin

if L_index eq -1 then L_index = info.L_number_of_images - 1

info.L_current_image_number = L_index
widget_control,info.L_widget_image_sequence_slider,set_value = info.L_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.L_current_image_number,info.L_background_image_number,info.L_difference_imaging, $
    info.L_list_of_image_data,info.L_image_saturation_value,info.L_coronagraph_image_object,info.L_border_image_object, $
    info.CME_matches_image_L_Image_number,info.L_current_background_color, $
    info.background_color,info.L_current_text_color,info.color_stereo_B,info.L_cme_outline,info.L_cme_MATCH_outline, $
    info.L_widget_outline_matches_image,info.CME_matches_image_L_CME_outline, $
    info.L_ut_string_object,info.L_list_of_full_time_strings,info.L_title_object,info.L_Window,info.L_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color

swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline
  if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.L_cme_outline -> setProperty, hide = 1
   endif else begin
   info.L_cme_outline -> setProperty, hide = 0 
   endelse 


endif else begin

info.L_Window->Draw, info.L_both_views
info.L_window->erase, color=info.background_color_stereo_B
info.L_plot->SetProperty, color=info.color_stereo_B

endelse

endif


if info.C_number_of_images gt 0 and info.clicked_C eq 1 then begin

; ECDP
;if (info.which_mode eq 1) then begin
;  info.C_cme_outline -> setProperty, hide = 1
;endif else begin
;  info.C_cme_outline -> setProperty, hide = 0
;endelse
; ECDP  
info.C_cme_MATCH_outline-> setProperty, hide = 1
info.C2_cme_MATCH_outline-> setProperty, hide = 1

info.currently_showing_lasco_C2_or_C3 = 'C3'
swpc_cat_calculate_c2_c3_lead_edge, info
;print, 'Here!!! =', info.C_leading_edge_in_rs
widget_control,info.widget_show_C2_or_C3,set_value='Show LASCO C2'

C_julian = (info.C_list_of_datetime_Julian).toarray()
C_index = (where(this_julian-C_julian lt 0.0))[0]

if C_index gt 0 then begin
   if abs(C_julian[C_index - 1] - this_julian) lt abs(C_julian[C_index] - this_julian) then C_index --
endif

if abs(C_julian[C_index] - this_julian) lt (1./48.) then begin

if C_index eq -1 then C_index = info.C_number_of_images - 1
info.C_current_image_number = C_index
widget_control, info.C_widget_image_sequence_slider,set_slider_max = n_elements(info.C_list_of_datetime_Julian)
widget_control,info.C_widget_image_sequence_slider,set_value = info.C_current_image_number + 1

info.C_title_object -> setproperty, strings = 'SOHO LASCO C3'

swpc_cat_REDRAW_THE_IMAGE, $
    info.C_current_image_number,info.C_background_image_number,info.C_difference_imaging, $
    info.C_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c3,info.C_cme_outline,info.C_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C_CME_outline, $
    info.C_ut_string_object,info.C_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color

swpc_cat_Calculate_Earth_B_Angle,(info.C_list_of_datetime_Julian)[0],B_angle_degrees
info.C_HEEQ_coords[1] = B_angle_degrees

info.C_coords_offset_model -> SetProperty,transform = info.initial_transform
info.C_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.C_coords_offset_model->rotate,[1,0,0], info.C_HEEQ_coords[1]
info.C_coords_offset_model_copy->rotate,[1,0,0], info.C_HEEQ_coords[1]

info.C_telescope_field_of_view_Rs = (256. * ((info.C_list_of_pixel_scales)[0] / (info.C_list_of_image_scaling_factors)[0])) / (info.C_list_of_rsuns)[0]

info.C_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]
info.C_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]

the_day = long((info.C_list_of_datetime_Julian)[0])
i_day = where(the_day lt info.Julian_day_for_Earth_pos)
i_day = i_day[0]

info.C_second_view -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99  ; 0.99 factor is for L1 as opposed to Earth.
info.C_second_view_copy -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99

info.C_Fixed_Sun_Satellite_model -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
info.C_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model_copy -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
;info.C_cme_outline -> setProperty, hide = 0

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
;print, 'WHich mode = ', info.which_mode
   if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.C_cme_outline -> setProperty, hide = 1
   endif else begin
   info.C_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP

info.C_Window->Draw, info.C_both_views

;set_timeline_highlight_block, info.C_plot, info.C_number_of_images, C_index, info.color_C3, info.cme_outline_color

endif else begin

info.C_window->erase, color=info.background_color_lasco
info.C_plot->SetProperty, color = info.color_C3

endelse

endif

if info.C2_number_of_images gt 0 and info.clicked_C2 eq 1 then begin

; ECDP
;if (info.which_mode eq 1) then begin
;  info.C_cme_outline -> setProperty, hide = 1
;endif else begin
;  info.C_cme_outline -> setProperty, hide = 0
;endelse
; ECDP  
;info.C_cme_outline -> setProperty, hide = 1
info.C_cme_MATCH_outline-> setProperty, hide = 1
info.C2_cme_MATCH_outline-> setProperty, hide = 1

info.currently_showing_lasco_C2_or_C3 = 'C2'
swpc_cat_calculate_c2_c3_lead_edge, info
;print, 'Here!!!! =', info.C_leading_edge_in_rs
widget_control,info.widget_show_C2_or_C3,set_value='Show LASCO C3'

C2_julian = (info.C2_list_of_datetime_Julian).toarray()
C2_index = (where(this_julian-C2_julian lt 0.0))[0]
if C2_index gt 0 then begin
if abs(C2_julian[C2_index - 1] - this_julian) lt abs(C2_julian[C2_index] - this_julian) then C2_index --
endif

if abs(C2_julian[C2_index] - this_julian) lt (1./48.) then begin

if C2_index eq -1 then C2_index = info.C2_number_of_images - 1
info.C2_current_image_number = C2_index
widget_control, info.C_widget_image_sequence_slider,set_slider_max = n_elements(info.C2_list_of_datetime_Julian)
widget_control,info.C_widget_image_sequence_slider,set_value = info.C2_current_image_number + 1

info.C_title_object -> setproperty, strings = 'SOHO LASCO C2'

swpc_cat_REDRAW_THE_IMAGE, $
    info.C2_current_image_number,info.C2_background_image_number,info.C2_difference_imaging, $
    info.C2_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C2_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c2,info.C_cme_outline,info.C2_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C2_CME_outline, $
    info.C_ut_string_object,info.C2_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, info.C2_current_image_number, info.color_C2, info.cme_outline_color

swpc_cat_Calculate_Earth_B_Angle,(info.C2_list_of_datetime_Julian)[0],B_angle_degrees
info.C2_HEEQ_coords[1] = B_angle_degrees

info.C_coords_offset_model -> SetProperty,transform = info.initial_transform
info.C_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.C_coords_offset_model->rotate,[1,0,0], info.C2_HEEQ_coords[1]
info.C_coords_offset_model_copy->rotate,[1,0,0], info.C2_HEEQ_coords[1]

info.C2_telescope_field_of_view_Rs = (256. * ((info.C2_list_of_pixel_scales)[0] / (info.C2_list_of_image_scaling_factors)[0])) / (info.C2_list_of_rsuns)[0]

info.C_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.C2_telescope_field_of_view_Rs,2.0*info.C2_telescope_field_of_view_Rs]
info.C_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.C2_telescope_field_of_view_Rs,2.0*info.C2_telescope_field_of_view_Rs]

the_day = long((info.C2_list_of_datetime_Julian)[0])
i_day = where(the_day lt info.Julian_day_for_Earth_pos)
i_day = i_day[0]

info.C_second_view -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99  ; 0.99 factor is for L1 as opposed to Earth.
info.C_second_view_copy -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99

info.C_Fixed_Sun_Satellite_model -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model -> translate, info.C2_telescope_field_of_view_Rs , info.C2_telescope_field_of_view_Rs, 0.0
info.C_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model_copy -> translate, info.C2_telescope_field_of_view_Rs , info.C2_telescope_field_of_view_Rs, 0.0

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
   if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.C_cme_outline -> setProperty, hide = 1
   endif else begin
   info.C_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP
info.C_Window->Draw, info.C_both_views

;set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, C2_index, info.color_C2, info.cme_outline_color

endif else begin

info.C_window->erase, color=info.background_color_lasco
info.C2_plot->SetProperty, color = info.color_C2

endelse

endif


if info.R_number_of_images gt 0 and info.clicked_R eq 1 then begin

;  if (info.which_mode eq 1) then begin
;   info.R_cme_outline -> setProperty, hide = 1
;   endif else begin
;   info.R_cme_outline -> setProperty, hide = 0 
;   endelse 
; ECDP

R_julian = (info.R_list_of_datetime_Julian).toarray()
R_index = (where(this_julian-R_julian lt 0.0))[0]
if R_index gt 0 then begin 
if abs(R_julian[R_index - 1] - this_julian) lt abs(R_julian[R_index] - this_julian) then R_index --
endif

if abs(R_julian[R_index] - this_julian) lt (1./48.) then begin

if R_index eq -1 then R_index = info.R_number_of_images - 1
info.R_current_image_number = R_index
widget_control,info.R_widget_image_sequence_slider,set_value = info.R_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.R_current_image_number,info.R_background_image_number,info.R_difference_imaging, $
    info.R_list_of_image_data,info.R_image_saturation_value,info.R_coronagraph_image_object,info.R_border_image_object, $
    info.CME_matches_image_R_Image_number,info.R_current_background_color, $
    info.background_color,info.R_current_text_color,info.color_stereo_A,info.R_cme_outline,info.R_cme_MATCH_outline, $
    info.R_widget_outline_matches_image,info.CME_matches_image_R_CME_outline, $
    info.R_ut_string_object,info.R_list_of_full_time_strings,info.R_title_object,info.R_Window,info.R_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_a, info.cme_outline_color

swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline
if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
  info.R_cme_outline -> setProperty, hide = 1
endif else begin
  info.R_cme_outline -> setProperty, hide = 0
endelse
; ECDP  
endif else begin

info.R_Window->Draw, info.R_both_views
info.R_window->erase, color=info.background_color_stereo_A
info.R_plot->SetProperty, color=info.color_stereo_A

endelse

endif



endelse ; not right-click context menu



       END
       
       
   'MOTION': BEGIN ; Trackball events
   
   if info.timeline_left_mouse_button_being_pressed eq 1 then begin
   
if event.y gt 10 then begin
   
time_line_position = ((this_julian - info.start_julian)/(info.end_julian - info.start_julian) * 0.9) + 0.05
line_data=fltarr(2,2)
line_data[0,0] = time_line_position
line_data[1,0] = 0.
line_data[0,1] = time_line_position
line_data[1,1] = 1.
info.animation_current_time_marker->setproperty,data=line_data
info.images_timeline_window->Draw, info.images_timeline_view


   
   info.clicked_L = 1
   info.clicked_C = 1
   info.clicked_C2 = 0
   info.clicked_R = 1
   
if info.L_number_of_images gt 0 and info.clicked_L eq 1 then begin

L_julian = (info.L_list_of_datetime_Julian).toarray()
L_index = (where(this_julian-L_julian lt 0.0))[0]
if L_index gt 0 then begin
if abs(L_julian[L_index - 1] - this_julian) lt abs(L_julian[L_index] - this_julian) then L_index --
endif

if abs(L_julian[L_index] - this_julian) lt (1./48.) then begin

if L_index eq -1 then L_index = info.L_number_of_images - 1

info.L_current_image_number = L_index
widget_control,info.L_widget_image_sequence_slider,set_value = info.L_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.L_current_image_number,info.L_background_image_number,info.L_difference_imaging, $
    info.L_list_of_image_data,info.L_image_saturation_value,info.L_coronagraph_image_object,info.L_border_image_object, $
    info.CME_matches_image_L_Image_number,info.L_current_background_color, $
    info.background_color,info.L_current_text_color,info.color_stereo_B,info.L_cme_outline,info.L_cme_MATCH_outline, $
    info.L_widget_outline_matches_image,info.CME_matches_image_L_CME_outline, $
    info.L_ut_string_object,info.L_list_of_full_time_strings,info.L_title_object,info.L_Window,info.L_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.L_plot, info.L_number_of_images, info.L_current_image_number, info.color_stereo_B, info.cme_outline_color
; ECDP
swpc_cat_update_cme_outline,info.L_Window_copy,info.L_second_view_copy,info.L_cme_outline
   if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.L_cme_outline -> setProperty, hide = 1
   endif else begin
   info.L_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP
   
endif else begin


info.L_Window->Draw, info.L_both_views
info.L_window->erase, color=info.background_color_stereo_B
info.L_plot->SetProperty, color=info.color_stereo_B

endelse

endif





if info.C_number_of_images gt 0 and info.clicked_C eq 1 then begin

info.C_cme_outline -> setProperty, hide = 1
info.C_cme_MATCH_outline-> setProperty, hide = 1
info.C2_cme_MATCH_outline-> setProperty, hide = 1

info.currently_showing_lasco_C2_or_C3 = 'C3'
widget_control,info.widget_show_C2_or_C3,set_value='Show LASCO C2'

C_julian = (info.C_list_of_datetime_Julian).toarray()
C_index = (where(this_julian-C_julian lt 0.0))[0]
if C_index gt 0 then begin
if abs(C_julian[C_index - 1] - this_julian) lt abs(C_julian[C_index] - this_julian) then C_index --
endif

if abs(C_julian[C_index] - this_julian) lt (1./48.) then begin

if C_index eq -1 then C_index = info.C_number_of_images - 1
info.C_current_image_number = C_index
widget_control, info.C_widget_image_sequence_slider,set_slider_max = n_elements(info.C_list_of_datetime_Julian)
widget_control,info.C_widget_image_sequence_slider,set_value = info.C_current_image_number + 1

info.C_title_object -> setproperty, strings = 'SOHO LASCO C3'

swpc_cat_REDRAW_THE_IMAGE, $
    info.C_current_image_number,info.C_background_image_number,info.C_difference_imaging, $
    info.C_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c3,info.C_cme_outline,info.C_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C_CME_outline, $
    info.C_ut_string_object,info.C_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.C_plot, info.C_number_of_images, info.C_current_image_number, info.color_C3, info.cme_outline_color

swpc_cat_Calculate_Earth_B_Angle,(info.C_list_of_datetime_Julian)[0],B_angle_degrees
info.C_HEEQ_coords[1] = B_angle_degrees

info.C_coords_offset_model -> SetProperty,transform = info.initial_transform
info.C_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.C_coords_offset_model->rotate,[1,0,0], info.C_HEEQ_coords[1]
info.C_coords_offset_model_copy->rotate,[1,0,0], info.C_HEEQ_coords[1]

info.C_telescope_field_of_view_Rs = (256. * ((info.C_list_of_pixel_scales)[0] / (info.C_list_of_image_scaling_factors)[0])) / (info.C_list_of_rsuns)[0]

info.C_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]
info.C_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.C_telescope_field_of_view_Rs,2.0*info.C_telescope_field_of_view_Rs]

the_day = long((info.C_list_of_datetime_Julian)[0])
i_day = where(the_day lt info.Julian_day_for_Earth_pos)
i_day = i_day[0]

info.C_second_view -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99  ; 0.99 factor is for L1 as opposed to Earth.
info.C_second_view_copy -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99

info.C_Fixed_Sun_Satellite_model -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
info.C_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model_copy -> translate, info.C_telescope_field_of_view_Rs , info.C_telescope_field_of_view_Rs, 0.0
;info.C_cme_outline -> setProperty, hide = 0

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
   if (info.which_mode eq 1) then begin

;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
;    print, info.c_mode_name, string(info.c_mode_name)
;   stop 
   info.C_cme_outline -> setProperty, hide = 1
   endif else begin
   info.C_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP
;info.C_Window->Draw, info.C_both_views

;set_timeline_highlight_block, info.C_plot, info.C_number_of_images, C_index, info.color_C3, info.cme_outline_color

endif else begin
info.C_window->Draw, info.C_both_views
info.C_window->erase, color=info.background_color_lasco
info.C_plot->SetProperty, color = info.color_C3

endelse

endif


; ECDP
if info.C2_number_of_images gt 0 and info.clicked_C2 eq 1 then begin

info.C_cme_outline -> setProperty, hide = 1
info.C_cme_MATCH_outline-> setProperty, hide = 1
info.C2_cme_MATCH_outline-> setProperty, hide = 1

info.currently_showing_lasco_C2_or_C3 = 'C2'
widget_control,info.widget_show_C2_or_C3,set_value='Show LASCO C3'

C2_julian = (info.C2_list_of_datetime_Julian).toarray()
C2_index = (where(this_julian-C2_julian lt 0.0))[0]
if C2_index gt 0 then begin
if abs(C2_julian[C2_index - 1] - this_julian) lt abs(C2_julian[C2_index] - this_julian) then C2_index --
endif

if abs(C2_julian[C2_index] - this_julian) lt (1./48.) then begin

if C2_index eq -1 then C2_index = info.C2_number_of_images - 1
info.C2_current_image_number = C2_index
widget_control, info.C_widget_image_sequence_slider,set_slider_max = n_elements(info.C2_list_of_datetime_Julian)
widget_control,info.C_widget_image_sequence_slider,set_value = info.C2_current_image_number + 1

info.C_title_object -> setproperty, strings = 'SOHO LASCO C2'

swpc_cat_REDRAW_THE_IMAGE, $
    info.C2_current_image_number,info.C2_background_image_number,info.C2_difference_imaging, $
    info.C2_list_of_image_data,info.C_image_saturation_value,info.C_coronagraph_image_object,info.C_border_image_object, $
    info.CME_matches_image_C2_Image_number,info.C_current_background_color, $
    info.background_color,info.C_current_text_color,info.color_c2,info.C_cme_outline,info.C2_cme_MATCH_outline, $
    info.C_widget_outline_matches_image,info.CME_matches_image_C2_CME_outline, $
    info.C_ut_string_object,info.C2_list_of_full_time_strings,info.C_title_object,info.C_Window,info.C_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, info.C2_current_image_number, info.color_C2, info.cme_outline_color

swpc_cat_Calculate_Earth_B_Angle,(info.C2_list_of_datetime_Julian)[0],B_angle_degrees
info.C2_HEEQ_coords[1] = B_angle_degrees

info.C_coords_offset_model -> SetProperty,transform = info.initial_transform
info.C_coords_offset_model_copy -> SetProperty,transform = info.initial_transform

info.C_coords_offset_model->rotate,[1,0,0], info.C2_HEEQ_coords[1]
info.C_coords_offset_model_copy->rotate,[1,0,0], info.C2_HEEQ_coords[1]

info.C2_telescope_field_of_view_Rs = (256. * ((info.C2_list_of_pixel_scales)[0] / (info.C2_list_of_image_scaling_factors)[0])) / (info.C2_list_of_rsuns)[0]

info.C_second_view->SetProperty, Viewplane_Rect=[0,0,2.0*info.C2_telescope_field_of_view_Rs,2.0*info.C2_telescope_field_of_view_Rs]
info.C_second_view_copy->SetProperty, Viewplane_Rect=[0,0,2.0*info.C2_telescope_field_of_view_Rs,2.0*info.C2_telescope_field_of_view_Rs]

the_day = long((info.C2_list_of_datetime_Julian)[0])
i_day = where(the_day lt info.Julian_day_for_Earth_pos)
i_day = i_day[0]

info.C_second_view -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99  ; 0.99 factor is for L1 as opposed to Earth.
info.C_second_view_copy -> setproperty, eye = 215. * info.Earth_pos_AU[i_day] * 0.99

info.C_Fixed_Sun_Satellite_model -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model -> translate, info.C2_telescope_field_of_view_Rs , info.C2_telescope_field_of_view_Rs, 0.0
info.C_Fixed_Sun_Satellite_model_copy -> SetProperty,transform = info.C_Fixed_Sun_Satellite_model_original_transform
info.C_Fixed_Sun_Satellite_model_copy -> translate, info.C2_telescope_field_of_view_Rs , info.C2_telescope_field_of_view_Rs, 0.0

swpc_cat_update_cme_outline,info.C_Window_copy,info.C_second_view_copy,info.C_cme_outline
info.C_Window->Draw, info.C_both_views

;set_timeline_highlight_block, info.C2_plot, info.C2_number_of_images, C2_index, info.color_C2, info.cme_outline_color

endif else begin

info.C_Window->Draw, info.C_both_views
info.C_window->erase, color=info.background_color_lasco
info.C2_plot->SetProperty, color = info.color_C2

endelse

endif
; ECDP






if info.R_number_of_images gt 0 and info.clicked_R eq 1 then begin

R_julian = (info.R_list_of_datetime_Julian).toarray()
R_index = (where(this_julian-R_julian lt 0.0))[0]
if R_index gt 0 then begin 
if abs(R_julian[R_index - 1] - this_julian) lt abs(R_julian[R_index] - this_julian) then R_index --
endif

if abs(R_julian[R_index] - this_julian) lt (1./48.) then begin

if R_index eq -1 then R_index = info.R_number_of_images - 1
info.R_current_image_number = R_index
widget_control,info.R_widget_image_sequence_slider,set_value = info.R_current_image_number + 1

swpc_cat_REDRAW_THE_IMAGE, $
    info.R_current_image_number,info.R_background_image_number,info.R_difference_imaging, $
    info.R_list_of_image_data,info.R_image_saturation_value,info.R_coronagraph_image_object,info.R_border_image_object, $
    info.CME_matches_image_R_Image_number,info.R_current_background_color, $
    info.background_color,info.R_current_text_color,info.color_stereo_A,info.R_cme_outline,info.R_cme_MATCH_outline, $
    info.R_widget_outline_matches_image,info.CME_matches_image_R_CME_outline, $
    info.R_ut_string_object,info.R_list_of_full_time_strings,info.R_title_object,info.R_Window,info.R_both_views,0,0, info.i_log_scale

swpc_cat_set_timeline_highlight_block, info.R_plot, info.R_number_of_images, info.R_current_image_number, info.color_stereo_a, info.cme_outline_color

; ECDP
swpc_cat_update_cme_outline,info.R_Window_copy,info.R_second_view_copy,info.R_cme_outline
   if (info.which_mode eq 1) then begin
;if ((info.which_mode eq 1) or (string(info.c_mode_name) eq 'SELECT MODE')) then begin
   info.R_cme_outline -> setProperty, hide = 1
   endif else begin
   info.R_cme_outline -> setProperty, hide = 0 
   endelse 
; ECDP
;info.R_Window->Draw, info.R_both_views

endif else begin
info.R_Window->Draw, info.R_both_views
info.R_window->erase, color=info.background_color_stereo_A
info.R_plot->SetProperty, color=info.color_stereo_A

endelse

endif




endif else begin  ; if event.y lt 10

; determine whether we are moving the start or end markers
 
       myLoc = [event.x, event.y]    
       oObjArr = info.images_timeline_window->Select(info.images_timeline_view, myLoc, dimensions=[100,20])
       nSel = N_ELEMENTS(oObjArr)
       if nsel gt 1 then begin
       iobject = intarr(nSel)
       
       FOR i=0, nSel-1 DO BEGIN
          oObjArr[i]->GetProperty, NAME=name
          if name eq 'animation_start_time_marker' then begin

animation_start_position = ((this_julian - info.start_julian)/(info.end_julian - info.start_julian) * 0.9) + 0.05
line_data=fltarr(2,2)
line_data[0,0] = animation_start_position
line_data[1,0] = 0.
line_data[0,1] = animation_start_position
line_data[1,1] = 1.
info.animation_start_time_marker->setproperty,data=line_data

xpos = animation_start_position
info.animation_start_time_marker_handle_data=[[xpos,0.],[xpos+0.005,0.],[xpos+0.005,0.08],[xpos,0.08]]
info.animation_start_time_marker_handle->setproperty,data=info.animation_start_time_marker_handle_data

info.images_timeline_window->Draw, info.images_timeline_view

CASE info.which_window_to_animate OF
   0 : BEGIN
array_of_julians = info.L_list_of_datetime_Julian->toarray()
   ENDCASE
   1 : BEGIN
array_of_julians = info.C_list_of_datetime_Julian->toarray()
   ENDCASE
   2 : BEGIN
array_of_julians = info.R_list_of_datetime_Julian->toarray()
   ENDCASE
   
ENDCASE

l = where(array_of_julians gt this_julian)
info.anim_start_frame = l[0]
info.anim_start_julian = this_julian

goto, jumphere

          endif
          if name eq 'animation_end_time_marker' then begin

animation_end_position = ((this_julian - info.start_julian)/(info.end_julian - info.start_julian) * 0.9) + 0.05
line_data=fltarr(2,2)
line_data[0,0] = animation_end_position
line_data[1,0] = 0.
line_data[0,1] = animation_end_position
line_data[1,1] = 1.
info.animation_end_time_marker->setproperty,data=line_data

xpos = animation_end_position
info.animation_end_time_marker_handle_data=[[xpos,0.],[xpos-0.005,0.],[xpos-0.005,0.08],[xpos,0.08]]
info.animation_end_time_marker_handle->setproperty,data=info.animation_end_time_marker_handle_data

info.images_timeline_window->Draw, info.images_timeline_view

CASE info.which_window_to_animate OF

   0 : BEGIN
array_of_julians = info.L_list_of_datetime_Julian->toarray()
   ENDCASE
   1 : BEGIN
array_of_julians = info.C_list_of_datetime_Julian->toarray()
   ENDCASE
   2 : BEGIN
array_of_julians = info.R_list_of_datetime_Julian->toarray()
   ENDCASE
   
ENDCASE

l = where(array_of_julians lt this_julian)
info.anim_end_frame = l[-1]
info.anim_end_julian = this_julian

goto, jumphere

          endif
       ENDFOR
       
jumphere :

       endif



endelse


   
;   print, ' MOTION ', info.timeline_left_mouse_button_being_pressed ;event.x, event.y
;print, 'Here2'
;print, 'MOTION', event.x, event.y

endif   ; if left mouse button currently clicked
   


END

   ELSE:


ENDCASE

    ; Draw the view.

info.images_timeline_window->Draw, info.images_timeline_view

    ;Put the info structure back.
    
endif



Widget_Control, event.top, Set_UValue=info, /No_Copy
END

PRO swpc_cat_calculate_c2_c3_lead_edge, info
compile_opt idl2
;Widget_Control, event.top, Get_UValue=info, /No_Copy

; ECDP experiment with this
          case info.currently_showing_lasco_C2_or_C3 of

          'C3' : begin
            
          C_leading_edge_float = (info.C_now_dist * (info.C_list_of_pixel_scales)[0] / (info.C_list_of_image_scaling_factors)[0]) / (info.C_list_of_rsuns)[0]  

          endcase
          
          
          'C2' : begin
            
          C_leading_edge_float = (info.C_now_dist * (info.C2_list_of_pixel_scales)[0] / (info.C2_list_of_image_scaling_factors)[0]) / (info.C2_list_of_rsuns)[0]  

          endcase
                   
          endcase  
; ECDP

  info.C_leading_edge_in_rs = C_leading_edge_float      
  C_lead_edge_string = string(C_leading_edge_float,format='(f5.1)')

;      if info2.leading_edge_in_rs lt 10 then lead_edge_string = ' ' + lead_edge_string
;         ;Displays new value of the leading edge near the yellow leading edge line
  info.C_leading_edge_string->SetProperty,strings=C_lead_edge_string


;Widget_Control, event.top, Set_UValue=info, /No_Copy
end

; ECDP
PRO select_mode, event

compile_opt idl2

Widget_Control, event.top, Get_UValue=info, /No_Copy

Widget_Control, event.id, Get_UValue=which_mode


CASE which_mode OF

   'NORMAL': BEGIN
   
      info.which_mode = 0
      info.C_clock_angle_model -> SetProperty, Hide = 1
      info.C_leading_edge_model -> SetProperty, Hide = 1
      info.C_cme_outline->SetProperty, Hide = 0
      info.L_cme_outline->SetProperty, Hide = 0
      info.R_cme_outline->SetProperty, Hide = 0
      widget_control,info.widget_latitude_slider,/sensitive
      widget_control,info.widget_longitude_slider,/sensitive
      widget_control,info.widget_angular_width_slider,/sensitive
      widget_control,info.widget_radial_distance_slider,/sensitive
;      info.C_image_model->setproperty, strings = 'NORMAL (SOHO + STEREO)'

      if (info.compact eq 0) then begin
      widget_control, info.c_mode_name, set_value = 'NORMAL'
;      info.C_Window->Draw, info.C_both_views
      widget_control, info.time_plot_title_widget, set_value = 'CME LEADING EDGE (Rs) vs TIME (UT)' 
      endif
      
      widget_control, info.C_widget_outline_matches_image, sensitive = 1
      widget_control, info.C_widget_pos_outline_matches_image, sensitive = 0
      widget_control, info.C_widget_representative_image, sensitive = 1
      widget_control, info.L_widget_representative_image, sensitive = 1
      widget_control, info.R_widget_representative_image, sensitive = 1
      widget_control, info.L_widget_outline_matches_image, sensitive = 1
      widget_control, info.R_widget_outline_matches_image, sensitive = 1
      
;      if (info.C_the_action eq 0) then begin
;      	Widget_Control, info.draw_C, Draw_Motion_Events=1
;      endif

      info.C_Window->Draw, info.C_both_views

;      info.C_clock_angle_model -> GetProperty, renderer

      ENDCASE
 
      
   'POS': BEGIN
   
      info.which_mode = 1

;      if (info.C_the_action eq 0) then begin
;      	Widget_Control, info.draw_C, Draw_Motion_Events=1
;      endif	

      info.C_clock_angle_model -> SetProperty, Hide = 0
      info.C_leading_edge_model -> SetProperty, Hide = 0
      info.C_cme_outline->SetProperty, Hide = 1
      info.L_cme_outline->SetProperty, Hide = 1
      info.R_cme_outline->SetProperty, Hide = 1
      widget_control,info.widget_latitude_slider, sensitive = 0
      widget_control,info.widget_longitude_slider, sensitive = 0
      widget_control,info.widget_angular_width_slider, sensitive = 0
      widget_control,info.widget_radial_distance_slider, sensitive = 0
      
      if (info.compact eq 0) then begin
      widget_control, info.time_plot_title_widget, set_value = 'POS LEADING EDGE (Rs) vs TIME (UT)' 
      widget_control, info.c_mode_name, set_value = 'PLANE OF SKY'
      endif
      
      widget_control, info.C_widget_pos_outline_matches_image, sensitive = 1
      widget_control, info.C_widget_outline_matches_image, sensitive = 0
      widget_control, info.C_widget_representative_image, sensitive = 0
      widget_control, info.L_widget_representative_image, sensitive = 0
      widget_control, info.R_widget_representative_image, sensitive = 0
      widget_control, info.L_widget_outline_matches_image, sensitive = 0
      widget_control, info.R_widget_outline_matches_image, sensitive = 0
      info.C_Window->Draw, info.C_both_views
       
      ENDCASE
      
;   'FULL HALO': BEGIN
   
;      info.which_mode = 2
;      info.c_hide = 1
;      info.C_clock_angle_model -> SetProperty, Hide = info.c_hide
;      info.C_leading_edge_model -> SetProperty, Hide = info.c_hide
;      widget_control, info.time_plot_title_widget, set_value = 'FULL HALO CME' 
;      widget_control, info.draw_plot_window, set_value = 'hello3'
;      info.plot_velocity_title -> SetProperty, Hide = info.c_hide
;      print, 'C_hide =', info.c_hide 
;      info2.background_image_number = info2.current_image_number

;      ENDCASE

ENDCASE

info.C_Window->Draw, info.C_both_views

Widget_Control, event.top, Set_UValue=info, /No_Copy

END



PRO SWPC_CAT

compile_opt idl2

; a variable to make the window always come up at center of the main monitor 
; ECDP
i_center_on_main_monitor = 1
; ECDP
; 
; new compact variable builds a smaller application to fit on the 
; ops laptop screen (res = 1600 * 900)

compact = 1
compact_frame = 1
if compact eq 1 then compact_frame = 0

swpc_cat_1point0_production_build = 1
; do we want to export line plots for velocity and POS velocity (0 = no, 1 = yes)
output_line_plots = 0

!except = 0
retain = 2
allow_Bernouli = 0
allow_manual_data_point = 1
show_image_line_plot = 1
allow_show_3D_view = 1
show_cone_Z_axis = 1
calculate_individual_velocities_for_each_telescope = 1
debug_mode = 1
output_matched_line_data_in_txt_file = 0
output_matched_pos_line_data_in_txt_file = 0


if swpc_cat_1point0_production_build eq 1 then begin
  allow_Bernouli = 0
  allow_manual_data_point = 1
  show_image_line_plot = 1
  allow_show_3D_view = 0
  show_cone_Z_axis = 0
  calculate_individual_velocities_for_each_telescope = 0
  debug_mode = 0
  output_matched_line_data_in_txt_file = 0
endif



color_red = [138,226,70]
color_yellow = [220,220,0]
color_fill = [50,50,50]
background_color = [35,35,35]
text_color = [153,153,153]
cme_outline_color = [255,255,0]
color_stereo_B = [150,150,255]
color_stereo_A = [255,100,100]
color_c3 = [100,255,100]
;color_c2 = [0,100,0]
color_c2 = [100,220,220]

color_cme_orange = [220,120,0]

L_current_background_color = background_color
L_current_text_color = color_stereo_B

C_current_background_color = background_color
C_current_text_color = color_c3

R_current_background_color = background_color
R_current_text_color = color_stereo_A

; set show_copy_windows to 1 in order to see the black/white mattes that are used
; to create the lemniscate outlines. Normally this is set to 0 and the matte process
; is done in an invisible buffer window. 
show_copy_windows = 0

;get directory path associated with this IDL procedure
source_path = swpc_cat_SourcePath()

;get input values that can change, from a configuration file
;swpc_cat_get_config,source_path,sep,telescope_array,image_in_folder_array, $
;           image_in_root,export_location_root, $
;           max_interval_in_days,OPS_or_VnV, $
;           input_file

;get input values that can change, from a configuration file
; ECDP
;swpc_cat_get_config,source_path,sep,telescope_array,image_in_folder_array, $
;           image_in_root,export_location_root,preferences_folder, $
;           max_interval_in_days,OPS_or_VnV, $
;           input_file
swpc_cat_get_config,source_path,sep,telescope_array,image_in_folder_array, $
           image_in_root,export_location_root,preferences_folder, $
           max_interval_in_days,OPS_or_VnV, $
           earth_position_filename, $
           input_file



           
Earth_coords_file = source_path + sep + 'swpc_cat_data' + sep + earth_position_filename
;print, Earth_coords_file
one_line = ''
openr,lun,Earth_coords_file,/get_lun
readf,lun,one_line ; header line before actual numbers
Julian_day_for_Earth_pos = list()
Earth_pos_AU = list()
Earth_pos_HG_LAT_deg = list()
while ~ eof(lun) do begin
readf, lun, one_line
result = strsplit(one_line,/extract)
;YYYY DOY     AU    HG_LAT  HG_LON  HGI_LON
yyyy = long(result[0])
doy = long(result[1])
Julian_day_for_Earth_pos.add, JULDAY(1, 1, yyyy) + doy - 1L
Earth_pos_AU.add, float(result[2])
Earth_pos_HG_LAT_deg.add, float(result[3])
endwhile
close, lun
Julian_day_for_Earth_pos = Julian_day_for_Earth_pos.toarray()
Earth_pos_AU = Earth_pos_AU.toarray()
Earth_pos_HG_LAT_deg = Earth_pos_HG_LAT_deg.toarray()

L_number_of_images = 0
C_number_of_images = 0
C2_number_of_images = 0
R_number_of_images = 0

Xsize = 512
Ysize = 512

L_both_views = OBJ_NEW('IDLgrViewGroup')
L_image_view = OBJ_NEW('IDLgrView', Color=background_color, Viewplane_Rect=[0.0,0.0,Xsize,Ysize])
L_both_views -> add, L_image_view
L_image_model = OBJ_NEW('IDLgrModel') 
L_image_view -> Add, L_image_model
L_image_model->SetProperty, hide = 0

L_clock_angle_model = OBJ_NEW('IDLgrModel') 
L_clock_angle_model_location = [0.,0.]
L_clock_angle_model -> translate,255,255,0.0
L_clock_angle_model_location = [255,255]
L_clock_angle_model->GetProperty, transform = L_clock_angle_model_original_transform
; Marker for the 'clock angle' line
L_clock_angle_marker_line_x = [0.,0.]
L_clock_angle_marker_line_y = [0.,ysize]
L_clock_angle_marker_line = Obj_New("IDLgrPolyline", L_clock_angle_marker_line_x , L_clock_angle_marker_line_y, $
name = 'L_clock_angle', Color=[255,0,0], thick=1)
L_clock_angle_model -> Add, L_clock_angle_marker_line
L_image_view -> Add, L_clock_angle_model

; Model and Marker for the 'leading edge' line

;leading_edge_string_value = ' 0'
L_leading_edge_model = OBJ_NEW('IDLgrModel')

;marker for leading edge line:
L_leading_edge_marker_line_x = [-xsize/20.,+xsize/20.]
L_leading_edge_marker_line_y = [ysize/7.,ysize/7.]
L_leading_edge_marker_line = Obj_New("IDLgrPolyline", L_leading_edge_marker_line_x , L_leading_edge_marker_line_y, $
name = 'leading_edge', Color=[255,255,0], thick=2)
L_leading_edge_model -> Add, L_leading_edge_marker_line
L_leading_edge_model -> SetProperty, Hide = 0
;Displays the value of the leading edge near the leading edge line
;leading_edge_string = obj_new("idlgrtext", strings = leading_edge_string_value,color=[255,255,0], $
;font=Courier12pt, /onglass, locations = [+xsize/20.,(ysize/7.) + 10],hide=0)
;leading_edge_model -> Add, leading_edge_string
L_clock_angle_model -> Add, L_leading_edge_model


L_image_lineplot_model = OBJ_NEW('IDLgrModel')
x = [0,1]
y = [0,1]
xrange = [0 , 1]
yrange = [0 , 1]
L_image_lineplot = OBJ_NEW('IDLgrPlot', x, y, color=[255,0,0])
position_L_image_lineplot = [10, 10, 502, 100]
L_image_xs = swpc_cat_FSC_Normalize(xrange, Position=[position_L_image_lineplot[0], position_L_image_lineplot[2]])
L_image_ys = swpc_cat_FSC_Normalize(yrange, Position=[position_L_image_lineplot[1], position_L_image_lineplot[3]])
L_image_lineplot->SetProperty, XCoord_Conv=L_image_xs, YCoord_Conv=L_image_ys
L_image_lineplot_model -> Add, L_image_lineplot
L_image_view -> Add, L_image_lineplot_model

;if show_image_line_plot eq 0 then begin
L_image_lineplot_model -> SetProperty, hide = 1
L_clock_angle_model -> SetProperty, hide = 1
;endif

C_both_views = OBJ_NEW('IDLgrViewGroup')
C_image_view = OBJ_NEW('IDLgrView', Color=background_color, Viewplane_Rect=[0.0,0.0,Xsize,Ysize])
C_both_views -> add, C_image_view
C_image_model = OBJ_NEW('IDLgrModel') 
C_image_view -> Add, C_image_model
C_image_model->SetProperty, hide = 0

C_clock_angle_model = OBJ_NEW('IDLgrModel') 
C_clock_angle_model_location = [0.,0.]
C_clock_angle_model -> translate,255,255,0.0
C_clock_angle_model_location = [255,255]
C_clock_angle_model->GetProperty, transform = C_clock_angle_model_original_transform
; Marker for the 'clock angle' line
C_clock_angle_marker_line_x = [0.,0.]
C_clock_angle_marker_line_y = [0.,ysize/2.5]
C_clock_angle_marker_line = Obj_New("IDLgrPolyline", C_clock_angle_marker_line_x , C_clock_angle_marker_line_y, $
name = 'C_clock_angle', Color=[255,0,0], thick=1)
C_clock_angle_model -> Add, C_clock_angle_marker_line
C_clock_angle_string_value = ' 0'

;Displays the clock angle degrees near the clock angle line
C_clock_angle_string = obj_new("idlgrtext", strings = C_clock_angle_string_value,color=[255,0,0], $
font=Courier12pt, /onglass, locations = [3,(ysize/3.) + 23],fill_background=1,fill_color=color_fill, alpha_channel=0.7,hide=0)
;C_clock_angle_string = obj_new("idlgrtext", strings = C_clock_angle_string_value,color=[255,0,0], $
;font=Courier12pt, /onglass,fill_background=1,fill_color=color_fill, alpha_channel=0.7,hide=0)

C_clock_angle_model -> Add, C_clock_angle_string
C_clock_angle_model -> SetProperty, Hide = 1
C_image_view -> Add, C_clock_angle_model

; ECDP
;leading_edge_string_value = ' 0'
C_leading_edge_model = OBJ_NEW('IDLgrModel')

;marker for leading edge line:
C_leading_edge_marker_line_x = [-xsize/20.,+xsize/20.]
C_leading_edge_marker_line_y = [ysize/7.,ysize/7.]
C_leading_edge_marker_line = Obj_New("IDLgrPolyline", C_leading_edge_marker_line_x , C_leading_edge_marker_line_y, $
name = 'leading_edge', Color=[255,255,0], thick=2,alpha_channel = 0.7)
; include the old leading edge lines


C_leading_edge_model -> Add, C_leading_edge_marker_line
;C_leading_edge_model -> SetProperty, Hide = 0
C_leading_edge_model -> SetProperty, Hide=1
;Displays the value of the leading edge near the leading edge line
C_leading_edge_string = obj_new("idlgrtext", strings = leading_edge_string_value,color=[255,255,0], $
font=Courier12pt, /onglass, locations = [+xsize/20.,(ysize/7.) + 10],fill_background=1,fill_color=color_fill, hide=0)
C_leading_edge_model -> Add, C_leading_edge_string
C_clock_angle_model -> Add, C_leading_edge_model
;Main_Model -> Add, C_clock_angle_model
; ECDP

C_image_lineplot_model = OBJ_NEW('IDLgrModel')
x = [0,1]
y = [0,1]
xrange = [0 , 1]
yrange = [0 , 1]
C_image_lineplot = OBJ_NEW('IDLgrPlot', x, y, color=[255,0,0])
position_C_image_lineplot = [10, 10, 502, 100]
C_image_xs = swpc_cat_FSC_Normalize(xrange, Position=[position_C_image_lineplot[0], position_C_image_lineplot[2]])
C_image_ys = swpc_cat_FSC_Normalize(yrange, Position=[position_C_image_lineplot[1], position_C_image_lineplot[3]])
C_image_lineplot->SetProperty, XCoord_Conv=C_image_xs, YCoord_Conv=C_image_ys
C_image_lineplot_model -> Add, C_image_lineplot
C_image_view -> Add, C_image_lineplot_model

;if show_image_line_plot eq 0 then begin
C_image_lineplot_model -> SetProperty, hide = 1
C_clock_angle_model -> SetProperty, hide = 1
;endif

R_both_views = OBJ_NEW('IDLgrViewGroup')
R_image_view = OBJ_NEW('IDLgrView', Color=background_color, Viewplane_Rect=[0.0,0.0,Xsize,Ysize])
R_both_views -> add, R_image_view
R_image_model = OBJ_NEW('IDLgrModel') 
R_image_view -> Add, R_image_model
R_image_model->SetProperty, hide = 0

R_clock_angle_model = OBJ_NEW('IDLgrModel') 
R_clock_angle_model_location = [0.,0.]
R_clock_angle_model -> translate,255,255,0.0
R_clock_angle_model_location = [255,255]
R_clock_angle_model->GetProperty, transform = R_clock_angle_model_original_transform
; Marker for the 'clock angle' line
R_clock_angle_marker_line_x = [0.,0.]
R_clock_angle_marker_line_y = [0.,ysize/2.5]
R_clock_angle_marker_line = Obj_New("IDLgrPolyline", R_clock_angle_marker_line_x , R_clock_angle_marker_line_y, $
name = 'R_clock_angle', Color=[255,0,0], thick=1)
R_clock_angle_model -> Add, R_clock_angle_marker_line
R_image_view -> Add, R_clock_angle_model


R_image_lineplot_model = OBJ_NEW('IDLgrModel')
x = [0,1]
y = [0,1]
xrange = [0 , 1]
yrange = [0 , 1]
R_image_lineplot = OBJ_NEW('IDLgrPlot', x, y, color=[255,0,0])
position_R_image_lineplot = [10, 10, 502, 100]
R_image_xs = swpc_cat_FSC_Normalize(xrange, Position=[position_R_image_lineplot[0], position_R_image_lineplot[2]])
R_image_ys = swpc_cat_FSC_Normalize(yrange, Position=[position_R_image_lineplot[1], position_R_image_lineplot[3]])
R_image_lineplot->SetProperty, XCoord_Conv=R_image_xs, YCoord_Conv=R_image_ys
R_image_lineplot_model -> Add, R_image_lineplot
R_image_view -> Add, R_image_lineplot_model

;if show_image_line_plot eq 0 then begin
R_image_lineplot_model -> SetProperty, hide = 1
R_clock_angle_model -> SetProperty, hide = 1
;endif



images_timeline_view = OBJ_NEW('IDLgrView', Color=background_color, Viewplane_Rect=[0.0,0.0,1.0,1.05])
images_timeline_model = OBJ_NEW('IDLgrModel') 
images_timeline_view -> Add, images_timeline_model


; ECDP ()
;app_data_path = getenv('APPDATA')
;swpc_cat_preferences_file = app_data_path + sep + 'swpc_cat_prefs'

;openr,lun,swpc_cat_preferences_file, error = err, /get_lun
;if err eq 0 then begin
;   readf,lun, x_pos , y_pos
;   free_lun,lun
;endif else begin
;   x_pos = 20 & y_pos = 75
;endelse

;openr,lun,preferences_folder+sep+'last_window_position.txt', error = err, /get_lun
;if err eq 0 then begin
;   readf,lun, x_pos , y_pos
;   free_lun,lun
;endif else begin
   x_pos = 20 & y_pos = 75
;endelse

tlb_position = [x_pos,y_pos]
; ECDP



; ECDP
tlb = Widget_Base(Title='CAT (CME Analysis Tool)', row=1, xoffset = x_pos , yoffset = y_pos, /tlb_move_events, MBar=menubase)
; 

if i_center_on_main_monitor eq 1 then cgCenterTLB, tlb, -200, -200
;Widget_Control, tlb, /Realize
; ECDP

horizontal_base=widget_base(tlb,row=1)
;vertical_base=widget_base(horizontal_base,row=1)
vertical_base=widget_base(horizontal_base,column=1)
swpc_cat_base = widget_base(vertical_base,row=1)
L_view_base = widget_base(swpc_cat_base,column=1)
C_view_base = widget_base(swpc_cat_base,column=1)
R_view_base = widget_base(swpc_cat_base,column=1)
swpc_cat_base2 = widget_base(vertical_base,row=1)
horizontal_base2 = widget_base(vertical_base,row=1,/align_left)
vertical_base2 = widget_base(horizontal_base2,column=1)

if (compact eq 0) then begin
  c_mode_name = WIDGET_LABEL(C_view_base, value = '', /align_center, frame=0, xsize = 150, ysize = 20, font='Helvetica*BOLD*18')
; so that they align
  l_empty_mode_name = WIDGET_LABEL(L_view_base, value = '', /align_center, frame=0, xsize = 100, ysize = 20)
  r_empty_mode_name = WIDGET_LABEL(r_view_base, value = '', /align_center, frame=0, xsize = 100, ysize = 20)
;mode_name = WIDGET_LABEL(tlb, value = ' Stuff ',frame=1, /align_center, font='Helvetica*ITALIC*16', xoffset = 255)
endif else begin
  c_mode_name = 0
  l_empty_mode_name = 0
  r_empty_mode_name = 0
endelse  

draw_L = Widget_Draw(L_view_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, $
   Expose_Events=1, Event_Pro='swpc_cat_Window_button_click_events', Button_Events=1,/motion_events,render = 0) 

;draw_C = Widget_Draw(C_view_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, $
;   Expose_Events=1, Event_Pro='swpc_cat_Window_button_click_events', Button_Events=1,/motion_events,render = 1)
draw_C = Widget_Draw(C_view_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, $
   Expose_Events=1, Event_Pro='swpc_cat_Window_button_click_events', Button_Events=1,/motion_events,render = 0)
   
draw_R = Widget_Draw(R_view_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, $
   Expose_Events=1, Event_Pro='swpc_cat_Window_button_click_events', Button_Events=1,/motion_events,render = 0) 
 
;draw_L = Widget_Draw(L_view_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, $
;   Expose_Events=1, Event_Pro='swpc_cat_Window_button_click_events', Button_Events=1,motion_events=0,render = 0) 

;draw_C = Widget_Draw(C_view_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, $
;   Expose_Events=1, Event_Pro='swpc_cat_Window_button_click_events', Button_Events=1,motion_events=0,render = 0)
;   
;draw_R = Widget_Draw(R_view_base, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, $
;   Expose_Events=1, Event_Pro='swpc_cat_Window_button_click_events', Button_Events=1,motion_events=0,render = 0)   


if show_copy_windows eq 1 then begin

draw_L_copy = Widget_Draw(swpc_cat_base2, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, render = 0)  
draw_C_copy = Widget_Draw(swpc_cat_base2, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, render = 0)
draw_R_copy = Widget_Draw(swpc_cat_base2, XSize=Xsize, YSize=Ysize, Graphics_Level=2, retain=retain, render = 0)

endif




L_contextBase = WIDGET_BASE(draw_L, /CONTEXT_MENU, UNAME = 'L_drawContext',sensitive=0)
L_widget_outline_matches_image = Widget_Button(L_contextBase, Value='CME Matches Image',Event_Pro='swpc_cat_cme_match_or_unmatch')
widget_image_L = Widget_button(L_contextBase,Value='Image Differencing',/Menu)
difference_image = Widget_Button(widget_image_L, Value='None',Uvalue='L Normal',Event_Pro='swpc_cat_diff_format')
difference_image = Widget_Button(widget_image_L, Value='Running Difference',Uvalue='L Diff: Running',Event_Pro='swpc_cat_diff_format')
difference_image = Widget_Button(widget_image_L, Value='Set Current Image as Background',Uvalue='L Diff: Set Current Image as Background',Event_Pro='swpc_cat_diff_format')
L_widget_representative_image = Widget_Button(L_contextBase, Value='Set As Representative Image',Event_Pro='swpc_cat_define_representative_image')
L_widget_remove_this_image = Widget_Button(L_contextBase, Value='Remove This Image',Event_Pro='swpc_cat_remove_this_image',sensitive=0)
;L_widget_show_line_plot = Widget_Button(L_contextBase, Value='Show line plot',Event_Pro='swpc_cat_show_line_plot')

C_contextBase = WIDGET_BASE(draw_C, /CONTEXT_MENU, UNAME = 'C_drawContext',sensitive=0)
;if (which_mode eq 0) then begin

; ECDP
C_widget_pos_outline_matches_image = Widget_Button(C_contextBase, Value='POS CME Matches Image',Event_Pro='swpc_cat_pos_cme_match_or_unmatch')
; ECDP
C_widget_outline_matches_image = Widget_Button(C_contextBase, Value='CME Matches Image',Event_Pro='swpc_cat_cme_match_or_unmatch')
;endif else begin
;C_widget_outline_matches_image = Widget_Button(C_contextBase, Value='CME Matches Image',Event_Pro='swpc_cat_cme_match_or_unmatch', sensitive=0)
;endelse  
widget_image_C = Widget_button(C_contextBase,Value='Image Differencing',/Menu)
difference_image = Widget_Button(widget_image_C, Value='None',Uvalue='C Normal',Event_Pro='swpc_cat_diff_format')
difference_image = Widget_Button(widget_image_C, Value='Running Difference',Uvalue='C Diff: Running',Event_Pro='swpc_cat_diff_format')
difference_image = Widget_Button(widget_image_C, Value='Set Current Image as Background',Uvalue='C Diff: Set Current Image as Background',Event_Pro='swpc_cat_diff_format')
widget_show_C2_or_C3 = Widget_Button(C_contextBase, Value='Show LASCO C2',Event_Pro='swpc_cat_Show_C2_or_C3')
C_widget_representative_image = Widget_Button(C_contextBase, Value='Set As Representative Image',Event_Pro='swpc_cat_define_representative_image')
C_widget_remove_this_image = Widget_Button(C_contextBase, Value='Remove This Image',Event_Pro='swpc_cat_remove_this_image',sensitive=0)
;C_widget_show_line_plot = Widget_Button(C_contextBase, Value='Show line plot',Event_Pro='swpc_cat_show_line_plot')


R_contextBase = WIDGET_BASE(draw_R, /CONTEXT_MENU, UNAME = 'R_drawContext',sensitive=0)
R_widget_outline_matches_image = Widget_Button(R_contextBase, Value='CME Matches Image',Event_Pro='swpc_cat_cme_match_or_unmatch')
widget_image_R = Widget_button(R_contextBase,Value='Image Differencing',/Menu)
difference_image = Widget_Button(widget_image_R, Value='None',Uvalue='R Normal',Event_Pro='swpc_cat_diff_format')
difference_image = Widget_Button(widget_image_R, Value='Running Difference',Uvalue='R Diff: Running',Event_Pro='swpc_cat_diff_format')
difference_image = Widget_Button(widget_image_R, Value='Set Current Image as Background',Uvalue='R Diff: Set Current Image as Background',Event_Pro='swpc_cat_diff_format')
R_widget_representative_image = Widget_Button(R_contextBase, Value='Set As Representative Image',Event_Pro='swpc_cat_define_representative_image')
R_widget_remove_this_image = Widget_Button(R_contextBase, Value='Remove This Image',Event_Pro='swpc_cat_remove_this_image',sensitive=0)
;R_widget_show_line_plot = Widget_Button(R_contextBase, Value='Show line plot',Event_Pro='swpc_cat_show_line_plot')

; Plot A is the 'available images timeline'.....

xsize_timeline = (3.*Xsize)+20
ysize_timeline = Ysize/8
if compact eq 1 then ysize_timeline = Ysize/10

draw_available_images_timeline = Widget_Draw(swpc_cat_base2, XSize=xsize_timeline, YSize=ysize_timeline, Graphics_Level=2, retain=retain, $
   Expose_Events=1, Event_Pro='swpc_cat_images_timeline_window_button_click_events', Button_Events=1,/motion_events,render = 0)

x=[0.,1.]
y=[0.,1.]
position_timeline = [0.05, 0.3, 0.95, 0.9]
;if compact eq 1 then position_timeline = [0.05, 0.36, 0.95, 0.98]
if compact eq 1 then position_timeline = [0.05, 0.42, 0.95, 0.98]
xtitle_A = ''
ytitle_A = ''
plottitle_A = ''
psym = 0
symsize = 1.
exact = [0,0]
thisSymbol = Obj_New('IDLgrSymbol', psym, Color=[0, 255, 255])
;helvetica10pt = Obj_New('IDLgrFont', 'Helvetica', Size=10)
;helvetica12pt = Obj_New('IDLgrFont', 'Helvetica', Size=12)
;helvetica14pt = Obj_New('IDLgrFont', 'Helvetica', Size=14)
;courier10pt = obj_new('IDLgrFont','Courier', Size=10)
;courier12pt = obj_new('IDLgrFont','Courier', Size=12)
;courier14pt = obj_new('IDLgrFont','Courier', Size=14)
;monospace_symbol14pt = obj_new('IDLgrFont','Monospace Symbol', Size=14)

helvetica10pt = Obj_New('IDLgrFont', 'Helvetica', Size=12)
helvetica12pt = Obj_New('IDLgrFont', 'Helvetica', Size=14)
helvetica14pt = Obj_New('IDLgrFont', 'Helvetica', Size=16)
courier10pt = obj_new('IDLgrFont','Courier', Size=12)
courier12pt = obj_new('IDLgrFont','Courier', Size=14)
courier14pt = obj_new('IDLgrFont','Courier', Size=16)
monospace_symbol14pt = obj_new('IDLgrFont','Monospace Symbol', Size=16)

thisSymbol_L = obj_new("IDLgrsymbol",data=6)
thisSymbol_C = obj_new("IDLgrsymbol",data=6)
thisSymbol_C2 = obj_new("IDLgrsymbol",data=6)
thisSymbol_R = obj_new("IDLgrsymbol",data=6)

L_plot = Obj_New("IDLgrPLOT", x, y, Symbol=thisSymbol_L, Thick=1 , linestyle = 6)
C_plot = Obj_New("IDLgrPLOT", x, y, Symbol=thisSymbol_C, Thick=1, linestyle = 6)
R_plot = Obj_New("IDLgrPLOT", x, y, Symbol=thisSymbol_R, Thick=1, linestyle = 6)
C2_plot = Obj_New("IDLgrPLOT", x, y, Symbol=thisSymbol_C2, Thick=1, linestyle = 6)

xaxis_images_timeline = Obj_New("IDLgrAxis", 0, Color=[35,35,35], Ticklen=0.01, Minor=0, $
                                Location=[1000, position_timeline[1] ,0], Exact=exact[0],/use_text_color)
xaxis_images_timeline->GetProperty, Ticktext=xAxisText
xAxisText->SetProperty, Recompute_Dimensions=2    ; ensure that when axes are changed, the text size is recomputed properly.
xAxisText->SetProperty, Font=helvetica10pt
xAxisText->SetProperty, color=[255,255,255]

line_data = fltarr(2,2)
animation_start_time_marker = obj_new("idlgrpolyline",line_data,color=[0,255,0],name='animation_start_time_marker')
animation_end_time_marker = obj_new("idlgrpolyline",line_data,color=[255,0,0],name='animation_end_time_marker')
animation_start_time_marker_handle_data = fltarr(2,4)
animation_start_time_marker_handle = obj_new("idlgrpolygon",animation_start_time_marker_handle_data,color=[0,255,0],name='animation_start_time_marker_handle')
animation_end_time_marker_handle_data = fltarr(2,4)
animation_end_time_marker_handle = obj_new("idlgrpolygon",animation_end_time_marker_handle_data,color=[255,0,0],name='animation_end_time_marker_handle')
animation_current_time_marker = obj_new("idlgrpolyline",line_data,color=[255,255,0],name='animation_current_time_marker')

images_timeline_model->Add, animation_start_time_marker
images_timeline_model->Add, animation_end_time_marker
images_timeline_model->Add, animation_start_time_marker_handle
images_timeline_model->Add, animation_end_time_marker_handle
images_timeline_model->Add, animation_current_time_marker

images_timeline_model->Add, L_plot
images_timeline_model->Add, C_plot
images_timeline_model->Add, C2_plot
images_timeline_model->Add, R_plot
images_timeline_model->Add, xaxis_images_timeline

timeline_contextBase = WIDGET_BASE(draw_available_images_timeline, /CONTEXT_MENU, UNAME = 'timeline_contextBase',sensitive=0)
timeline_set_animation_start = Widget_Button(timeline_contextBase, Value='Set Start',Event_Pro='swpc_cat_set_animation_start')
timeline_set_animation_end = Widget_Button(timeline_contextBase, Value='Set End',Event_Pro='swpc_cat_set_animation_end')

; figure out the current time so we can provide some decent defaults for the coming
; date/time text boxes.............

current_julian_time = systime(/julian,/utc)
CALDAT, current_julian_time, current_Month, current_Day, current_Year, current_Hour, current_Minute

julian_this_hour = JULDAY(current_Month, current_Day, current_Year, current_Hour) 
julian_this_hour_ending = julian_this_hour + (1.d/24.d)
julian_this_hour_starting = julian_this_hour_ending - 1.d

CALDAT, julian_this_hour_ending, current_Month_ending, current_Day_ending, current_Year_ending, current_Hour_ending
CALDAT, julian_this_hour_starting, current_Month_starting, current_Day_starting, current_Year_starting, current_Hour_starting

current_year_starting = strcompress(string(current_year_starting),/remove_all)
current_month_starting = strcompress(string(current_month_starting),/remove_all)
current_day_starting = strcompress(string(current_day_starting),/remove_all)
current_hour_starting = strcompress(string(current_hour_starting),/remove_all)
current_minute_starting = ""

current_year_ending = strcompress(string(current_year_ending),/remove_all)
current_month_ending = strcompress(string(current_month_ending),/remove_all)
current_day_ending = strcompress(string(current_day_ending),/remove_all)
current_hour_ending = strcompress(string(current_hour_ending),/remove_all)
current_minute_ending = ""

load_images_controls_outer_base = widget_base(horizontal_base2,column=1,frame=compact_frame)

if compact eq 0 then Result = WIDGET_LABEL(load_images_controls_outer_base, /ALIGN_CENTER, value = ' START / END TIMES ',frame=2)

labelbase = Widget_Base(load_images_controls_outer_base, Row=6)
row1base = Widget_Base(labelbase, Row=1)
row2base = Widget_Base(labelbase, Row=1)
row3base = Widget_Base(labelbase, Row=1)
row4base = Widget_Base(labelbase, Row=1)

start_label = Widget_Label(row1base, Value='Start [Y M D H M]')

start_yearID = swpc_cat_FSC_InputField(row2base, Value=current_Year_starting, /IntegerValue, digits = 4 , xsize = 4 , title='',event_pro = 'swpc_cat_datebox_event')
start_monthID = swpc_cat_FSC_InputField(row2base, Value=current_Month_starting, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')
start_dayID = swpc_cat_FSC_InputField(row2base, Value=current_Day_starting, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')
start_hourID = swpc_cat_FSC_InputField(row2base, Value=current_Hour_starting, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')
start_minuteID = swpc_cat_FSC_InputField(row2base, Value=current_Minute_starting, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')

end_label = Widget_Label(row3base, Value='End  [Y M D H M]')

widget_plus_12h_button = Widget_Button(row3base, Value='+12h',Event_Pro='swpc_cat_make_end_date_Xhours_after_start')
widget_plus_24h_button = Widget_Button(row3base, Value='+24h',Event_Pro='swpc_cat_make_end_date_Xhours_after_start')

end_yearID = swpc_cat_FSC_InputField(row4base, Value=current_Year_ending, /IntegerValue, digits = 4 , xsize = 4 , title='',event_pro = 'swpc_cat_datebox_event')
end_monthID = swpc_cat_FSC_InputField(row4base, Value=current_Month_ending, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')
end_dayID = swpc_cat_FSC_InputField(row4base, Value=current_Day_ending, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')
end_hourID = swpc_cat_FSC_InputField(row4base, Value=current_Hour_ending, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')
end_minuteID = swpc_cat_FSC_InputField(row4base, Value=current_Minute_ending, /IntegerValue, digits = 2 , xsize = 2 , title='',event_pro = 'swpc_cat_datebox_event')

buttonBase = Widget_Base(load_images_controls_outer_base, Row=1, /align_center)
acceptID = Widget_Button(buttonBase, Value='Load Images',event_pro = 'swpc_cat_set_start_and_end_time')

; Mode button below load button
filer3 = widget_base(load_images_controls_outer_base, row=1, /align_center, yoffset = 30)
;filer3 = Widget_Button(buttonbase, Value='Mode', /MENU)
which_mode = ''
widget_select_mode2 = Widget_Button(filer3, Value='Select Mode',Event_Pro='select_mode',/Menu,sensitive = 1, xoffset = 30, yoffset =30, font='Helvetica*ITALIC*15', xsize = 100, ysize = 30)
select_mode2 = Widget_Button(widget_select_mode2, Value='NORMAL',Uvalue='NORMAL',Event_Pro='select_mode')
select_mode2 = Widget_Button(widget_select_mode2, Value='Plane of Sky (POS)',Uvalue='POS',Event_Pro='select_mode')
; 


date_array = [current_year_starting,current_month_starting,current_day_starting,current_hour_starting,current_minute_starting, $
              current_year_ending,current_month_ending,current_day_ending,current_hour_ending,current_minute_ending]       


animation_controls_outer_base = widget_base(horizontal_base2,column=1,frame=compact_frame,sensitive=0)
if compact eq 0 then Result = WIDGET_LABEL(animation_controls_outer_base, /ALIGN_CENTER, value = ' ANIMATION CONTROLS ',frame=2)

values = ['L', 'C', 'R'] 
select_anims_base = WIDGET_BASE(animation_controls_outer_base,/ROW,/align_center) 
widget_which_window_to_animate = CW_BGROUP(select_anims_base, values, /row, /EXCLUSIVE, /FRAME, SET_VALUE=1, event_funct='swpc_cat_which_window_to_animate')


animation_controls_base=widget_base(animation_controls_outer_base,column=1)
playbuttonBase = Widget_Base(animation_controls_base, Row=1, /align_center)
play_pause_button = Widget_Button(playbuttonBase, Value='Play ',Event_Pro='swpc_cat_xmovie_event')
animation_controls_base2=widget_base(animation_controls_base,row=1)
widget_animation_delay = WIDGET_SLIDER(animation_controls_base2, Max=100, Min=0, /Suppress_Value, $
   Title='Speed', Value=50, Event_Pro='swpc_cat_XMovie_Delay_Events')

widget_anim_back_and_forth = cw_bgroup(animation_controls_base2, 'Altern8',/nonexclusive,Event_Funct='swpc_cat_xmovie_back_and_forth')
;animation_set_start_end_base=widget_base(animation_controls_outer_base,row=1,/align_center)
;set_anim_start = Widget_Button(animation_set_start_end_base, Value='Set Start',Event_Pro='set_animation_start')
;set_anim_end = Widget_Button(animation_set_start_end_base, Value='Set End',Event_Pro='set_animation_end')


slider_size = 200

image_controls_base = widget_base(horizontal_base2,/column,frame=compact_frame,sensitive=0)
if compact eq 0 then Result = WIDGET_LABEL(image_controls_base, /ALIGN_CENTER, value = ' IMAGE ADJUST ',frame=2)
image_controls_tab = widget_tab(image_controls_base,/align_center, event_pro = 'swpc_cat_change_image_tab')

L_controls = widget_base(image_controls_tab,column=1,title='     L     ',/align_left)

   L_widget_botSlider = Widget_Slider(L_controls, Value=0, Min=0, $
      Max=100, XSize=slider_size,Event_Pro='swpc_cat_XColors_Bottom_Slider', $
      Title='Stretch Bottom', Drag=1,/Suppress_Value,sensitive = 0, font='Helvetica*ITALIC*12')
   L_widget_topSlider = Widget_Slider(L_controls, Value=100, Min=0, $
      Max=100, XSize=slider_size, Event_Pro='swpc_cat_XColors_Top_Slider', $
      Title='Stretch Top', Drag=1, /Suppress_Value,sensitive = 0, font='Helvetica*ITALIC*12')
   L_widget_gammaSlider = Widget_Slider(L_controls, Value=10, Min=1, Max=30, $
      Drag=1, XSize=slider_size, /Suppress_Value, Event_Pro='swpc_cat_XColors_Gamma_Slider', $
      Title='Gamma Correction',sensitive = 0, font='Helvetica*ITALIC*12')
   L_widget_saturationSlider = Widget_Slider(L_controls, Value=200, Min=1, Max=4000, $
      Drag=1, XSize=slider_size, Event_Pro='swpc_cat_change_the_image_saturation', $
      Title='Image Saturation Value', $
      /suppress_value,sensitive = 0, font='Helvetica*ITALIC*12')
      
L_reset_controls = widget_base(L_controls,row=1,/align_center)      
L_reset = Widget_Button(L_reset_controls, Value='Reset', Event_Pro='swpc_cat_Reset_L', font='Helvetica*ITALIC*12')
      
L_copy_controls = widget_base(L_controls,row=1,/align_center)      
copy_L_to_C = Widget_Button(L_copy_controls, Value='Copy to C ->', Event_Pro='swpc_cat_copy_L_to_C', font='Helvetica*ITALIC*12')
copy_L_to_R = Widget_Button(L_copy_controls, Value='Copy to R ->', Event_Pro='swpc_cat_copy_L_to_R', font='Helvetica*ITALIC*12')

C_controls = widget_base(image_controls_tab,column=1,title='     C     ',/align_center)

   C_widget_botSlider = Widget_Slider(C_controls, Value=0, Min=0, $
      Max=100, XSize=slider_size,Event_Pro='swpc_cat_XColors_Bottom_Slider', $
      Title='Stretch Bottom', Drag=1,/Suppress_Value,sensitive = 0, font='Helvetica*ITALIC*12')
   C_widget_topSlider = Widget_Slider(C_controls, Value=100, Min=0, $
      Max=100, XSize=slider_size, Event_Pro='swpc_cat_XColors_Top_Slider', $
      Title='Stretch Top', Drag=1, /Suppress_Value,sensitive = 0, font='Helvetica*ITALIC*12')
   C_widget_gammaSlider = Widget_Slider(C_controls, Value=10, Min=1, Max=30, $
      Drag=1, XSize=slider_size, /Suppress_Value, Event_Pro='swpc_cat_XColors_Gamma_Slider', $
      Title='Gamma Correction',sensitive = 0, font='Helvetica*ITALIC*12')
   C_widget_saturationSlider = Widget_Slider(C_controls, Value=50, Min=1, Max=4000, $
      Drag=1, XSize=slider_size, Event_Pro='swpc_cat_change_the_image_saturation', $
      Title='image saturation value', $
      /suppress_value,sensitive = 0, font='Helvetica*ITALIC*12')
      
C_reset_controls = widget_base(C_controls,row=1,/align_center)      
C_reset = Widget_Button(C_reset_controls, Value='Reset', Event_Pro='swpc_cat_Reset_C', font='Helvetica*ITALIC*12')

C_copy_controls = widget_base(C_controls,row=1,/align_center)      
copy_C_to_L = Widget_Button(C_copy_controls, Value='<- Copy to L',Event_Pro='swpc_cat_copy_C_to_L', font='Helvetica*ITALIC*12')
copy_C_to_R = Widget_Button(C_copy_controls, Value='Copy to R ->',Event_Pro='swpc_cat_copy_C_to_R', font='Helvetica*ITALIC*12')


R_controls = widget_base(image_controls_tab,column=1,title='     R     ',/align_right)
         
   R_widget_botSlider = Widget_Slider(R_controls, Value=0, Min=0, $
      Max=100, XSize=slider_size,Event_Pro='swpc_cat_XColors_Bottom_Slider', $
      Title='Stretch Bottom', Drag=1,/Suppress_Value,sensitive = 0, font='Helvetica*ITALIC*12')
   R_widget_topSlider = Widget_Slider(R_controls, Value=100, Min=0, $
      Max=100, XSize=slider_size, Event_Pro='swpc_cat_XColors_Top_Slider', $
      Title='Stretch Top', Drag=1, /Suppress_Value,sensitive = 0, font='Helvetica*ITALIC*12')
   R_widget_gammaSlider = Widget_Slider(R_controls, Value=10, Min=1, Max=30, $
      Drag=1, XSize=slider_size, Suppress_Value=1, Event_Pro='swpc_cat_XColors_Gamma_Slider', $
      Title='Gamma Correction',sensitive = 0, font='Helvetica*ITALIC*12')
   R_widget_saturationSlider = Widget_Slider(R_controls, Value=200, Min=1, Max=4000, $
      Drag=1, XSize=slider_size, Event_Pro='swpc_cat_change_the_image_saturation', $
      Title='image saturation value', $
      suppress_value=1,sensitive = 0, font='Helvetica*ITALIC*12')
      
R_reset_controls = widget_base(R_controls,row=1,/align_center)      
R_reset = Widget_Button(R_reset_controls, Value='Reset', Event_Pro='swpc_cat_Reset_R', font='Helvetica*ITALIC*12')

R_copy_controls = widget_base(R_controls,row=1,/align_center)       
copy_R_to_L = Widget_Button(R_copy_controls, Value='<- Copy to L',Event_Pro='swpc_cat_copy_R_to_L', font='Helvetica*ITALIC*12')
copy_R_to_C = Widget_Button(R_copy_controls, Value='<- Copy to C',Event_Pro='swpc_cat_copy_R_to_C', font='Helvetica*ITALIC*12')

if allow_show_3D_view eq 1 then begin
widget_show_3d_view = Widget_Button(image_controls_base, Value='Show 3D View',Event_Pro='swpc_cat_show_3d_view')
endif
Three_D_view_tlb_ID = 0L
Three_D_view_event_widget = 0L



cme_controls_base = widget_base(horizontal_base2,/column,frame=compact_frame,sensitive=0)
if compact eq 0 then Result = WIDGET_LABEL(cme_controls_base, /ALIGN_CENTER, value = ' CME CONTROLS ',frame=2)
cme_parameters_base = widget_base(cme_controls_base,column=1,frame=1,/align_center)
                               
widget_latitude_slider = widget_slider(cme_parameters_base, Value=0, $
                              event_pro='swpc_cat_change_latitude', $ 
                              minimum = -900, maximum = 900, scroll = 1, $
                              drag = 1, title = 'Latitude', $
                              suppress_value=1,sensitive=0, font='Helvetica*ITALIC*16')
                              
widget_longitude_slider = widget_slider(cme_parameters_base, Value=0, $
                                event_pro='swpc_cat_change_longitude', $
                                minimum = -1800, maximum = 1800, scroll = 1, $
                                drag = 1, title = 'Longitude', $
                                suppress_value=1,sensitive=0, font='Helvetica*ITALIC*16')                
                
                               
widget_angular_width_slider = widget_slider(cme_parameters_base, Value=90, $
                              event_pro='swpc_cat_change_lemniscate_angular_width', $ 
                              minimum = 20, maximum = 160, scroll = 1, $
                              drag = 1, title = 'Angular Width (cone)', $
                              suppress_value=1,sensitive=0, font='Helvetica*ITALIC*16')
                              
widget_radial_distance_slider = widget_slider(cme_parameters_base, Value=80, $
                                event_pro='swpc_cat_change_lemniscate_radial_distance', $
                                minimum = 10, maximum = 800, scroll = 1, $
                                drag = 1, title = 'Radial Distance (dist)', $
                                suppress_value=1,sensitive=0, font='Helvetica*ITALIC*16')

Ysize_cme_params = 120                                
if compact eq 1 then Ysize_cme_params = 90

draw_info_ID2 = Widget_Draw(cme_parameters_base, XSize=250, YSize=Ysize_cme_params, $
                Graphics_Level=2, retain=retain, render = 0)
      
if allow_Bernouli eq 1 then begin                
widget_which_style_lemniscate = cw_bgroup(cme_parameters_base, 'Bernouli',/nonexclusive,Event_Funct='swpc_cat_which_style_lemniscate')
endif
                
;   CME leading edge vs time plot....
plot_View = OBJ_NEW('IDLgrView', Color=background_color, Viewplane_Rect=[0.0,0.0,1.0,1.05])
plot_Model = OBJ_NEW('IDLgrModel') 
plot_View -> Add, plot_Model
  
plot_window_base = widget_base(horizontal_base2,/column,frame=compact_frame,sensitive=0)

; ECDP                        
;if compact eq 0 then time_plot_title_widget = WIDGET_LABEL(plot_window_base, /ALIGN_CENTER, value = ' CME LEADING EDGE (Rs) vs TIME (UT) ',frame=0)
if (compact eq 0) then begin 
  time_plot_title_widget = WIDGET_LABEL(plot_window_base, /ALIGN_CENTER, value = ' CME LEADING EDGE (Rs) vs TIME (UT) ',frame=0)
endif else begin
  time_plot_title_widget = 0 
endelse  
; ECDP


ysize_line_plot = 284.444
if compact eq 1 then ysize_line_plot = 234.444
draw_plot_window = Widget_Draw(plot_window_base, XSize=Xsize * 0.8, YSize=ysize_line_plot, Graphics_Level=2, retain=retain, $
   Expose_Events=1, Event_Pro='swpc_cat_plot_window_button_click_events', Button_Events=1,motion_events=0,render = 0)
;draw_plot_window = Widget_Draw(plot_window_base, XSize=Xsize * 0.8, YSize=ysize_line_plot, Graphics_Level=2, retain=retain, $
;   Expose_Events=1, Event_Pro='swpc_cat_plot_window_button_click_events', Button_Events=1,/motion_events,render = 0)


x=[0.,1.]
y=[0.,1.]
position_B = [0.1, 0.1, 0.95, 0.95]
xtitle_B = '' ;'Date/Time'
ytitle_B = '' ;'Leading Edge (Rs)'
psym = 0
symsize = 1.
exact = [0,0]

swpc_cat_create_plot_objects,position_B,xTitle_B,yTitle_B, $
                    LE_plot_matched_CMEs,thisPlot_B2,xAxis1_B,xAxis2_B,yAxis1_B,yAxis2_B, $
                    thisSymbol,symsize,exact, $
                    helvetica10pt,helvetica14pt,x,y, $
                    thisPlot2_c3, thisPlot2_c2, thisPlot2_stereo_a, thisPlot2_stereo_b
                    
thisPlot2_c3->SetProperty, color=color_c3
thisPlot2_c2->SetProperty, color=color_c2
thisPlot2_stereo_a->SetProperty, color=color_stereo_a
thisPlot2_stereo_b->SetProperty, color=color_stereo_b

plot_Model->Add, LE_plot_matched_CMEs
plot_Model->Add, thisPlot_B2
plot_Model->Add, thisPlot2_c3
plot_Model->Add, thisPlot2_c2
plot_Model->Add, thisPlot2_stereo_a
plot_Model->Add, thisPlot2_stereo_b
plot_Model->Add, xAxis1_B
plot_Model->Add, xAxis2_B
plot_Model->Add, yAxis1_B
plot_Model->Add, yAxis2_B


if allow_manual_data_point eq 1 then begin
plot_window_contextBase = WIDGET_BASE(draw_plot_window, /CONTEXT_MENU, UNAME = 'plot_window_contextBase',sensitive=1)
plot_window_manually_add_point = Widget_Button(plot_window_contextBase, Value='Manually Add Point',Event_Pro='swpc_cat_manually_add_or_remove_point')
endif else begin
plot_window_manually_add_point = 0 ; dummy to satisfy info structure
endelse


; end of CME leading edge vs time plot....

enlil_info_base = widget_base(horizontal_base2,/column,frame=compact_frame,sensitive=0)
if (compact eq 0) then Result = WIDGET_LABEL(enlil_info_base, /ALIGN_CENTER, value = ' ENLIL PARAMETERS ',frame=2)

y_enlil_info = 200
if (compact eq 1) then y_enlil_info = 170

;draw_enlil_info_window = Widget_Draw(enlil_info_base, XSize=Xsize * 0.4, YSize=y_enlil_info, Graphics_Level=2, retain=retain)
; ECDP, for some reason, this does not work with retain = 2, but does work with retain = 0 or 1.
draw_enlil_info_window = Widget_Draw(enlil_info_base, XSize=Xsize * 0.47, YSize=y_enlil_info, Graphics_Level=2, retain=retain, render=0)

if (compact eq 0) then begin
; ECDP
widget_calculate_pos_velocity = Widget_Button(enlil_info_base, Value='Calculate POS Velocity',Event_Pro='swpc_cat_calculate_pos_velocity',sensitive=0)
widget_use_emp_cone = Widget_Button(enlil_info_base, Value='Copy Empirical Angle to CME controls',Event_Pro='swpc_cat_set_opening_angle',sensitive=0)
; ECDP
widget_calculate_velocity = Widget_Button(enlil_info_base, Value='Calculate Velocity',Event_Pro='swpc_cat_calculate_velocity',sensitive=0)
widget_export_cone = Widget_Button(enlil_info_base, Value='Export Analysis',Event_Pro='swpc_cat_export_cme_analysis_data',sensitive=0)
widget_reset_cme_analysis = Widget_Button(enlil_info_base, Value='Reset Analysis',Event_Pro='swpc_cat_reset_cme_analysis',sensitive=1)

endif else begin
  
  enlil_info_buttons_base = widget_base(enlil_info_base, row = 1)
  widget_export_cone = Widget_Button(enlil_info_base, Value='Export Analysis',Event_Pro='swpc_cat_export_cme_analysis_data',sensitive=0)
  ; ECDP
  widget_calculate_pos_velocity = Widget_Button(enlil_info_buttons_base, Value='Calc POS Velo',Event_Pro='swpc_cat_calculate_pos_velocity',sensitive=0)
  widget_use_emp_cone = Widget_Button(enlil_info_buttons_base, Value='Copy Emp Angle',Event_Pro='swpc_cat_set_opening_angle',sensitive=0)
; ECDP
  widget_calculate_velocity = Widget_Button(enlil_info_buttons_base, Value='Calc Velo',Event_Pro='swpc_cat_calculate_velocity',sensitive=0)
;  widget_export_cone = Widget_Button(enlil_info_buttons_base, Value='Export Analysis',Event_Pro='swpc_cat_export_cme_analysis_data',sensitive=0)
  widget_reset_cme_analysis = 1
 ; widget_reset_cme_analysis = Widget_Button(enlil_info_base, Value='Reset Analysis',Event_Pro='swpc_cat_reset_cme_analysis',sensitive=1)
  
  
endelse


enlil_info_View = OBJ_NEW('IDLgrView', Color=background_color, Viewplane_Rect=[0.0,0.0,1.0,1.05])
;enlil_info_View = OBJ_NEW('IDLgrView', Color=[0,0,255], Viewplane_Rect=[0.0,0.0,1.0,1.05])
enlil_info_Model = OBJ_NEW('IDLgrModel') 
enlil_info_View -> Add, enlil_info_Model


L_widget_image_sequence_slider = widget_slider(L_view_base, Value=1, $
                            event_pro='swpc_cat_change_L_C_or_R_image_with_slider',  $
                            minimum = 1, maximum = 2, scroll = 1, drag = 1, $
                            sensitive=0, $
                            suppress_value=1) ;,Title='Scroll STEREO B')
                            
C_widget_image_sequence_slider = widget_slider(C_view_base, Value=1, $
                            event_pro='swpc_cat_change_L_C_or_R_image_with_slider',  $
                            minimum = 1, maximum = 2, scroll = 1, drag = 1, $
                            sensitive=0, $
                            suppress_value=1) ; ,Title='Scroll LASCO')
                            
R_widget_image_sequence_slider = widget_slider(R_view_base, Value=1, $
                            event_pro='swpc_cat_change_L_C_or_R_image_with_slider',  $
                            minimum = 1, maximum = 2, scroll = 1, drag = 1, $
                            sensitive=0, $
                            suppress_value=1) ; ,Title='Scroll STEREO A')
      


Widget_Control, tlb, /Realize
Widget_Control, draw_L, Get_Value=L_Window
Widget_Control, draw_C, Get_Value=C_Window
Widget_Control, draw_R, Get_Value=R_Window
if show_copy_windows eq 1 then begin
Widget_Control, draw_L_copy, Get_Value=L_Window_copy
Widget_Control, draw_C_copy, Get_Value=C_Window_copy
Widget_Control, draw_R_copy, Get_Value=R_Window_copy
endif else begin
L_Window_copy = OBJ_NEW( 'IDLgrBuffer' )
L_Window_copy -> SetProperty, DIMENSIONS = [512,512]
C_Window_copy = OBJ_NEW( 'IDLgrBuffer' )
C_Window_copy -> SetProperty, DIMENSIONS = [512,512]
R_Window_copy = OBJ_NEW( 'IDLgrBuffer' )
R_Window_copy -> SetProperty, DIMENSIONS = [512,512]
endelse
Widget_Control, draw_available_images_timeline, Get_Value=images_timeline_window
Widget_Control, draw_plot_window, Get_Value=plot_Window
Widget_Control, draw_enlil_info_window, Get_Value=enlil_info_Window
Widget_Control, draw_info_ID2, Get_Value=cme_info_Window


;
;  End of defining all of the GUI widget stuff.........
;
L_image_color_palette = OBJ_NEW('IDLgrPalette')
L_image_color_palette -> loadct,0
C_image_color_palette = OBJ_NEW('IDLgrPalette')
C_image_color_palette -> loadct,0
R_image_color_palette = OBJ_NEW('IDLgrPalette')
R_image_color_palette -> loadct,0

;  Process the input diff data array - which contains
;  all of the diff images in a single array.
;  Define an object array (all_diff_image_objects)
;  that will hold all of the images.....

files_directory = ''
filename = ''

L_coronagraph_image_object = OBJ_NEW('IDLgrImage')
L_coronagraph_image_object -> SetProperty, location = [0 , 0]
L_coronagraph_image_object -> setproperty, palette = L_image_color_palette
L_coronagraph_image_object -> setproperty, alpha_channel = 1.0
L_coronagraph_image_object -> setproperty, blend_function = [3,4]

C_coronagraph_image_object = OBJ_NEW('IDLgrImage')
C_coronagraph_image_object -> SetProperty, location = [0 , 0]
C_coronagraph_image_object -> setproperty, palette = C_image_color_palette
C_coronagraph_image_object -> setproperty, alpha_channel = 1.0
C_coronagraph_image_object -> setproperty, blend_function = [3,4]

R_coronagraph_image_object = OBJ_NEW('IDLgrImage')
R_coronagraph_image_object -> SetProperty, location = [0 , 0]
R_coronagraph_image_object -> setproperty, palette = R_image_color_palette
R_coronagraph_image_object -> setproperty, alpha_channel = 1.0
R_coronagraph_image_object -> setproperty, blend_function = [3,4]

L_border_image_object = OBJ_NEW('IDLgrImage')
L_border_image_object -> SetProperty, location = [0 , 0]
;L_border_image_object -> setproperty, alpha_channel = 0.3
L_border_image_object -> setproperty, blend_function = [3,4]
border_data = bytarr(4,512,512)
border_data[0,*,0:18] = background_color[0]
border_data[1,*,0:18] = background_color[1]
border_data[2,*,0:18] = background_color[2]
border_data[3,*,0:18] = 160
border_data[0,*,493:511] = background_color[0]
border_data[1,*,493:511] = background_color[1]
border_data[2,*,493:511] = background_color[2]
border_data[3,*,493:511] = 160
border_data[0,0:18,*] = background_color[0]
border_data[1,0:18,*] = background_color[1]
border_data[2,0:18,*] = background_color[2]
border_data[3,0:18,*] = 160
border_data[0,493:511,*] = background_color[0]
border_data[1,493:511,*] = background_color[1]
border_data[2,493:511,*] = background_color[2]
border_data[3,493:511,*] = 160
L_border_image_object -> setproperty, data = border_data

C_border_image_object = OBJ_NEW('IDLgrImage')
C_border_image_object -> SetProperty, location = [0 , 0]
C_border_image_object -> setproperty, blend_function = [3,4]
C_border_image_object -> setproperty, data = border_data

R_border_image_object = OBJ_NEW('IDLgrImage')
R_border_image_object -> SetProperty, location = [0 , 0]
R_border_image_object -> setproperty, blend_function = [3,4]
R_border_image_object -> setproperty, data = border_data

L_image_model -> Add, L_coronagraph_image_object
L_image_model -> Add, L_border_image_object
C_image_model -> Add, C_coronagraph_image_object
C_image_model -> Add, C_border_image_object
R_image_model -> Add, R_coronagraph_image_object
R_image_model -> Add, R_border_image_object

L_coronagraph_image_object->SetProperty, hide = 0
C_coronagraph_image_object->SetProperty, hide = 0
R_coronagraph_image_object->SetProperty, hide = 0

center_of_SunX = 100.
center_of_SunY = 100.

;  variables to keep track of the current image we are viewing
;  and the current radial distance for the cme outline

current_R_in_Solar_radii = 0.0

;  The calculated CME radial velocity in km/s....

radial_velocity = 0.
POS_velocity = 0.
;  The calculated time of the CME at the Enlil boundary (21.5 Rs)
;  as both a double precision Julian date and the yyyymmddhhmmss
;  array of integers

time_at_Enlil_boundary_Julian = 0.d
time_at_Enlil_boundary = intarr(6)

initial_cme_outline_associated_with_image_number = -1

L_difference_imaging = 'running'
C_difference_imaging = 'running'
C2_difference_imaging = 'running'
R_difference_imaging = 'running'

L_background_image_number = 0
C_background_image_number = 0
C2_background_image_number = 0
R_background_image_number = 0

diff_red_levels = intarr(256)
diff_green_levels = intarr(256)
diff_blue_levels = intarr(256)

palette_0 = OBJ_NEW('IDLgrPalette')    ; , rval, gval, bval)
palette_0 -> loadct,0
   wcolors = !D.Table_Size
   IF N_Elements(wcolors) EQ 1 THEN BEGIN
      IF wcolors[0] EQ 1 THEN wcolors = [!D.Table_Size, 0] ELSE wcolors = [wcolors, 0]
   ENDIF
   wcolors[0] = (wcolors[1] + wcolors[0]) < (!D.Table_Size - wcolors[1])


xmin = 0.d
xmax = 0.d
ymin = 0.
ymax = 0.
xtickinterval = 0.25


latitude_text = obj_new("idlgrtext", strings = 'Lat :', $
                color=[255,255,255], locations = [0.03,0.76],hide=1)
latitude_text -> setproperty, font = courier12pt

longitude_text = obj_new("idlgrtext", strings = 'Lon :', $
                color=[255,255,255], locations = [0.03,0.66],hide=1)
longitude_text -> setproperty, font = courier12pt
                
cone_angle_text = obj_new("idlgrtext", strings = '1/2 Angle :', $
                color=[255,255,255], locations = [0.03,0.56],hide=1)
cone_angle_text -> setproperty, font = courier12pt                
               
velocity_text = obj_new("idlgrtext", strings = 'Velocity  :', $
                color=[255,255,255], locations = [0.03,0.46],hide=1)
velocity_text -> setproperty, font = courier12pt                

; ECDP          
pos_velocity_text = obj_new("idlgrtext", strings = 'POS Velocity  :', $
                color=[255,255,255], locations = [0.03,0.36],hide=1)
pos_velocity_text -> setproperty, font = courier12pt       

emp_cone_angle_text = obj_new("idlgrtext", strings = 'Empirical Cone Angle  :', $
                color=[255,255,255], locations = [0.03,0.26],hide=1)
emp_cone_angle_text -> setproperty, font = courier12pt  
; 

             
T21_5_text = obj_new("idlgrtext", strings = 'T21.5 :', $
                color=[255,255,255], locations = [0.03,0.86],hide=1)
T21_5_text -> setproperty, font = courier12pt                 
                
T21_5_string = ''

;cone_tool_version = 'CME Analysis Tool (CAT) v1.0'
cone_tool_version = 'CME Analysis Tool (CAT) v1.5'

enlil_info_model -> add, latitude_text
enlil_info_model -> add, longitude_text
enlil_info_model -> add, cone_angle_text
enlil_info_model -> add, velocity_text
; ECDP
enlil_info_model -> add, pos_velocity_text
enlil_info_model -> add, emp_cone_angle_text
; ECDP
enlil_info_model -> add, T21_5_text

velocity_text_c3 = obj_new("idlgrtext", strings = '', $
                color=color_c3, locations = [0.03,0.16],hide=0, font = courier12pt)
velocity_text_c2 = obj_new("idlgrtext", strings = '', $
                color=color_c2, locations = [0.03,0.06],hide=0, font = courier12pt)
velocity_text_stereo_a = obj_new("idlgrtext", strings = '', $
                color=color_stereo_a, locations = [0.55,0.16],hide=0, font = courier12pt)
velocity_text_stereo_b = obj_new("idlgrtext", strings = '', $
                color=color_stereo_b, locations = [0.55,0.06],hide=0, font = courier12pt)
                
enlil_info_model -> add, velocity_text_c3
enlil_info_model -> add, velocity_text_c2
enlil_info_model -> add, velocity_text_stereo_a
enlil_info_model -> add, velocity_text_stereo_b

which_telescope = ''
L_telescope_field_of_view_Rs = 15.
C_telescope_field_of_view_Rs = 30.
C2_telescope_field_of_view_Rs = 6.
R_telescope_field_of_view_Rs = 15.
image_size_pixels = 512.
Solar_radius_km = 695500.

B_angle_degrees = -999.

L_image_saturation_value = 200.
C_image_saturation_value = 50.
R_image_saturation_value = 200.
palette_top_normal = 100.
palette_bottom_normal = 0.
palette_gamma_normal = 10.
palette_top_diff = 100.
palette_bottom_diff = 0.
palette_gamma_diff = 10.

L_ut_string_object = OBJ_NEW('IDLgrText','')
L_ut_string_object -> setproperty, location = [255, 2]
L_ut_string_object -> setproperty, alignment = 0.5
L_ut_string_object -> setproperty, color = [100,100,255]
L_ut_string_object -> setproperty, font = courier12pt
L_image_model -> Add, L_ut_string_object

C_ut_string_object = OBJ_NEW('IDLgrText','')
C_ut_string_object -> setproperty, location = [255, 2]
C_ut_string_object -> setproperty, alignment = 0.5
C_ut_string_object -> setproperty, color = [100,255,100]
C_ut_string_object -> setproperty, font = courier12pt
C_image_model -> Add, C_ut_string_object

R_ut_string_object = OBJ_NEW('IDLgrText','')
R_ut_string_object -> setproperty, location = [255, 2]
R_ut_string_object -> setproperty, alignment = 0.5
R_ut_string_object -> setproperty, color = [255,100,100]
R_ut_string_object -> setproperty, font = courier12pt
R_image_model -> Add, R_ut_string_object

L_title_object = OBJ_NEW('IDLgrText','STEREO B COR2')
L_title_object -> setproperty, location = [255, 495]
L_title_object -> setproperty, alignment = 0.5
L_title_object -> setproperty, color = [100,100,255]
L_title_object -> setproperty, font = courier12pt
L_image_model -> Add, L_title_object

C_title_object = OBJ_NEW('IDLgrText','SOHO LASCO C3')
C_title_object -> setproperty, location = [255, 495]
C_title_object -> setproperty, alignment = 0.5
C_title_object -> setproperty, color = [100,255,100]
C_title_object -> setproperty, font = courier12pt
C_image_model -> Add, C_title_object

R_title_object = OBJ_NEW('IDLgrText','STEREO A COR2')
R_title_object -> setproperty, location = [255, 495]
R_title_object -> setproperty, alignment = 0.5
R_title_object -> setproperty, color = [255,100,100]
R_title_object -> setproperty, font = courier12pt
R_image_model -> Add, R_title_object


; All the 3D view stuff.......


Sun_earth_distance = 215.0
draw_window_size = 512
angular_width_lemniscate = 90.
radial_distance_lemniscate = 8.
zclip_near = 100.
zclip_far = - zclip_near 
zclip = [zclip_near , zclip_far]

L_second_view = OBJ_NEW('IDLgrView', Color=[100,0,0], $
                      Viewplane_Rect=[0,0,2.0*L_telescope_field_of_view_Rs,2.0*L_telescope_field_of_view_Rs])
L_second_view -> setproperty, transparent = 1
L_second_view -> setproperty, projection = 2
L_second_view_copy = OBJ_NEW('IDLgrView', Color=[0,0,0], $
                      Viewplane_Rect=[0,0,2.0*L_telescope_field_of_view_Rs,2.0*L_telescope_field_of_view_Rs])
L_second_view_copy -> setproperty, projection = 2
L_both_views -> add, L_second_view
L_second_view -> setproperty, eye = Sun_earth_distance
L_second_view -> setproperty, zclip = zclip
L_second_view_copy -> setproperty, eye = Sun_earth_distance
L_second_view_copy -> setproperty, zclip = zclip
L_Fixed_Sun_Satellite_model = OBJ_NEW('IDLgrModel')
L_Fixed_Sun_Satellite_model_copy = OBJ_NEW('IDLgrModel')
L_Fixed_Sun_Satellite_model -> getproperty, transform = L_Fixed_Sun_Satellite_model_original_transform

L_coords_offset_model = OBJ_NEW('IDLgrModel')
L_coords_offset_model_copy = OBJ_NEW('IDLgrModel')
L_Fixed_Sun_Satellite_model -> Add, L_coords_offset_model
L_Fixed_Sun_Satellite_model_copy -> Add, L_coords_offset_model_copy

L_second_view -> Add, L_Fixed_Sun_Satellite_model
L_second_view_copy -> Add, L_Fixed_Sun_Satellite_model_copy
L_trackball = OBJ_NEW('Trackball', [draw_window_size / 2. , draw_window_size / 2.], draw_window_size)
L_Rotating_model_Transform = dblarr(4,4)

C_second_view = OBJ_NEW('IDLgrView', Color=[100,0,0], $
                      Viewplane_Rect=[0,0,2.0*C_telescope_field_of_view_Rs,2.0*C_telescope_field_of_view_Rs])
C_second_view -> setproperty, transparent = 1
C_second_view -> setproperty, projection = 2
C_second_view_copy = OBJ_NEW('IDLgrView', Color=[0,0,0], $
                      Viewplane_Rect=[0,0,2.0*C_telescope_field_of_view_Rs,2.0*C_telescope_field_of_view_Rs])
C_second_view_copy -> setproperty, projection = 2
C_both_views -> add, C_second_view
C_second_view -> setproperty, eye = Sun_earth_distance
C_second_view -> setproperty, zclip = zclip
C_second_view_copy -> setproperty, eye = Sun_earth_distance
C_second_view_copy -> setproperty, zclip = zclip
C_Fixed_Sun_Satellite_model = OBJ_NEW('IDLgrModel')
C_Fixed_Sun_Satellite_model_copy = OBJ_NEW('IDLgrModel')
C_Fixed_Sun_Satellite_model -> getproperty, transform = C_Fixed_Sun_Satellite_model_original_transform


C_coords_offset_model = OBJ_NEW('IDLgrModel')
C_coords_offset_model_copy = OBJ_NEW('IDLgrModel')
C_Fixed_Sun_Satellite_model -> Add, C_coords_offset_model
C_Fixed_Sun_Satellite_model_copy -> Add, C_coords_offset_model_copy

C_second_view -> Add, C_Fixed_Sun_Satellite_model
C_second_view_copy -> Add, C_Fixed_Sun_Satellite_model_copy
C_trackball = OBJ_NEW('Trackball', [draw_window_size / 2. , draw_window_size / 2.], draw_window_size)
C_Rotating_model_Transform = dblarr(4,4)

R_second_view = OBJ_NEW('IDLgrView', Color=[100,0,0], $
                      Viewplane_Rect=[0,0,2.0*R_telescope_field_of_view_Rs,2.0*R_telescope_field_of_view_Rs])
R_second_view -> setproperty, transparent = 1
R_second_view -> setproperty, projection = 2
R_second_view_copy = OBJ_NEW('IDLgrView', Color=[0,0,0], $
                      Viewplane_Rect=[0,0,2.0*R_telescope_field_of_view_Rs,2.0*R_telescope_field_of_view_Rs])
R_second_view_copy -> setproperty, projection = 2
R_both_views -> add, R_second_view
R_second_view -> setproperty, eye = Sun_earth_distance
R_second_view -> setproperty, zclip = zclip
R_second_view_copy -> setproperty, eye = Sun_earth_distance
R_second_view_copy -> setproperty, zclip = zclip
R_Fixed_Sun_Satellite_model = OBJ_NEW('IDLgrModel')
R_Fixed_Sun_Satellite_model_copy = OBJ_NEW('IDLgrModel')
R_Fixed_Sun_Satellite_model -> getproperty, transform = R_Fixed_Sun_Satellite_model_original_transform


R_coords_offset_model = OBJ_NEW('IDLgrModel')
R_coords_offset_model_copy = OBJ_NEW('IDLgrModel')
R_Fixed_Sun_Satellite_model -> Add, R_coords_offset_model
R_Fixed_Sun_Satellite_model_copy -> Add, R_coords_offset_model_copy

R_second_view -> Add, R_Fixed_Sun_Satellite_model
R_second_view_copy -> Add, R_Fixed_Sun_Satellite_model_copy
R_trackball = OBJ_NEW('Trackball', [draw_window_size / 2. , draw_window_size / 2.], draw_window_size)
R_Rotating_model_Transform = dblarr(4,4)


cme_info_view = OBJ_NEW('IDLgrView', Color=background_color, Viewplane_Rect=[0.0,0.0,200,150])

cme_info_model = OBJ_NEW('IDLgrModel') 

cme_info_view -> Add, cme_info_model
cme_info_model->setProperty, Hide = 0


shift_up = 4

info_title = 'CME Params'

info_title_object = OBJ_NEW('IDLgrText',info_title)
info_title_object -> setproperty, location = [100, 110 + shift_up]
info_title_object -> setproperty, alignment = 0.5
info_title_object -> setproperty, color = color_cme_orange
info_title_object -> setproperty, font = courier12pt
cme_info_model -> Add, info_title_object



lat_string = 'lat  :'
lon_string = 'lon  :'

lat_string_object = OBJ_NEW('IDLgrText',lat_string)
lat_string_object -> setproperty, location = [100, 83 + shift_up]
lat_string_object -> setproperty, alignment = 0.5
lat_string_object -> setproperty, color = color_cme_orange
lat_string_object -> setproperty, font = courier10pt
cme_info_model -> Add, lat_string_object

lon_string_object = OBJ_NEW('IDLgrText',lon_string)
lon_string_object -> setproperty, location = [100, 55 + shift_up]
lon_string_object -> setproperty, alignment = 0.5
lon_string_object -> setproperty, color = color_cme_orange
lon_string_object -> setproperty, font = courier10pt
cme_info_model -> Add, lon_string_object


angular_width_string = 'cone :'
radial_distance_string = 'dist :'

angular_width_string_object = OBJ_NEW('IDLgrText',angular_width_string)
angular_width_string_object -> setproperty, location = [100, 29 + shift_up]
angular_width_string_object -> setproperty, alignment = 0.5
angular_width_string_object -> setproperty, color = color_cme_orange
angular_width_string_object -> setproperty, font = courier10pt
cme_info_model -> Add, angular_width_string_object

radial_distance_string_object = OBJ_NEW('IDLgrText',radial_distance_string)
radial_distance_string_object -> setproperty, location = [100, 3 + shift_up]
radial_distance_string_object -> setproperty, alignment = 0.5
radial_distance_string_object -> setproperty, color = color_cme_orange
radial_distance_string_object -> setproperty, font = courier10pt
cme_info_model -> Add, radial_distance_string_object

initial_radial_distance_lemniscate = 8.0
initial_angular_width_lemniscate = 90.
lemnsicate_style = 0

swpc_cat_define_cme_lemniscate, initial_radial_distance_lemniscate, initial_angular_width_lemniscate,lemnsicate_style, fitted_cme_info


sun_radius = 1.

L_cme_model = obj_new("idlgrmodel")

L_cme_model->GetProperty, transform = L_initial_transform

;print, ' L_initial_transform ',L_initial_transform

L_cme_model_copy = obj_new("idlgrmodel")
L_cme_fitted_surf = obj_new("idlgrpolygon",sun_radius*fitted_cme_info.vertices $
                         ,polygons=fitted_cme_info.connections               $
                         ,color=[255,127,0],alpha_channel=0.5, hide = 1)
L_cme_fitted_surf_copy = obj_new("idlgrpolygon",sun_radius*fitted_cme_info.vertices $
                         ,polygons=fitted_cme_info.connections               $
                         ,color=[255,255,255],alpha_channel=1.0)
L_cme_model->add,L_cme_fitted_surf
L_cme_model_copy->add,L_cme_fitted_surf_copy

L_coords_offset_model -> Add, L_cme_model
L_coords_offset_model_copy -> Add, L_cme_model_copy

C_cme_model = obj_new("idlgrmodel")
C_cme_model_copy = obj_new("idlgrmodel")
C_cme_fitted_surf = obj_new("idlgrpolygon",sun_radius*fitted_cme_info.vertices $
                         ,polygons=fitted_cme_info.connections               $
                         ,color=[255,127,0],alpha_channel=0.5, hide = 1)
C_cme_fitted_surf_copy = obj_new("idlgrpolygon",sun_radius*fitted_cme_info.vertices $
                         ,polygons=fitted_cme_info.connections               $
                         ,color=[255,255,255],alpha_channel=1.0)
C_cme_model->add,C_cme_fitted_surf
C_cme_model_copy->add,C_cme_fitted_surf_copy

C_coords_offset_model -> Add, C_cme_model
C_coords_offset_model_copy -> Add, C_cme_model_copy

R_cme_model = obj_new("idlgrmodel")
R_cme_model_copy = obj_new("idlgrmodel")
R_cme_fitted_surf = obj_new("idlgrpolygon",sun_radius*fitted_cme_info.vertices $
                         ,polygons=fitted_cme_info.connections               $
                         ,color=[255,127,0],alpha_channel=0.5, hide = 1)
R_cme_fitted_surf_copy = obj_new("idlgrpolygon",sun_radius*fitted_cme_info.vertices $
                         ,polygons=fitted_cme_info.connections               $
                         ,color=[255,255,255],alpha_channel=1.0)
R_cme_model->add,R_cme_fitted_surf
R_cme_model_copy->add,R_cme_fitted_surf_copy

R_coords_offset_model -> Add, R_cme_model
R_coords_offset_model_copy -> Add, R_cme_model_copy



cone_X_axis_data = fltarr(3,2)
cone_X_axis_data[0,1] = radial_distance_lemniscate
L_cone_X_axis = OBJ_NEW('IDLgrpolyline',data = cone_X_axis_data,color=[255,100,100])
C_cone_X_axis = OBJ_NEW('IDLgrpolyline',data = cone_X_axis_data,color=[255,100,100])
R_cone_X_axis = OBJ_NEW('IDLgrpolyline',data = cone_X_axis_data,color=[255,100,100])
L_CME_model-> Add, L_cone_X_axis
C_CME_model-> Add, C_cone_X_axis
R_CME_model-> Add, R_cone_X_axis

cone_Y_axis_data = fltarr(3,2)
cone_Y_axis_data[1,1] = radial_distance_lemniscate
L_cone_Y_axis = OBJ_NEW('IDLgrpolyline',data = cone_Y_axis_data,color=[100,255,100])
C_cone_Y_axis = OBJ_NEW('IDLgrpolyline',data = cone_Y_axis_data,color=[100,255,100])
R_cone_Y_axis = OBJ_NEW('IDLgrpolyline',data = cone_Y_axis_data,color=[100,255,100])
L_CME_model-> Add, L_cone_Y_axis
C_CME_model-> Add, C_cone_Y_axis
R_CME_model-> Add, R_cone_Y_axis

cone_Z_axis_data = fltarr(3,2)
cone_Z_axis_data[2,1] = radial_distance_lemniscate * 3.
L_cone_Z_axis = OBJ_NEW('IDLgrpolyline',data = cone_Z_axis_data,color=[100,100,255])
C_cone_Z_axis = OBJ_NEW('IDLgrpolyline',data = cone_Z_axis_data,color=[100,100,255])
R_cone_Z_axis = OBJ_NEW('IDLgrpolyline',data = cone_Z_axis_data,color=[100,100,255])
L_CME_model-> Add, L_cone_Z_axis
C_CME_model-> Add, C_cone_Z_axis
R_CME_model-> Add, R_cone_Z_axis

L_cone_X_axis-> SetProperty, hide = 1
C_cone_X_axis-> SetProperty, hide = 1
R_cone_X_axis-> SetProperty, hide = 1
L_cone_Y_axis-> SetProperty, hide = 1
C_cone_Y_axis-> SetProperty, hide = 1
R_cone_Y_axis-> SetProperty, hide = 1
if show_cone_Z_axis eq 1 then begin
L_cone_Z_axis-> SetProperty, hide = 0
C_cone_Z_axis-> SetProperty, hide = 0
R_cone_Z_axis-> SetProperty, hide = 0
endif else begin
L_cone_Z_axis-> SetProperty, hide = 1
C_cone_Z_axis-> SetProperty, hide = 1
R_cone_Z_axis-> SetProperty, hide = 1
endelse

lem_data = fltarr(2,3)
lem_data[0,0] = 200
lem_data[0,1] = 250
lem_data[0,2] = 275
lem_data[1,0] = 200
lem_data[1,1] = 250
lem_data[1,2] = 300
L_cme_outline = obj_new("idlgrroi",data=lem_data, color = cme_outline_color,style=2,hide = 1)
L_image_model-> Add, L_cme_outline
C_cme_outline = obj_new("idlgrroi",data=lem_data, color = cme_outline_color,style=2,hide = 1)
C_image_model-> Add, C_cme_outline
R_cme_outline = obj_new("idlgrroi",data=lem_data, color = cme_outline_color,style=2,hide = 1)
R_image_model-> Add, R_cme_outline


L_cme_MATCH_outline = obj_new("idlgrroi",data=lem_data, color = color_stereo_B,style=2,hide = 1)
L_image_model-> Add, L_cme_MATCH_outline
C_cme_MATCH_outline = obj_new("idlgrroi",data=lem_data, color = color_c3,style=2,hide = 1)
C_image_model-> Add, C_cme_MATCH_outline
C2_cme_MATCH_outline = obj_new("idlgrroi",data=lem_data, color = color_c2,style=2,hide = 1)
C_image_model-> Add, C2_cme_MATCH_outline
R_cme_MATCH_outline = obj_new("idlgrroi",data=lem_data, color = color_stereo_A,style=2,hide = 1)
R_image_model-> Add, R_cme_MATCH_outline

;lem_region = obj_new('IDLanROI')







;  End of the new 3D stuff.....





cone_radius = 0.

rotation = 0.
pixel_scale = 0.

L_telescope_code = ''
C_telescope_code = ''
C2_telescope_code = ''
R_telescope_code = ''
date_array_int = intarr(10)
start_date = intarr(5)
end_date = intarr(5)
start_julian = -1.0D
end_julian = -1.0D
animation_start_julian = -1.0D
animation_end_julian = -1.0D

rotation = 0.
rsun = 950.   ;set initial value, solar radius from L1 in arcseconds
;pixel_scale = 1.
rotate_3D_cme =  1   ; trackball won't work.
;be sure to use reals, unless you only want integers...
scaling_factor = 1.0
;rotation = 0.

L_list_of_image_names = list()
L_list_of_image_data = list()
L_list_of_datetime_strings = list()
L_list_of_datetime_Julian = list()
L_list_of_full_time_strings = list()
L_list_of_image_exposure_times = list()
L_list_of_image_offsets = list()
L_list_of_image_scaling_factors = list()

C_list_of_image_names = list()
C_list_of_image_data = list()
C_list_of_datetime_strings = list()
C_list_of_datetime_Julian = list()
C_list_of_full_time_strings = list()
C_list_of_image_exposure_times = list()
C_list_of_image_offsets = list()
C_list_of_image_scaling_factors = list()

C2_list_of_image_names = list()
C2_list_of_image_data = list()
C2_list_of_datetime_strings = list()
C2_list_of_datetime_Julian = list()
C2_list_of_full_time_strings = list()
C2_list_of_image_exposure_times = list()
C2_list_of_image_offsets = list()
C2_list_of_image_scaling_factors = list()

R_list_of_image_names = list()
R_list_of_image_data = list()
R_list_of_datetime_strings = list()
R_list_of_datetime_Julian = list()
R_list_of_full_time_strings = list()
R_list_of_image_exposure_times = list()
R_list_of_image_offsets = list()
R_list_of_image_scaling_factors = list()

L_list_of_HEEQ_coords = list()
C_list_of_HEEQ_coords = list()
C2_list_of_HEEQ_coords = list()
R_list_of_HEEQ_coords = list()

L_HEEQ_coords = [0.,0.]
C_HEEQ_coords = [0.,0.]
C2_HEEQ_coords = [0.,0.]
R_HEEQ_coords = [0.,0.]
master_list = list()
L_indexes = list()
C_indexes = list()
C2_indexes = list()
R_indexes = list()
master_list_size = 0
L_current_image_number = 0
C_current_image_number = 0
C2_current_image_number = 0
R_current_image_number = 0

initial_transform = dblarr(4,4)
initial_transform[0,0] = 1.
initial_transform[1,1] = 1.
initial_transform[2,2] = 1.
initial_transform[3,3] = 1.

latitude_degrees = 0.
longitude_degrees = 0.

CME_matches_image_Julian = list()
CME_matches_image_DateTime_string = list()
CME_matches_image_telescope = list()
CME_matches_image_Rs_leading_edge = list()
CME_matches_image_Image_number = list()
CME_matches_image_CME_outline = list()
CME_matches_image_L_Image_number = list()
CME_matches_image_C_Image_number = list()
CME_matches_image_C2_Image_number = list()
CME_matches_image_R_Image_number = list()
CME_matches_image_L_CME_outline = list()
CME_matches_image_C_CME_outline = list()
CME_matches_image_C2_CME_outline = list()
CME_matches_image_R_CME_outline = list()

; ECDP, lists for POS mode
POS_Julian = list()
POS_leading_edge = list()
POS_matches_image_DateTime_string = list()
POS_matches_image_telescope = list()

L_list_of_pixel_scales = list()
C_list_of_pixel_scales = list()
C2_list_of_pixel_scales = list()
R_list_of_pixel_scales = list()

L_list_of_rsuns = list()
C_list_of_rsuns = list()
C2_list_of_rsuns = list()
R_list_of_rsuns = list()

L_list_of_Sun_satellite_distances = list()
C_list_of_Sun_satellite_distances = list()
C2_list_of_Sun_satellite_distances = list()
R_list_of_Sun_satellite_distances = list()

anim_start_frame = 1
anim_end_frame = 1
anim_start_julian = 0.d
anim_end_julian = 0.d
i_anim_back_and_forth = 0
i_move = 1

which_window_to_animate = 1

clicked_L = 0
clicked_C = 0
clicked_R = 0
clicked_C2 = 0

currently_showing_lasco_C2_or_C3 = 'C3'

xSymbolSize_timeline = 0.
ySymbolSize_timeline = 0.

representative_image_data = bytarr(512,512,3)
representative_image_has_been_defined = 0

images_are_loaded = 0

timeline_left_mouse_button_being_pressed = 0

background_color_lasco = [175,175,175]
background_color_stereo_A = [175,175,175]
background_color_stereo_B = [175,175,175]

animation_frames = 0
animate_L = 0
animate_C = 1
animate_R = 0

timeline_normalized_x = 0.

lemniscate_style = 0

manual_point = list()
use_manual_point = 0

i_log_scale = 0

; all to do with the red line....

L_the_action = 0
L_latest_click_X = 0
L_latest_click_Y = 0
L_previous_click_X = 0
L_previous_click_Y = 0
L_click_and_drag  = 0
L_clock_angle_degrees = 0
L_rotate_x = 0
L_rotate_y = 0

C_the_action = 0
C_latest_click_X = 0
C_latest_click_Y = 0
C_previous_click_X = 0
C_previous_click_Y = 0
C_click_and_drag  = 0
C_clock_angle_degrees = 0
C_rotate_x = 0
C_rotate_y = 0

half_cone_angle_emp = 0.

;which_telescope = ''
;telescope_field_of_view_Rs = 30.
;image_size_pixels = 512.
;conversion_factor_Rs_per_pixel = telescope_field_of_view_Rs * 2. / image_size_pixels

rsun = 950.   ;set initial value, solar radius from L1 in arcseconds
pixel_scale = 1.
;clock_angle_degrees = 0
C_leading_edge_in_rs = 0.
;telescope_code = ''

;image_type = 'fits'
;image_fits_header = ptr_new(/allocate_heap)
;image_saturation_value = 50.
;C_the_action = 0

;full_halo = 0
; ECDP


R_the_action = 0
R_latest_click_X = 0
R_latest_click_Y = 0
R_previous_click_X = 0
R_previous_click_Y = 0
R_click_and_drag  = 0
R_clock_angle_degrees = 0
R_rotate_x = 0
R_rotate_y = 0

pressed_the_right_button = 0

position_image_lineplot = [29, 29, 482, 100]

c_now_dist = 0.
;c_now_x = 0.
;c_now_y = 0.

mode_name = ''

info = $
       { L_Window:L_Window, $ 
         L_Window_copy:L_Window_copy, $ 
         L_both_views:L_both_views, $     
         L_image_view:L_image_view, $     
         L_second_view:L_second_view, $
         L_second_view_copy:L_second_view_copy, $
         C_Window:C_Window, $ 
         C_Window_copy:C_Window_copy, $ 
         C_both_views:C_both_views, $     
         C_image_view:C_image_view, $     
         C_second_view:C_second_view, $
         C_second_view_copy:C_second_view_copy, $
         R_Window:R_Window, $ 
         R_Window_copy:R_Window_copy, $ 
         R_both_views:R_both_views, $     
         R_image_view:R_image_view, $     
         R_second_view:R_second_view, $
         R_second_view_copy:R_second_view_copy, $
         draw_L:draw_L, $
         draw_C:draw_C, $
         draw_R:draw_R, $  
         cme_info_Window:cme_info_Window, $
         cme_info_view:cme_info_view, $
         cone_radius:cone_radius, $
         widget_angular_width_slider:widget_angular_width_slider, $
         widget_radial_distance_slider:widget_radial_distance_slider, $
         widget_select_mode2:widget_select_mode2, $
         L_trackball:L_trackball, $ 
         L_Fixed_Sun_Satellite_model:L_Fixed_Sun_Satellite_model, $
         L_Fixed_Sun_Satellite_model_copy:L_Fixed_Sun_Satellite_model_copy, $
         L_Fixed_Sun_Satellite_model_original_transform:L_Fixed_Sun_Satellite_model_original_transform, $        
         L_Rotating_model_Transform:L_Rotating_model_Transform, $
         C_trackball:C_trackball, $ 
         C_Fixed_Sun_Satellite_model:C_Fixed_Sun_Satellite_model, $
         C_Fixed_Sun_Satellite_model_copy:C_Fixed_Sun_Satellite_model_copy, $
         C_Fixed_Sun_Satellite_model_original_transform:C_Fixed_Sun_Satellite_model_original_transform, $        
         C_Rotating_model_Transform:C_Rotating_model_Transform, $
         R_trackball:R_trackball, $ 
         R_Fixed_Sun_Satellite_model:R_Fixed_Sun_Satellite_model, $
         R_Fixed_Sun_Satellite_model_copy:R_Fixed_Sun_Satellite_model_copy, $
         R_Fixed_Sun_Satellite_model_original_transform:R_Fixed_Sun_Satellite_model_original_transform, $        
         R_Rotating_model_Transform:R_Rotating_model_Transform, $
         L_cone_Z_axis:L_cone_Z_axis, $
         C_cone_Z_axis:C_cone_Z_axis, $
         R_cone_Z_axis:R_cone_Z_axis, $
         lat_string:lat_string, $
         lon_string:lon_string, $
         lat_string_object:lat_string_object,$
         lon_string_object:lon_string_object, $
         angular_width_string:angular_width_string, $
         angular_width_string_object:angular_width_string_object, $
         radial_distance_string:radial_distance_string, $
         radial_distance_string_object:radial_distance_string_object, $            
         images_timeline_window:images_timeline_window, $      
         images_timeline_view:images_timeline_view, $            
         images_timeline_model:images_timeline_model, $
         which_mode:which_mode, $            
         draw_available_images_timeline:draw_available_images_timeline, $
         position_timeline: position_timeline, $
         L_plot: L_plot, $
         C_plot: C_plot, $
         C2_plot: C2_plot, $
         R_plot: R_plot, $
;         mode_array:mode_array, $
         xaxis_images_timeline  : xaxis_images_timeline, $
         position_B: position_B, $
         LE_plot_matched_CMEs:LE_plot_matched_CMEs, $
         thisPlot_B2:thisPlot_B2, $
         thisPlot2_c3:thisPlot2_c3, $
         thisPlot2_c2:thisPlot2_c2, $
         thisPlot2_stereo_a:thisPlot2_stereo_a, $
         thisPlot2_stereo_b:thisPlot2_stereo_b, $       
         xAxis1_B  : xAxis1_B, $
         xAxis2_B  : xAxis2_B, $
         yAxis1_B  : yAxis1_B, $
         yAxis2_B  : yAxis2_B, $
         xmin:xmin, $
         xmax:xmax, $
         ymin:ymin, $
         ymax:ymax, $
         xtickinterval:xtickinterval, $
         latitude_text:latitude_text, $
         longitude_text:longitude_text, $
         cone_angle_text:cone_angle_text, $
         velocity_text:velocity_text, $
         T21_5_text:T21_5_text, $
         T21_5_string:T21_5_string, $
         cone_tool_version:cone_tool_version, $
         tlb_position:tlb_position, $           
         files_directory:files_directory, $
         L_number_of_images:L_number_of_images, $
         C_number_of_images:C_number_of_images, $
         C2_number_of_images:C2_number_of_images, $
         R_number_of_images:R_number_of_images, $
         L_image_saturation_value:L_image_saturation_value, $
         C_image_saturation_value:C_image_saturation_value, $
         R_image_saturation_value:R_image_saturation_value, $
         palette_top_normal:palette_top_normal, $
         palette_bottom_normal:palette_bottom_normal, $
         palette_gamma_normal:palette_gamma_normal, $
         palette_top_diff:palette_top_diff, $
         palette_bottom_diff:palette_bottom_diff, $
         palette_gamma_diff:palette_gamma_diff, $
         L_coronagraph_image_object:L_coronagraph_image_object, $
         C_coronagraph_image_object:C_coronagraph_image_object, $
         R_coronagraph_image_object:R_coronagraph_image_object, $
         L_border_image_object:L_border_image_object, $
         C_border_image_object:C_border_image_object, $
         R_border_image_object:R_border_image_object, $                  
         L_image_color_palette:L_image_color_palette, $
         C_image_color_palette:C_image_color_palette, $
         R_image_color_palette:R_image_color_palette, $
         widget_calculate_velocity:widget_calculate_velocity, $
         widget_calculate_pos_velocity:widget_calculate_pos_velocity, $
         widget_export_cone:widget_export_cone, $
         L_widget_botSlider:L_widget_botSlider, $
         L_widget_topSlider:L_widget_topSlider, $
         L_widget_gammaSlider:L_widget_gammaSlider, $
         L_widget_saturationSlider:L_widget_saturationSlider, $
         C_widget_botSlider:C_widget_botSlider, $
         C_widget_topSlider:C_widget_topSlider, $
         C_widget_gammaSlider:C_widget_gammaSlider, $
         C_widget_saturationSlider:C_widget_saturationSlider, $
         R_widget_botSlider:R_widget_botSlider, $
         R_widget_topSlider:R_widget_topSlider, $
         R_widget_gammaSlider:R_widget_gammaSlider, $
         R_widget_saturationSlider:R_widget_saturationSlider, $
         center_of_sunX:center_of_sunX, $
         center_of_sunY:center_of_sunY, $
         current_R_in_Solar_radii:current_R_in_Solar_radii, $
         radial_velocity:radial_velocity, $
         time_at_Enlil_boundary_Julian:time_at_Enlil_boundary_Julian, $
         time_at_Enlil_boundary:time_at_Enlil_boundary, $
         initial_cme_outline_associated_with_image_number:initial_cme_outline_associated_with_image_number, $
         L_difference_imaging:L_difference_imaging, $
         C_difference_imaging:C_difference_imaging, $
         C2_difference_imaging:C2_difference_imaging, $
         R_difference_imaging:R_difference_imaging, $
         L_background_image_number:L_background_image_number, $
         C_background_image_number:C_background_image_number, $
         C2_background_image_number:C2_background_image_number, $
         R_background_image_number:R_background_image_number, $
         diff_red_levels:diff_red_levels, $
         diff_green_levels:diff_green_levels, $
         diff_blue_levels:diff_blue_levels, $
         wcolors:wcolors, $
         palette_0:palette_0, $
         which_telescope:which_telescope, $
         L_telescope_field_of_view_Rs:L_telescope_field_of_view_Rs, $
         C_telescope_field_of_view_Rs:C_telescope_field_of_view_Rs, $
         C2_telescope_field_of_view_Rs:C2_telescope_field_of_view_Rs, $
         R_telescope_field_of_view_Rs:R_telescope_field_of_view_Rs, $
         image_size_pixels:image_size_pixels, $
         Solar_radius_km:Solar_radius_km, $
         B_angle_degrees:B_angle_degrees, $
         L_cme_fitted_surf:L_cme_fitted_surf, $
         L_cme_fitted_surf_copy:L_cme_fitted_surf_copy, $
         C_cme_fitted_surf:C_cme_fitted_surf, $
         C_cme_fitted_surf_copy:C_cme_fitted_surf_copy, $
         R_cme_fitted_surf:R_cme_fitted_surf, $
         R_cme_fitted_surf_copy:R_cme_fitted_surf_copy, $
         angular_width_lemniscate:angular_width_lemniscate, $
         radial_distance_lemniscate:radial_distance_lemniscate, $
         L_cme_outline:L_cme_outline, $
         C_cme_outline:C_cme_outline, $
         R_cme_outline:R_cme_outline, $
         L_cme_MATCH_outline:L_cme_MATCH_outline, $
         C_cme_MATCH_outline:C_cme_MATCH_outline, $
         C2_cme_MATCH_outline:C2_cme_MATCH_outline, $
         R_cme_MATCH_outline:R_cme_MATCH_outline, $
;         lem_region:lem_region, $
         L_telescope_code:L_telescope_code, $
         C_telescope_code:C_telescope_code, $
         C2_telescope_code:C2_telescope_code, $
         R_telescope_code:R_telescope_code, $
         date_array_int:date_array_int, $
         start_date:start_date, $
         end_date:end_date, $
         start_julian:start_julian, $
         end_julian:end_julian, $
         telescope_array:telescope_array, $
         image_in_folder_array:image_in_folder_array, $
         image_in_root:image_in_root, $
         export_location_root:export_location_root, $
         max_interval_in_days:max_interval_in_days, $
         sep:sep, $
         rotation:rotation, $
         rsun:rsun, $
         pixel_scale:pixel_scale, $
         L_list_of_image_names:L_list_of_image_names, $
         L_list_of_image_data:L_list_of_image_data, $
         L_list_of_datetime_strings:L_list_of_datetime_strings, $
         L_list_of_datetime_Julian:L_list_of_datetime_Julian, $
         L_list_of_full_time_strings:L_list_of_full_time_strings, $
         L_list_of_image_exposure_times:L_list_of_image_exposure_times, $
         L_list_of_image_offsets:L_list_of_image_offsets, $
         L_list_of_image_scaling_factors:L_list_of_image_scaling_factors, $
         C_list_of_image_names:C_list_of_image_names, $
         C_list_of_image_data:C_list_of_image_data, $
         C_list_of_datetime_strings:C_list_of_datetime_strings, $
         C_list_of_datetime_Julian:C_list_of_datetime_Julian, $
         C_list_of_full_time_strings:C_list_of_full_time_strings, $
         C_list_of_image_exposure_times:C_list_of_image_exposure_times, $
         C_list_of_image_offsets:C_list_of_image_offsets, $
         C_list_of_image_scaling_factors:C_list_of_image_scaling_factors, $
         C2_list_of_image_names:C2_list_of_image_names, $
         C2_list_of_image_data:C2_list_of_image_data, $
         C2_list_of_datetime_strings:C2_list_of_datetime_strings, $
         C2_list_of_datetime_Julian:C2_list_of_datetime_Julian, $
         C2_list_of_full_time_strings:C2_list_of_full_time_strings, $
         C2_list_of_image_exposure_times:C2_list_of_image_exposure_times, $
         C2_list_of_image_offsets:C2_list_of_image_offsets, $
         C2_list_of_image_scaling_factors:C2_list_of_image_scaling_factors, $
         R_list_of_image_names:R_list_of_image_names, $
         R_list_of_image_data:R_list_of_image_data, $
         R_list_of_datetime_strings:R_list_of_datetime_strings, $
         R_list_of_datetime_Julian:R_list_of_datetime_Julian, $
         R_list_of_full_time_strings:R_list_of_full_time_strings, $
         R_list_of_image_exposure_times:R_list_of_image_exposure_times, $
         R_list_of_image_offsets:R_list_of_image_offsets, $
         R_list_of_image_scaling_factors:R_list_of_image_scaling_factors, $
         rotate_3D_cme:rotate_3D_cme, $
         cme_info_model:cme_info_model, $
;         widget_set_start_and_end_time:widget_set_start_and_end_time, $
         master_list:master_list, $
         L_indexes:L_indexes, $
         C_indexes:C_indexes, $
         C2_indexes:C2_indexes, $
         R_indexes:R_indexes, $
         master_list_size:master_list_size, $
         L_current_image_number:L_current_image_number, $
         C_current_image_number:C_current_image_number, $
         C2_current_image_number:C2_current_image_number, $
         R_current_image_number:R_current_image_number, $
         L_list_of_HEEQ_coords:L_list_of_HEEQ_coords, $
         C_list_of_HEEQ_coords:C_list_of_HEEQ_coords, $
         C2_list_of_HEEQ_coords:C2_list_of_HEEQ_coords, $
         R_list_of_HEEQ_coords:R_list_of_HEEQ_coords, $
         L_list_of_pixel_scales:L_list_of_pixel_scales, $
         C_list_of_pixel_scales:C_list_of_pixel_scales, $
         C2_list_of_pixel_scales:C2_list_of_pixel_scales, $
         R_list_of_pixel_scales:R_list_of_pixel_scales, $
         L_list_of_rsuns:L_list_of_rsuns, $
         C_list_of_rsuns:C_list_of_rsuns, $
         C2_list_of_rsuns:C2_list_of_rsuns, $
         R_list_of_rsuns:R_list_of_rsuns, $
         L_list_of_Sun_satellite_distances:L_list_of_Sun_satellite_distances, $
         C_list_of_Sun_satellite_distances:C_list_of_Sun_satellite_distances, $
         C2_list_of_Sun_satellite_distances:C2_list_of_Sun_satellite_distances, $
         R_list_of_Sun_satellite_distances:R_list_of_Sun_satellite_distances, $
         L_HEEQ_coords:L_HEEQ_coords, $
         C_HEEQ_coords:C_HEEQ_coords, $
         C2_HEEQ_coords:C2_HEEQ_coords, $
         R_HEEQ_coords:R_HEEQ_coords, $
         L_coords_offset_model:L_coords_offset_model, $
         L_coords_offset_model_copy:L_coords_offset_model_copy, $
         C_coords_offset_model:C_coords_offset_model, $
         C_coords_offset_model_copy:C_coords_offset_model_copy, $
         R_coords_offset_model:R_coords_offset_model, $
         R_coords_offset_model_copy:R_coords_offset_model_copy, $
         initial_transform:initial_transform, $
         latitude_degrees:latitude_degrees, $
         longitude_degrees:longitude_degrees, $
         L_cme_model:L_cme_model, $
         C_cme_model:C_cme_model, $
         R_cme_model:R_cme_model, $
         L_cme_model_copy:L_cme_model_copy, $
         C_cme_model_copy:C_cme_model_copy, $
         R_cme_model_copy:R_cme_model_copy, $
         L_ut_string_object:L_ut_string_object, $
         C_ut_string_object:C_ut_string_object, $
         R_ut_string_object:R_ut_string_object, $
         xsize:xsize, $
         ysize:ysize, $
         L_widget_image_sequence_slider:L_widget_image_sequence_slider, $
         C_widget_image_sequence_slider:C_widget_image_sequence_slider, $
         R_widget_image_sequence_slider:R_widget_image_sequence_slider, $
         L_widget_outline_matches_image:L_widget_outline_matches_image, $
         C_widget_outline_matches_image:C_widget_outline_matches_image, $
         C_widget_pos_outline_matches_image:C_widget_pos_outline_matches_image, $
         R_widget_outline_matches_image:R_widget_outline_matches_image, $
         widget_latitude_slider:widget_latitude_slider, $
         widget_longitude_slider:widget_longitude_slider, $
         CME_matches_image_Julian:CME_matches_image_Julian, $
         CME_matches_image_DateTime_string:CME_matches_image_DateTime_string, $
         CME_matches_image_telescope:CME_matches_image_telescope, $
         CME_matches_image_Rs_leading_edge:CME_matches_image_Rs_leading_edge, $
         CME_matches_image_Image_number:CME_matches_image_Image_number, $
         CME_matches_image_CME_outline:CME_matches_image_CME_outline, $
         C_leading_edge_string:C_leading_edge_string, $
         plot_window:plot_window, $
         plot_view:plot_view, $
         Julian_day_for_Earth_pos:Julian_day_for_Earth_pos, $
         Earth_pos_AU:Earth_pos_AU, $
         Earth_pos_HG_LAT_deg:Earth_pos_HG_LAT_deg, $
         source_path:source_path, $
         play_pause_button:play_pause_button, $
         delay:0.25, $
         stopflag:0, $
         thisFrame:0, $
         anim_start_frame:anim_start_frame, $
         anim_end_frame:anim_end_frame, $
         i_anim_back_and_forth:i_anim_back_and_forth, $
         i_move:i_move, $
         which_window_to_animate:which_window_to_animate, $
         enlil_info_View:enlil_info_View, $
         enlil_info_window:enlil_info_window, $
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
         acceptID:acceptID, $
         date_array:date_array, $
         clicked_L:clicked_L, $
         clicked_C:clicked_C, $
         clicked_C2:clicked_C2, $
         clicked_R:clicked_R, $
         currently_showing_lasco_C2_or_C3:currently_showing_lasco_C2_or_C3, $
         L_title_object:L_title_object, $
         C_title_object:C_title_object, $
         R_title_object:R_title_object, $
         color_stereo_B:color_stereo_B, $
         color_stereo_A:color_stereo_A, $
         color_c3:color_c3, $
         color_c2:color_c2, $
         cme_outline_color:cme_outline_color, $
         background_color:background_color, $
         xsize_timeline:xsize_timeline, $
         ysize_timeline:ysize_timeline, $
         CME_matches_image_L_Image_number:CME_matches_image_L_Image_number, $
         CME_matches_image_C_Image_number:CME_matches_image_C_Image_number, $
         CME_matches_image_C2_Image_number:CME_matches_image_C2_Image_number, $
         CME_matches_image_R_Image_number:CME_matches_image_R_Image_number, $
         CME_matches_image_L_CME_outline:CME_matches_image_L_CME_outline, $
         CME_matches_image_C_CME_outline:CME_matches_image_C_CME_outline, $
         CME_matches_image_C2_CME_outline:CME_matches_image_C2_CME_outline, $
         CME_matches_image_R_CME_outline:CME_matches_image_R_CME_outline, $
         L_current_background_color:L_current_background_color, $
         L_current_text_color:L_current_text_color, $
         C_current_background_color:C_current_background_color, $
         C_current_text_color:C_current_text_color, $
         R_current_background_color:R_current_background_color, $
         R_current_text_color:R_current_text_color, $
         widget_show_C2_or_C3:widget_show_C2_or_C3, $
         xSymbolSize_timeline:xSymbolSize_timeline, $
         ySymbolSize_timeline:ySymbolSize_timeline, $
         widget_reset_cme_analysis:widget_reset_cme_analysis, $
         L_widget_representative_image:L_widget_representative_image, $
         C_widget_representative_image:C_widget_representative_image, $
         R_widget_representative_image:R_widget_representative_image, $
         representative_image_data:representative_image_data, $
         representative_image_has_been_defined:representative_image_has_been_defined, $
         OPS_or_VnV:OPS_or_VnV, $
         images_are_loaded:images_are_loaded, $
         animation_controls_outer_base:animation_controls_outer_base, $
         image_controls_base:image_controls_base, $
         cme_controls_base:cme_controls_base, $
         plot_window_base:plot_window_base, $
         enlil_info_base:enlil_info_base, $
         timeline_left_mouse_button_being_pressed:timeline_left_mouse_button_being_pressed, $
         background_color_lasco:background_color_lasco, $
         background_color_stereo_A:background_color_stereo_A, $
         background_color_stereo_B:background_color_stereo_B, $
         animation_frames:animation_frames, $
         animate_L:animate_L, $
         animate_C:animate_C, $
         animate_R:animate_R, $
         animation_start_time_marker:animation_start_time_marker, $
         animation_end_time_marker:animation_end_time_marker, $
         animation_start_time_marker_handle:animation_start_time_marker_handle, $
         animation_end_time_marker_handle:animation_end_time_marker_handle, $
         animation_start_time_marker_handle_data:animation_start_time_marker_handle_data, $
         animation_end_time_marker_handle_data:animation_end_time_marker_handle_data, $
         animation_current_time_marker:animation_current_time_marker, $
         animation_start_julian:animation_start_julian, $
         animation_end_julian:animation_end_julian, $
         timeline_normalized_x:timeline_normalized_x, $
         anim_start_julian:anim_start_julian, $
         anim_end_julian:anim_end_julian, $
         widget_animation_delay:widget_animation_delay, $
         widget_which_window_to_animate:widget_which_window_to_animate, $
         widget_anim_back_and_forth:widget_anim_back_and_forth, $
         L_widget_remove_this_image:L_widget_remove_this_image, $
         C_widget_remove_this_image:C_widget_remove_this_image, $
         R_widget_remove_this_image:R_widget_remove_this_image, $
         widget_plus_12h_button:widget_plus_12h_button, $
         widget_plus_24h_button:widget_plus_24h_button, $
         lemniscate_style:lemniscate_style, $
         manual_point:manual_point, $
         use_manual_point:use_manual_point, $
         plot_window_manually_add_point:plot_window_manually_add_point, $
         allow_manual_data_point:allow_manual_data_point, $
         i_log_scale:i_log_scale, $
         L_the_action:L_the_action, $
         L_latest_click_X:L_latest_click_X, $
         L_latest_click_Y:L_latest_click_Y, $
         L_previous_click_X:L_previous_click_X, $
         L_previous_click_Y:L_previous_click_X, $
         L_click_and_drag:L_click_and_drag, $
         L_clock_angle_degrees:L_clock_angle_degrees, $
         L_clock_angle_model:L_clock_angle_model, $
         L_rotate_x:L_rotate_x, $
         L_rotate_y:L_rotate_y, $
         L_image_lineplot:L_image_lineplot, $
         C_the_action:C_the_action, $
         C_latest_click_X:C_latest_click_X, $
         C_latest_click_Y:C_latest_click_Y, $
         C_previous_click_X:C_previous_click_X, $
         C_previous_click_Y:C_previous_click_X, $
         C_click_and_drag:C_click_and_drag, $
         C_clock_angle_degrees:C_clock_angle_degrees, $
         C_clock_angle_model:C_clock_angle_model, $
         C_rotate_x:C_rotate_x, $
         C_rotate_y:C_rotate_y, $
         C_clock_angle_marker_line:C_clock_angle_marker_line, $
         C_image_lineplot:C_image_lineplot, $
         R_the_action:R_the_action, $
         R_latest_click_X:R_latest_click_X, $
         R_latest_click_Y:R_latest_click_Y, $
         R_previous_click_X:R_previous_click_X, $
         R_previous_click_Y:R_previous_click_X, $
         R_click_and_drag:R_click_and_drag, $
         R_clock_angle_degrees:R_clock_angle_degrees, $
         R_clock_angle_model:R_clock_angle_model, $
         R_rotate_x:R_rotate_x, $
         R_rotate_y:R_rotate_y, $
         R_image_lineplot:R_image_lineplot, $
         pressed_the_right_button:pressed_the_right_button, $
         position_image_lineplot:position_image_lineplot, $
         show_image_line_plot:show_image_line_plot, $
         Three_D_view_tlb_ID:Three_D_view_tlb_ID, $
         allow_show_3D_view:allow_show_3D_view, $
         Three_D_view_event_widget:Three_D_view_event_widget, $
;         L_widget_show_line_plot:L_widget_show_line_plot, $
;         C_widget_show_line_plot:C_widget_show_line_plot, $
;         R_widget_show_line_plot:R_widget_show_line_plot, $
;         L_image_lineplot_model:L_image_lineplot_model, $
;         C_image_lineplot_model:C_image_lineplot_model, $
;         R_image_lineplot_model:R_image_lineplot_model, $
         velocity_text_c3:velocity_text_c3, $
         velocity_text_c2:velocity_text_c2, $
         velocity_text_stereo_a:velocity_text_stereo_a, $
         velocity_text_stereo_b:velocity_text_stereo_b, $
         calculate_individual_velocities_for_each_telescope:calculate_individual_velocities_for_each_telescope, $
         debug_mode:debug_mode, $
;         swpc_cat_preferences_file:swpc_cat_preferences_file, $
         C_leading_edge_model:C_leading_edge_model, $
         C_leading_edge_in_rs:C_leading_edge_in_rs, $
         scaling_factor:scaling_factor, $
         C_clock_angle_string:C_clock_angle_string, $
         POS_leading_edge:POS_leading_edge, $
         POS_Julian:POS_Julian, $
         POS_velocity:POS_velocity, $
         POS_velocity_text:POS_velocity_text, $
         half_cone_angle_emp:half_cone_angle_emp, $
         emp_cone_angle_text:emp_cone_angle_text, $
         time_plot_title_widget:time_plot_title_widget, $
         C_now_dist:C_now_dist, $ 
         widget_use_emp_cone:widget_use_emp_cone, $
         mode_name:mode_name, $
         c_mode_name:c_mode_name, $
;         c_now_x:c_now_x, $
;         c_now_y:c_now_y, $
;        mode_title_object:mode_title_object, $
;        mode_title_view:mode_title_view, $
;        mode_title_model:mode_title_model, $
         compact:compact, $
	 preferences_folder:preferences_folder, $
         output_line_plots:output_line_plots, $
         POS_matches_image_telescope:POS_matches_image_telescope, $
         output_matched_pos_line_data_in_txt_file:output_matched_pos_line_data_in_txt_file, $
         POS_matches_image_DateTime_string:POS_matches_image_DateTime_string, $
         C_clock_angle_model_original_transform:C_clock_angle_model_original_transform, $
         output_matched_line_data_in_txt_file:output_matched_line_data_in_txt_file}
 
         
Widget_Control, tlb, Set_UValue=info, /No_Copy

XManager, 'CAT', tlb, Cleanup='swpc_cat_Cleanup',  /No_Block, $
   Event_Handler='swpc_cat_top_level_base_events', Group_Leader=groupLeader
END

