
; Compute CESM/LENS TRACER budget terms

;*****************************************************************************;
load "eoshift.ncl" 

;*****************************************************************************;
undef("tracer_budget_vol3d")
function tracer_budget_vol3d(fname:string, nlon:integer, nlat:integer, nz:integer)
;
;  Returns global 3D volume array: vold3(nz,ny,nx)
;  NOTE: does not include SSH variations
;
local vol3d, fnames, file_ind, f, dz, TAREA, KMT, MASK_2D
begin
   ; initialize vol3d 
   vol3d = new((/ nz, nlat, nlon /), "double")
   ;vol3d!0 = "nz_G"
   vol3d!0 = "z_t"
   ;vol3d!1 = "nlat_G"
   vol3d!1 = "nlat"
   ;vol3d!2 = "nlon_G"
   vol3d!2 = "nlon"
   vol3d@long_name = "Tcell volume"
   vol3d@units = "cm^3"
   vol3d = 0.0

   file_ind = 0
   f = addfile(fname, "r")
   TAREA = f->TAREA
   KMT = f->KMT
   dz = f->dz

   vol3d = 0.0
   do k = 0, nz-1
      MASK_2D = where(k .lt. KMT, 1.0, 0.0)
      if (sum(MASK_2D) .gt. 0.0) then
         vol3d(k,:,:) = (/ (dz(k) * TAREA * MASK_2D) /)
      end if
   end do
  return(vol3d)
end

;*****************************************************************************;
undef("tracer_budget_var_3d_zint_map")
function tracer_budget_var_3d_zint_map(var_name:string, fname:string, tlo:integer, thi:integer,\ 
            nlon:integer, nlat:integer, klo:integer, khi:integer,vol3d[*][*][*]:double,\
            annual_cycle:logical, ssh_flag:logical, fssh:string)
local var_zint_map, FIELD, dz, workvol, MASK_2D, VOL_WEIGHT
begin
    f = addfile(fname,"r")
    time = f->time(tlo:thi)
    workvol = vol3d(klo:khi,:,:)
    KMT = f->KMT
    MASK_2D = where(KMT .gt. 0, 1.0, 0.0)

    FIELD1 = f->$var_name$(tlo:thi,klo:khi,:,:)
    if annual_cycle then
        FIELD2 = rmMonAnnCycLLLT(FIELD1(z_t|:, nlat|:, nlon|:, time|:))
        FIELD = FIELD2(time|:, z_t|:, nlat|:, nlon|:)
        delete(FIELD2)
    else 
       FIELD = FIELD1
    end if
    delete(FIELD1)
   
    ; initialize var_zint_map
    var_zint_map = new((/ dimsizes(time), nlat, nlon /), "double")
    var_zint_map!0 = "time"
    var_zint_map&time = time
    var_zint_map!1 = "nlat"
    var_zint_map!2 = "nlon"
    var_zint_map@long_name = var_name+" vertical integral"
    var_zint_map = 0.0
   
    var_zint_map@units = FIELD@units + "cm^3"
    dz = f->dz
    zlo = f->z_w(klo)
    zhi = f->z_w(khi+1)

    var_zint_map@k_range = tostring(klo)+" - "+tostring(khi)
    var_zint_map@depth_range = sprintf("%6.1f",zlo/100.)+"m -"+sprintf("%6.1f",zhi/100.)+"m"
      
    VOL_WEIGHT = conform(FIELD, workvol, (/1,2,3/))
    ; add SSH volume variation
    if (ssh_flag.and.klo.eq.0) then
        f2 = addfile(fssh,"r")
        SSH = f2->SSH(tlo:thi,:,:)
        TAREA = f2->TAREA
        VOL_WEIGHT(:,0,:,:) = VOL_WEIGHT(:,0,:,:) + SSH*conform(SSH, TAREA * MASK_2D, (/ 1, 2 /))
        var_zint_map@long_name = var_name+" vertical integral + SSH vol variation"
        delete(f2)
    end if
    var_zint_map = (/ dim_sum_n(FIELD * VOL_WEIGHT,1) /)
    return(var_zint_map)
end

;*****************************************************************************;
undef("tracer_budget_tend_appr")
function tracer_budget_tend_appr(time_bound[*][*]:double, var_zint[*][*][*]:double, annual_cycle:logical)
;
; Computes approximate TRACER budget tendency given vertically-integrated POP
; TRACER based on differencing successive monthly means
; NOTE: Assumes monthly POP output with timestamp at end-of-month
;       rather than mid-month; assumes time has dimension "days".
;
local date,timec,deltat
begin
    secperday = 60.*60.*24.
    var_zint_tend = var_zint
    var_zint_tend@units = var_zint@units + "/s"

    nt0 = dimsizes(time_bound)
    nt = nt0(0)
    deltat = (time_bound(:,1)-time_bound(:,0))*secperday
    if annual_cycle then
        work1 = rmMonAnnCycTLL(var_zint)
    else
        work1 = var_zint
    end if
    work1(0:nt-2,:,:) = 0.5*(work1(0:nt-2,:,:)+work1(1:nt-1,:,:))	; end-of-month TRACER
    work1(nt-1,:,:) = work1@_FillValue

    work2 = var_zint
    work2(1:nt-2,:,:) = (work1(1:nt-2,:,:) - work1(0:nt-3,:,:))
    work2(0,:,:) = work2@_FillValue
    work2(nt-1,:,:) = work2@_FillValue

    var_zint_tend = (/ work2/conform(work2,deltat,(/0/)) /)
   return(var_zint_tend)
end

;*****************************************************************************;
undef("get_tb")
function get_tb(fname:string)
local fnames, f, tb
begin
   f = addfile(fname, "r")
   tb = f->time_bound
   return(tb)
end

;*****************************************************************************;
undef("tracer_budget_lat_adv_resolved") 
function tracer_budget_lat_adv_resolved(fname1:string, fname2:string, \
   nlon:integer, nlat:integer, klo:integer, khi:integer, tlo:integer, thi:integer, vol3d[*][*][*]:double)
;
; horizontal divergence: advection
;
local var_lat_adv_res_map, \
      FIELD1e, FIELD1w, FIELD2n, FIELD2s, MASK_3D, volc,vols,volw,work, work1, \
      time, work2, work3, work4, dz
begin
   work = vol3d				            ; Tcell_i,j
   work1 = (/ eoshift(work,1,-1) /)		; Tcell_i-1,j
   work2 = (/ eoshift(work,2,-1) /)		; Tcell_i,j-1
   volc = work(klo:khi,:,:)
   volw = work1(klo:khi,:,:)
   vols = work2(klo:khi,:,:)
   delete(work)
   delete(work1)
   delete(work2)

   ; initialize :
   var_name1 = "UET"
   var_name2 = "VNT"
   f1 = addfile(fname1,"r")
   f2 = addfile(fname2,"r")
   time = f1->time(tlo:thi)
   var_lat_adv_res_map = new((/ dimsizes(time), nlat, nlon /), "double")
   var_lat_adv_res_map!0 = "time"
   var_lat_adv_res_map&time = time
   var_lat_adv_res_map!1 = "nlat"
   var_lat_adv_res_map!2 = "nlon"
   var_lat_adv_res_map@long_name = "lateral advective flux (resolved)"
   var_lat_adv_res_map = 0.0

   ; advective flux variable name
   var_lat_adv_res_map@description = "Int_z{-Div[<"+var_name1+">, <"+var_name2+">]}"
   dz = f1->dz
   zlo = f1->z_w(klo)
   zhi = f1->z_w(khi+1)
   var_lat_adv_res_map@k_range = tostring(klo)+" - "+tostring(khi)
   var_lat_adv_res_map@depth_range = sprintf("%6.1f",zlo/100.)+"m -"+sprintf("%6.1f",zhi/100.)+"m"

   FIELD1e = f1->$var_name1$(tlo:thi, klo:khi, :, :)
   FIELD2n = f2->$var_name2$(tlo:thi, klo:khi, :, :)
      if (isatt(FIELD1e, "scale_factor")) then
         FIELD1e = (/ FIELD1e@scale_factor * FIELD1e /)
         delete(FIELD1e@scale_factor)
      end if
      if (isatt(FIELD2n, "scale_factor")) then
         FIELD2n = (/ FIELD2n@scale_factor * FIELD2n /)
         delete(FIELD2n@scale_factor)
      end if

   FIELD1w = eoshift(FIELD1e,1,-1) 
   FIELD2s = eoshift(FIELD2n,2,-1) 
   work1 = todouble(FIELD1e)
     work1 = (/ work1*conform(work1,volc,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
   work2 = todouble(FIELD1w)
     work2 = (/ work2*conform(work2,volw,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
   work3 = todouble(FIELD2n)	
     work3 = (/ work3*conform(work3,volc,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
   work4 = todouble(FIELD2s)	
     work4 = (/ work4*conform(work4,vols,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
     
   work1 = (work2-work1)+(work4-work3)
   MASK_3D = where(volc.gt.0,1,0)
   work1 = (/ work1*conform(work1,MASK_3D,(/1,2,3/)) /)
   var_lat_adv_res_map = (/ dim_sum_n(work1,1) /)

   return(var_lat_adv_res_map)
end

;*****************************************************************************;
undef("traver_budget_hmix")
function tracer_budget_hmix(fname1:string, fname2:string, \
   nlon:integer, nlat:integer, klo:integer, khi:integer, tlo:integer, thi:integer, vol3d[*][*][*]:double)
;
; horizontal divergence: mixing
;
local var_lat_adv_res_map, \
      FIELD1e, FIELD1w, FIELD2n, FIELD2s, MASK_3D, volc,vols,volw,work, work1, \
      time, work2, work3, work4, dz
begin
   work = vol3d				            ; Tcell_i,j
   work1 = (/ eoshift(work,1,-1) /)		; Tcell_i-1,j
   work2 = (/ eoshift(work,2,-1) /)		; Tcell_i,j-1
   volc = work(klo:khi,:,:)
   volw = work1(klo:khi,:,:)
   vols = work2(klo:khi,:,:)
   delete(work)
   delete(work1)
   delete(work2)

   ; initialize :
   var_name1 = "HDIFE_TEMP"
   var_name2 = "HDIFN_TEMP"
   f1 = addfile(fname1,"r")
   f2 = addfile(fname2,"r")
   time = f1->time(tlo:thi)
   var_lat_mix_res_map = new((/ dimsizes(time), nlat, nlon /), "double")
   var_lat_mix_res_map!0 = "time"
   var_lat_mix_res_map&time = time
   var_lat_mix_res_map!1 = "nlat"
   var_lat_mix_res_map!2 = "nlon"
   var_lat_mix_res_map@long_name = "lateral diffusive flux"
   var_lat_mix_res_map = 0.0

   ; mixective flux variable name
   var_lat_mix_res_map@description = "Int_z{-Div[<"+var_name1+">, <"+var_name2+">]}"
   dz = f1->dz
   zlo = f1->z_w(klo)
   zhi = f1->z_w(khi+1)
   var_lat_mix_res_map@k_range = tostring(klo)+" - "+tostring(khi)
   var_lat_mix_res_map@depth_range = sprintf("%6.1f",zlo/100.)+"m -"+sprintf("%6.1f",zhi/100.)+"m"

   FIELD1e = f1->$var_name1$(tlo:thi, klo:khi, :, :)
   FIELD2n = f2->$var_name2$(tlo:thi, klo:khi, :, :)
      if (isatt(FIELD1e, "scale_factor")) then
         FIELD1e = (/ FIELD1e@scale_factor * FIELD1e /)
         delete(FIELD1e@scale_factor)
      end if
      if (isatt(FIELD2n, "scale_factor")) then
         FIELD2n = (/ FIELD2n@scale_factor * FIELD2n /)
         delete(FIELD2n@scale_factor)
      end if

   FIELD1w = eoshift(FIELD1e,1,-1) 
   FIELD2s = eoshift(FIELD2n,2,-1) 
   work1 = todouble(FIELD1e)
     work1 = (/ work1*conform(work1,volc,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
   work2 = todouble(FIELD1w)
     work2 = (/ work2*conform(work2,volw,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
   work3 = todouble(FIELD2n)	
     work3 = (/ work3*conform(work3,volc,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
   work4 = todouble(FIELD2s)	
     work4 = (/ work4*conform(work4,vols,(/1,2,3/)) /)	; e.g., "degC cm^3/s"
     
   work1 = (work1-work2)+(work3-work4)
   MASK_3D = where(volc.gt.0,1,0)
   work1 = (/ work1*conform(work1,MASK_3D,(/1,2,3/)) /)
   var_lat_mix_res_map = (/ dim_sum_n(work1,1) /)

   return(var_lat_mix_res_map)
end

;*****************************************************************************;
undef("tracer_budget_vert_adv_resolved")
function tracer_budget_vert_adv_resolved (fname1:string,nlon:integer, nlat:integer,\
    nz:integer, klo:integer, khi:integer, tlo:integer, thi:integer, vol3d[*][*][*]:double)
;
begin
    nx = nlon
    ny = nlat
   
    volc = vol3d				; Tcell_i,j
    
    ; k-slice or k-full
    if (khi.eq.nz) then
      khip1 = khi
    else
      khip1 = khi+1
    end if

    f = addfile(fname1, "r")
    time = f->time(tlo:thi)

    ; initialize :
    var_vert_adv_res_map = new((/ dimsizes(time), ny, nx /), "double")
    var_vert_adv_res_map!0 = "time"
    var_vert_adv_res_map&time = time
    var_vert_adv_res_map!1 = "nlat"
    var_vert_adv_res_map!2 = "nlon"
    var_vert_adv_res_map@long_name = "vertical advective flux (resolved)"
    var_vert_adv_res_map = 0.0

    ; advective flux variable name
    var_name = "WTT"
    var_vert_adv_res_map@description = "Int_z{-d[<"+var_name+">]/dz}"
      
    ;var_vert_adv_res_map@units = tracer_budget_var_3d_int_units(f->$var_name$@units)
    dz = f->dz
    zlo = f->z_w(klo)
    zhi = f->z_w(khi+1)
    var_vert_adv_res_map@k_range = tostring(klo)+" - "+tostring(khi)
    var_vert_adv_res_map@depth_range = sprintf("%6.1f",zlo/100.)+"m -"+sprintf("%6.1f",zhi/100.)+"m"
    

    FIELD = f->$var_name$(tlo:thi, :, :, :)
    if (isatt(FIELD, "scale_factor")) then
       FIELD = (/ FIELD@scale_factor * FIELD /)
       delete(FIELD@scale_factor)
    end if

;;  Since WTT has _FillValue masking:
    FIELD = where(ismissing(FIELD),0., FIELD)

    work1 = todouble(FIELD(:,klo,:,:))
      work1 = (/ work1*conform(work1,volc(klo,:,:),(/1,2/)) /)	; e.g., "degC cm^3/s"
    work2 = todouble(FIELD(:,khip1,:,:))
    if (khi.lt.nz) then
        work2 = (/ work2*conform(work2,volc(khip1,:,:),(/1,2/)) /)	; e.g., "degC cm^3/s"
    else
        work2 = (/ work2*0.0 /)	
    end if 
    var_vert_adv_res_map(:,:,:) = (/ work2 - work1 /)		; vertical convergence
   return(var_vert_adv_res_map)
end


;*****************************************************************************;
undef("tracer_budget_dia_vmix")
function tracer_budget_dia_vmix(fname1:string,nlon:integer, nlat:integer,\
    nz:integer, klo:integer, khi:integer, tlo:integer, thi:integer)
;
;  Computes vertical integral of diabatic vertical mixing (DIA_IMPVF_), ie. KPP
;
begin

    f = addfile(fname1, "r")
    time = f->time(tlo:thi)
    TAREA = f->TAREA
    KMT = f->KMT

    ; initialize :
    var_vert_mix_map = new((/ dimsizes(time), nlat, nlon /), "double")
    var_vert_mix_map!0 = "time"
    var_vert_mix_map&time = time
    var_vert_mix_map!1 = "nlat"
    var_vert_mix_map!2 = "nlon"
    var_vert_mix_map@long_name = "vertical (diabatic) mixing flux"
    var_vert_mix_map = 0.0

    ; variable name
    var_name = "DIA_IMPVF_TEMP"
    var_vert_mix_map@description = "Int_z{-d[<"+var_name+">]/dz}"

    ;var_vert_mix_map@units = tracer_budget_var_2d_int_units(f->$var_name$@units)
    dz = f->dz
    zlo = f->z_w(klo)
    zhi = f->z_w(khi+1)
    var_vert_mix_map@k_range = tostring(klo)+" - "+tostring(khi)
    var_vert_mix_map@depth_range = sprintf("%6.1f",zlo/100.)+"m -"+sprintf("%6.1f",zhi/100.)+"m"

    FIELD = f->$var_name$(tlo:thi, :, :, :)
    if (isatt(FIELD, "scale_factor")) then
        FIELD = (/ FIELD@scale_factor * FIELD /)
        delete(FIELD@scale_factor)
    end if

;; Do this since DIA_IMPVF_TEMP is masked with _FillValue
    FIELD = where(ismissing(FIELD),0.,FIELD)
    if (klo.eq.0) then
        FIELD_TOP = FIELD(:,klo,:,:)
        FIELD_TOP = 0.0					; zero diffusive flux across sea surface
    else
        FIELD_TOP = FIELD(:,klo-1,:,:)
    end if
    FIELD_BOT = FIELD(:,khi,:,:)

    MASK_2D = where(klo .lt. KMT, 1.0, 0.0)
    work1 = FIELD_TOP*conform(FIELD_TOP,TAREA*MASK_2D,(/1,2/))      ; e.g., "degC cm^3/s"
    MASK_2D = where(khi .lt. KMT, 1.0, 0.0)
    work2 = FIELD_BOT*conform(FIELD_BOT,TAREA*MASK_2D,(/1,2/))      ; e.g., "degC cm^3/s"

    var_vert_mix_map(:,:,:) = (/ -(work2 - work1) /)	      
    var_vert_mix_map = where(ismissing(var_vert_mix_map),0.0,var_vert_mix_map)

    return(var_vert_mix_map)
end

;*****************************************************************************;
undef("tracer_budget_adi_vmix")
function tracer_budget_adi_vmix(fname1:string,nlon:integer, nlat:integer,\
    nz:integer, klo:integer, khi:integer, tlo:integer, thi:integer, vol3d[*][*][*]:double)
;
;  Computes vertical integral of adiabatic vertical mixing (HDIFB_), ie. GM+Submeso
;
local var_vert_mix_map, fnames, i0, i1, nx, lhi, llo, khip1, file_ind, f, tlen_file, \
      FIELD, FIELD_TOP, FIELD_BOT, work, work1, work2, MASK_2D, dz
begin
    
    f = addfile(fname1,"r")
    time = f->time(tlo:thi)
    volc = vol3d            ; Tcell_i,j
    if (klo.eq.0) then
        klom1 = klo
    else
        klom1 = klo-1
    end if

   ; initialize :
   var_vert_mix_map = new((/ dimsizes(time), nlat, nlon /), "double")
   var_vert_mix_map!0 = "time"
   var_vert_mix_map&time = time
   var_vert_mix_map!1 = "nlat"
   var_vert_mix_map!2 = "nlon"
   var_vert_mix_map@long_name = "vertical (adiabatic) mixing flux"
   var_vert_mix_map = 0.0

    ; variable name
    var_name = "HDIFB_TEMP"
    var_vert_mix_map@description = "Int_z{-d[<"+var_name+">]/dz}"


    ;var_vert_mix_map@units = tracer_budget_var_3d_int_units(f->$var_name$@units)
    dz = f->dz
    zlo = f->z_w(klo)
    zhi = f->z_w(khi+1)
    var_vert_mix_map@k_range = tostring(klo)+" - "+tostring(khi)
    var_vert_mix_map@depth_range = sprintf("%6.1f",zlo/100.)+"m -"+sprintf("%6.1f",zhi/100.)+"m"

    FIELD = f->$var_name$(tlo:thi, :, :, :)
    if (isatt(FIELD, "scale_factor")) then
        FIELD = (/ FIELD@scale_factor * FIELD /)
        delete(FIELD@scale_factor)
    end if

;; Do this since HDIFB_TEMP is masked with _FillValue
    FIELD = where(ismissing(FIELD),0.,FIELD)

    if (klo.eq.0) then
        FIELD_TOP = FIELD(:,klo,:,:)
        FIELD_TOP = 0.0					; zero flux across sea surface
    else
        FIELD_TOP = FIELD(:,klom1,:,:)
    end if
    FIELD_BOT = FIELD(:,khi,:,:)

    work1 = FIELD_TOP*conform(FIELD_TOP,volc(klom1,:,:),(/1,2/))    ; e.g., "degC cm^3/s"
    work2 = FIELD_BOT*conform(FIELD_BOT,volc(khi,:,:),(/1,2/))      ; e.g., "degC cm^3/s"

    var_vert_mix_map(:,:,:) = (/  -(work2 - work1) /)	      

    var_vert_mix_map = where(ismissing(var_vert_mix_map),0.0,var_vert_mix_map)
    return(var_vert_mix_map)
end
;*****************************************************************************;

