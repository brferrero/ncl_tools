;
; calculate local steric sea level from in-situ density annual means 
; https://journals.ametsoc.org/doi/full/10.1175/2010JCLI3533.1
;

begin

;*****************************************************************************;
  ; I/O
  ;varin = "RHO"
  ;dirin = "/chuva/wk2/bruno/cesm_lens_slr/B20TRC5CNBDRD/"+varin+"/"
  ;fin = "b.e11.B20TRC5CNBDRD.f09_g16.008.pop.h."+varin+".yr.192001-200512.nc"
  ;filein = dirin + fin
  ;
  ;varout = "ssl"
  ;dirout = "./"
  ;fout = "b.e11.B20TRC5CNBDRD.f09_g16.008.pop.h."+varout+".yr.192001-200512.nc"
  ;fileout = dirout + fout
  ;refyear = 1920

;*****************************************************************************;
  ; reading variables
  in = addfile(filein,"r")
  region_mask = in->REGION_MASK
  tarea = in->TAREA
  dz = in->dz
  ; 4d var
  var0 = in->$varin$    ;(0:10,:,:,:) ; for testing
  ; var0(time,depth,nlat,nlon)
  sz = dimsizes(var0)
  tlong = in->TLONG
  tlat = in->TLAT
  ;rho_0
  rho_zero = 1.026
  time = var0&time
  yyyy = (cd_calendar (time, -1)/100)
  tlong = in->TLONG
  tlat = in->TLAT
  z_t = in->z_t

;*****************************************************************************;
  ; masking areas 
  mask_land = where(region_mask.eq.0,tarea@_FillValue,1.0)
  mask_land@_FillValue = tarea@_FillValue
  
  ; areacello (double)
  areacello = tarea*mask_land
  copy_VarMeta(tarea,areacello)
  
  ; volume mask (double)
  vmask = where(var0(0,:,:,:).eq.var0@_FillValue,tarea@_FillValue,1)
  vmask@_FillValue = tarea@_FillValue
  dv = vmask
  volcello = dv
  
  ; cell volume
  do i=0,sz(1)-1
    dv(i,:,:) = vmask(i,:,:)*areacello*dz(i)    ; cell volume 
    volcello(i,:,:) = vmask(i,:,:)*dz(i)        ; depth
  end do
  copy_VarMeta(var0(0,:,:,:),dv)
  copy_VarMeta(var0(0,:,:,:),volcello)
  dv@long_name = "Ocean Grid-Cell Volume"
  dv@units = "cm3"
  volcello@long_name = dz@long_name
  volcello@units = dz@units
  ; total Volume & Area
  V = sum(dv)
  A = sum(areacello)
 
;*****************************************************************************;
  var1 = new(sz,double)
  var2 = new(sz,double)
  ;loop t dimension
  do i=0,sz(0)-1 ; time   
    var1(i,:,:,:) = var0(i,:,:,:)*dv
    var2(i,:,:,:) = volcello - (volcello*(var0(i,:,:,:)/rho_zero))
  end do
  
  ;global mean rho
  gmrho = dim_sum_n(var1,(/1,2,3/))/V
  
  gmssl = gmrho
  va = V/A
  
  ; global mean steric sea level rise
  gmssl = (va)*(1 - (gmrho/rho_zero))
  
  ; gmssl relative to refyear
  ind_refyear = get1Dindex(yyyy,refyear)
  gmssl = gmssl - gmssl(ind_refyear)
  
  ; local mean steric sea level relative to refyear
  ssl_ref = var2(ind_refyear,:,:,:)
  ssl = tofloat(var2)
  do i=0,sz(0)-1 ; time   
    ssl(i,:,:,:) =  tofloat(var2(i,:,:,:) - ssl_ref)
  end do
  
  ; metadata
  copy_VarMeta(var0(:,0,0,0),gmssl)
  copy_VarMeta(var0,ssl)
  gmssl@long_name = "Global mean steric sea level"
  gmssl@units = "cm"
  ssl@long_name = "Relative steric sea level rise"
  ssl@units = "cm"

;*****************************************************************************;
  ; write to netcdf

  system("rm "+fileout)
  setfileoption("nc","Format","NetCDF4Classic")
  ; setfileoption("nc","CompressionLevel",3)    ; slow
  out = addfile(fileout,"c")
  out->ssl = ssl
  out->gmssl = gmssl
  out->TLONG = tlong
  out->TLAT = tlat
  out->dz = volcello
  out->areacello = areacello
end
