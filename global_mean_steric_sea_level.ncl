;
; global mean steric sea level from annual means 
; as in Fig. 3/4. from Griffies et. al, 2014
; https://journals.ametsoc.org/doi/full/10.1175/2010JCLI3533.1
;

begin

;*****************************************************************************;
  ; I/O
  varin = "RHO"
  dirin = "/chuva/wk2/bruno/cesm_lens_slr/B20TRC5CNBDRD/"+varin+"/"
  fin = "b.e11.B20TRC5CNBDRD.f09_g16.001.pop.h."+varin+".yr.185001-200512.nc"
  filein = dirin + fin
  ;
  varout = "gmssl"
  dirout = "./"
  fout = "b.e11.B20TRC5CNBDRD.f09_g16.008.pop.h."+varout+".yr.185001-200512.nc"
  fileout = dirout + fout
  basin="full" 

;*****************************************************************************;
  ; reading variables
  in = addfile(filein,"r")
  region_mask = in->REGION_MASK
  tarea = in->TAREA
  dz = in->dz
  ; 4d var
  var0 = in->$varin$    ;(0:3,:,:,:) ; for testing
  ; var0(time,depth,nlat,nlon)
  sz = dimsizes(var0)
  tlong = in->TLONG
  tlat = in->TLAT
  ;rho_0
  rho_zero = 1.026
  time = var0&time
  yyyy = (cd_calendar (time, -1)/100)

;*****************************************************************************;
  ; masking areas 
  mask_land = where(region_mask.eq.0,tarea@_FillValue,1.0)
  mask_land@_FillValue = tarea@_FillValue
  
  ; Southern Ocean - 1 ; Pacific - 2 ; Indian - 3; Atlantic - 6
  if basin.eq."all" then
    mask_area2d = where((region_mask.eq.1).or.(region_mask.eq.2).or.(region_mask.eq.3).or.(region_mask.eq.6)\
                    ,1.0,tarea@_FillValue)
    print("all basins")
  else if basin.eq."so" then
    mask_area2d = where(region_mask.eq.1,1.0,tarea@_FillValue)
    print("southern ocean")
  else if basin.eq."pac" then
    mask_area2d = where(region_mask.eq.2,1.0,tarea@_FillValue)
    print("pacific ocean")
  else if basin.eq."ind" then
    mask_area2d = where(region_mask.eq.3,1.0,tarea@_FillValue)
    print("indian ocean")
  else if basin.eq."atl" then
    mask_area2d = where(region_mask.eq.6,1.0,tarea@_FillValue)
    print("atlantic ocean")
  else if basin.eq."full" then
    mask_area2d = mask_land
    print("full ocean area")
  end if
  end if
  end if
  end if
  end if
  end if

  mask_area2d@_FillValue = tarea@_FillValue

  ; areacello (double)
  areacello = tarea*mask_land*mask_area2d
  copy_VarMeta(tarea,areacello)
  
  ; volume mask (double)
  vmask = where(var0(0,:,:,:).eq.var0@_FillValue,tarea@_FillValue,1)
  vmask@_FillValue = tarea@_FillValue
  dv = vmask
  
  ; cell volume
  do i=0,sz(1)-1
    dv(i,:,:) = vmask(i,:,:)*areacello*dz(i) 
  end do
  copy_VarMeta(var0(0,:,:,:),dv)
  dv@long_name = "Ocean Grid-Cell Volume"
  dv@units = "cm3"
  ; total Volume & Area
  V = sum(dv)
  A = sum(areacello)
 
;*****************************************************************************;
  var1 = new(sz,double)
  ;loop t dimension
  do i=0,sz(0)-1 ; time   
    var1(i,:,:,:) = var0(i,:,:,:)*dv
  end do
  
  ;global mean rho
  gmrho = dim_sum_n(var1,(/1,2,3/))/V
  
  dh = gmrho
  va = V/A
  
  ; global mean steric sea level rise
  dh = (va)*(1 - (gmrho/rho_zero))
  
  ; gmssl relative to 1920
  ind_1920 = get1Dindex(yyyy,1920)
  dh = dh - dh(ind_1920)
  copy_VarMeta(var0(:,0,0,0),dh)
  
  dh@long_name = "Global mean steric sea level"
  dh@units = "cm"
  
;*****************************************************************************;
  ; write to netcdf

  print(fout)
  system("rm "+fout)
  out=addfile(fout,"c")
  out->$varout$ = dh
end
