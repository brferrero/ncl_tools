;
; Volume-weighted mean for any 4D variable
; run with: run_volume_wgt_mean.sh
;

begin

;*****************************************************************************;
  ; I/O
  ;varin = "RHO"
  ;dirin = "/chuva/db2/CESM-LENS/fully_coupled/mon/pop/"+varin+"/B20TRC5CNBDRD/"
  ;fin = "b.e11.B20TRC5CNBDRD.f09_g16.008.pop.h."+varin+".192001-200512.nc"
  ;filein = dirin + fin
  ;
  ;varout = "gmrho"
  ;dirout = "./"
  ;fout = "b.e11.B20TRC5CNBDRD.f09_g16.008.pop.h."+varout+".yr.192001-200512.nc"
  ;fileout = dirout + fout
  ;basin="full" 

;*****************************************************************************;
  ; reading variables
  in = addfile(filein,"r")
  region_mask = in->REGION_MASK
  tarea = in->TAREA
  dz = in->dz
  ; 4d var
  var0 = in->$varin$(0:3,:,:,:) ; for testing
  ; var0(time,depth,nlat,nlon)
  sz = dimsizes(var0)

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
  
  ; cell thickness
  do i=0,sz(1)-1
    dv(i,:,:) = vmask(i,:,:)*areacello*dz(i) 
  end do
  copy_VarMeta(var0(0,:,:,:),dv)
  dv@long_name = "Ocean Grid-Cell Volume"
  dv@units = "cm3"
  V = sum(dv)
  
;*****************************************************************************;
  ; volume wgt mean
  dv_4d = conform_dims(sz,dv,(/1,2,3/))
  ; var*dV
  var1 = var0*dv_4d
  ; volume sum
  gmvar = dim_sum_n_Wrap(var1,(/1,2,3/))/V
  copy_VarMeta(var0(:,0,0,0),gmvar)
  
;*****************************************************************************;
  ; write to netcdf
  system("rm -f " + fileout)
  out = addfile (fileout,"c")
  out->$varout$ = gmvar
  ;out->area = areacello
end
