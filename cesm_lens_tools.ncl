; FUNCTIONS
;
; some tools to handle CESM/LENS data
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
;*****************************************************************************; 
;-get files
undef ("get_filemember")
function get_filemember (basedir, scenario, freq, realm, ens_member, varname, single)
local fnstr, dirin, f
begin
    fnstr = systemfunc("printf %03d\\n "+ tostring(ens_member))
    dirin = basedir + "/"+scenario+"/" + realm +"/" + freq + "/" +varname + "/"
    if single.eq."xaer" then
        f = "b.e11."+scenario+".f09_g16.xaer."+fnstr+ ".pop.h."+varname+".yr.*.nc"
    elseif single.eq."xghg" then
        f = "b.e11."+scenario+".f09_g16.xghg."+fnstr+ ".pop.h."+varname+".yr.*.nc"
    else
        f = "b.e11."+scenario+".f09_g16."+fnstr+ ".pop.h."+varname+".yr.*.nc"
    end if
    return (systemfunc("ls "+dirin+f))
end

;*****************************************************************************;
;-get time index (freq: year)
undef ("get_timeindex")
function get_timeindex (time, year)
local yyyy
begin
    yyyy = cd_calendar(time, -1 )/100
    return (closest_val(year,yyyy))
end

;*****************************************************************************; 
;-calculate regression coeficient with sig995
undef ("linear_trend")
function linear_trend (time, x)
local rc
begin
    rc = regCoef(time,x(nlat|:,nlon|:,time|:))                        
    rc@long_name = x@long_name + " Linear Trend"
    rc@units         = x@units+"/year"
    copy_VarCoords(x(0,:,:), rc)   ; copy lat,lon coords
    rc@coordinates = "TLAT TLONG"
    delete_VarAtts(rc,(/"nptxy","rstd","yintercept","tval"/))
    return rc
end

;*****************************************************************************;
; remap to rectilinear grid
undef ("pop2rect")
function pop2rect(xpop,resolution)
begin
    if resolution.eq.0 then
    ; 360x180
        return (PopLatLon(xpop,"gx1v6","1x1d","bilin","da","100716"))
    elseif resolution.eq.1 then
    ; 192x288
        return (PopLatLon(xpop,"gx1v6","fv0.9x1.25","bilin","da","090309"))
    else
    ; 96x144
        return (PopLatLon(xpop,"gx1v6","fv1.9x2.5","bilin","da","090206"))
    end if
end

;*****************************************************************************; 
; devolve uma mascara o a bacia de interesse
undef ("get_pop_ocean_basin") 
function get_pop_ocean_basin (basin,remap)
local fgrid, in, region,basin_area
begin
    fgrid="/wk2/bruno/cesm_lens_slr/analysis/grid_fx_ocn/cesm_fx_grid_ocn.nc"
    in = addfile(fgrid,"r")
    region = in->region_mask

    if basin.eq."atlantic" then
        basin_area =  where(region.eq.6,1.0,default_fillvalue("float"))
    elseif basin.eq."pacific" then
        basin_area =  where(region.eq.2,1.0,default_fillvalue("float"))
    elseif basin.eq."southern" then
        basin_area =  where(region.eq.1,1.0,default_fillvalue("float"))
    elseif basin.eq."indian" then
        basin_area =  where(region.eq.3,1.0,default_fillvalue("float"))
    elseif basin.eq."oceans" then
        basin_area =  where(region.eq.0,1.0,default_fillvalue("float"))
    else
        basin_area = region
    end if
    basin_area@_FillValue = default_fillvalue("float")
    return basin_area
end

;*****************************************************************************;
; recebe TAREA e REGION_MASK e devolve area somente para oceano
undef ("get_ocean_area") 
function get_ocean_area (tarea, region)
local ocean_area
begin
    ocean_area = where(region.eq.0 , default_fillvalue("double"), tarea)
    return (ocean_area)
end

;*****************************************************************************;
; recebe OCEAN_AREA, e VAR3D e devolve o global area mean
undef ("global_area_mean") 
function global_area_mean (var3d, ocean_area)
local wgt, wgt3d, x
begin
    wgt = ocean_area/sum(ocean_area)
    wgt3d = conform_dims(dimsizes(var3d),wgt,(/1,2/))
    x = var3d*wgt3d
    copy_VarMeta(var3d,x)
    return (dim_sum_n_Wrap(x,(/1,2/)))
end

;*****************************************************************************;
undef ("get_vol3d") 
function get_vol3d (ocean_area, dz, kmt)
local vol3d, sz, k
begin
    sz = dimsizes(ocean_area)
    nz = dimsizes(dz)
    vol3d = new((/ nz, sz(0), sz(1) /), "double")
    vol3d!0 = "nz_G"
    vol3d!1 = "nlat_G"
    vol3d!2 = "nlon_G"
    vol3d@long_name = "Tcell volume"
    vol3d@units = "cm^3"
    vol3d = 0.
    do k = 0, nz-1
        MASK_2D = where(k .lt. kmt, 1.0, 0.0)
        ;MASK_2D = where(k .lt. kmt, 1.0, default_fillvalue("double"))
        if (sum(MASK_2D) .gt. 0.0) then
            vol3d(k,:,:) = (/ (dz(k) * ocean_area * MASK_2D) /)
        end if
   end do
  return(vol3d)
end

;*****************************************************************************;
; recebe VOL3D e VAR4D e devolve global volume mean
undef ("global_volume_mean") 
function global_volume_mean (var4d, vol3d)
local wgt, wgt4d, x, V
begin
    V = sum(vol3d)
    wgt = vol3d/V
    ;print(sum(wgt))
    wgt4d = conform_dims(dimsizes(var4d),wgt,(/1,2,3/))
    x = var4d*wgt4d
    copy_VarMeta(var4d,x)
    return (dim_sum_n_Wrap(x,(/1,2,3/)))
end

;*****************************************************************************;
