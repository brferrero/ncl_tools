load "tracer_budget_tools.ncl"

;*****************************************************************************;

begin
    nlon = 320
    nlat = 384
    nz = 60
    klo = 0
    khi = 25
    tlo = 491
    thi = 610
    dir_root = "/chuva/db2/CESM-LENS/download/budget/"
    fout = "temp_budget.nc"
    tracer_name = "TEMP"

;*****************************************************************************;
; temp_tend
    vname = "TEMP"
    fname = "/chuva/db2/CESM-LENS/fully_coupled/mon/pop/"+vname+"/B20TRC5CNBDRD/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    vssh = "SSH" 
    fssh = "/chuva/db2/CESM-LENS/fully_coupled/mon/pop/"+vssh+"/B20TRC5CNBDRD/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vssh+".192001-200512.nc"
    time_bound = get_tb(fname)
    print("Computing vol3d")
    vol3d = tracer_budget_vol3d(fname, nlon, nlat, nz)
    print("Computing temp_zint_map")
    temp_zint_map = tracer_budget_var_3d_zint_map(vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False, True, fssh)
    print("Computing temp_tend")
    temp_tend = tracer_budget_tend_appr(time_bound(tlo:thi,:), temp_zint_map, False)
    
    system("rm -f " + fout)
    out = addfile(fout,"c")
    out->vol3d = vol3d
    out->temp_zint_map = temp_zint_map
    out->temp_tend = temp_tend
    out->time_bound = time_bound(tlo:thi,:)

;*****************************************************************************;
; temp_lat_adv_res
    fname1 = dir_root + "UET/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.UET.192001-200512.nc"
    fname2 = dir_root + "VNT/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.VNT.192001-200512.nc"
    print("Computing temp_lat_adv_res")
    temp_lat_adv_res = tracer_budget_lat_adv_resolved(fname1, fname2,nlon, nlat, klo, khi, tlo, thi, vol3d)
    out->temp_lat_adv_res = temp_lat_adv_res

;*****************************************************************************;
; temp_lat_mix
    fname1 = dir_root + "HDIFE_TEMP/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.HDIFE_TEMP.192001-200512.nc"
    fname2 = dir_root + "HDIFN_TEMP/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h.HDIFN_TEMP.192001-200512.nc"
    print("Computing temp_lat_mix")
    temp_lat_mix = tracer_budget_hmix(fname1, fname2, nlon, nlat, klo, khi, tlo, thi, vol3d)
    out->temp_lat_mix = temp_lat_mix
;*****************************************************************************;
; temp_vert_adv_res
    vname = "WTT"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Computing temp_vert_adv_res")
    temp_vert_adv_res = tracer_budget_vert_adv_resolved(fname, nlon, nlat, nz, klo, khi, tlo, thi, vol3d)
    out->temp_vert_adv_res = temp_vert_adv_res

;*****************************************************************************;
; temp_dia_vmix
    vname = "DIA_IMPVF_TEMP" 
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Computing temp_dia_vmix")
    temp_dia_vmix = tracer_budget_dia_vmix(fname, nlon, nlat, nz, klo, khi, tlo, thi)
    out->temp_dia_vmix = temp_dia_vmix

;*****************************************************************************;
; temp_adi_vmix
    vname = "HDIFB_TEMP"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Computing temp_adi_vmix")
    temp_adi_vmix = tracer_budget_adi_vmix(fname, nlon, nlat, nz, klo, khi, tlo, thi, vol3d)
    out->temp_adi_vmix = temp_adi_vmix

;*****************************************************************************;
; temp_kpp_src
    vname = "KPP_SRC_TEMP"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Computing kpp_src_temp_integral")
    temp_kpp_src = tracer_budget_var_3d_zint_map(vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False, False, fssh)
    out->temp_kpp_src = temp_kpp_src

;*****************************************************************************;
; sfluxes
; get surface flux related terms
; NOTE: right now, is hardwired for tracer="TEMP"
   
    vname = "QFLUX"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_qflux = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)

    vname = "SHF"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_shf = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)
    
    vname = "SENH_F"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_senh_f = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)
    
    vname = "LWDN_F"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_lwdn_f = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)
    
    vname = "LWUP_F"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_lwup_f = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)
    
    vname = "MELTH_F"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_melth_f = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)

    vname = "SHF_QSW"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_shf_qsw = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)

    vname = "EVAP_F"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_evap_f = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)
    
    vname = "SNOW_F"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_snow_f = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)

    vname = "IOFF_F"
    fname = dir_root+vname+"/b.e11.B20TRC5CNBDRD.f09_g16.004.pop.h."+vname+".192001-200512.nc"
    print("Compute " + vname)
    temp_ioff_f = tracer_budget_sflux(tracer_name, vname,fname,tlo,thi,nlon, nlat, klo, khi, vol3d, False,False,fssh)

;*****************************************************************************;
; net shf tendency
    print("computing net shf tendency")
    temp_shf_tend = temp_lat_adv_res 
    temp_shf_tend = (/ temp_qflux +  temp_senh_f + temp_lwdn_f + temp_lwup_f + \
	temp_melth_f + temp_shf_qsw + temp_evap_f + temp_snow_f + temp_ioff_f + temp_kpp_src /)
    out->temp_shf_tend = temp_shf_tend

;*****************************************************************************;
; anomalies
    print("Computing climatological anomalies")
    out->temp_tend_anom = rmMonAnnCycTLL(temp_tend)
    out->temp_lat_adv_res_anom = rmMonAnnCycTLL(temp_lat_adv_res)
    out->temp_lat_mix_anom = rmMonAnnCycTLL(temp_lat_mix)
    out->temp_vert_adv_res_anom = rmMonAnnCycTLL(temp_vert_adv_res)
    out->temp_dia_vmix_anom = rmMonAnnCycTLL(temp_dia_vmix)
    out->temp_adi_vmix_anom = rmMonAnnCycTLL(temp_adi_vmix)
    out->temp_shf_tend_anom = rmMonAnnCycTLL(temp_shf_tend)
;*****************************************************************************;
end
