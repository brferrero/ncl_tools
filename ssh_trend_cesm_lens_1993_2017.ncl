;
; ssh linear trend for CESM/LENS 
; individual member trend (panel) and trend ensemble mean
;

begin
;*****************************************************************************; 
  filein = new(2,string)
;-[ensemble | 35, nlat | 384, nlon |320]
  ssh_trend = new((/35,384,320/),float)
  ssh_trend!0 = "ensemble"
  do fn = 1,35
    fnstr = systemfunc("printf %03d\\n "+ tostring(fn))
    scenario = "B20TRC5CNBDRD"
    dirin = "/chuva/wk2/bruno/cesm_lens_slr/"+scenario+"/SSH/"
    f = "b.e11."+scenario+".f09_g16."+fnstr+".pop.h.SSH.yr.*.nc" 
    filein(0) = systemfunc("ls "+dirin+f)
;--  
    scenario = "BRCP85C5CNBDRD"
    dirin = "/chuva/wk2/bruno/cesm_lens_slr/"+scenario+"/SSH/"
    f = "b.e11."+scenario+".f09_g16."+fnstr+".pop.h.SSH.yr.*.nc" 
    filein(1) = systemfunc("ls "+dirin+f)
    ;print(filein)
    in = addfiles(filein,"r")
    ListSetType (in, "cat") 
    time =in[:]->time
    tsz = dimsizes(time)
    i2017 = tsz - (2100 - 2017) - 1
    i1993 = i2017 - (2017-1993)
    yyyy = cd_calendar(time(i1993:i2017), -1 )/100
    ;print(yyyy)
;-read 25yrs ssh    
    ssh = in[:]->SSH(i1993:i2017,:,:)
    ;printVarSummary(ssh)

;*****************************************************************************; 
;-calculate regression coeficient with sig995
    rc = regCoef(yyyy,ssh(nlat|:,nlon|:,time|:))                        
     
    rc@long_name = "SSH Trend"
    rc@units         = ssh@units+"/year"
    copy_VarCoords(ssh(0,:,:), rc)   ; copy lat,lon coords
    rc@coordinates = "TLAT TLONG"
    delete_VarAtts(rc,(/"nptxy","rstd","yintercept","tval"/))
    ssh_trend(fn-1,:,:) = rc
    delete(rc)
    delete(ssh)
    delete(time)
  end do
  ; read coords
  tlong = in[0]->TLONG
  tlat = in[0]->TLAT
  
;*****************************************************************************; 
  ; prepare to plot
  ssh_trend@lon2d  = tlong
  ssh_trend@lat2d  = tlat
  ssh_trend = ssh_trend*10;*nyrs               
  ssh_trend@units     = "mm/year"
  ens_trend = dim_avg_n_Wrap(ssh_trend,0)
  ;printVarSummary(ens_trend)
  figname = "ssh_trend_cesm_lens_1993_2017"  

;*****************************************************************************; 
  wks  = gsn_open_wks("ps",figname)         ; send graphics to PNG file
   
  ; minimalist plot 
  res                       = True
  res@gsnMaximize           = True             ; make large
  res@gsnAddCyclic          = True
  res@cnFillOn              = True             ; turn on color
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@cnLinesOn             = False            ; turn off contour lines
  res@cnLineLabelsOn        = False            ; turn off contour line labels
  res@cnFillMode            = "RasterFill"

  res@cnLevelSelectionMode  = "ManualLevels"   ; set manual contour levels
  res@cnMinLevelValF        = -3.0            ; set min contour level
  res@cnMaxLevelValF        =  3.0            ; set max contour level
  res@cnLevelSpacingF       =  0.1            ; set contour interval
  res@lbLabelStride         =  10
  res@cnFillPalette         = "BlueWhiteOrangeRed"
  res@lbBoxSeparatorLinesOn = False
  res@lbBoxEndCapStyle      = "TriangleBothEnds" 
  res@lbBoxLinesOn          = False

  res@mpFillOn              = True            ; default background gray
  res@mpCenterLonF          = -60

  res@mpProjection          = "Robinson"
  res@mpPerimOn             = False
  res@mpLandFillColor       = "grey40"
  res@tiMainString          = "CESM/LENS ssh trend: 25yr ensemble mean"
  res@gsnLeftString = "Linear trend"
  
;*****************************************************************************; 
  sres = True
  sres@gsnDraw               = False
  sres@gsnFrame              = False
  sres@cnLevelSelectionMode = "ExplicitLevels" ; 
  sres@cnLevels             = 0;(/-3,0,3/)             
  sres@cnFillOn             = False     
  sres@cnLineLabelsOn       = False
  sres@cnLinesOn            = True   
  sres@cnLineThicknessF     = 2.5
  sres@cnLineDashPattern    = 1
  sres@cnInfoLabelOn        = False   ; turn off contour label

;*****************************************************************************; 
;-Panel plot: True 
  resP = True
  if resP.eq.False then
    res@gsnDraw = True
    res@gsnFrame = True
    if sres then
        sres@gsnDraw = True
        sres@gsnFrame = True
        plot = gsn_csm_contour_map_overlay(wks,ens_trend,ens_trend,res,sres)
    else
        plot = gsn_csm_contour_map(wks,ens_trend,res)
    end if
  else
      print("Panel plot")
      res@tiMainString = ""
      res@gsnCenterString = ""
      res@gsnRightString = ""
      res@gsnLeftString = ""
      res@lbLabelBarOn = False
      plot = new(35,graphic)
      do i=0,8
        fnstr = systemfunc("printf %03d\\n "+ tostring(i+1))
        res@gsnCenterString = "member #"+fnstr
        plot(i) = gsn_csm_contour_map_overlay(wks,ssh_trend(i,:,:),ssh_trend(i,:,:),res,sres)   
      end do
      resP@gsnPanelMainString   =  "CESM/LENS ssh trend: 35 ensemble members"
      resP@gsnPanelLabelBar     = True                ; add common colorbar
      resP@lbLabelFontHeightF   = 0.007               ; make labels smaller
      resP@lbLabelStride        =  10
      resP@lbBoxSeparatorLinesOn= False
      resP@lbBoxLinesOn         = False
      gsn_panel(wks,plot,(/3,3/),resP) 
  end if
  delete(wks)
  system("convert -density 144 -rotate 270 -trim +repage -border 8 -bordercolor white "+figname+".ps "+figname+ ".png")
end
