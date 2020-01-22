#########################################################################                                  
:FileType          rvt ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using weighted model options                                                            
#------------------------------------------------------------------------

# meteorological forcings
:Gauge
  :Latitude    54.09639
  :Longitude -122.67972
  :Elevation  606.0
  #
  # energy limited  (100% rain+snow = original; dryness index PET/P = 0.734646456)
  :RedirectToFile data_obs/Salmon-River-Near-Prince-George_meteo_daily_1.00_energy-limited.rvt
  #
  # water limited   ( 33% snow+rain; dryness index PET/P = 2.203961407)
  #:RedirectToFile data_obs/Salmon-River-Near-Prince-George_meteo_daily_0.33_water-limited.rvt  
  #
  # not limited     ( 75% snow+rain; dryness index PET/P = 0.979528608)
  #:RedirectToFile data_obs/Salmon-River-Near-Prince-George_meteo_daily_0.75_not-limited.rvt
:EndGauge

# observed streamflow
:RedirectToFile data_obs/Salmon-River-Near-Prince-George_Qobs_daily.rvt
