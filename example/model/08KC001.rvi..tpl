#########################################################################                                  
:FileType          rvi ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using weighted model options                                                        
#------------------------------------------------------------------------
#
:StartDate               1989-01-01 00:00:00 # 1954-01-01 00:00:00 
:EndDate                 2010-12-31 00:00:00  
:EvaluationTime          1991-01-01 00:00:00                                                        
# :Duration              20819                                                                                   
:TimeStep                1.0                                                                                   
:Method                  ORDERED_SERIES 

:PotentialMeltMethod     POTMELT_HMETS
:RainSnowFraction        RAINSNOW_DATA
:SWRadiationMethod       SW_RAD_NONE         # no radiation is faster
:Evaporation             PET_DATA
:CatchmentRoute          ROUTE_DUMP
:Routing                 ROUTE_NONE 
:SoilModel               SOIL_TWO_LAYER

:Alias DELAYED_RUNOFF CONVOLUTION[1] 

:HydrologicProcesses
  :Precipitation   RAVEN_DEFAULT                         ATMOS_PRECIP   MULTIPLE 
  :ProcessGroup #infiltration group
                :Infiltration    INF_HMETS               PONDED_WATER   MULTIPLE 
                :Infiltration    INF_VIC_ARNO            PONDED_WATER   MULTIPLE 
                :Infiltration    INF_HBV                 PONDED_WATER   MULTIPLE 
  :EndProcessGroup {weights[process1]}
                  :Overflow      OVERFLOW_RAVEN          SOIL[0]        DELAYED_RUNOFF
  :ProcessGroup #quickflow group
                :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[0]        SURFACE_WATER   # interflow, really
                :Baseflow        BASE_VIC                SOIL[0]        SURFACE_WATER
                :Baseflow        BASE_TOPMODEL           SOIL[0]        SURFACE_WATER
  :EndProcessGroup {weights[process2]}
  :Percolation                   PERC_LINEAR             SOIL[0]        SOIL[1]         # recharge
    :Overflow                    OVERFLOW_RAVEN          SOIL[1]        DELAYED_RUNOFF
  :ProcessGroup #evaporation group
                :SoilEvaporation SOILEVAP_ALL            SOIL[0]        ATMOSPHERE      # AET
                :SoilEvaporation SOILEVAP_TOPMODEL       SOIL[0]        ATMOSPHERE      # AET
  :EndProcessGroup {weights[process3]}
  :Convolve                      CONVOL_GAMMA            CONVOLUTION[0] SURFACE_WATER   # 'surface runoff'
  :Convolve                      CONVOL_GAMMA_2          DELAYED_RUNOFF SURFACE_WATER   # 'delayed runoff'
  :ProcessGroup #quickflow group
                :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[1]        SURFACE_WATER
                :Baseflow        BASE_POWER_LAW          SOIL[1]        SURFACE_WATER
  :EndProcessGroup {weights[process4]} 
  :ProcessGroup #snow balance group
                :SnowBalance     SNOBAL_HMETS            MULTIPLE       MULTIPLE
                :SnowBalance     SNOBAL_SIMPLE_MELT      SNOW           PONDED_WATER
                :SnowBalance     SNOBAL_HBV              MULTIPLE       MULTIPLE
                #:SnowBalance     SNOBAL_GAWSER           MULTIPLE       MULTIPLE
  :EndProcessGroup {weights[process5]} 
:EndHydrologicProcesses

#:CreateRVPTemplate

#---------------------------------------------------------
# Output Options
#
# :WriteForcingFunctions
# :WriteNetcdfFormat

# Accumulated Infiltration volume
:CustomOutput DAILY AVERAGE Between:PONDED_WATER.And.SOIL[0] BY_BASIN

# Actual snow on ground over whole basin 
:CustomOutput DAILY AVERAGE SNOW BY_BASIN

# Accumulated Snow melt volume
:CustomOutput DAILY AVERAGE From:SNOW BY_BASIN

:EvaluationMetrics NASH_SUTCLIFFE RMSE
:SilentMode
:DontWriteWatershedStorage
#
