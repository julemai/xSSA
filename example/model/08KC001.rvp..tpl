#########################################################################                                  
:FileType          rvp ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using weighted model options                                                              
#------------------------------------------------------------------------                                 
#

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x25" and "par_x14" and "par_x10" wouldn't be detectable)
#    para_sum_x24_x25 = {dpar[sum_x24_x25]} =  para_x24 + para_x25 = {par[x24]} + {par[x25]}
#    para_sum_x13_x14 = {dpar[sum_x13_x14]} =  para_x13 + para_x14 = {par[x13]} + {par[x14]}
#    para_sum_x09_x10 = {dpar[sum_x09_x10]} =  para_x09 + para_x10 = {par[x09]} + {par[x10]}
#    para_pow_x04     = {dpar[pow_x04]}     =  10^(para_x04)       = 10^{par[x04]}
#    para_pow_x11     = {dpar[pow_x11]}     =  10^(para_x11)       = 10^{par[x11]}

#-----------------------------------------------------------------
# Soil Classes
#-----------------------------------------------------------------
:SoilClasses
  :Attributes,
  :Units,
  TOPSOIL,
  PHREATIC,    
:EndSoilClasses

#-----------------------------------------------------------------
# Land Use Classes
#-----------------------------------------------------------------
:LandUseClasses, 
  :Attributes,        IMPERM,    FOREST_COV, 
       :Units,          frac,          frac, 
       FOREST,           0.0,           1.0,    
:EndLandUseClasses

#-----------------------------------------------------------------
# Vegetation Classes
#-----------------------------------------------------------------
:VegetationClasses, 
  :Attributes,        MAX_HT,       MAX_LAI, MAX_LEAF_COND, 
       :Units,             m,          none,      mm_per_s, 
       FOREST,             4,             5,             5,     
:EndVegetationClasses

#-----------------------------------------------------------------
# Soil Profiles
#-----------------------------------------------------------------
:SoilProfiles
         LAKE, 0
         ROCK, 0
  DEFAULT_P, 2, TOPSOIL, {par[x29]}, PHREATIC, {par[x30]}, 
# DEFAULT_P, 2, TOPSOIL,      x(29), PHREATIC,      x(30), 
:EndSoilProfiles

#-----------------------------------------------------------------
# Terrain Classes
#-----------------------------------------------------------------
:TerrainClasses
  :Attributes,        hillslope_len, drainage_dens,            lambda,       
       :Units,                   ??,            ??,                ??
    DEFAULT_T,                  1.0,           1.0,        {par[x07]}    
#                                                     TOPMODEL_LAMBDA x(7)  
:EndTerrainClasses

#-----------------------------------------------------------------
# Global Parameters
#-----------------------------------------------------------------
:GlobalParameter         SNOW_SWI_MIN {par[x13]}            # x(13)    
:GlobalParameter         SNOW_SWI_MAX {dpar[sum_x13_x14]}   # x(13)+x(14) 
:GlobalParameter     SWI_REDUCT_COEFF {par[x15]}            # x(15)
:GlobalParameter             SNOW_SWI {par[x19]}            # x(19)
#:GlobalParameter      TOC_MULTIPLIER 1.0                   # 

#-----------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------
:SoilParameterList
  :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION,  BASEFLOW_COEFF,      B_EXP,   HBV_BETA, MAX_BASEFLOW_RATE, BASEFLOW_N,      FIELD_CAPACITY,   SAT_WILT,
       :Units,               -,             1/d,               -,             1/d             
      TOPSOIL,             1.0,      {par[x28]},      {par[x08]}, {dpar[pow_x04]}, {par[x02]}, {par[x03]},        {par[x05]}, {par[x06]}, {dpar[sum_x09_x10]}, {par[x09]},
     PHREATIC,             1.0,             0.0,             0.0, {dpar[pow_x11]},        0.0,        0.0,               0.0, {par[x12]},                 0.0,        0.0,   
 #    TOPSOIL,             1.0,           x(28),           x(08),           x(04),      x(02),      x(03),             x(05),      x(06),         x(09)+x(10),      x(09),
 #   PHREATIC,             1.0,             0.0,             0.0,           x(11),        0.0,        0.0,               0.0,      x(12),                 0.0,        0.0,
:EndSoilParameterList

#-----------------------------------------------------------------
# Land Use Parameters
#-----------------------------------------------------------------
:LandUseParameterList
  :Parameters, MIN_MELT_FACTOR,     MAX_MELT_FACTOR,    DD_MELT_TEMP,  DD_AGGRADATION, REFREEZE_FACTOR, REFREEZE_EXP, DD_REFREEZE_TEMP, HMETS_RUNOFF_COEFF,
       :Units,          mm/d/C,              mm/d/C,               C,            1/mm,          mm/d/C,            -,                C,                  -,
    [DEFAULT],      {par[x24]}, {dpar[sum_x24_x25]},      {par[x26]},      {par[x27]},      {par[x18]},   {par[x17]},       {par[x16]},         {par[x01]},
#                        x(24),         x(24)+x(25),           x(26),           x(27),           x(18),        x(17),            x(16),              x(01),        
:EndLandUseParameterList
:LandUseParameterList
  :Parameters,   GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,    FOREST_SPARSENESS,
       :Units,             -,               -,               -,               -,                    -,
    [DEFAULT],    {par[x20]},       {par[x21]},      {par[x22]},      {par[x23]},                 0.0,
    #                  x(20),           x(21),           x(22),           x(23),                  0.0,
:EndLandUseParameterList

#-----------------------------------------------------------------
# Vegetation Parameters
#-----------------------------------------------------------------
:VegetationParameterList
  :Parameters,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,    SAI_HT_RATIO
       :Units,               -,               -,               -
    [DEFAULT],             0.0,             0.0,             0.0
:EndVegetationParameterList

:SeasonalRelativeLAI
  FOREST, 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
:EndSeasonalRelativeLAI
:SeasonalRelativeHeight
  FOREST, 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
:EndSeasonalRelativeHeight
