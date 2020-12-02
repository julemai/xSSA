#!/usr/bin/env python

# Copyright 2019-2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#
# run with:
#     run sa_model_multiple_processes.py 

#!/usr/bin/env python
from __future__ import print_function
from __future__ import with_statement

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../../lib')

import numpy as np
import copy
import sobol
import sobol_index
import pickle
import json
import msgpack
import msgpack_numpy as m
import netCDF4 as nc4
from autostring import astr
import PieShareDistribution as psd
from collections import OrderedDict

# force all msgpack serialization and deserialization routines (and other packages that use them) to become numpy-aware
m.patch()

__all__ = ['sa_model_multiple_processes']

"""
Perform Sensitivity Analysis for models with multiple process options. 
Derives Sensitivity index estimates for:
- each parameter          (p1, p2, p3, p4, ...)
- each process option     (A1, A2, ..., B1, B2, ...)
- each process            (A, B, C)

History
-------
Written,  JM, Jun 2019
"""

def sa_model_multiple_processes(paras_per_option, para_ranges, model_function, basin_prop, constants=None, nsets=None, budget=None,
                                    save_pkl=None,
                                    save_json=None,
                                    save_msgpack=None,
                                    save_nc4=None):
    """
        This function that estimates the Sobol' sensitivity estimates for models with mutiple process options. 
        The options and the parameters of those options are given in a nested list 'paras_per_option'. 
        Further, the range of each parameter needs to be given and a function that returns model outputs 
        when a set of parameters and weights are given. The weights are used to weight all the process option 
        outputs. Hence, the returned model output is a weighted model output. The sampling of all weights and 
        parameters is done internally in this method. Sobol' sequences are used for this purpose.

        Definition
        ----------
        def sa_model_multiple_processes(paras_per_option, para_ranges, model_function, constants=None, nsets=None)


        Input               Format               Description
        -----               -----                -----------
        paras_per_option    list of lists        lists the parameters of each process option
                                                 example: process A has option A1 using {x1,x2} and option A2 using {x1}
                                                          process B has option B1 using {x3,x4} and option B2 using {}
                                                          process C has option C1 using {} and option C2 using {x5,x6} and C3 using {x5}
                                                          --> paras_per_option = [[[0,1],[0]],   # process A
                                                                                  [[2,3],[]],    # process B
                                                                                  [[],[4,5],[4]] # process C
        para_ranges         list of lists        lists lower and upper bound for each parameter
                                                 example: para_ranges = [[ 0.0,3.0], #x1
                                                                         [ 1.0,3.0], #x2
                                                                         [-1.0,0.0], #x3
                                                                         [ 0.0,3.0], #x4
                                                                         [ 0.0,3.0], #x5
                                                                         [ 5.0,9.0]] #x6
        model_function      function             - function dictionary that can contain multiple weighted model output of all model options (scalar or 1D output), e.g.,
                                                       { 'Q':        np.array([2.5,3.5,...,4.5]),
                                                         'NSE':      0.6,
                                                         'baseflow': np.array([20.5,23.5,...,24.5])
                                                       }
                                                 - the internal sampling will make sure that all weights are between 0 and 1 
                                                   and sum up to 1 for each process
                                                 - interface must look like: 
                                                   model_function(set_of_parameters, set_of_weights, constants=constants, run_id=None) 
                                                 - weights are given as nested list of list (similar to 'paras_per_option')
                                                 example: 
                                                       def model_function(pp,ww,constant=None,run_id=None):
                                                            # process A
                                                            proc_a = ww[0][0] * (pp[0]**2+pp[1]) + ww[0][1] * (sin(pp[0]))                  
                                                            # process B
                                                            proc_b = ww[1][0] * (pp[2]**4+pp[3]**2) + ww[1][1] * (7.0)                      
                                                            # process C
                                                            proc_c = ww[2][0] * (9.81) + ww[2][1] * (pp[4]+cos(pp[5])) + ww[2][2] * pp[4]   
                                                            # model output: this is totally fake...
                                                            model['out'] = proc_a * proc_b + proc_c
                                                            return model
        basin_prop          dictionary           basin properties
                                                 example:
                                                       {'area_km2':     2303.95,
                                                        'elevation_m':  250.31,
                                                        'forest_frac':  0.9063,
                                                        'id':           '01013500',
                                                        'lat_deg':      47.23739,
                                                        'lon_deg':      -68.58264,
                                                        'name':         'Fish River near Fort Kent, Maine',
                                                        'slope_m_km-1': 21.64152}
        constants           list                 optional: list of constants that 'model_function' might need
                                                 default: None
        nsets               integer              optional: number of reference parameter sets
                                                 default: 1000
        budget              integer              optional: total number of model runs allowed for analysis; overwrites nsets
                                                 budget = nsets / ((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2))
                                                 default: None
        save_pkl            string               filename to save parameter sets and respective model outputs in pickle file (slow but small but requires lots of RAM)
                                                 default: None (nothing saved to file)
        save_json           string               filename to save parameter sets and respective model outputs in JSON file (fast but large files but requires not a lot of RAM)
                                                 default: None (nothing saved to file)
        save_msgpack        string               filename to save parameter sets and respective model outputs in MessagePack file  (fast and small and does not require lots of RAM)
                                                 default: None (nothing saved to file)
        save_nc4            string               filename to save parameter sets and respective model outputs in NetCDF4 file  (fast and small and does not require lots of RAM)
                                                 default: None (nothing saved to file)
        

        Output          Format        Description
        -----           -----         -----------
        sobol_indexes   dict          if model output is scalar:
                                          sobol_indexes['paras'][0]           ... main  Sobol' index of all parameters and weights (dims: nparas+nrand)
                                          sobol_indexes['paras'][1]           ... total Sobol' index of all parameters and weights (dims: nparas+nrand)
                                          sobol_indexes['process_options'][0] ... main  Sobol' index of all process options and weights (dims: nprocopts+nrand)
                                          sobol_indexes['process_options'][1] ... total Sobol' index of all process options and weights (dims: nprocopts+nrand)
                                          sobol_indexes['processes'][0]       ... main  Sobol' index of all processes (weights included) (dims: nprocesses)
                                          sobol_indexes['processes'][1]       ... total Sobol' index of all processes (weights included) (dims: nprocesses)
                                      if model output is 1D:
                                          sobol_indexes['paras'][0]           ...               main  Sobol' index of all parameters and weights (dimes: ntime, nparas+nrand)
                                          sobol_indexes['paras'][1]           ...               total Sobol' index of all parameters and weights (dimes: ntime, nparas+nrand)
                                          sobol_indexes['paras'][2]           ...          mean main  Sobol' index of all parameters and weights (dimes: ntime, nparas+nrand)
                                          sobol_indexes['paras'][3]           ...          mean total Sobol' index of all parameters and weights (dimes: ntime, nparas+nrand)
                                          sobol_indexes['paras'][4]           ... weighted mean main  Sobol' index of all parameters and weights (dimes: ntime, nparas+nrand)
                                          sobol_indexes['paras'][5]           ... weighted mean total Sobol' index of all parameters and weights (dimes: ntime, nparas+nrand)
                                          sobol_indexes['process_options'][0] ...               main  Sobol' index of all process options and weights (dimes: ntime, nprocopts+nrand)
                                          sobol_indexes['process_options'][1] ...               total Sobol' index of all process options and weights (dimes: ntime, nprocopts+nrand)
                                          sobol_indexes['process_options'][2] ...          mean main  Sobol' index of all process options and weights (dimes: ntime, nprocopts+nrand)
                                          sobol_indexes['process_options'][3] ...          mean total Sobol' index of all process options and weights (dimes: ntime, nprocopts+nrand)
                                          sobol_indexes['process_options'][4] ... weighted mean main  Sobol' index of all process options and weights (dimes: ntime, nprocopts+nrand)
                                          sobol_indexes['process_options'][5] ... weighted mean total Sobol' index of all process options and weights (dimes: ntime, nprocopts+nrand)
                                          sobol_indexes['processes'][0]       ...               main  Sobol' index of all processes (weights included) (dimes: ntime, nprocesses)
                                          sobol_indexes['processes'][1]       ...               total Sobol' index of all processes (weights included) (dimes: ntime, nprocesses)
                                          sobol_indexes['processes'][2]       ...          mean main  Sobol' index of all processes (weights included) (dimes: ntime, nprocesses)
                                          sobol_indexes['processes'][3]       ...          mean total Sobol' index of all processes (weights included) (dimes: ntime, nprocesses)
                                          sobol_indexes['processes'][4]       ... weighted mean main  Sobol' index of all processes (weights included) (dimes: ntime, nprocesses)
                                          sobol_indexes['processes'][5]       ... weighted mean total Sobol' index of all processes (weights included) (dimes: ntime, nprocesses)

        Description
        -----------
        


        Restrictions
        ------------
        Parameters can only be uniformly distributed in a range [a,b]. 
        No Gaussian distribution etc possible yet.

        Examples
        --------

        >>> import numpy as np
        >>> nsets = 1000

        --------------------------------------------------
        Simple setup
        --------------------------------------------------

        >>> # list of parameters that go into each option (numbering starts with 0)
        >>> # (a) simple setup
        >>> paras_per_option = [ 
        ...       [[0], []],             # parameters of process options A1 and A2
        ...       [[1], [2], [3,4]],     # parameters of process options B1, B2, and B3
        ...       [[5], [6]]             # parameters of process options A1 and A2
        ...     ]
        >>> para_ranges = [ 
        ...       [-np.pi,np.pi],      # parameter range of x1
        ...       [-np.pi,np.pi],      # parameter range of x2
        ...       [-np.pi,np.pi],      # parameter range of x3
        ...       [-np.pi,np.pi],      # parameter range of x4
        ...       [-np.pi,np.pi],      # parameter range of x5
        ...       [-np.pi,np.pi],      # parameter range of x6
        ...       [-np.pi,np.pi]       # parameter range of x7
        ...     ]  
        >>> basin_prop = {'area_km2':     2303.95,
        ...               'elevation_m':  250.31,
        ...               'forest_frac':  0.9063,
        ...               'id':           '01013500',
        ...               'lat_deg':      47.23739,
        ...               'lon_deg':      -68.58264,
        ...               'name':         'Fish River near Fort Kent, Maine',
        ...               'slope_m_km-1': 21.64152} 
        >>> def model_function(paras, weights, basin_prop, constants=None, run_id=None):
        ...     # input:
        ...     #     paras     ... list of model parameters scaled to their range;
        ...     #                   values for all N model parameters have to be provided
        ...     #                   example:
        ...     #                        [ x1, x2, x3, x4, .... ]
        ...     #     weights   ... list of lists of weights to weight options of each process;
        ...     #                   each list of the lists need to sum up to 1.0;
        ...     #                   each sublist is the N_i weights for the N_i process options of process i;
        ...     #                   example:
        ...     #                        [ [w_a1, w_a2, ...], [w_b1, w_b2, w_b3, ...], [w_c1, w_c2, ...], ... ]
        ...     #     constants ... optional list of constants that are same for all models;
        ...     #                   like parameters a and b in Ishigami-Homma function
        ...     #                   example:
        ...     #                        [2.0, 1.0]
        ...     # output:
        ...     #     model output
        ...     #     example:
        ...     #           { 'Q':        np.array([2.5,3.5,...,4.5]),
        ...     #             'NSE':      0.6,
        ...     #             'baseflow': np.array([20.5,23.5,...,24.5])
        ...     #           }
        ...
        ...     # check that provided number of weights is correct:
        ...     # --> one weight per option per process
        ...     if ( [len(ilist) for ilist in weights] != [2,3,2] ):
        ...         print("Number of weights: ",[len(ilist) for ilist in weights])
        ...         raise ValueError("sa_model_multiple_processes: model_function: provided number of weights must be [2,3,2]")
        ...     # check if sum up to 1.0:
        ...     if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
        ...         print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
        ...         raise ValueError("sa_model_multiple_processes: model_function: sum of weights must be 1.0 for all processes")
        ...     # check if weights <= 1.0:
        ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
        ...         print("Weights: ",weights)
        ...         raise ValueError("sa_model_multiple_processes: model_function: weights must be all less or equal 1.0")
        ...     # check if weights >= 0.0:
        ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
        ...         print("Weights: ",weights)
        ...         raise ValueError("sa_model_multiple_processes: model_function: weights must be all greater or equal 0.0")
        ...     # check if number of parameters is correct:
        ...     if (len(paras) != 7):
        ...         print("Number of parameters: ",len(paras))
        ...         raise ValueError("sa_model_multiple_processes: model_function: provided number of parameters must be 7")
        ...        
        ...     out = 0.0
        ...
        ...     if constants is None:
        ...         aa = 2.0
        ...         bb = 1.0
        ...     else:
        ...         aa = constants[0]
        ...         bb = constants[1]
        ...
        ...     # ---------------
        ...     # simple model
        ...     # ---------------
        ...        
        ...     # process A
        ...     out += ( weights[0][0] * np.sin(paras[0]) +              # A1
        ...              weights[0][1] * 1.0 )                           # A2
        ...     # process B
        ...     out *= ( weights[1][0] * (1.0 + bb * paras[1]**4) +      # B1
        ...              weights[1][1] * (1.0 + bb * paras[2]**2) +      # B2
        ...              weights[1][2] * (paras[3] + bb * paras[4]) )    # B3
        ...     # process C
        ...     out += ( weights[2][0] * (aa * np.sin(paras[5])**2) +    # C1
        ...              weights[2][1] * (1.0 + bb * paras[6]**4) )      # C2
        ...
        ...     model = {}
        ...     model['result_0d'] = out
        ...     model['result_1d'] = np.array( [ out for itime in range(4) ] )
        ...
        ...     return model

        >>> # this is calling the actual tool
        >>> sobol_indexes = sa_model_multiple_processes(paras_per_option, para_ranges, model_function, basin_prop, nsets=nsets)

        >>> # printing
        >>> print("parameter sensitivities:      S_xi  = ",astr(sobol_indexes['paras']['si']['result_0d'],prec=5))
        parameter sensitivities:      S_xi  =  [' 0.02295' ' 0.02361' '-0.00099' ' 0.00048' ' 0.00056' ' 0.00061' ' 0.49380' '-0.00491' ' 0.02107' ' 0.00436' ' 0.09056']
        >>> print("parameter sensitivities:      ST_xi = ",astr(sobol_indexes['paras']['sti']['result_0d'],prec=5))
        parameter sensitivities:      ST_xi =  ['0.06706' '0.16565' '0.00204' '0.00075' '0.00074' '0.00047' '0.64679' '0.08102' '0.11487' '0.00187' '0.24561']
        >>> print("process option sensitivities: S_Ai  = ",astr(sobol_indexes['process_options']['si']['result_0d'],prec=5))
        process option sensitivities: S_Ai  =  [' 0.02295' ' 0.00000' ' 0.02361' '-0.00099' ' 0.00105' ' 0.00061' ' 0.49380' '-0.00491' ' 0.02107' ' 0.00436' ' 0.09056']
        >>> print("process option sensitivities: ST_Ai = ",astr(sobol_indexes['process_options']['sti']['result_0d'],prec=5))
        process option sensitivities: ST_Ai =  ['0.06706' '0.00000' '0.16565' '0.00204' '0.00153' '0.00047' '0.64679' '0.08102' '0.11487' '0.00187' '0.24561']
        >>> print("process sensitivities:        S_A   = ",astr(sobol_indexes['processes']['si']['result_0d'],prec=5))
        process sensitivities:        S_A   =  ['0.03576' '0.05918' '0.76269']
        >>> print("process sensitivities:        ST_A  = ",astr(sobol_indexes['processes']['sti']['result_0d'],prec=5))
        process sensitivities:        ST_A  =  ['0.13269' '0.20039' '0.74577']
        >>> print("process sensitivities:        ST_A_msti  = ",astr(sobol_indexes['processes']['msti']['result_1d'],prec=5))
        process sensitivities:        ST_A_msti  =  ['0.13269' '0.20039' '0.74577']
        >>> print("process sensitivities:        ST_A_wsti  = ",astr(sobol_indexes['processes']['wsti']['result_1d'],prec=5))
        process sensitivities:        ST_A_wsti  =  ['0.13269' '0.20039' '0.74577']

        # --------------------------------------------------
        # Realistic setup
        # --------------------------------------------------

        # >>> # list of parameters that go into each option (numbering starts with 0)
        # >>> # (a) simple setup
        # >>> paras_per_option = [
        # ...       [[0], [0,1]],             # parameters of process options A1 and A2
        # ...       [[1], [2], [3,4]],        # parameters of process options B1, B2, and B3
        # ...       [[5], [2,6]]              # parameters of process options A1 and A2
        # ...     ]
        # >>> para_ranges = [
        # ...       [-np.pi,np.pi],      # parameter range of x1
        # ...       [-np.pi,np.pi],      # parameter range of x2
        # ...       [-np.pi,np.pi],      # parameter range of x3
        # ...       [-np.pi,np.pi],      # parameter range of x4
        # ...       [-np.pi,np.pi],      # parameter range of x5
        # ...       [-np.pi,np.pi],      # parameter range of x6
        # ...       [-np.pi,np.pi]       # parameter range of x7
        # ...     ]
        # >>> basin_prop = {'area_km2':     2303.95,
        # ...               'elevation_m':  250.31,
        # ...               'forest_frac':  0.9063,
        # ...               'id':           '01013500',
        # ...               'lat_deg':      47.23739,
        # ...               'lon_deg':      -68.58264,
        # ...               'name':         'Fish River near Fort Kent, Maine',
        # ...               'slope_m_km-1': 21.64152}
        # >>> def model_function(paras, weights, constants=None, run_id=None):
        # ...     # input:
        # ...     #     paras     ... list of model parameters scaled to their range;
        # ...     #                   values for all N model parameters have to be provided
        # ...     #                   example:
        # ...     #                        [ x1, x2, x3, x4, .... ]
        # ...     #     weights   ... list of lists of weights to weight options of each process;
        # ...     #                   each list of the lists need to sum up to 1.0;
        # ...     #                   each sublist is the N_i weights for the N_i process options of process i;
        # ...     #                   example:
        # ...     #                        [ [w_a1, w_a2, ...], [w_b1, w_b2, w_b3, ...], [w_c1, w_c2, ...], ... ]
        # ...     #     constants ... optional list of constants that are same for all models;
        # ...     #                   like parameters a and b in Ishigami-Homma function
        # ...     #                   example:
        # ...     #                        [2.0, 1.0]
        # ...     # output:
        # ...     #     model output
        # ...     #     example:
        # ...     #           { 'Q':        np.array([2.5,3.5,...,4.5]),
        # ...     #             'NSE':      0.6,
        # ...     #             'baseflow': np.array([20.5,23.5,...,24.5])
        # ...     #           }
        # ...
        # ...     # check that provided number of weights is correct:
        # ...     # --> one weight per option per process
        # ...     if ( [len(ilist) for ilist in weights] != [2,3,2] ):
        # ...         print("Number of weights: ",[len(ilist) for ilist in weights])
        # ...         raise ValueError("sa_model_multiple_processes: model_function: provided number of weights must be [2,3,2]")
        # ...     # check if sum up to 1.0:
        # ...     if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
        # ...         print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
        # ...         raise ValueError("sa_model_multiple_processes: model_function: sum of weights must be 1.0 for all processes")
        # ...     # check if weights <= 1.0:
        # ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
        # ...         print("Weights: ",weights)
        # ...         raise ValueError("sa_model_multiple_processes: model_function: weights must be all less or equal 1.0")
        # ...     # check if weights >= 0.0:
        # ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
        # ...         print("Weights: ",weights)
        # ...         raise ValueError("sa_model_multiple_processes: model_function: weights must be all greater or equal 0.0")
        # ...     # check if number of parameters is correct:
        # ...     if (len(paras) != 7):
        # ...         print("Number of parameters: ",len(paras))
        # ...         raise ValueError("sa_model_multiple_processes: model_function: provided number of parameters must be 7")
        # ...
        # ...     out = 0.0
        # ...
        # ...     if constants is None:
        # ...         aa = 2.0
        # ...         bb = 1.0
        # ...     else:
        # ...         aa = constants[0]
        # ...         bb = constants[1]
        # ...
        # ...     # ---------------
        # ...     # realistic model
        # ...     # ---------------
        # ...
        # ...     # process D
        # ...     out += ( weights[0][0] * np.sin(paras[0]) +                            # D1
        # ...              weights[0][1] * (paras[0]+paras[1]**2) )                      # D2
        # ...     # process E
        # ...     out *= ( weights[1][0] * (1.0 + bb * paras[1]**4) +                    # E1
        # ...              weights[1][1] * (1.0 + bb * paras[2]**2) +                    # E2
        # ...              weights[1][2] * (paras[3] + bb * paras[4]) )                  # E3
        # ...     # process F
        # ...     out += ( weights[2][0] * (aa * np.sin(paras[5])**2) +                  # F1
        # ...              weights[2][1] * (1.0 + bb * paras[6]**4) + paras[2]**2 )      # F2
        # ...
        # ...     model = {}
        # ...     model['result'] = out
        # ...
        # ...     return model

        # >>> # this is calling the actual tool
        # >>> sobol_indexes = sa_model_multiple_processes(paras_per_option, para_ranges, model_function, basin_prop, nsets=nsets)

        # >>> # printing
        # >>> print("parameter sensitivities:      S_xi  = ",astr(sobol_indexes['paras']['si']['result'],prec=5))
        # parameter sensitivities:      S_xi  =  [' 0.02671' ' 0.40059' ' 0.00444' ' 0.00133' ' 0.00112' '-0.00009' ' 0.04133' ' 0.03133' ' 0.07610' ' 0.00051' ' 0.00949']
        # >>> print("parameter sensitivities:      ST_xi = ",astr(sobol_indexes['paras']['sti']['result'],prec=5))
        # parameter sensitivities:      ST_xi =  ['0.05547' '0.83153' '0.00687' '0.00093' '0.00104' '0.00004' '0.05498' '0.26319' '0.44192' '0.00257' '0.02088']
        # >>> print("process option sensitivities: S_Ai  = ",astr(sobol_indexes['process_options']['si']['result'],prec=5))
        # process option sensitivities: S_Ai  =  [' 0.02671' ' 0.42217' ' 0.40059' ' 0.00444' ' 0.00245' '-0.00009' ' 0.04577' ' 0.03133' ' 0.07610' ' 0.00051' ' 0.00949']
        # >>> print("process option sensitivities: ST_Ai = ",astr(sobol_indexes['process_options']['sti']['result'],prec=5))
        # process option sensitivities: ST_Ai =  ['0.05547' '0.85117' '0.83153' '0.00687' '0.00202' '0.00004' '0.06311' '0.26319' '0.44192' '0.00257' '0.02088']
        # >>> print("process sensitivities:        S_A   = ",astr(sobol_indexes['processes']['si']['result'],prec=5))
        # process sensitivities:        S_A   =  ['0.61152' '0.57158' '0.07124']
        # >>> print("process sensitivities:        ST_A  = ",astr(sobol_indexes['processes']['sti']['result'],prec=5))
        # process sensitivities:        ST_A  =  ['0.93249' '0.86664' '0.07131']


        
        License
        -------
        This file is part of the "SA for Models with Multiple Processes" Python package.

        The "SA for Models with Multiple Processes" Python package is free software: you 
        can redistribute it and/or modify it under the terms of the GNU Lesser General 
        Public License as published by the Free Software Foundation, either version 3 of 
        the License, or (at your option) any later version.

        The "SA for Models with Multiple Processes" Python package is distributed in the 
        hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
        warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the PieShareDistribution project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2019 Juliane Mai - juliane.mai@uwaterloo.ca


        History
        -------
        Written,  Juliane Mai, June 2019
    """

    # initialize return variable
    sobol_indexes = OrderedDict()
    ntime         = OrderedDict()

    para_ranges = np.array(para_ranges)

    # number of parameters
    nparas   = np.shape(para_ranges)[0]
    # number of options per process
    noptions = np.array([ len(oo) for oo in paras_per_option ])
    # number of processes
    nprocess = np.shape(noptions)[0]
    # number of weights required
    nweights = np.sum(noptions-1)

    if ( nsets is None ):
        nsets = 1000
        
    if not( budget is None ):
        # overwrite nsets if budget is given
        nsets = budget / ((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2))
        if nsets <= 0:
            print("Minimal budget:     ",((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2)))
            print("Recommended budget: ",((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2))*1000)
            raise ValueError("sa_model_multiple_processes: Budget is too small!")

    if not( (save_pkl is None) and (save_json is None) and (save_msgpack is None) ):
        save_to_file = OrderedDict()

    if not(save_nc4 is None):
        
        nc4_out = nc4.Dataset(save_nc4, "w", format="NETCDF4")
        # create dimensions
        dim_nsets    = nc4_out.createDimension("nsets",    None)                 # unlimited
        dim_nparas   = nc4_out.createDimension("nparas",   nparas+nweights)
        dim_noptions = nc4_out.createDimension("noptions", np.sum(noptions)+nweights)
        dim_nprocess = nc4_out.createDimension("nprocesses", nprocess)
        nc4_out.close()
        

    # (A) Sampling parameters and weights in unit interval using Sobol' sequences
    sobol_sets = sobol.i4_sobol_generate((nparas+nweights)*2,nsets,40000)
    sobol_sets = np.transpose(sobol_sets)

    # (B) Scaling of parameters using given 'para_ranges'
    block_a_paras  = copy.deepcopy(sobol_sets[:,0:nparas])
    block_a_paras *= (para_ranges[:,1]-para_ranges[:,0])
    block_a_paras += para_ranges[:,0]
    
    block_b_paras  = copy.deepcopy(sobol_sets[:,nparas+nweights:2*nparas+nweights])
    block_b_paras *= (para_ranges[:,1]-para_ranges[:,0])
    block_b_paras += para_ranges[:,0]

    # (C) Weights still in unit range (will be scaled later; otherwise weights might not sum up to one when Ci blocks are created)
    block_a_weights = copy.deepcopy(sobol_sets[:,nparas:nparas+nweights])
    block_b_weights = copy.deepcopy(sobol_sets[:,2*nparas+nweights:2*nparas+2*nweights])

    def random_to_weights_to_nested(rnd,noptions):
        # convert 1D random numbers first to weights (without remainder) and then to nested list
        # e.g:     random  =
        #      --> weights = [ 0.4,      0.5,0.1,      0.3]   with noptions [2,3,2]
        #      --> nested  = [[0.4,0.6],[0.5,0.1,0.4],[0.3,0.7]]
        nprocess = len(noptions)
        nsets    = np.shape(rnd)[0]
        
        start = 0
        csum  = np.cumsum(noptions-1)
        weights = np.ones([np.shape(rnd)[0],np.shape(rnd)[1]]) * -9999.0
        # convert random numbers of [0,1] to weights using PieShareDistribution function
        for iprocess in range(nprocess):
            weights[:,start:csum[iprocess]] = psd.PieShareDistribution(nsets,noptions[iprocess],remainder=False,randomnumbers=rnd[:,start:csum[iprocess]])
            weights[:,start:csum[iprocess]] = psd.PieShareDistribution(nsets,noptions[iprocess],remainder=False,randomnumbers=rnd[:,start:csum[iprocess]])
            start = csum[iprocess]

        csum_tmp = np.append(np.cumsum(noptions-1),0)    
        weights_to_nested = [ [ np.append(np.round(weights[iset,csum_tmp[io-1]:csum_tmp[io]],6),
                                          1.0-np.sum(np.round(weights[iset,csum_tmp[io-1]:csum_tmp[io]],6)))
                              for io,oo in enumerate(noptions) ]
                            for iset in range(nsets) ]
        return weights_to_nested

    def weights_to_nested(weights,noptions):
        # convert 1D weights (without remainder) to nested list
        # e.g: weights = [ 0.4,      0.5,0.1,      0.3]   with noptions [2,3,2]
        #      -->       [[0.4,0.6],[0.5,0.1,0.4],[0.3,0.7]]
        nsets    = np.shape(weights)[0]
        csum_tmp = np.append(np.cumsum(noptions-1),0)
        weights_to_nested = [ [ np.append(weights[iset,csum_tmp[io-1]:csum_tmp[io]],
                                          1.0-np.sum(weights[iset,csum_tmp[io-1]:csum_tmp[io]]))
                              for io,oo in enumerate(noptions) ]
                            for iset in range(nsets) ]
        return weights_to_nested

    # (D) f_A and f_B
    block_a_weights_nested = random_to_weights_to_nested(block_a_weights,noptions)
    # f_a = np.array([ model_function(block_a_paras[iset], block_a_weights_nested[iset], basin_prop, constants=constants, run_id=basin_prop['id']+"_a_set_"+str(iset)) for iset in range(nsets) ])
    f_a = np.array([ model_function(block_a_paras[iset], block_a_weights_nested[iset], basin_prop, constants=constants, run_id=basin_prop['id']+"_a_set") for iset in range(nsets) ])
    # print('block_a_paras   = ',block_a_paras)
    # print('block_a_weights = ',block_a_weights)

    # convert list of dicts into dict of lists:
    #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
    keys = f_a[0].keys()
    tmp = OrderedDict()
    for ikey in keys:
        tmp_key = []
        for iset in range(nsets):
            tmp_key.append( f_a[iset][ikey] )
            
        tmp[ikey] = np.array(tmp_key)
    f_a = tmp
    # print('keys = ',keys)
    # print('f_a = ',f_a)
    
    # if vector of model outputs f_a has shape (nsets,ntime) --> must be (ntime, nsets)
    for ikey in keys:
        if (len(np.shape(f_a[ikey])) == 2):
            f_a[ikey] = np.transpose(f_a[ikey])

    for ikey in keys:
        if (len(np.shape(f_a[ikey])) == 2):
            ntime[ikey] = np.shape(f_a[ikey])[0]
        else:
            ntime[ikey] = 0
    # print('ntime = ',ntime)
    # store to save in pickle later
    if not( (save_pkl is None) and (save_json is None) and (save_msgpack is None) ):
        save_to_file["ntime"] = ntime

    # ----------------------------
    # NetCDF
    #    create groups for each model output (only for sets A because groups dont exist yet)
    # ----------------------------
    if not(save_nc4 is None):
        
        nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
        for ikey in keys:
            grp = nc4_out.createGroup(ikey)
            if ntime[ikey] > 0:
                grp_time_dim = grp.createDimension("ntime", ntime[ikey])
        nc4_out.close()

    # ----------------------------
    # NetCDF
    #    create variables now that time dimension size is known
    # ----------------------------
    if not(save_nc4 is None):
        
        nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
        # create groups (only for sets A because groups dont exist yet)
        grp_a = {} ; grp_b = {}; grp_c_paras = {}; grp_c_options = {}; grp_c_processes = {}
        for ikey in keys:
            if ntime[ikey] > 0:
                grp_var = nc4_out.createVariable(ikey+'/f_a'            , "f4",("ntime","nsets",),               zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_b'            , "f4",("ntime","nsets",),               zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_c_paras'      , "f4",("ntime","nparas",     "nsets",), zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_c_options'    , "f4",("ntime","noptions",   "nsets",), zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_c_processes'  , "f4",("ntime","nprocesses", "nsets",), zlib=True)

                grp_var = nc4_out.createVariable(ikey+'/si_paras'       , "f4",("ntime","nparas",),              zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/sti_paras'      , "f4",("ntime","nparas",),              zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/msi_paras'      , "f4",("nparas",),                      zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/msti_paras'     , "f4",("nparas",),                      zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/wsi_paras'      , "f4",("nparas",),                      zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/wsti_paras'     , "f4",("nparas",),                      zlib=True)

                grp_var = nc4_out.createVariable(ikey+'/si_options'     , "f4",("ntime","noptions",),            zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/sti_options'    , "f4",("ntime","noptions",),            zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/msi_options'    , "f4",("noptions",),                    zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/msti_options'   , "f4",("noptions",),                    zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/wsi_options'    , "f4",("noptions",),                    zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/wsti_options'   , "f4",("noptions",),                    zlib=True)

                grp_var = nc4_out.createVariable(ikey+'/si_processes'   , "f4",("ntime","nprocesses",),          zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/sti_processes'  , "f4",("ntime","nprocesses",),          zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/msi_processes'  , "f4",("nprocesses",),                  zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/msti_processes' , "f4",("nprocesses",),                  zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/wsi_processes'  , "f4",("nprocesses",),                  zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/wsti_processes' , "f4",("nprocesses",),                  zlib=True)
            else:
                grp_var = nc4_out.createVariable(ikey+'/f_a'            , "f4",("nsets",),                       zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_b'            , "f4",("nsets",),                       zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_c_paras'      , "f4",("nparas",     "nsets",),         zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_c_options'    , "f4",("noptions",   "nsets",),         zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/f_c_processes'  , "f4",("nprocesses", "nsets",),         zlib=True)
                
                grp_var = nc4_out.createVariable(ikey+'/si_paras'       , "f4",("nparas",),                      zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/sti_paras'      , "f4",("nparas",),                      zlib=True)

                grp_var = nc4_out.createVariable(ikey+'/si_options'     , "f4",("noptions",),                    zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/sti_options'    , "f4",("noptions",),                    zlib=True)

                grp_var = nc4_out.createVariable(ikey+'/si_processes'   , "f4",("nprocesses",),                  zlib=True)
                grp_var = nc4_out.createVariable(ikey+'/sti_processes'  , "f4",("nprocesses",),                  zlib=True)
        nc4_out.close()

    # ----------------------------
    # NetCDF
    #    save model output A in NetCDF
    # ----------------------------
    if not(save_nc4 is None):
        nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
        for ikey in keys:
            nc4_out.groups[ikey].variables['f_a'][:] = f_a[ikey]
        nc4_out.close()

    block_b_weights_nested = random_to_weights_to_nested(block_b_weights,noptions)
    # f_b = np.array([ model_function(block_b_paras[iset], block_b_weights_nested[iset], basin_prop, constants=constants, run_id=basin_prop['id']+"_b_set_"+str(iset)) for iset in range(nsets) ])
    f_b = np.array([ model_function(block_b_paras[iset], block_b_weights_nested[iset], basin_prop, constants=constants, run_id=basin_prop['id']+"_b_set") for iset in range(nsets) ])
    # print('block_b_paras   = ',block_b_paras)
    # print('block_b_weights = ',block_b_weights)

    # convert list of dicts into dict of lists:
    #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
    # keys = f_b[0].keys()
    tmp = OrderedDict()
    for ikey in keys:
        tmp_key = []
        for iset in range(nsets):
            tmp_key.append( f_b[iset][ikey] )
            
        tmp[ikey] = np.array(tmp_key)
    f_b = tmp
    # print('f_b = ',f_b)
    
    # if vector of model outputs f_b has shape (nsets,ntime) --> must be (ntime, nsets)
    for ikey in keys:
        if (len(np.shape(f_b[ikey])) == 2):
            f_b[ikey] = np.transpose(f_b[ikey])

    # ----------------------------
    # NetCDF
    #    save model output A in NetCDF
    # ----------------------------
    if not(save_nc4 is None):
        nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
        for ikey in keys:
            nc4_out.groups[ikey].variables['f_b'][:] = f_b[ikey]
        nc4_out.close()

    # store to save in pickle later
    if not( (save_pkl is None) and (save_json is None) and (save_msgpack is None) ):
        save_to_file["block_a_paras"]          = copy.deepcopy(block_a_paras)
        save_to_file["block_b_paras"]          = copy.deepcopy(block_b_paras)
        save_to_file["block_a_weights_nested"] = copy.deepcopy(block_a_weights_nested)
        save_to_file["block_b_weights_nested"] = copy.deepcopy(block_b_weights_nested)
        save_to_file["f_a"]                    = copy.deepcopy(f_a)
        save_to_file["f_b"]                    = copy.deepcopy(f_b)
        
    # (1) parameter sensitivities:
    #     main effects:  [S_x1,  S_x2,  ..., S_w1,  S_w2,  ...]
    #     total effects: [ST_x1, ST_x2, ..., ST_w1, ST_w2, ...]
    col_changes_Ci = [ [ii] for ii in range(nparas+nweights) ]  # list of columns to change at ones for Ci
    f_c = OrderedDict()
    for ikey in keys:
        if (ntime[ikey] == 0):
            f_c[ikey] = np.ones([len(col_changes_Ci),nsets]) * -9999.
        else:
            # print('nsets = ',nsets)
            tttmp = np.ones([ntime[ikey],len(col_changes_Ci),nsets]) * -9999.
            f_c[ikey] = tttmp

    block_c_paras          = [ [] for icol in col_changes_Ci ]
    block_c_weights        = [ [] for icol in col_changes_Ci ]
    block_c_weights_nested = [ [] for icol in col_changes_Ci ]
    for iicol, icol in enumerate(col_changes_Ci):
        # (1a) create Ci/j blocks for all parameters xi and all weights wj
        block_c_paras[iicol]    = copy.deepcopy(block_a_paras)
        block_c_weights[iicol]  = copy.deepcopy(block_a_weights)
        for ipara in icol:
            if ipara < nparas:
                block_c_paras[iicol][:,ipara] = copy.deepcopy(block_b_paras[:,ipara])
            else:
                block_c_weights[iicol][:,ipara-nparas] = copy.deepcopy(block_b_weights[:,ipara-nparas])
        # (1b) run model for all Ci/j blocks
        block_c_weights_nested[iicol] = random_to_weights_to_nested(block_c_weights[iicol],noptions)
        # f_c_tmp = np.array([ model_function(block_c_paras[iicol][iset],
        #                                         block_c_weights_nested[iicol][iset],
        #                                         basin_prop,
        #                                         constants=constants,
        #                                         run_id=basin_prop['id']+"_c_para_"+str(iicol)+"_set_"+str(iset)) for iset in range(nsets) ])
        f_c_tmp = np.array([ model_function(block_c_paras[iicol][iset],
                                                block_c_weights_nested[iicol][iset],
                                                basin_prop,
                                                constants=constants,
                                                run_id=basin_prop['id']+"_c_set") for iset in range(nsets) ])

        # convert list of dicts into dict of lists:
        #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
        # keys = f_c_tmp[0].keys()
        tmp = OrderedDict()
        for ikey in keys:
            tmp_key = []
            for iset in range(nsets):
                tmp_key.append( f_c_tmp[iset][ikey] )
            
            tmp[ikey] = np.array(tmp_key)
        f_c_tmp = tmp
        # print('f_c_tmp = ',f_c_tmp)

        for ikey in keys:
            if (ntime[ikey] == 0):
                f_c[ikey][iicol,:] = f_c_tmp[ikey]
            else:
                # if vector of model outputs f_c_tmp has shape (nsets,ntime) --> must be (ntime, nsets)
                f_c[ikey][:,iicol,:] = np.transpose(f_c_tmp[ikey])

        # ----------------------------
        # NetCDF
        #    save model output C_paras in NetCDF
        # ----------------------------
        if not(save_nc4 is None):
            nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
            for ikey in keys:
                if (ntime[ikey] == 0):
                    nc4_out.groups[ikey].variables['f_c_paras'][iicol:iicol+1,:] = f_c[ikey][iicol:iicol+1,:]
                else:
                    nc4_out.groups[ikey].variables['f_c_paras'][:,iicol:iicol+1,:] = f_c[ikey][:,iicol:iicol+1,:]
            nc4_out.close()

    # store to save in pickle later
    if not( (save_pkl is None) and (save_json is None) and (save_msgpack is None) ):
        save_to_file["block_c_paras"]                = copy.deepcopy(block_c_paras)
        save_to_file["block_c_weights_nested_paras"] = copy.deepcopy(block_c_weights_nested)
        save_to_file["f_c_paras"]                    = copy.deepcopy(f_c)
        
    # (1c) calculate Sobol' indexes
    si   = OrderedDict()
    sti  = OrderedDict()
    msi  = OrderedDict()
    msti = OrderedDict()
    wsi  = OrderedDict()
    wsti = OrderedDict()
    for ikey in keys:
        if (ntime[ikey] == 0):
            si[ikey], sti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                          yb=f_b[ikey],
                                                          yc=f_c[ikey],
                                                          si=True,
                                                          sti=True,
                                                          method='Mai1999')
            msi[ikey]  = None
            msti[ikey] = None
            wsi[ikey]  = None
            wsti[ikey] = None

        else:
            si[ikey], sti[ikey], msi[ikey], msti[ikey], wsi[ikey], wsti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                                                                        yb=f_b[ikey],
                                                                                                        yc=f_c[ikey],
                                                                                                        si=True,
                                                                                                        sti=True,
                                                                                                        mean=True,
                                                                                                        wmean=True,
                                                                                                        method='Mai1999')

    tmp = OrderedDict()
    tmp['si']   = si
    tmp['sti']  = sti
    tmp['msi']  = msi
    tmp['msti'] = msti
    tmp['wsi']  = wsi
    tmp['wsti'] = wsti
    sobol_indexes['paras'] = tmp

    # ----------------------------
    # NetCDF
    #    save sensitivity indexes 
    # ----------------------------    
    if not(save_nc4 is None):
        nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
        for ikey in keys:
            nc4_out.groups[ikey].variables['si_paras'][:]  = si[ikey]
            nc4_out.groups[ikey].variables['sti_paras'][:] = sti[ikey]
            if (ntime[ikey] > 0):
                nc4_out.groups[ikey].variables['msi_paras'][:]  = msi[ikey]
                nc4_out.groups[ikey].variables['msti_paras'][:] = msti[ikey]
                nc4_out.groups[ikey].variables['wsi_paras'][:]  = wsi[ikey]
                nc4_out.groups[ikey].variables['wsti_paras'][:] = wsti[ikey]
        nc4_out.close()

    # print("si['"+ikey+"']   = ", si[ikey])
    # print("sti['"+ikey+"']  = ", sti[ikey])
    # print("msi['"+ikey+"']  = ", msi[ikey])
    # print("msti['"+ikey+"'] = ", msti[ikey])
    # print("wsi['"+ikey+"']  = ", wsi[ikey])
    # print("wsti['"+ikey+"'] = ", wsti[ikey])

    #print("shape si['"+ikey+"']   = ",np.shape(si[ikey]))
    #print("shape sti['"+ikey+"']  = ",np.shape(sti[ikey]))
    #print("shape msi['"+ikey+"']  = ",np.shape(msi[ikey]))
    #print("shape msti['"+ikey+"'] = ",np.shape(msti[ikey]))
    #print("shape wsi['"+ikey+"']  = ",np.shape(wsi[ikey]))
    #print("shape wsti['"+ikey+"'] = ",np.shape(wsti[ikey]))


    

    # (2) process option sensitivities:
    #     main effects:  [S_A1,  S_A2,  ..., S_w1,  S_w2,  ...]
    #     total effects: [ST_A1, ST_A2, ..., ST_w1, ST_w2, ...]
    col_changes_Ci = [ ii for ioption in paras_per_option for ii in ioption ] + [ [ii] for ii in range(nparas,nparas+nweights) ]  # list of columns to change at ones for Ci
    f_c = OrderedDict()
    for ikey in keys:
        if (ntime[ikey] == 0):
            f_c[ikey] = np.ones([len(col_changes_Ci),nsets]) * -9999.
        else:
            f_c[ikey] = np.ones([ntime[ikey],len(col_changes_Ci),nsets]) * -9999.
        
    block_c_paras          = [ [] for icol in col_changes_Ci ]
    block_c_weights        = [ [] for icol in col_changes_Ci ]
    block_c_weights_nested = [ [] for icol in col_changes_Ci ]
    for iicol, icol in enumerate(col_changes_Ci):
        # (2a) create Ci/j blocks for all options Ai and weights wj
        block_c_paras[iicol]    = copy.deepcopy(block_a_paras)
        block_c_weights[iicol]  = copy.deepcopy(block_a_weights)
        for ipara in icol:
            if ipara < nparas:
                block_c_paras[iicol][:,ipara] = copy.deepcopy(block_b_paras[:,ipara])
            else:
                block_c_weights[iicol][:,ipara-nparas] = copy.deepcopy(block_b_weights[:,ipara-nparas])
        # (2b) run model for all Ci/j blocks
        block_c_weights_nested[iicol] = random_to_weights_to_nested(block_c_weights[iicol],noptions)
        # f_c_tmp = np.array([ model_function(block_c_paras[iicol][iset],
        #                                         block_c_weights_nested[iicol][iset],
        #                                         basin_prop,
        #                                         constants=constants,
        #                                         run_id=basin_prop['id']+"_c_para_"+str(iicol)+"_set_"+str(iset)) for iset in range(nsets) ])
        f_c_tmp = np.array([ model_function(block_c_paras[iicol][iset],
                                                block_c_weights_nested[iicol][iset],
                                                basin_prop,
                                                constants=constants,
                                                run_id=basin_prop['id']+"_c_set") for iset in range(nsets) ])

        # convert list of dicts into dict of lists:
        #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
        # keys = f_c_tmp[0].keys()
        tmp = OrderedDict()
        for ikey in keys:
            tmp_key = []
            for iset in range(nsets):
                tmp_key.append( f_c_tmp[iset][ikey] )
            
            tmp[ikey] = np.array(tmp_key)
        f_c_tmp = tmp
        # print('f_c_tmp = ',f_c_tmp)

        for ikey in keys:
            if (ntime[ikey] == 0):
                f_c[ikey][iicol,:] = f_c_tmp[ikey]
            else:
                # if vector of model outputs f_c_tmp has shape (nsets,ntime) --> must be (ntime, nsets)
                f_c[ikey][:,iicol,:] = np.transpose(f_c_tmp[ikey])

        # ----------------------------
        # NetCDF
        #    save model output C_options in NetCDF
        # ----------------------------
        if not(save_nc4 is None):
            nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
            for ikey in keys:
                if (ntime[ikey] == 0):
                    nc4_out.groups[ikey].variables['f_c_options'][iicol:iicol+1,:] = f_c[ikey][iicol:iicol+1,:]
                else:
                    nc4_out.groups[ikey].variables['f_c_options'][:,iicol:iicol+1,:] = f_c[ikey][:,iicol:iicol+1,:]
            nc4_out.close()

    # print("f_c = ",f_c)

    # store to save in pickle later
    if not( (save_pkl is None) and (save_json is None) and (save_msgpack is None) ):
        save_to_file["block_c_process_options"]                = copy.deepcopy(block_c_paras)
        save_to_file["block_c_weights_nested_process_options"] = copy.deepcopy(block_c_weights_nested)
        save_to_file["f_c_process_options"]                    = copy.deepcopy(f_c)
        
    # (2c) calculate Sobol' indexes
    si   = OrderedDict()
    sti  = OrderedDict()
    msi  = OrderedDict()
    msti = OrderedDict()
    wsi  = OrderedDict()
    wsti = OrderedDict()
    for ikey in keys:
        if (ntime[ikey] == 0):
            si[ikey], sti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                          yb=f_b[ikey],
                                                          yc=f_c[ikey],
                                                          si=True,
                                                          sti=True,
                                                          method='Mai1999')
            msi[ikey]  = None
            msti[ikey] = None
            wsi[ikey]  = None
            wsti[ikey] = None
        else:
            si[ikey], sti[ikey], msi[ikey], msti[ikey], wsi[ikey], wsti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                                                                        yb=f_b[ikey],
                                                                                                        yc=f_c[ikey],
                                                                                                        si=True,
                                                                                                        sti=True,
                                                                                                        mean=True,
                                                                                                        wmean=True,
                                                                                                        method='Mai1999')
    tmp = OrderedDict()
    tmp['si']   = si
    tmp['sti']  = sti
    tmp['msi']  = msi
    tmp['msti'] = msti
    tmp['wsi']  = wsi
    tmp['wsti'] = wsti
    sobol_indexes['process_options'] = tmp

    # ----------------------------
    # NetCDF
    #    save sensitivity indexes 
    # ----------------------------
    if not(save_nc4 is None):
        nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
        for ikey in keys:
            nc4_out.groups[ikey].variables['si_options'][:]  = si[ikey]
            nc4_out.groups[ikey].variables['sti_options'][:] = sti[ikey]
            if (ntime[ikey] > 0):
                nc4_out.groups[ikey].variables['msi_options'][:]  = msi[ikey]
                nc4_out.groups[ikey].variables['msti_options'][:] = msti[ikey]
                nc4_out.groups[ikey].variables['wsi_options'][:]  = wsi[ikey]
                nc4_out.groups[ikey].variables['wsti_options'][:] = wsti[ikey]
        nc4_out.close()

    # print("si['"+ikey+"']   = ", si[ikey])
    # print("sti['"+ikey+"']  = ", sti[ikey])
    # print("msi['"+ikey+"']  = ", msi[ikey])
    # print("msti['"+ikey+"'] = ", msti[ikey])
    # print("wsi['"+ikey+"']  = ", wsi[ikey])
    # print("wsti['"+ikey+"'] = ", wsti[ikey])

    #print("shape si['"+ikey+"']   = ",np.shape(si[ikey]))
    #print("shape sti['"+ikey+"']  = ",np.shape(sti[ikey]))
    #print("shape msi['"+ikey+"']  = ",np.shape(msi[ikey]))
    #print("shape msti['"+ikey+"'] = ",np.shape(msti[ikey]))
    #print("shape wsi['"+ikey+"']  = ",np.shape(wsi[ikey]))
    #print("shape wsti['"+ikey+"'] = ",np.shape(wsti[ikey]))

    

    # (3) process sensitivities:
    #     main effects:  [S_A,   S_B,   ...]
    #     total effects: [ST_A,  ST_B,  ...]
    csum_tmp = np.append(np.cumsum(noptions-1),0) 
    col_changes_Ci = [ [ ipara for ilist in ioption for ipara in ilist ] +  # list of parameters in this process
                       list(range(nparas,nparas+nweights)[csum_tmp[iioption-1]:csum_tmp[iioption]]) for iioption,ioption in enumerate(paras_per_option) # list of weights for this process
                     ]  # list of columns to change at ones for Ci
    f_c = OrderedDict()
    for ikey in keys:
        if (ntime[ikey] == 0):
            f_c[ikey] = np.ones([len(col_changes_Ci),nsets]) * -9999.
        else:
            f_c[ikey] = np.ones([ntime[ikey],len(col_changes_Ci),nsets]) * -9999.
        
    block_c_paras          = [ [] for icol in col_changes_Ci ]
    block_c_weights        = [ [] for icol in col_changes_Ci ]
    block_c_weights_nested = [ [] for icol in col_changes_Ci ]
    for iicol, icol in enumerate(col_changes_Ci):
        # (3a) create Ci blocks for all processes A, B, C, ...
        block_c_paras[iicol]    = copy.deepcopy(block_a_paras)
        block_c_weights[iicol]  = copy.deepcopy(block_a_weights)
        for ipara in icol:
            if ipara < nparas:
                block_c_paras[iicol][:,ipara] = copy.deepcopy(block_b_paras[:,ipara])
            else:
                block_c_weights[iicol][:,ipara-nparas] = copy.deepcopy(block_b_weights[:,ipara-nparas])
        # (3b) run model for all Ci blocks
        block_c_weights_nested[iicol] = random_to_weights_to_nested(block_c_weights[iicol],noptions)
        # f_c_tmp = np.array([ model_function(block_c_paras[iicol][iset],
        #                                         block_c_weights_nested[iicol][iset],
        #                                         basin_prop,
        #                                         constants=constants,
        #                                         run_id=basin_prop['id']+"_c_para_"+str(iicol)+"_set_"+str(iset)) for iset in range(nsets) ])
        f_c_tmp = np.array([ model_function(block_c_paras[iicol][iset],
                                                block_c_weights_nested[iicol][iset],
                                                basin_prop,
                                                constants=constants,
                                                run_id=basin_prop['id']+"_c_set") for iset in range(nsets) ])

        # convert list of dicts into dict of lists:
        #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
        # keys = f_c_tmp[0].keys()
        tmp = OrderedDict()
        for ikey in keys:
            tmp_key = []
            for iset in range(nsets):
                tmp_key.append( f_c_tmp[iset][ikey] )
            
            tmp[ikey] = np.array(tmp_key)
        f_c_tmp = tmp
        # print('f_c_tmp = ',f_c_tmp)

        for ikey in keys:
            if (ntime[ikey] == 0):
                f_c[ikey][iicol,:] = f_c_tmp[ikey]
            else:
                # if vector of model outputs f_c_tmp has shape (nsets,ntime) --> must be (ntime, nsets)
                f_c[ikey][:,iicol,:] = np.transpose(f_c_tmp[ikey])

        # ----------------------------
        # NetCDF
        #    save model output C_options in NetCDF
        # ----------------------------
        if not(save_nc4 is None):
            nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
            for ikey in keys:
                if (ntime[ikey] == 0):
                    nc4_out.groups[ikey].variables['f_c_processes'][iicol:iicol+1,:] = f_c[ikey][iicol:iicol+1,:]
                else:
                    nc4_out.groups[ikey].variables['f_c_processes'][:,iicol:iicol+1,:] = f_c[ikey][:,iicol:iicol+1,:]
            nc4_out.close()
                

    # print("f_c = ",f_c)

    # store to save in pickle later
    if not( (save_pkl is None) and (save_json is None) and (save_msgpack is None) ):
        save_to_file["block_c_processes"]                = copy.deepcopy(block_c_paras)
        save_to_file["block_c_weights_nested_processes"] = copy.deepcopy(block_c_weights_nested)
        save_to_file["f_c_processes"]                    = copy.deepcopy(f_c)

    # (3c) calculate Sobol' indexes
    si   = OrderedDict()
    sti  = OrderedDict()
    msi  = OrderedDict()
    msti = OrderedDict()
    wsi  = OrderedDict()
    wsti = OrderedDict()
    for ikey in keys:
        if (ntime[ikey] == 0):
            si[ikey], sti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                          yb=f_b[ikey],
                                                          yc=f_c[ikey],
                                                          si=True,
                                                          sti=True,
                                                          method='Mai1999')
            msi[ikey]  = None
            msti[ikey] = None
            wsi[ikey]  = None
            wsti[ikey] = None
        else:
            si[ikey], sti[ikey], msi[ikey], msti[ikey], wsi[ikey], wsti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                                                                        yb=f_b[ikey],
                                                                                                        yc=f_c[ikey],
                                                                                                        si=True,
                                                                                                        sti=True,
                                                                                                        mean=True,
                                                                                                        wmean=True,
                                                                                                        method='Mai1999')
    tmp = OrderedDict()
    tmp['si']   = si
    tmp['sti']  = sti
    tmp['msi']  = msi
    tmp['msti'] = msti
    tmp['wsi']  = wsi
    tmp['wsti'] = wsti
    sobol_indexes['processes'] = tmp

    # ----------------------------
    # NetCDF
    #    save sensitivity indexes 
    # ----------------------------
    if not(save_nc4 is None):
        nc4_out = nc4.Dataset(save_nc4, "a", format="NETCDF4")
        for ikey in keys:
            nc4_out.groups[ikey].variables['si_processes'][:]  = si[ikey]
            nc4_out.groups[ikey].variables['sti_processes'][:] = sti[ikey]
            if (ntime[ikey] > 0):
                nc4_out.groups[ikey].variables['msi_processes'][:]  = msi[ikey]
                nc4_out.groups[ikey].variables['msti_processes'][:] = msti[ikey]
                nc4_out.groups[ikey].variables['wsi_processes'][:]  = wsi[ikey]
                nc4_out.groups[ikey].variables['wsti_processes'][:] = wsti[ikey]
        nc4_out.close()

    # print("si['"+ikey+"']   = ", si[ikey])
    # print("sti['"+ikey+"']  = ", sti[ikey])
    # print("msi['"+ikey+"']  = ", msi[ikey])
    # print("msti['"+ikey+"'] = ", msti[ikey])
    # print("wsi['"+ikey+"']  = ", wsi[ikey])
    # print("wsti['"+ikey+"'] = ", wsti[ikey])

    #print("shape si['"+ikey+"']   = ",np.shape(si[ikey]))
    #print("shape sti['"+ikey+"']  = ",np.shape(sti[ikey]))
    #print("shape msi['"+ikey+"']  = ",np.shape(msi[ikey]))
    #print("shape msti['"+ikey+"'] = ",np.shape(msti[ikey]))
    #print("shape wsi['"+ikey+"']  = ",np.shape(wsi[ikey]))
    #print("shape wsti['"+ikey+"'] = ",np.shape(wsti[ikey]))



    # store to save in pickle later
    if not( (save_pkl is None) and (save_json is None) and (save_msgpack is None) ):
        save_to_file["sobol_indexes"] = sobol_indexes
        
    # save to pickle
    if not(save_pkl is None):
        pickle.dump( save_to_file, open( save_pkl, "wb" ) )

    # enable serialization of numpy arrays to be stored in JSON
    # found on: https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
    class NumpyEncoder(json.JSONEncoder):
        """ Special json encoder for numpy types """
        def default(self, obj):
            if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                np.int16, np.int32, np.int64, np.uint8,
                np.uint16, np.uint32, np.uint64)):
                return int(obj)
            elif isinstance(obj, (np.float_, np.float16, np.float32, 
                np.float64)):
                return float(obj)
            elif isinstance(obj,(np.ndarray,)): #### This is the fix
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)

    # save to json
    if not(save_json is None):
        # save everything
        json_file_handle = open(save_json, 'w')
        try:
            json.dump( json.dumps(save_to_file, cls=NumpyEncoder), json_file_handle )
        finally:
            json_file_handle.close()


        # save only Sobol' indexes
        save_json_si = '/'.join(save_json.split('/')[:-1])+'/sensitivity_'+'_'.join(save_json.split('/')[-1].split('_')[1:])
        json_file_handle = open(save_json_si, 'w')
        try:
            json.dump( json.dumps(save_to_file['sobol_indexes'], cls=NumpyEncoder), json_file_handle )
        finally:
            json_file_handle.close()

    if not(save_msgpack is None):

        # save everything
        with open(save_msgpack, 'wb') as msgpack_file_handle:
            msgpack.pack(save_to_file, msgpack_file_handle)

        # save only Sobol' indexes
        save_msgpack_si = '/'.join(save_msgpack.split('/')[:-1])+'/sensitivity_'+'_'.join(save_msgpack.split('/')[-1].split('_')[1:])
        with open(save_msgpack_si, 'wb') as msgpack_file_handle:
            msgpack.pack(save_to_file['sobol_indexes'], msgpack_file_handle)

            
    # # read PICKLE with:
    # import pickle
    # setup = pickle.load( open( <save_pkl>, "rb" ) )
    # setup.keys()

    # iset=3, ipara=18
    #     setup["f_a"][:,3]           --> model_function_raven(setup["block_a_paras"][3], setup["block_a_weights_nested"][3])
    #     setup["f_b"][:,3]           --> model_function_raven(setup["block_b_paras"][3], setup["block_b_weights_nested"][3])
    #     setup["f_c_paras"][:,18,3]  --> model_function_raven(setup["block_c_paras"][18][3], setup["block_c_weights_nested_paras"][18][3])

    # failed runs:
    #

    # in case of NetCDF write token file to indicate processing is finished (since outputs are continueously appended) 
    if not(save_nc4 is None):
        token_file = '.'.join(save_nc4.split('.')[:-1])+'.token'
        ff = open(token_file,"w")
        ff.write('Done.')
        ff.close()

    # Done.
    return sobol_indexes


if __name__ == '__main__':
    
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    
    # nsets = 1000

    # # list of parameters that go into each option (numbering starts with 0)
    # # (a) simple setup
    # paras_per_option = [ 
    #     [[0], []],             # parameters of process options A1 and A2
    #     [[1], [2], [3,4]],     # parameters of process options B1, B2, and B3
    #     [[5], [6]]             # parameters of process options A1 and A2
    # ]

    # para_ranges = [ 
    #     [-np.pi,np.pi],      # parameter range of x1
    #     [-np.pi,np.pi],      # parameter range of x2
    #     [-np.pi,np.pi],      # parameter range of x3
    #     [-np.pi,np.pi],      # parameter range of x4
    #     [-np.pi,np.pi],      # parameter range of x5
    #     [-np.pi,np.pi],      # parameter range of x6
    #     [-np.pi,np.pi]       # parameter range of x7
    # ]

    # basin_prop = {}

    # def model_function(paras, weights, basin_prop, constants=None, run_id=None):
    #     # input:
    #     #     paras     ... list of model parameters scaled to their range;
    #     #                   values for all N model parameters have to be provided
    #     #                   example:
    #     #                        [ x1, x2, x3, x4, .... ]
    #     #     weights   ... list of lists of weights to weight options of each process;
    #     #                   each list of the lists need to sum up to 1.0;
    #     #                   each sublist is the N_i weights for the N_i process options of process i;
    #     #                   example:
    #     #                        [ [w_a1, w_a2, ...], [w_b1, w_b2, w_b3, ...], [w_c1, w_c2, ...], ... ]
    #     #     constants ... optional list of constants that are same for all models;
    #     #                   like parameters a and b in Ishigami-Homma function
    #     #                   example:
    #     #                        [2.0, 1.0]
    #     # output:
    #     #     model output
    #     #     example:
    #     #           { 'Q':        np.array([2.5,3.5,...,4.5]),
    #     #             'NSE':      0.6,
    #     #             'baseflow': np.array([20.5,23.5,...,24.5])
    #     #           }

    #     # check that provided number of weights is correct:
    #     # --> one weight per option per process
    #     if ( [len(ilist) for ilist in weights] != [2,3,2] ):
    #         print("Number of weights: ",[len(ilist) for ilist in weights])
    #         raise ValueError("sa_model_multiple_processes: model_function: provided number of weights must be [2,3,2]")
    #     # check if sum up to 1.0:
    #     if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
    #         print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
    #         raise ValueError("sa_model_multiple_processes: model_function: sum of weights must be 1.0 for all processes")
    #     # check if weights <= 1.0:
    #     if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
    #         print("Weights: ",weights)
    #         raise ValueError("sa_model_multiple_processes: model_function: weights must be all less or equal 1.0")
    #     # check if weights >= 0.0:
    #     if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
    #         print("Weights: ",weights)
    #         raise ValueError("sa_model_multiple_processes: model_function: weights must be all greater or equal 0.0")
    #     # check if number of parameters is correct:
    #     if (len(paras) != 7):
    #         print("Number of parameters: ",len(paras))
    #         raise ValueError("sa_model_multiple_processes: model_function: provided number of parameters must be 7")
      
    #     out = 0.0

    #     if constants is None:
    #         aa = 2.0
    #         bb = 1.0
    #     else:
    #         aa = constants[0]
    #         bb = constants[1]

    #     # ---------------
    #     # simple model
    #     # ---------------
      
    #     # process A
    #     out += ( weights[0][0] * np.sin(paras[0]) +              # A1
    #              weights[0][1] * 1.0 )                           # A2
    #     # process B
    #     out *= ( weights[1][0] * (1.0 + bb * paras[1]**4) +      # B1
    #              weights[1][1] * (1.0 + bb * paras[2]**2) +      # B2
    #              weights[1][2] * (paras[3] + bb * paras[4]) )    # B3
    #     # process C
    #     out += ( weights[2][0] * (aa * np.sin(paras[5])**2) +    # C1
    #              weights[2][1] * (1.0 + bb * paras[6]**4) )      # C2

    #     model = {}
    #     model['result'] = out

    #     return model

    # # this is calling the actual tool
    # sobol_indexes = sa_model_multiple_processes(paras_per_option, para_ranges, model_function, basin_prop, constants=None, nsets=nsets)

    # # printing
    # print("---------------------------------------")
    # print("SIMPLE SETUP")
    # print("---------------------------------------")
    # print("parameter sensitivities: ")
    # print("    S_xi  = ",astr(sobol_indexes['paras']['si']['result'],prec=5))
    # print("    ST_xi = ",astr(sobol_indexes['paras']['sti']['result'],prec=5))
    # print("")
    # print("process option sensitivities: ")
    # print("    S_Ai  = ",astr(sobol_indexes['process_options']['si']['result'],prec=5))
    # print("    ST_Ai = ",astr(sobol_indexes['process_options']['sti']['result'],prec=5))
    # print("")
    # print("process sensitivities: ")
    # print("    S_A   = ",astr(sobol_indexes['processes']['si']['result'],prec=5))
    # print("    ST_A  = ",astr(sobol_indexes['processes']['sti']['result'],prec=5))
    # print("")


    # # list of parameters that go into each option (numbering starts with 0)
    # # (b) realistic setup
    # paras_per_option = [ 
    #     [[0], [0,1]],          # parameters of process options A1 and A2
    #     [[1], [2], [3,4]],     # parameters of process options B1, B2, and B3
    #     [[5], [2,6]]           # parameters of process options A1 and A2
    # ]

    # para_ranges = [ 
    #     [-np.pi,np.pi],      # parameter range of x1
    #     [-np.pi,np.pi],      # parameter range of x2
    #     [-np.pi,np.pi],      # parameter range of x3
    #     [-np.pi,np.pi],      # parameter range of x4
    #     [-np.pi,np.pi],      # parameter range of x5
    #     [-np.pi,np.pi],      # parameter range of x6
    #     [-np.pi,np.pi]       # parameter range of x7
    # ]

    # basin_prop = {}

    # def model_function(paras, weights, basin_prop, constants=None, run_id=None):
    #     # input:
    #     #     paras     ... list of model parameters scaled to their range;
    #     #                   values for all N model parameters have to be provided
    #     #                   example:
    #     #                        [ x1, x2, x3, x4, .... ]
    #     #     weights   ... list of lists of weights to weight options of each process;
    #     #                   each list of the lists need to sum up to 1.0;
    #     #                   each sublist is the N_i weights for the N_i process options of process i;
    #     #                   example:
    #     #                        [ [w_a1, w_a2, ...], [w_b1, w_b2, w_b3, ...], [w_c1, w_c2, ...], ... ]
    #     #     constants ... optional list of constants that are same for all models;
    #     #                   like parameters a and b in Ishigami-Homma function
    #     #                   example:
    #     #                        [2.0, 1.0]
    #     # output:
    #     #     model output
    #     #     example:
    #     #           { 'Q':        np.array([2.5,3.5,...,4.5]),
    #     #             'NSE':      0.6,
    #     #             'baseflow': np.array([20.5,23.5,...,24.5])
    #     #           }

    #     # check that provided number of weights is correct:
    #     # --> one weight per option per process
    #     if ( [len(ilist) for ilist in weights] != [2,3,2] ):
    #         print("Number of weights: ",[len(ilist) for ilist in weights])
    #         raise ValueError("sa_model_multiple_processes: model_function: provided number of weights must be [2,3,2]")
    #     # check if sum up to 1.0:
    #     if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
    #         print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
    #         raise ValueError("sa_model_multiple_processes: model_function: sum of weights must be 1.0 for all processes")
    #     # check if weights <= 1.0:
    #     if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
    #         print("Weights: ",weights)
    #         raise ValueError("sa_model_multiple_processes: model_function: weights must be all less or equal 1.0")
    #     # check if weights >= 0.0:
    #     if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
    #         print("Weights: ",weights)
    #         raise ValueError("sa_model_multiple_processes: model_function: weights must be all greater or equal 0.0")
    #     # check if number of parameters is correct:
    #     if (len(paras) != 7):
    #         print("Number of parameters: ",len(paras))
    #         raise ValueError("sa_model_multiple_processes: model_function: provided number of parameters must be 7")
      
    #     out = 0.0

    #     if constants is None:
    #         aa = 2.0
    #         bb = 1.0
    #     else:
    #         aa = constants[0]
    #         bb = constants[1]


    #     # ---------------
    #     # realistic model
    #     # ---------------
      
    #     # process D
    #     out += ( weights[0][0] * np.sin(paras[0]) +                            # D1
    #              weights[0][1] * (paras[0]+paras[1]**2) )                      # D2
    #     # process E
    #     out *= ( weights[1][0] * (1.0 + bb * paras[1]**4) +                    # E1
    #              weights[1][1] * (1.0 + bb * paras[2]**2) +                    # E2
    #              weights[1][2] * (paras[3] + bb * paras[4]) )                  # E3
    #     # process F
    #     out += ( weights[2][0] * (aa * np.sin(paras[5])**2) +                  # F1          
    #              weights[2][1] * (1.0 + bb * paras[6]**4) + paras[2]**2 )      # F2

    #     model = {}
    #     model['result'] = out

    #     return model

    # # this is calling the actual tool
    # sobol_indexes = sa_model_multiple_processes(paras_per_option, para_ranges, model_function, basin_prop, constants=None, nsets=nsets)

    # # printing
    # print("---------------------------------------")
    # print("REALISTIC SETUP")
    # print("---------------------------------------")
    # print("parameter sensitivities: ")
    # print("    S_xi  = ",astr(sobol_indexes['paras']['si']['result'],prec=5))
    # print("    ST_xi = ",astr(sobol_indexes['paras']['sti']['result'],prec=5))
    # print("")
    # print("process option sensitivities: ")
    # print("    S_Ai  = ",astr(sobol_indexes['process_options']['si']['result'],prec=5))
    # print("    ST_Ai = ",astr(sobol_indexes['process_options']['sti']['result'],prec=5))
    # print("")
    # print("process sensitivities: ")
    # print("    S_A   = ",astr(sobol_indexes['processes']['si']['result'],prec=5))
    # print("    ST_A  = ",astr(sobol_indexes['processes']['sti']['result'],prec=5))
    # print("")

