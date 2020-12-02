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
#     run sa_model_multiple_processes_DVM.py -p sa_model_multiple_processes_DVM -t pdf

#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import copy
import sobol
import sobol_index
import pickle
from autostring import astr
import PieShareDistribution as psd
from collections import OrderedDict

__all__ = ['sa_model_multiple_processes_DVM']

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

def sa_model_multiple_processes_DVM(paras_per_option, para_ranges, model_function, constants=None, nsets=None, budget=None, save=None, nsets_parasets=128):
    """
        This function that estimates the Sobol' sensitivity estimates for models with mutiple process options. 
        The options and the parameters of those options are given in a nested list 'paras_per_option'. 
        Further, the range of each parameter needs to be given and a function that returns model outputs 
        when a set of parameters and weights are given. The weights are used to weight all the process option 
        outputs. Hence, the returned model output is a weighted model output. The sampling of all weights and 
        parameters is done internally in this method. Sobol' sequences are used for this purpose.

        Definition
        ----------
        def sa_model_multiple_processes_DVM(paras_per_option, para_ranges, model_function, constants=None, nsets=None)


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
        constants           list                 optional: list of constants that 'model_function' might need
                                                 default: None
        nsets               integer              optional: number of reference parameter sets
                                                 default: 1024
        budget              integer              optional: total number of model runs allowed for analysis; overwrites nsets
                                                 budget = nsets / ((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2))
                                                 default: None
        save                string               filename to save parameter sets and respective model outputs in pickle file
                                                 default: None (nothing saved to file)
        nsets_parasets      integer              optional: number of pre-sampled parameter sets per process
                                                 default: 128
        

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
        >>> nsets = 1024
        >>> nsets_parasets = 128

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
        >>> def model_function(paras, weights, constants=None, run_id=None):
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
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: provided number of weights must be [2,3,2]")
        ...     # check if sum up to 1.0:
        ...     if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
        ...         print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: sum of weights must be 1.0 for all processes")
        ...     # check if weights <= 1.0:
        ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
        ...         print("Weights: ",weights)
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: weights must be all less or equal 1.0")
        ...     # check if weights >= 0.0:
        ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
        ...         print("Weights: ",weights)
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: weights must be all greater or equal 0.0")
        ...     # check if number of parameters is correct:
        ...     if (len(paras) != 7):
        ...         print("Number of parameters: ",len(paras))
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: provided number of parameters must be 7")
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
        ...     model['result'] = out
        ...
        ...     return model

        >>> # this is calling the actual tool
        >>> sobol_indexes = sa_model_multiple_processes_DVM(paras_per_option, para_ranges, model_function, constants=None, nsets=nsets, nsets_parasets=nsets_parasets)

        >>> # printing
        >>> print("process sensitivities:        S_A   = ",astr(sobol_indexes['processes']['si']['result'],prec=5))
        process sensitivities:        S_A   =  ['0.06200' '0.09874' '0.68941']
        >>> print("process sensitivities:        ST_A  = ",astr(sobol_indexes['processes']['sti']['result'],prec=5))
        process sensitivities:        ST_A  =  ['0.17770' '0.23122' '0.70231']

        --------------------------------------------------
        Realistic setup
        --------------------------------------------------

        >>> # list of parameters that go into each option (numbering starts with 0)
        >>> # (a) simple setup
        >>> paras_per_option = [ 
        ...       [[0], [0,1]],             # parameters of process options A1 and A2
        ...       [[1], [2], [3,4]],        # parameters of process options B1, B2, and B3
        ...       [[5], [2,6]]              # parameters of process options A1 and A2
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
        >>> def model_function(paras, weights, constants=None, run_id=None):
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
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: provided number of weights must be [2,3,2]")
        ...     # check if sum up to 1.0:
        ...     if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
        ...         print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: sum of weights must be 1.0 for all processes")
        ...     # check if weights <= 1.0:
        ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
        ...         print("Weights: ",weights)
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: weights must be all less or equal 1.0")
        ...     # check if weights >= 0.0:
        ...     if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
        ...         print("Weights: ",weights)
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: weights must be all greater or equal 0.0")
        ...     # check if number of parameters is correct:
        ...     if (len(paras) != 7):
        ...         print("Number of parameters: ",len(paras))
        ...         raise ValueError("sa_model_multiple_processes_DVM: model_function: provided number of parameters must be 7")
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
        ...     # realistic model
        ...     # ---------------
        ...
        ...     # process D
        ...     out += ( weights[0][0] * np.sin(paras[0]) +                            # D1
        ...              weights[0][1] * (paras[0]+paras[1]**2) )                      # D2
        ...     # process E
        ...     out *= ( weights[1][0] * (1.0 + bb * paras[1]**4) +                    # E1
        ...              weights[1][1] * (1.0 + bb * paras[2]**2) +                    # E2
        ...              weights[1][2] * (paras[3] + bb * paras[4]) )                  # E3
        ...     # process F
        ...     out += ( weights[2][0] * (aa * np.sin(paras[5])**2) +                  # F1          
        ...              weights[2][1] * (1.0 + bb * paras[6]**4) + paras[2]**2 )      # F2
        ...
        ...     model = {}
        ...     model['result'] = out
        ...
        ...     return model

        >>> # this is calling the actual tool
        >>> sobol_indexes = sa_model_multiple_processes_DVM(paras_per_option, para_ranges, model_function, constants=None, nsets=nsets, nsets_parasets=nsets_parasets)

        >>> # printing
        >>> print("process sensitivities:        S_A   = ",astr(sobol_indexes['processes']['si']['result'],prec=5))
        process sensitivities:        S_A   =  ['0.07335' '0.55102' '0.04634']
        >>> print("process sensitivities:        ST_A  = ",astr(sobol_indexes['processes']['sti']['result'],prec=5))
        process sensitivities:        ST_A  =  ['0.34189' '0.87986' '0.06088']


        
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

        Copyright 2020 Juliane Mai - juliane.mai@uwaterloo.ca


        History
        -------
        Written,  Juliane Mai, April 2020
    """

    # initialize return variable
    sobol_indexes = OrderedDict()

    para_ranges = np.array(para_ranges)

    # number of parameters
    nparas   = np.shape(para_ranges)[0]
    # number of options per process
    noptions = np.array([ len(oo) for oo in paras_per_option ])
    # number of processes
    nprocess = np.shape(noptions)[0]
    # number of weights required
    nweights = np.sum(noptions)

    if ( nsets is None ):
        nsets = 1000

    if ( nsets_parasets is None ):
        nsets_parasets = 128
        
    if not( budget is None ):
        # overwrite nsets if budget is given
        nsets = budget / ((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2))
        if nsets <= 0:
            print("Minimal budget:     ",((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2)))
            print("Recommended budget: ",((nprocess+2)+(nweights+np.sum(noptions)+2)+(nparas+nweights+2))*1000)
            raise ValueError("sa_model_multiple_processes_DVM: Budget is too small!")

    if not(save is None):
        save_to_pickle = {}


    # sample parametersets for each process
    nsets_per_process = nsets_parasets  # DVM used 128 for observations, crop parameters, soil parameters, weather, and 8 for model structure
    para_sets = [ [] for iprocess in np.arange(nprocess) ]
    paras_per_process      = [ list(np.unique([item for sublist in iprocess for item in sublist])) for iprocess in paras_per_option ]  # [[0, 1], [1, 2, 3, 4], [2, 5, 6]]
    paras_per_process_flat = [item for sublist in paras_per_process for item in sublist]         # [[0, 1, 1, 2, 3, 4, 2, 5, 6]
    nparas_per_process     = [ len(iprocess)   for iprocess in paras_per_process ]               # [2,4,3]
    nweights_per_process   = [ len(iprocess)   for iprocess in paras_per_option ]                # [2,3,2]
    for iprocess in np.arange(nprocess):

        # sampling
        # (last column is a dummy to convert rn into weight later)
        para_sets[iprocess] = sobol.i4_sobol_generate((nparas_per_process[iprocess]+nweights_per_process[iprocess]),nsets_per_process,40000) 
        para_sets[iprocess] = np.transpose(para_sets[iprocess])
        para_sets[iprocess][:,-1] = -9999.0

        # scaling
        for iipara,ipara in enumerate(paras_per_process[iprocess]):
            para_sets[iprocess][:,iipara] *= (para_ranges[ipara,1]-para_ranges[ipara,0])
            para_sets[iprocess][:,iipara] += para_ranges[ipara,0]

        # convert weight random numbers in weights
        para_sets[iprocess][:,nparas_per_process[iprocess]:] = psd.PieShareDistribution(nsets_per_process,nweights_per_process[iprocess],remainder=True,randomnumbers=para_sets[iprocess][:,nparas_per_process[iprocess]:nparas_per_process[iprocess]+nweights_per_process[iprocess]-1])
        

    # (A) Sampling integer number that indicates parameter set ID for each process in unit interval using Sobol' sequences
    sobol_sets = sobol.i4_sobol_generate((nprocess)*2,nsets,40000)
    sobol_sets = np.transpose(sobol_sets)

    # (B) Convert into integer
    sobol_sets_int = np.array([ [ int(ival) for ival in sobol_sets[iset]*nsets_per_process ] for iset in np.arange(np.shape(sobol_sets)[0]) ])

    # (C) Construct block A, B and Ci (for processes)
    block_a  = copy.deepcopy(sobol_sets_int[:,0*nprocess:1*nprocess])    
    block_b  = copy.deepcopy(sobol_sets_int[:,1*nprocess:2*nprocess])
    block_c  = [[] for iprocess in np.arange(nprocess)]
    for iprocess in np.arange(nprocess):
        block_c[iprocess]             = copy.deepcopy(block_a)
        block_c[iprocess][:,iprocess] = copy.deepcopy(block_b[:,iprocess])

    # --------------------------------------------
    # Derive model output for A sets
    # --------------------------------------------
        
    f_a = []
    for iset in np.arange(nsets):

        paraweightset_proc = [ para_sets[iprocess][block_a[iset,iprocess]] for iprocess in np.arange(nprocess) ]      # parameter values and weights for each process
        paraset_proc_flat  = [ item for isublist,sublist in enumerate(paraweightset_proc) for item in sublist[0:nparas_per_process[isublist]] ]
        paraset            = [ -9999.0 for ipara in np.arange(nparas) ]   # parameter set combined in right order
        weightset          = [ sublist[nparas_per_process[isublist]:nparas_per_process[isublist]+nweights_per_process[isublist]] for isublist,sublist in enumerate(paraweightset_proc) ] # weight    set combined in right order
        for iipara,ipara in enumerate(paras_per_process_flat):
            # MAJOR PROBLEM is that parameters that appear in multiple processes get overwritten
            paraset[ipara] = paraset_proc_flat[iipara]
            
        f_a.append( model_function(paraset, weightset, constants=constants, run_id="a_set_"+str(iset)) )

    # convert list of dicts into dict of lists:
    #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
    keys = f_a[0].keys()
    tmp = {}
    for ikey in keys:
        tmp_key = []
        for iset in np.arange(nsets):
            tmp_key.append( f_a[iset][ikey] )
            
        tmp[ikey] = np.array(tmp_key)
    f_a = tmp
    # print('f_a = ',f_a)

    # --------------------------------------------
    # Derive model output for B sets
    # --------------------------------------------

    f_b = []
    for iset in np.arange(nsets):

        paraweightset_proc = [ para_sets[iprocess][block_b[iset,iprocess]] for iprocess in np.arange(nprocess) ]      # parameter values and weights for each process
        paraset_proc_flat  = [ item for isublist,sublist in enumerate(paraweightset_proc) for item in sublist[0:nparas_per_process[isublist]] ]
        paraset            = [ -9999.0 for ipara in np.arange(nparas) ]   # parameter set combined in right order
        weightset          = [ sublist[nparas_per_process[isublist]:nparas_per_process[isublist]+nweights_per_process[isublist]] for isublist,sublist in enumerate(paraweightset_proc) ] # weight    set combined in right order
        for iipara,ipara in enumerate(paras_per_process_flat):
            # MAJOR PROBLEM is that parameters that appear in multiple processes get overwritten
            paraset[ipara] = paraset_proc_flat[iipara]
            
        f_b.append( model_function(paraset, weightset, constants=constants, run_id="b_set_"+str(iset)) )
        
    # convert list of dicts into dict of lists:
    #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
    keys = f_b[0].keys()
    tmp = {}
    for ikey in keys:
        tmp_key = []
        for iset in np.arange(nsets):
            tmp_key.append( f_b[iset][ikey] )
            
        tmp[ikey] = np.array(tmp_key)
    f_b = tmp

    # --------------------------------------------
    # Derive model output for Ci sets
    # --------------------------------------------

    ntime = {}
    for ikey in keys:
        if (len(np.shape(f_b[ikey])) == 2):
            ntime[ikey] = np.shape(f_b[ikey])[0]
        else:
            ntime[ikey] = 0
        
    # (1) parameter sensitivities:
    #     main effects:  [S_x1,  S_x2,  ..., S_w1,  S_w2,  ...]
    #     total effects: [ST_x1, ST_x2, ..., ST_w1, ST_w2, ...]
    f_c = {}
    for ikey in keys:
        if (ntime[ikey] == 0):
            f_c[ikey] = np.ones([nprocess,nsets]) * -9999.
        else:
            print('nsets = ',nsets)
            tttmp = np.ones([ntime[ikey],nprocess,nsets]) * -9999.
            f_c[ikey] = tttmp

    for iprocess in np.arange(nprocess):

        f_c_tmp = []
        for iset in np.arange(nsets):

            paraweightset_proc = [ para_sets[iiprocess][block_c[iprocess][iset,iiprocess]] for iiprocess in np.arange(nprocess) ]      # parameter values and weights for each process
            paraset_proc_flat  = [ item for isublist,sublist in enumerate(paraweightset_proc) for item in sublist[0:nparas_per_process[isublist]] ]
            paraset            = [ -9999.0 for ipara in np.arange(nparas) ]   # parameter set combined in right order
            weightset          = [ sublist[nparas_per_process[isublist]:nparas_per_process[isublist]+nweights_per_process[isublist]] for isublist,sublist in enumerate(paraweightset_proc) ] # weight    set combined in right order
            for iipara,ipara in enumerate(paras_per_process_flat):
                # MAJOR PROBLEM is that parameters that appear in multiple processes get overwritten
                paraset[ipara] = paraset_proc_flat[iipara]
                
            f_c_tmp.append( model_function(paraset, weightset, constants=constants,
                                                run_id="c_process_"+str(iprocess)+"_set_"+str(iset)) )
        f_c_tmp = np.array(f_c_tmp)

        # convert list of dicts into dict of lists:
        #       [{'result_1':1.0, 'result_2:2.0'}, {'result_1':3.0, 'result_2:4.0'}, ...] --> [{'result_1':[1.0,3.0,...],'result_2':[2.0,4.0,...]}]
        keys = f_c_tmp[0].keys()
        tmp = {}
        for ikey in keys:
            tmp_key = []
            for iset in np.arange(nsets):
                tmp_key.append( f_c_tmp[iset][ikey] )
            
            tmp[ikey] = np.array(tmp_key)
        f_c_tmp = tmp
        # print('f_c_tmp = ',f_c_tmp)

        for ikey in keys:
            if (ntime[ikey] == 0):
                f_c[ikey][iprocess,:] = f_c_tmp[ikey]
            else:
                # if vector of model outputs f_c_tmp has shape (nsets,ntime) --> must be (ntime, nsets)
                f_c[ikey][:,iprocess,:] = np.transpose(f_c_tmp[ikey])

    # print("f_c = ",f_c)
        
    # (1c) calculate Sobol' indexes
    si   = {}
    sti  = {}
    msi  = {}
    msti = {}
    wsi  = {}
    wsti = {}
    for ikey in keys:
        if (ntime[ikey] == 0):
            si[ikey], sti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                          yb=f_b[ikey],
                                                          yc=f_c[ikey],
                                                          si=True,
                                                          sti=True,
                                                          method='Mai1999')
        else:
            si[ikey], sti[ikey], msi[ikey], msti[ikey], wsi[ikey], wsti[ikey] = sobol_index.sobol_index(ya=f_a[ikey],
                                                                                                        yb=f_b[ikey],
                                                                                                        yc=f_c[ikey],
                                                                                                        si=True,
                                                                                                        sti=True,
                                                                                                        mean=True,
                                                                                                        wmean=True,
                                                                                                        method='Mai1999')

    for ikey in keys:
        if (ntime[ikey] == 0):
            tmp = {}
            tmp['si']  = si
            tmp['sti'] = sti
            sobol_indexes['processes'] = tmp
        else:
            tmp = {}
            tmp['si']   = si
            tmp['sti']  = sti
            tmp['msi']  = msi
            tmp['msti'] = msti
            tmp['wsi']  = wsi
            tmp['wsti'] = wsti
            sobol_indexes['processes'] = tmp
            
            print("si['"+ikey+"']   = ", si[ikey])
            print("sti['"+ikey+"']  = ", sti[ikey])
            print("msi['"+ikey+"']  = ", msi[ikey])
            print("msti['"+ikey+"'] = ", msti[ikey])
            print("wsi['"+ikey+"']  = ", wsi[ikey])
            print("wsti['"+ikey+"'] = ", wsti[ikey])

    # store to save in pickle later
    if not(save is None):
        save_to_pickle["sobol_indexes"] = sobol_indexes
        
    # save to pickle
    if not(save is None):
        pickle.dump( save_to_pickle, open( save, "wb" ) )

        
    # Done.
    return sobol_indexes


if __name__ == '__main__':
    
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

