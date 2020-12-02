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
#     python __run_benchmark_DVM.py 

from __future__ import print_function

"""
Benchmark example to test Sensitivity Analysis for models with multiple process options 
using the DVM method introduced in:

Baroni, G., & Tarantola, S. (2014). 
A General Probabilistic Framework for uncertainty and global sensitivity analysis of 
deterministic models: A hydrological case study. 
Environmental Modelling & Software, 51, 26-34. 
http://doi.org/10.1016/j.envsoft.2013.09.022

History
-------
Written,  JM, Apr 2020
"""

# -------------------------------------------------------------------------
# Command line arguments - if script
#

# Comment|Uncomment - Begin
if __name__ == '__main__':

    import argparse
    import numpy as np
    import pickle

    do_disjoint_example_DVM     = True     # SA of disjoint-parameter    theoretical model - DVM
    do_shared_example_DVM       = True     # SA of shared-parameter theoretical model - DVM

    nsets          = [50,100,200,500,1000,2000,5000,10000,20000,50000,100000]      # number of Sobol sequences                        # DVM used 1024
    nsets_parasets = [32, 64, 128, 256, 512, 1024]                                 # number of pre-sampled parameter sets per process # DVM used 128
    nboot          = 1                                                             # Set to 1 for single run of SI and STI calculation
    
    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''Benchmark example to test Sensitivity Analysis for models with multiple process options.''')

    parser.add_argument('-b', '--nboot', action='store',
                        default=nboot, dest='nboot', metavar='nboot',
                        help='Number of bootstrap samples (default: nboot=10).')
    parser.add_argument('-n', '--nsets', action='store',
                        default=nsets, dest='nsets', metavar='nsets',
                        help='Number of Sobol reference samples (default: nsets=10).')
    parser.add_argument('-m', '--nsets_parasets', action='store',
                        default=nsets_parasets, dest='nsets_parasets', metavar='nsets_parasets',
                        help='Number of pre-sampled parameter sets per process (default: nsets_parasets=10).')

    args           = parser.parse_args()
    nboot          = np.int(args.nboot)
    nsets          = [ np.int(aa) for aa in args.nsets ]
    nsets_parasets = [ np.int(aa) for aa in args.nsets_parasets ]

    del parser, args
    # Comment|Uncomment - End

    # -----------------------
    # add subolder scripts/lib to search path
    # -----------------------
    import sys
    import os 
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(dir_path+'/../../lib')


    # -------------------------------------------------------------------------
    # Function definition - if function
    #

    import numpy as np
    import sa_model_multiple_processes_DVM
    import copy
    import sobol
    import time
    from   autostring      import astr     # in lib/
    
    t1 = time.time()


    # -------------------------------------------------------------------------
    # Setup Calculations
    #
        
    theo_si_sti_disjoint = {   # individual models with all possible process option combinations
                      'A1-B1-C1': np.array([[0.38302,0.00000,np.nan,np.nan,np.nan,0.00091,np.nan],[0.99909,0.61606,np.nan,np.nan,np.nan,0.00091,np.nan]]), 
                      'A1-B1-C2': np.array([[0.17167,0.00000,np.nan,np.nan,np.nan,np.nan,0.55222],[0.44778,0.27611,np.nan,np.nan,np.nan,np.nan,0.55222]]), 
                      'A1-B2-C1': np.array([[0.65581,np.nan,0.00000,np.nan,np.nan,0.03564,np.nan],[0.96436,np.nan,0.30856,np.nan,np.nan,0.03564,np.nan]]), 
                      'A1-B2-C2': np.array([[0.01337,np.nan,0.00000,np.nan,np.nan,np.nan,0.98034],[0.01966,np.nan,0.00629,np.nan,np.nan,np.nan,0.98034]]), 
                      'A1-B3-C1': np.array([[0.00000,np.nan,np.nan,0.00000,0.00000,0.13193,np.nan],[0.86807,np.nan,np.nan,0.43403,0.43403,0.13193,np.nan]]), 
                      'A1-B3-C2': np.array([[0.00000,np.nan,np.nan,0.00000,0.00000,np.nan,0.99515],[0.00485,np.nan,np.nan,0.00243,0.00243,np.nan,0.99515]]), 
                      'A2-B1-C1': np.array([[np.nan,0.99926,np.nan,np.nan,np.nan,0.00074,np.nan],[np.nan,0.99926,np.nan,np.nan,np.nan,0.00074,np.nan]]), 
                      'A2-B1-C2': np.array([[np.nan,0.50000,np.nan,np.nan,np.nan,np.nan,0.50000],[np.nan,0.50000,np.nan,np.nan,np.nan,np.nan,0.50000]]), 
                      'A2-B2-C1': np.array([[np.nan,np.nan,0.94541,np.nan,np.nan,0.05459,np.nan],[np.nan,np.nan,0.94541,np.nan,np.nan,0.05459,np.nan]]), 
                      'A2-B2-C2': np.array([[np.nan,np.nan,0.01267,np.nan,np.nan,np.nan,0.98733],[np.nan,np.nan,0.01267,np.nan,np.nan,np.nan,0.98733]]), 
                      'A2-B3-C1': np.array([[np.nan,np.nan,np.nan,0.46469,0.46469,0.07062,np.nan],[np.nan,np.nan,np.nan,0.46469,0.46469,0.07062,np.nan]]), 
                      'A2-B3-C2': np.array([[np.nan,np.nan,np.nan,0.00483,0.00483,np.nan,0.99034],[np.nan,np.nan,np.nan,0.00483,0.00483,np.nan,0.99034]]),
                      
                      # all process options weighted to one model: 3 processes (a,b,c)
                      'wa_A-wb_B-wc_C_process': np.array([[0.04999,0.09857,0.75286],
                                                          [0.14856,0.19715,0.75286]]),
                      # --> mae error = 0.0149428348938     (N=    1,000)
                      # --> mae error = 0.00307764826832    (N=   10,000)  
                      # --> mae error = 0.000469681942577   (N=  100,000)
                      # --> mae error = 2.5615005695e-05    (N=1,000,000)
                  }

    theo_si_sti_shared = {   # individual models with all possible process option combinations
                      'A1-B1-C1': np.array([[0.38300,0.00000,np.nan,np.nan,np.nan,0.00090,np.nan],[0.99910,0.61610,np.nan,np.nan,np.nan,0.00090,np.nan]]), 
                      'A1-B1-C2': np.array([[0.17050,0.00000,0.00700,np.nan,np.nan,np.nan,0.54830],[0.44460,0.27420,0.00700,np.nan,np.nan,np.nan,0.54830]]), 
                      'A1-B2-C1': np.array([[0.65580,np.nan,0.00000,np.nan,np.nan,0.03560,np.nan],[0.96440,np.nan,0.30860,np.nan,np.nan,0.03560,np.nan]]), 
                      'A1-B2-C2': np.array([[0.01320,np.nan,0.01240,np.nan,np.nan,np.nan,0.96820],[0.01940,np.nan,0.01860,np.nan,np.nan,np.nan,0.96820]]), 
                      'A1-B3-C1': np.array([[0.00000,np.nan,np.nan,0.00000,0.00000,0.13190,np.nan],[0.86810,np.nan,np.nan,0.43400,0.43400,0.13190,np.nan]]), 
                      'A1-B3-C2': np.array([[0.00000,np.nan,0.01260,0.00000,0.00000,np.nan,0.98260],[0.00480,np.nan,0.01260,0.00240,0.00240,np.nan,0.98260]]), 
                      'A2-B1-C1': np.array([[0.02420,0.93690,np.nan,np.nan,np.nan,0.00000,np.nan],[0.06310,0.97580,np.nan,np.nan,np.nan,0.00000,np.nan]]), 
                      'A2-B1-C2': np.array([[0.02390,0.92580,0.00010,np.nan,np.nan,np.nan,0.01170],[0.06240,0.96430,0.00010,np.nan,np.nan,np.nan,0.01170]]), 
                      'A2-B2-C1': np.array([[0.14500,0.38160,0.22440,np.nan,np.nan,0.00120,np.nan],[0.21320,0.56120,0.47220,np.nan,np.nan,0.00120,np.nan]]), 
                      'A2-B2-C2': np.array([[0.05230,0.13770,0.13770,np.nan,np.nan,np.nan,0.58300],[0.07690,0.20240,0.22710,np.nan,np.nan,np.nan,0.58300]]), 
                      'A2-B3-C1': np.array([[0.00000,0.00000,np.nan,0.23690,0.23690,0.00330,np.nan],[0.14400,0.37900,np.nan,0.49830,0.49830,0.00330,np.nan]]), 
                      'A2-B3-C2': np.array([[0.00000,0.00000,0.01040,0.04270,0.04270,np.nan,0.80980],[0.02600,0.06840,0.01040,0.08990,0.08990,np.nan,0.80980]]), 
                      #                      
                      # all process options weighted to one model: 3 processes (a,b,c)
                      'wa_A-wb_B-wc_C_process': np.array([[0.61237,0.60821,0.06485],
                                                          [0.87597,0.86055,0.06697]]), 
                      # --> mae error = 0.018469862885     (N=    1,000)
                      # --> mae error = 0.018132808782     (N=   10,000)  
                      # --> mae error = 0.00126200074934   (N=  100,000)
                      # --> mae error = 0.000984658078514  (N=1,000,000)
                  }

    # -------------------------------------------------------------------------
    # Set processes and options
    # -------------------------------------------------------------------------

    # list of parameters that go into each option (numbering starts with 0)
    options_paras_disjoint = [
        [[0],[]],              # parameters of process options A1 and A2       
        [[1],[2],[3,4]],       # parameters of process options B1, B2, and B3  
        [[5],[6]]              # parameters of process options A1 and A2       
        ]
        
    options_paras_shared = [ 
        [[0], [0,1]],          # parameters of process options D1 and D2
        [[1], [2], [3,4]],     # parameters of process options E1, E2, and E3
        [[5], [2,6]]           # parameters of process options F1 and F2
        ]

    para_ranges = [ 
        [-np.pi,np.pi],      # parameter range of x1
        [-np.pi,np.pi],      # parameter range of x2
        [-np.pi,np.pi],      # parameter range of x3
        [-np.pi,np.pi],      # parameter range of x4
        [-np.pi,np.pi],      # parameter range of x5
        [-np.pi,np.pi],      # parameter range of x6
        [-np.pi,np.pi]       # parameter range of x7
    ]

    def model_function_disjoint(paras, weights, constants=None, run_id=None):
        # input:
        #     paras     ... list of model parameters scaled to their range;
        #                   values for all N model parameters have to be provided
        #                   example:
        #                        [ x1, x2, x3, x4, .... ]
        #     weights   ... list of lists of weights to weight options of each process;
        #                   each list of the lists need to sum up to 1.0;
        #                   each sublist is the N_i weights for the N_i process options of process i;
        #                   example:
        #                        [ [w_a1, w_a2, ...], [w_b1, w_b2, w_b3, ...], [w_c1, w_c2, ...], ... ]
        #     constants ... optional list of constants that are same for all models;
        #                   like parameters a and b in Ishigami-Homma function
        #                   example:
        #                        [2.0, 1.0]
        #     run_id    ... optional name of this run (to, e.g., print or store in a file)
        #                   example:
        #                        run_aset_001
        # output:
        #     model output in dictionary
        #     example:
        #          model['out'] = 7.4

        # check that provided number of weights is correct:
        # --> one weight per option per process
        if ( [len(ilist) for ilist in weights] != [2,3,2] ):
            print("Number of weights: ",[len(ilist) for ilist in weights])
            raise ValueError("sa_model_multiple_processes: model_function: provided number of weights must be [2,3,2]")
        # check if sum up to 1.0:
        if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
            print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
            raise ValueError("sa_model_multiple_processes: model_function: sum of weights must be 1.0 for all processes")
        # check if weights <= 1.0:
        if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
            print("Weights: ",weights)
            raise ValueError("sa_model_multiple_processes: model_function: weights must be all less or equal 1.0")
        # check if weights >= 0.0:
        if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
            print("Weights: ",weights)
            raise ValueError("sa_model_multiple_processes: model_function: weights must be all greater or equal 0.0")
        # check if number of parameters is correct:
        if (len(paras) != 7):
            print("Number of parameters: ",len(paras))
            raise ValueError("sa_model_multiple_processes: model_function: provided number of parameters must be 7")
        
        model = {}

        if constants is None:
            aa = 2.0
            bb = 1.0
        else:
            aa = constants[0]
            bb = constants[1]

        # ---------------
        # disjoint-parameter model
        # ---------------
        out = 0.0
        
        # process A
        out += ( weights[0][0] * np.sin(paras[0]) +              # A1
                 weights[0][1] * 1.0 )                           # A2
        # process B
        out *= ( weights[1][0] * (1.0 + bb * paras[1]**4) +      # B1
                 weights[1][1] * (1.0 + bb * paras[2]**2) +      # B2
                 weights[1][2] * (paras[3] + bb * paras[4]) )    # B3
        # process C
        out += ( weights[2][0] * (aa * np.sin(paras[5])**2) +    # C1
                 weights[2][1] * (1.0 + bb * paras[6]**4) )      # C2

        model['out'] = out

        return model

    def model_function_shared(paras, weights, constants=None, run_id=None):
        # input:
        #     paras     ... list of model parameters scaled to their range;
        #                   values for all N model parameters have to be provided
        #                   example:
        #                        [ x1, x2, x3, x4, .... ]
        #     weights   ... list of lists of weights to weight options of each process;
        #                   each list of the lists need to sum up to 1.0;
        #                   each sublist is the N_i weights for the N_i process options of process i;
        #                   example:
        #                        [ [w_a1, w_a2, ...], [w_b1, w_b2, w_b3, ...], [w_c1, w_c2, ...], ... ]
        #     constants ... optional list of constants that are same for all models;
        #                   like parameters a and b in Ishigami-Homma function
        #                   example:
        #                        [2.0, 1.0]
        #     run_id    ... optional name of this run (to, e.g., print or store in a file)
        #                   example:
        #                        run_aset_001
        # output:
        #     model output in dictionary
        #     example:
        #          model['out'] = 7.4

        # check that provided number of weights is correct:
        # --> one weight per option per process
        if ( [len(ilist) for ilist in weights] != [2,3,2] ):
            print("Number of weights: ",[len(ilist) for ilist in weights])
            raise ValueError("sa_model_multiple_processes: model_function: provided number of weights must be [2,3,2]")
        # check if sum up to 1.0:
        if ( np.any(np.array([np.sum(ilist) for ilist in weights]) != 1.0) ):
            print("Sum of weights per process: ",[np.sum(ilist) for ilist in weights])
            raise ValueError("sa_model_multiple_processes: model_function: sum of weights must be 1.0 for all processes")
        # check if weights <= 1.0:
        if ( np.any(np.array([item for ilist in weights for item in ilist]) > 1.0) ):
            print("Weights: ",weights)
            raise ValueError("sa_model_multiple_processes: model_function: weights must be all less or equal 1.0")
        # check if weights >= 0.0:
        if ( np.any(np.array([item for ilist in weights for item in ilist]) < 0.0) ):
            print("Weights: ",weights)
            raise ValueError("sa_model_multiple_processes: model_function: weights must be all greater or equal 0.0")
        # check if number of parameters is correct:
        if (len(paras) != 7):
            print("Number of parameters: ",len(paras))
            raise ValueError("sa_model_multiple_processes: model_function: provided number of parameters must be 7")
        
        model = {}

        if constants is None:
            aa = 2.0
            bb = 1.0
        else:
            aa = constants[0]
            bb = constants[1]


        # ---------------
        # shared-parameter model
        # ---------------
        out = 0.0
        
        # process D
        out += ( weights[0][0] * np.sin(paras[0]) +                            # D1
                 weights[0][1] * (paras[0]+paras[1]**2) )                      # D2
        # process E
        out *= ( weights[1][0] * (1.0 + bb * paras[1]**4) +                    # E1
                 weights[1][1] * (1.0 + bb * paras[2]**2) +                    # E2
                 weights[1][2] * (paras[3] + bb * paras[4]) )                  # E3
        # process F
        out += ( weights[2][0] * (aa * np.sin(paras[5])**2) +                  # F1          
                 weights[2][1] * (1.0 + bb * paras[6]**4) + paras[2]**2 )      # F2

        model['out'] = out
        
        return model


    # ----------------------------------------------------
    # Sensitivity Analysis of 'disjoint-parameter setup' (DVM)
    # ----------------------------------------------------
    if do_disjoint_example_DVM:

        save_to_pickle = {}
        save_to_pickle["theo_si_sti_disjoint"] = theo_si_sti_disjoint

        sobol_indexes_disjoint = {}
        for iset in nsets:

            
            print('')
            print('SA of disjoint-parameter setup: nsets = ',iset)

            sobol_indexes_disjoint_tmp = {}
            for jset in nsets_parasets:

                print('   SA of disjoint-parameter setup: nsets_parasets = ',jset)

                sobol_indexes_disjoint_tmp[str(jset)] = sa_model_multiple_processes_DVM.sa_model_multiple_processes_DVM(options_paras_disjoint,
                                                                                                         para_ranges,
                                                                                                         model_function_disjoint,
                                                                                                         constants=None,
                                                                                                         nsets=iset,
                                                                                                         nsets_parasets=jset)
                # process sensitivities
                print("      si  numeric  (A,B,C) = ",astr(sobol_indexes_disjoint_tmp[str(jset)]['processes']['si']['out'],prec=5))
                print("      si  analytic (A,B,C) = ",astr(theo_si_sti_disjoint['wa_A-wb_B-wc_C_process'][0],prec=5))
                print("      sti numeric  (A,B,C) = ",astr(sobol_indexes_disjoint_tmp[str(jset)]['processes']['sti']['out'],prec=5))
                print("      sti analytic (A,B,C) = ",astr(theo_si_sti_disjoint['wa_A-wb_B-wc_C_process'][1],prec=5))
                print("      mae error    (A,B,C) = ",np.mean(np.array(
                    list(np.abs(sobol_indexes_disjoint_tmp[str(jset)]['processes']['si']['out']-theo_si_sti_disjoint['wa_A-wb_B-wc_C_process'][0]))+
                    list(np.abs(sobol_indexes_disjoint_tmp[str(jset)]['processes']['sti']['out']-theo_si_sti_disjoint['wa_A-wb_B-wc_C_process'][1])))))
            
            sobol_indexes_disjoint[str(iset)] = sobol_indexes_disjoint_tmp
            
        save_to_pickle['numerical_si_sti_disjoint'] = sobol_indexes_disjoint

        # save all to pickle
        filename = '../data_out/benchmark/results_disjoint-benchmark-model_DVM.pkl'
        pickle.dump( save_to_pickle, open( filename, "wb" ) )
        print("Wrote: ",filename)


    # ----------------------------------------------------
    # Sensitivity Analysis of 'shared-parameter setup' (DVM)
    # ----------------------------------------------------
    if do_shared_example_DVM:

        save_to_pickle = {}
        save_to_pickle["theo_si_sti_shared"] = theo_si_sti_shared

        sobol_indexes_shared = {}
        for iset in nsets:

            
            print('')
            print('SA of shared-parameter setup: nsets = ',iset)

            sobol_indexes_shared_tmp = {}
            for jset in nsets_parasets:

                print('')
                print('   SA of shared-parameter setup: nsets_parasets = ',jset)
        
                sobol_indexes_shared_tmp[str(jset)] = sa_model_multiple_processes_DVM.sa_model_multiple_processes_DVM(options_paras_shared,
                                                                                                         para_ranges,
                                                                                                         model_function_shared,
                                                                                                         constants=None,
                                                                                                         nsets=iset,
                                                                                                         nsets_parasets=jset)
                # process sensitivities
                print("      si  numeric  (A,B,C) = ",astr(sobol_indexes_shared_tmp[str(jset)]['processes']['si']['out'],prec=5))
                print("      si  analytic (A,B,C) = ",astr(theo_si_sti_shared['wa_A-wb_B-wc_C_process'][0],prec=5))
                print("      sti numeric  (A,B,C) = ",astr(sobol_indexes_shared_tmp[str(jset)]['processes']['sti']['out'],prec=5))
                print("      sti analytic (A,B,C) = ",astr(theo_si_sti_shared['wa_A-wb_B-wc_C_process'][1],prec=5))
                print("      mae error    (A,B,C) = ",np.mean(np.array(
                    list(np.abs(sobol_indexes_shared_tmp[str(jset)]['processes']['si']['out']-theo_si_sti_shared['wa_A-wb_B-wc_C_process'][0]))+
                    list(np.abs(sobol_indexes_shared_tmp[str(jset)]['processes']['sti']['out']-theo_si_sti_shared['wa_A-wb_B-wc_C_process'][1])))))
            
            sobol_indexes_shared[str(iset)] = sobol_indexes_shared_tmp
            
        save_to_pickle['numerical_si_sti_shared'] = sobol_indexes_shared

        # save all to pickle
        filename = '../data_out/benchmark/results_shared-benchmark-model_DVM.pkl'
        pickle.dump( save_to_pickle, open( filename, "wb" ) )
        print("Wrote: ",filename)


