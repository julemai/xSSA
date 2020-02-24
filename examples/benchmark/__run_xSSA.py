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
#     run __run_xSSA.py -n 1000

#!/usr/bin/env python
from __future__ import print_function

"""
Benchmark example to test Sensitivity Analysis for models with multiple process options. 
As an option, a file name can be specified to save all model outputs in a pickle file (-o). 
As another option, the user can specify the number of reference parameter sets for 
the Sobol' analysis (-n).

History
-------
Written,  JM, Jun 2019
"""

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../../lib')


# Comment|Uncomment - Begin
if __name__ == '__main__':

    import argparse
    import numpy as np
    import pickle
    import copy
    import collections
    
    import sobol                                      # in lib/  
    from   sobol_index     import sobol_index         # in lib/    
    from   fsread          import fsread              # in lib/
    from   autostring      import astr                # in lib/
    import sa_model_multiple_processes                # in lib/

    do_weighted_model    = True     # SA of weighted benchmark model 
    do_individual_models = True     # SA of every single model structure of benchamrk model

    nsets    = 1000          # number of reference parameter sets
    outfile  = None          # output file name stroing model runs, e.g., 'results_realistic-benchmark-model.pkl'
    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''Benchmark example to test Sensitivity Analysis for models with multiple process options.''')
    parser.add_argument('-o', '--outfile', action='store',
                        default=outfile, dest='outfile', metavar='outfile',
                        help='Optional output pickle file to store all model runs (default: None).')
    parser.add_argument('-n', '--nsets', action='store',
                        default=nsets, dest='nsets', metavar='nsets',
                        help='Optional number of reference parameter sets (default: 1000).')

    args     = parser.parse_args()
    outfile  = args.outfile
    nsets    = np.int(args.nsets)

    del parser, args
# Comment|Uncomment - End


# -------------------------------------------------------------------------
# Function definition - if function
#

    # Theoretical Sobol' sensitivity indexes (calculated with Mathematica based on integrals)
    theo_si_sti_realistic = {   # individual models with all possible process option combinations
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
                      # all process options weighted to one model: 7 parameters (x1,...,x7) + 4 weights
                      'wa_A-wb_B-wc_C_para': np.array([[0.02298,0.38064,0.00222,0.00023,0.00023,0.00003,0.03928,0.05163,0.05770,0.00043,0.01006],
                                                       [0.07533,0.77089,0.00450,0.00103,0.00103,0.00004,0.05237,0.26086,0.31784,0.00264,0.02333]]),
                      # --> mae error = 0.0129055976532     (N=    1,000)
                      # --> mae error = 0.00545848871962    (N=   10,000)    
                      # --> mae error = 0.000809206916118   (N=  100,000)
                      # --> mae error = 0.000532741670372   (N=1,000,000)                      
                      #
                      # all process options weighted to one model: 7 process options (a1,a2,b1,b2,b3,c1,c2) + 4 weights
                      'wa_A-wb_B-wc_C_option': np.array([[0.02298,0.42889,0.38064,0.00222,0.00046,0.00003,0.04149,0.05163,0.05770,0.00043,0.01006],
                                                         [0.07533,0.80441,0.77089,0.00450,0.00207,0.00004,0.05687,0.26086,0.31784,0.00264,0.02333]]),
                      # --> mae error = 0.0156002867833    (N=    1,000)
                      # --> mae error = 0.00949079962253   (N=   10,000)  
                      # --> mae error = 0.0011354491192    (N=  100,000)
                      # --> mae error = 0.000779724362322   (N=1,000,000)
                      
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
    # read parameter ranges from file
    paras = collections.OrderedDict()
    para_ranges = []
    for line in open('parameters.txt'):
        li=line.strip().split('#')[0].strip()
        if ( (li != '') and not(li.startswith('BeginParams')) and not(li.startswith('EndParams')) ):
            # print(li)
            para_name = li.split()[0]
            para_ini  = li.split()[1]
            para_min  = li.split()[2]
            para_max  = li.split()[3]

            paras[para_name]  = [para_ini,para_min,para_max]
            para_ranges      += [ [np.float(para_min),np.float(para_max)] ]

    # list of parameters that go into each option (numbering starts with 0)
    # e.g., options_paras_realistic = [ 
    #                                    [[0], [0,1]],          # parameters of process options A1 and A2
    #                                    [[1], [2], [3,4]],     # parameters of process options B1, B2, and B3
    #                                    [[5], [2,6]]           # parameters of process options C1 and C2
    #                                 ]
    options_paras_realistic = []
    for line in open('parameter-process-mapping.txt'):
        li=line.strip().split('#')[0].strip()
        if ( (li != '') and not(li.startswith('BeginParams')) and not(li.startswith('EndParams')) ):

            # [[par_x01], [par_x01, par_x02]]   -->   [['par_x01'], ['par_x01', 'par_x02']]
            li = li.replace('[','').replace(' ','').split('],')
            li = [ s.replace(']','').split(',') for s in li ]

            # find correct ID: [['par_x01'], ['par_x01', 'par_x02']] --> [[0],[0,1]]
            idx = [ [ paras.keys().index(item) for item in ilist ] for ilist in li ]
            options_paras_realistic += [idx]

    def model_function_realistic(paras, weights, constants=None, run_id=None):
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
        # realistic model
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
    # Sensitivity Analysis of 'realistic setup'
    # ----------------------------------------------------
    if do_weighted_model:

        save_to_pickle = {}
        save_to_pickle["theo_si_sti_realistic"] = theo_si_sti_realistic

        sobol_indexes_realistic = {}

        print('')
        print('SA of realistic setup: nsets = ',nsets)
        
        sobol_indexes_realistic[str(nsets)] = sa_model_multiple_processes.sa_model_multiple_processes(options_paras_realistic,
                                                                                                     para_ranges,
                                                                                                     model_function_realistic,
                                                                                                     constants=None,
                                                                                                     nsets=nsets)

        # parameter sensitivities
        print("")
        print("------------------------------")
        print("Sensitivities of parameters")
        print("------------------------------")
        print("   si  numeric  (x1,x2,x3,x4,x5,x6,x7,w1,..,w4) = ",astr(sobol_indexes_realistic[str(nsets)]['paras']['si']['out'],prec=5))
        print("   si  analytic (x1,x2,x3,x4,x5,x6,x7,w1,..,w4) = ",astr(theo_si_sti_realistic['wa_A-wb_B-wc_C_para'][0],prec=5))
        print("   sti numeric  (x1,x2,x3,x4,x5,x6,x7,w1,..,w4) = ",astr(sobol_indexes_realistic[str(nsets)]['paras']['sti']['out'],prec=5))
        print("   sti analytic (x1,x2,x3,x4,x5,x6,x7,w1,..,w4) = ",astr(theo_si_sti_realistic['wa_A-wb_B-wc_C_para'][1],prec=5))
        print("   mae error    (x1,x2,x3,x4,x5,x6,x7,w1,..,w4) = ",np.mean(np.array(
            list(np.abs(sobol_indexes_realistic[str(nsets)]['paras']['si']['out']-theo_si_sti_realistic['wa_A-wb_B-wc_C_para'][0]))+
            list(np.abs(sobol_indexes_realistic[str(nsets)]['paras']['sti']['out']-theo_si_sti_realistic['wa_A-wb_B-wc_C_para'][1])))))

        # process option sensitivities
        print("")
        print("------------------------------")
        print("Sensitivities of process options")
        print("------------------------------")
        print("   si  numeric  (A1,A2,B1,B2,B3,C1,C2,w1,..,w4) = ",astr(sobol_indexes_realistic[str(nsets)]['process_options']['si']['out'],prec=5))
        print("   si  analytic (A1,A2,B1,B2,B3,C1,C2,w1,..,w4) = ",astr(theo_si_sti_realistic['wa_A-wb_B-wc_C_option'][0],prec=5))
        print("   sti numeric  (A1,A2,B1,B2,B3,C1,C2,w1,..,w4) = ",astr(sobol_indexes_realistic[str(nsets)]['process_options']['sti']['out'],prec=5))
        print("   sti analytic (A1,A2,B1,B2,B3,C1,C2,w1,..,w4) = ",astr(theo_si_sti_realistic['wa_A-wb_B-wc_C_option'][1],prec=5))
        print("   mae error    (A1,A2,B1,B2,B3,C1,C2,w1,..,w4) = ",np.mean(np.array(
            list(np.abs(sobol_indexes_realistic[str(nsets)]['process_options']['si']['out']-theo_si_sti_realistic['wa_A-wb_B-wc_C_option'][0]))+
            list(np.abs(sobol_indexes_realistic[str(nsets)]['process_options']['sti']['out']-theo_si_sti_realistic['wa_A-wb_B-wc_C_option'][1])))))

        # process sensitivities
        print("")
        print("------------------------------")
        print("Sensitivities of processes")
        print("------------------------------")
        print("   si  numeric  (A,B,C) = ",astr(sobol_indexes_realistic[str(nsets)]['processes']['si']['out'],prec=5))
        print("   si  analytic (A,B,C) = ",astr(theo_si_sti_realistic['wa_A-wb_B-wc_C_process'][0],prec=5))
        print("   sti numeric  (A,B,C) = ",astr(sobol_indexes_realistic[str(nsets)]['processes']['sti']['out'],prec=5))
        print("   sti analytic (A,B,C) = ",astr(theo_si_sti_realistic['wa_A-wb_B-wc_C_process'][1],prec=5))
        print("   mae error    (A,B,C) = ",np.mean(np.array(
            list(np.abs(sobol_indexes_realistic[str(nsets)]['processes']['si']['out']-theo_si_sti_realistic['wa_A-wb_B-wc_C_process'][0]))+
            list(np.abs(sobol_indexes_realistic[str(nsets)]['processes']['sti']['out']-theo_si_sti_realistic['wa_A-wb_B-wc_C_process'][1])))))

        # process sensitivities
        noptions = np.sum(np.array([ len(oo)   for oo in options_paras_realistic ]))
        nweights = np.sum(np.array([ len(oo)-1 for oo in options_paras_realistic ]))
        nprocess = len(options_paras_realistic)
        nparas   = np.shape(para_ranges)[0]
        print("")
        print("------------------------------")
        print("Budget")
        print("------------------------------")
        print("    Parameters:      ",(nparas  +nweights+2)*nsets)
        print("    Process options: ",(noptions+nweights+2)*nsets)
        print("    Processes:       ",(nprocess         +2)*nsets)
        print("    ---------------------------")
        print("    TOTAL:           ",((nparas  +nweights+2)+(noptions+nweights+2)+(nprocess         +2))*nsets)
        print("")

        save_to_pickle['numerical_si_sti_realistic'] = sobol_indexes_realistic

        # save all to pickle
        if not( outfile is None ):
            filename = outfile
            pickle.dump( save_to_pickle, open( filename, "wb" ) )


    
    # --------------------
    # Global parameters
    # --------------------
    para_a = 2.0
    para_b = 1.0

    if do_individual_models: 
        print('--------------------------------------------')
        print('SA FOR INDIVIDUAL MODELS')
        print('--------------------------------------------')
    
        noptions = np.array([ len(oo) for oo in options_paras_realistic ])
        nprocess = len(options_paras_realistic)
        nparas   = np.shape(para_ranges)[0]
        zero_weight = [ list(np.zeros(noptions[iprocess])) for iprocess in range(nprocess) ]
        para_ranges = np.array(para_ranges)

        model_runs = 0
        for ioption_a in range(len(options_paras_realistic[0])):
            for ioption_b in range(len(options_paras_realistic[1])):
                for ioption_c in range(len(options_paras_realistic[2])):

                    model_str="A"+str(ioption_a+1)+"-B"+str(ioption_b+1)+"-C"+str(ioption_c+1)
                    print("Model: ",model_str)
                    print("SA:    each parameter")

                    paras_active_a = options_paras_realistic[0][ioption_a]
                    paras_active_b = options_paras_realistic[1][ioption_b]
                    paras_active_c = options_paras_realistic[2][ioption_c]
                    nnparas = len(np.unique(paras_active_a+paras_active_b+paras_active_c))

                    # --------------------
                    # Weights for current set of options to 1.0; others are 0.0
                    # --------------------
                    block_weights_nested = copy.deepcopy(zero_weight)
                    block_weights_nested[0][ioption_a] = 1.0
                    block_weights_nested[1][ioption_b] = 1.0
                    block_weights_nested[2][ioption_c] = 1.0

                    # --------------------
                    # get sobol sequences
                    # --------------------
                    sobol_sets = sobol.i4_sobol_generate(nparas*2,nsets,40000)
                    sobol_sets = np.transpose(sobol_sets)

                    # --------------------
                    # scale parameter sets A
                    # --------------------
                    block_a_paras  = copy.deepcopy(sobol_sets[:,0:nparas])
                    block_a_paras *= (para_ranges[:,1]-para_ranges[:,0])
                    block_a_paras += para_ranges[:,0]

                    # --------------------
                    # Run model for A-sets
                    # --------------------
                    f_a = np.array([ model_function_realistic(block_a_paras[iset], block_weights_nested, constants=[para_a,para_b])['out'] for iset in range(nsets) ])
                    model_runs += nsets

                    # --------------------
                    # scale parameter sets B
                    # --------------------
                    block_b_paras  = copy.deepcopy(sobol_sets[:,nparas:nparas*2])
                    block_b_paras *= (para_ranges[:,1]-para_ranges[:,0])
                    block_b_paras += para_ranges[:,0]

                    # --------------------
                    # Run model for B-sets
                    # --------------------
                    f_b = np.array([ model_function_realistic(block_b_paras[iset], block_weights_nested, constants=[para_a,para_b])['out'] for iset in range(nsets) ])
                    model_runs += nsets

                    # --------------------
                    # list of parameters active in this set of options
                    # --------------------
                    col_changes_Ci = np.unique(options_paras_realistic[0][ioption_a]+options_paras_realistic[1][ioption_b]+options_paras_realistic[2][ioption_c])
                    nnparas = len(col_changes_Ci)
                    
                    f_c = np.ones([nnparas,nsets])*-9999.0
                    for iipara,ipara in enumerate(col_changes_Ci):

                        # --------------------
                        # Create Ci-sets
                        # --------------------
                        block_c_paras          = copy.deepcopy(block_a_paras)
                        block_c_paras[:,ipara] = copy.deepcopy(block_b_paras[:,ipara])

                        # --------------------
                        # Run model for Ci-sets
                        # --------------------
                        f_c[iipara,:] = np.array([ model_function_realistic(block_c_paras[iset], block_weights_nested, constants=[para_a,para_b])['out'] for iset in range(nsets) ])
                        model_runs += nsets

                    # --------------------
                    # Calculate Sobol' index
                    # --------------------
                    si, sti = sobol_index(ya=f_a, yb=f_b, yc=f_c, si=True, sti=True, method='Mai1999')

                    # theoretical values
                    paras_str = ','.join([ 'x'+astr(ii) for ii in np.where(~np.isnan(theo_si_sti_realistic[model_str][0]))[0]+1 ])
                    si_theo  = theo_si_sti_realistic[model_str][0]
                    si_theo  = si_theo[~np.isnan(si_theo)]    # remove nan's
                    sti_theo = theo_si_sti_realistic[model_str][1]
                    sti_theo = sti_theo[~np.isnan(sti_theo)]  # remove nan's
                    print("   si  numeric  ("+paras_str+") = ",astr(si,prec=5))
                    print("   si  analytic ("+paras_str+") = ",astr(si_theo,prec=5))
                    print("   sti numeric  ("+paras_str+") = ",astr(sti,prec=5))
                    print("   sti analytic ("+paras_str+") = ",astr(sti_theo,prec=5))
                    print("   mae error    ("+paras_str+") = ",np.mean(np.array(list(np.abs(si_theo-si))+list(np.abs(sti_theo-sti)))))

        print("")
        print("------------------------------")
        print("Budget")
        print("------------------------------")
        print("    Parameters:      ",model_runs)

