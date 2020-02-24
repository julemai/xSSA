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

# You should have received a copy of the GNU Lesser General Public Licensefstop
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#
# run with:
#     # n=1000 takes about 10h
#     run __run_xSSA.py -n 2 -o output.pkl -i '08KC001'

#!/usr/bin/env python
from __future__ import print_function

"""
Raven example setup for Salmon River (BC, Canada) to derive sensitivity indexes for a model
with multiple process options. As an option, a file name can be specified to save all model 
outputs in a pickle file (-o). As another option, the user can specify the number of reference 
parameter sets for the Sobol' analysis (-n).

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
    import re
    import subprocess
    import shutil
    
    import sobol                                      # in lib/  
    from   sobol_index     import sobol_index         # in lib/
    from   fread           import fread               # in lib/
    from   fsread          import fsread              # in lib/
    from   autostring      import astr                # in lib/
    import sa_model_multiple_processes                # in lib/

    from   raven_templates import RVI, RVT, RVP, RVH, RVC
    from   raven_common    import writeString, makeDirectories
    from   pathlib2        import Path

    nsets       = 1000                # number of reference parameter sets
    outfile     = None                # output file name stroing model runs, e.g., 'results_raven-salmon.pkl'
    tmp_folder  = "/tmp/xSSA-test/"   # temporary folder to run model
    basin_id    = None                # Basin ID (for folder creation)
    
    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''Raven example setup for Salmon River (BC, Canada) to derive sensitivity indexes for a model with multiple process options.''')
    parser.add_argument('-o', '--outfile', action='store',
                        default=outfile, dest='outfile', metavar='outfile',
                        help='Optional output pickle file to store all model runs (default: None).')
    parser.add_argument('-n', '--nsets', action='store',
                        default=nsets, dest='nsets', metavar='nsets',
                        help='Optional number of reference parameter sets (default: 1000).')
    parser.add_argument('-t', '--tmp_folder', action='store',
                        default=tmp_folder, dest='tmp_folder', metavar='tmp_folder',
                        help='Temporary directory to run the model. (default: "/tmp/juletest/").')
    parser.add_argument('-i', '--basin_id', action='store',
                        default=basin_id, dest='basin_id', metavar='basin_id',
                        help='Basin ID of basins to analyze. Mandatory. (default: None).')

    args       = parser.parse_args()
    outfile    = args.outfile
    nsets      = np.int(args.nsets)
    tmp_folder = args.tmp_folder
    basin_id   = args.basin_id
    print('tmp_folder: ',tmp_folder)

    # Basin ID need to be specified
    if basin_id is None:
        raise ValueError('Basin ID (option -i) needs to given. ')

    del parser, args
# Comment|Uncomment - End


# -------------------------------------------------------------------------
# Function definition - if function

    # -------------------------------------------------------------------------
    # Basin properties (read from file)
    # -------------------------------------------------------------------------
    basin_prop = {}

    file_gauge_info = 'basin_physical_characteristics.txt'
    ff = open(file_gauge_info, "r")
    lines = ff.readlines()
    ff.close()

    found = False
    for ill,ll in enumerate(lines):
        if ill > 0:
            tmp = ll.strip().split(';')
            if (tmp[0] == basin_id):
                found = True
                # basin_id; basin_name; lat; lon; area_km2; elevation_m; slope_deg; forest_frac 
                
                basin_prop['id']           = str(basin_id)
                basin_prop['name']         = str(tmp[1].strip())
                basin_prop['lat_deg']      = np.float(tmp[2].strip())
                basin_prop['lon_deg']      = np.float(tmp[3].strip())
                basin_prop['area_km2']     = np.float(tmp[4].strip())
                basin_prop['elevation_m']  = np.float(tmp[5].strip())
                basin_prop['slope_deg']    = np.float(tmp[6].strip()) 
                basin_prop['forest_frac']  = np.float(tmp[7].strip())

    if not(found):
        raise ValueError('Basin ID not found in '+file_gauge_info)

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
    options_paras_raven = []
    for line in open('parameter-process-mapping.txt'):
        li=line.strip().split('#')[0].strip()
        if ( (li != '') and not(li.startswith('BeginParams')) and not(li.startswith('EndParams')) ):

            # [[par_x01], [par_x01, par_x02]]   -->   [['par_x01'], ['par_x01', 'par_x02']]
            li = ''.join(li.split())   # remove whitespaces
            li = li.replace('[','').replace(' ','').split('],')
            li = [ s.replace(']','').split(',') for s in li ]

            # find correct ID: [['par_x01'], ['par_x01', 'par_x02']] --> [[0],[0,1]]
            idx = [ [ paras.keys().index(item) for item in ilist if (item != '') ] for ilist in li ]
            options_paras_raven += [idx]

    def model_function_raven(paras, weights, basin_prop, constants=None, run_id=None, tmp_folder=tmp_folder):
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
        #     basin_prop ... dictionary of basin properties like elevation, latitude, longitude etc
        #                    example:
        #                         {'area_km2':     2303.95,
        #                          'elevation_m':  250.31,
        #                          'forest_frac':  0.9063,
        #                          'id':           '01013500',
        #                          'lat_deg':      47.23739,
        #                          'lon_deg':      -68.58264,
        #                          'name':         'Fish River near Fort Kent, Maine',
        #                          'slope_m_km-1': 21.64152}
        #     constants ... optional list of constants that are same for all models;
        #                   like parameters a and b in Ishigami-Homma function
        #                   example:
        #                        [2.0, 1.0]
        #     run_id    ... optional name of this run (to, e.g., print or store in a file)
        #                   example:
        #                        run_aset_001
        # output:
        #     (scalar) model output
        #     example:
        #          7.4

        if not(run_id is None):
            print("Run ID: ",run_id)

        # check that provided number of weights is correct:
        # --> one weight per option per process
        if ( [len(ilist) for ilist in weights] != [3,3,2,2,3,1,1,1,1] ):       # HMETS, SIMPLE_MELT, HBV
            print("Number of weights: ",[len(ilist) for ilist in weights])
            raise ValueError("sa_model_multiple_processes: model_function: provided number of weights must be [3,3,2,2,3,1,1,1,1]")
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
        if (len(paras) != 30):
            print("Number of parameters: ",len(paras))
            raise ValueError("sa_model_multiple_processes: model_function: provided number of parameters must be 30")

        # ---------------
        # derive some parameters
        # ---------------
        dict_dparas = {}
        dict_dparas['sum_x24_x25']  = paras[23]+paras[24]         # MAX_MELT_FACTOR > MIN_MELT_FACTOR
        dict_dparas['sum_x13_x14']  = paras[12]+paras[13]         # SNOW_SWI_MAX > SNOW_SWI_MIN
        dict_dparas['half_x29']     = paras[28] * 0.5 * 1000      # half the value but in [mm] not [m]
        dict_dparas['half_x30']     = paras[29] * 0.5 * 1000      # half the value but in [mm] not [m]
        dict_dparas['sum_x09_x10']  = paras[8]+paras[9]           # FIELD_CAPACITY > SAT_WILT
        dict_dparas['pow_x04']      = 10.0**(paras[3])            # BASEFLOW_COEFF TOPSOIL  = 10.0^x4
        dict_dparas['pow_x11']      = 10.0**(paras[10])           # BASEFLOW_COEFF PHREATIC = 10.0^x11
        
        # ---------------
        # paste all paras and weights into template files
        # ---------------
        #   ex.:    string = "parameter x01 = {par[x01]} and another parameter x02 = {par[x02]}"
        #           keys   = ['x01','x02']
        #           vals   = [1.0,3.0]
        #           string.format(par=dict(zip(keys,vals)))
        #
        #           --> 'parameter x01 = 1.0 and another parameter x02 = 3.0'
        #
        # to replace patterns: {par[x01]} by parameter value paras[0]
        #                      {par[x02]} by parameter value paras[1]
        #                      ...
        if len(paras) > 9 and len(paras) < 100:
            keys_paras = ["x{:02d}".format(ii)   for ii in range(1,len(paras)+1) ]
        elif len(paras) > 99 and len(paras) < 1000:
            keys_paras = ["x{:03d}"   for ii in range(1,len(paras)+1) ]
        elif len(paras) <= 9:
            keys_paras = ["x"+str(ii) for ii in range(1,len(paras)+1) ]
        else:
            raise ValueError("More than 999 parameters are not implemented yet!")
        vals_paras = paras
        dict_paras = dict(zip(keys_paras,vals_paras))

        # to replace patterns: {weights[process01]} by all weights of that process weights[0][0] weights[0][1] ...
        #                      {weights[process02]} by all weights of that process weights[1][0] weights[1][1] ...
        #                      ...
        if len(weights) > 9 and len(weights) < 100:
            keys_weights = ["process{:02d}".format(ii)   for ii in range(1,len(weights)+1) ]
        elif len(weights) > 99 and len(weights) < 1000:
            keys_weights = ["process{:03d}"   for ii in range(1,len(weights)+1) ]
        elif len(weights) <= 9:
            keys_weights = ["process"+str(ii) for ii in range(1,len(weights)+1) ]
        else:
            raise ValueError("More than 999 processes are not implemented yet!")
        vals_weights = [ ' '.join(map(str, iweight)) for iweight in weights ]    # [[0.2, 0.3, 0.4], [0.1, 0.9]] --> ['0.2 0.3 0.4', '0.1 0.9']
        dict_weights = dict(zip(keys_weights,vals_weights))

        # fill in to templates
        # templates need to have patterns:
        #         {par[x01]},  {par[x02]},                     ... for parameters
        #         {dpar[something]},  {dpar[somethingelse]},   ... for derived parameters
        #         {weights[process1]}, {weights[process2]},    ... for weights

        # ---------------
        # create a run folder
        # ---------------
        tmp_folder = tmp_folder+"/"+str(basin_id)+"/"+str(run_id) # "/tmp/juletest" #  TODO a generic folder name in /tmp
        raven_exe_name = "./model/Raven.exe"
        raven_obs_folder = "./model/data_obs"

        if os.path.exists(tmp_folder):
            shutil.rmtree(tmp_folder)

        # print("dict_paras   = ",dict_paras)
        # print("dict_dparas  = ",dict_dparas)
        # print("dict_weights = ",dict_weights)

        # all RAVEN setup files
        writeString( Path(tmp_folder,"raven_weighted_processes.rvi"), RVI.format(par=dict_paras,dpar=dict_dparas,weights=dict_weights,props=basin_prop) )
        writeString( Path(tmp_folder,"raven_weighted_processes.rvp"), RVP.format(par=dict_paras,dpar=dict_dparas,weights=dict_weights,props=basin_prop) )
        writeString( Path(tmp_folder,"raven_weighted_processes.rvh"), RVH.format(par=dict_paras,dpar=dict_dparas,weights=dict_weights,props=basin_prop) )
        writeString( Path(tmp_folder,"raven_weighted_processes.rvt"), RVT.format(par=dict_paras,dpar=dict_dparas,weights=dict_weights,props=basin_prop) )
        writeString( Path(tmp_folder,"raven_weighted_processes.rvc"), RVC.format(par=dict_paras,dpar=dict_dparas,weights=dict_weights,props=basin_prop) )        

        # link executable
        # print("path = ",str(Path(tmp_folder,os.path.basename(raven_exe_name))))
        # print("exists: ",os.path.exists(str(Path(tmp_folder,os.path.basename(raven_exe_name)))))
        if not(os.path.exists(str(Path(tmp_folder,os.path.basename(raven_exe_name))))):
            os.symlink(os.path.realpath(raven_exe_name), str(Path(tmp_folder,os.path.basename(raven_exe_name))))

        # link observations folder
        if not(os.path.exists(str(Path(tmp_folder,os.path.basename(raven_obs_folder))))):
            os.symlink(os.path.realpath(raven_obs_folder), str(Path(tmp_folder,os.path.basename(raven_obs_folder))))

        # create ouput folder
        # makeDirectories(Path(tmp_folder,"output"))
        out_folder = str(Path(tmp_folder,"output"))
        os.makedirs(out_folder)
        
        # ---------------
        # run the model with these input rv* files
        # ---------------        
        cmd = [str(Path(tmp_folder,os.path.basename(raven_exe_name))),str(Path(tmp_folder,"raven_weighted_processes")),"-o",str(Path(tmp_folder,"output"))+'/']

        #print("Raven run folder: ",tmp_folder)
        #print("raw cmd:          ",cmd)
        print("run cmd: ",' '.join(cmd))

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        print("")
        print("Raven standard output:")
        for line in process.stdout:
            print(">>> ",line.rstrip()) # rstrip removes trailing \n

        if not(os.path.exists(str(Path(tmp_folder,"output","Diagnostics.csv")))):            
            print("")
            print("ERROR: No Diagnostics.csv produced")
            print("")
            print("Raven error file content:")
            ff = open(str(Path(tmp_folder,"output","Raven_errors.txt")), "r")
            lines = ff.readlines()
            ff.close()
            for line in lines:
                print(">>> ",line.rstrip()) # rstrip removes trailing \n

            raise ValueError("ERROR: No Diagnostics.csv produced (scroll up to see content of error file)")
        # stop

        model = {}

        # ---------------
        # extract model output: Diagnostics: NSE
        # ---------------
        model['nse'] = 0.0
        ff = open(str(Path(tmp_folder,"output","Diagnostics.csv")), "r")
        lines = ff.readlines()
        ff.close()

        irow = 1  # starting with 0
        icol = 2  # starting with 0 and assuming "," as delimiter

        # print("lines: ",lines)
        model['nse'] = np.float((lines[irow]).split(',')[icol])
        print("NSE:     ",model['nse'])
        print("")
        
        # ---------------
        # extract model output: Hydrographs: simulated Q
        # ---------------
        model['Q']  = 0.0
        warmup = 2*365  # 1 # model timestep 1 day and want to skip 2 years  # first day 1991-01-01 00:00:00.00 (checked)
        model['Q']  = np.transpose(fread(str(Path(tmp_folder,"output","Hydrographs.csv")),skip=warmup+1,cskip=4,nc=1))[0]

        print("Q:              ",model['Q'][0:4],"...",model['Q'][-4:])
        print("Q_range:         [",np.min(model['Q']),",",np.max(model['Q']),"]")
        print("shape Q:        ",np.shape(model['Q']))
        print("")

        # # ---------------
        # # extract model output: BETWEEN_PONDED_WATER_AND_SOIL[0]_Daily_Average_BySubbasin.csv: accumulated infiltration volume
        # # ---------------
        # model['infiltration']  = 0.0
        # warmup = 2*365  # 1 # model timestep 1 day and want to skip 2 years  # first day 1990-12-31 00:00:00.00 (checked) But all timesteps are shifted by 1 day...
        # #
        # # accumulated infiltration volume
        # # model = np.transpose(fread(str(Path(tmp_folder,"output","BETWEEN_PONDED_WATER_AND_SOIL[0]_Daily_Average_BySubbasin.csv")),skip=warmup+1,cskip=2,nc=1))[0]
        # #
        # # de-accumulated infiltration volume
        # model['infiltration'] = np.transpose(fread(str(Path(tmp_folder,"output","BETWEEN_PONDED_WATER_AND_SOIL[0]_Daily_Average_BySubbasin.csv")),skip=warmup,cskip=2,nc=1))[0]
        # model['infiltration'] = np.diff(model['infiltration'])

        # print("Infiltration I: ",model['infiltration'][0:4],"...",model['infiltration'][-4:])
        # print("I_range:         [",np.min(model['infiltration']),",",np.max(model['infiltration']),"]")
        # print("shape I:        ",np.shape(model['infiltration']))
        # print("")

        # # ---------------
        # # extract model output: FROM_SNOW_Daily_Average_BySubbasin.csv: accumulated snow melt volume
        # # ---------------
        # model['snowmelt']  = 0.0
        # warmup = 2*365  # 1 # model timestep 1 day and want to skip 2 years  # first day 1990-12-31 00:00:00.00 (checked) But all timesteps are shifted by 1 day...
        # #
        # # de-accumulated snow melt volume
        # model['snowmelt'] = np.transpose(fread(str(Path(tmp_folder,"output","FROM_SNOW_Daily_Average_BySubbasin.csv")),skip=warmup,cskip=2,nc=1))[0]
        # model['snowmelt'] = np.diff(model['snowmelt'])

        # print("Snowmelt SM:    ",model['snowmelt'][0:4],"...",model['snowmelt'][-4:])
        # print("SM_range:        [",np.min(model['snowmelt']),",",np.max(model['snowmelt']),"]")
        # print("shape SM:       ",np.shape(model['snowmelt']))
        # print("")
        
        # # ---------------
        # # extract model output: SNOW_Daily_Average_BySubbasin.csv: Actual snow on ground over whole basin (not accumulated)
        # # ---------------
        # model['snowdepth']  = 0.0
        # warmup = 2*365  # 1 # model timestep 1 day and want to skip 2 years  # first day 1990-12-31 00:00:00.00 (checked) But all timesteps are shifted by 1 day...
        # #
        # # read data (no de-accumulation)
        # model['snowdepth'] = np.transpose(fread(str(Path(tmp_folder,"output","SNOW_Daily_Average_BySubbasin.csv")),skip=warmup+1,cskip=2,nc=1))[0]

        # print("Snowdepth SD:   ",model['snowdepth'][0:4],"...",model['snowdepth'][-4:])
        # print("SD_range:        [",np.min(model['snowdepth']),",",np.max(model['snowdepth']),"]")
        # print("shape SD:       ",np.shape(model['snowdepth']))
        # print("")

        return model


    
    
    # ----------------------------------------------------
    # Sensitivity Analysis of 'raven setup'
    # ----------------------------------------------------
    if True:
        sobol_indexes_raven = sa_model_multiple_processes.sa_model_multiple_processes(options_paras_raven,
                                                                                      para_ranges,
                                                                                      model_function_raven,
                                                                                      constants=None,
                                                                                      basin_prop=basin_prop,
                                                                                      nsets=nsets,
                                                                                      save_pkl=outfile)

        keys = sobol_indexes_raven['paras']['si'].keys()

        print("")
        for ikey in keys:

            print("")
            print("---------------------")
            print(ikey)
            print("---------------------")
            
            # scalar model output (such as NSE)
            if (len(np.shape(sobol_indexes_raven['paras']['si'][ikey])) == 1):
                # parameter sensitivities
                print("   si  numeric  (x1,x2,x3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['paras']['si'][ikey] ])
                print("   sti numeric  (x1,x2,x3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['paras']['sti'][ikey] ])
                
                # process option sensitivities
                print("   si  numeric  (a1,a2,a3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['process_options']['si'][ikey] ])
                print("   sti numeric  (a1,a2,a3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['process_options']['sti'][ikey] ])
                
                # process sensitivities
                print("   si  numeric  (A,B,C,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['processes']['si'][ikey] ])
                print("   sti numeric  (A,B,C,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['processes']['sti'][ikey] ])
                
            # 1D model output (such as discharge time series)
            if (len(np.shape(sobol_indexes_raven['paras']['si'][ikey])) == 2):
                # parameter sensitivities
                print("   mean si  numeric  (x1,x2,x3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['paras']['msi'][ikey] ])
                print("   mean sti numeric  (x1,x2,x3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['paras']['msti'][ikey] ])
                
                # process option sensitivities
                print("   mean si  numeric  (a1,a2,a3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['process_options']['msi'][ikey] ])
                print("   mean sti numeric  (a1,a2,a3,...,w1,w2,w3,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['process_options']['msti'][ikey] ])
                
                # process sensitivities
                print("   mean si  numeric  (A,B,C,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['processes']['msi'][ikey] ])
                print("   mean sti numeric  (A,B,C,...) = ",[ astr(pp,prec=5) if ~np.isnan(pp) else "nan" for pp in sobol_indexes_raven['processes']['msti'][ikey] ])


    if False:
        
        print('--------------------------------------------')
        print('SA FOR INDIVIDUAL MODELS')
        print('--------------------------------------------')
    
        noptions = np.array([ len(oo) for oo in options_paras_raven ])
        nprocess = len(options_paras_raven)

        model_runs = 0

        for ioption_m in range(len(options_paras_raven[0])):
            for ioption_n in range(len(options_paras_raven[1])):
                for ioption_o in range(len(options_paras_raven[2])):
                    for ioption_p in range(len(options_paras_raven[3])):
                        for ioption_q in range(len(options_paras_raven[4])):
                            for ioption_r in range(len(options_paras_raven[5])):
                                for ioption_s in range(len(options_paras_raven[6])):
                                    for ioption_t in range(len(options_paras_raven[7])):
                                        for ioption_u in range(len(options_paras_raven[8])):

                                            npara_m = len(options_paras_raven[0][ioption_m])
                                            npara_n = len(options_paras_raven[1][ioption_n])
                                            npara_o = len(options_paras_raven[2][ioption_o])
                                            npara_p = len(options_paras_raven[3][ioption_p])
                                            npara_q = len(options_paras_raven[4][ioption_q])
                                            npara_r = len(options_paras_raven[5][ioption_r])
                                            npara_s = len(options_paras_raven[6][ioption_s])
                                            npara_t = len(options_paras_raven[7][ioption_t])
                                            npara_u = len(options_paras_raven[8][ioption_u])

                                            nnparas = npara_m + npara_n + npara_o + npara_p + npara_q + npara_r + npara_s + npara_t + npara_u

                                            model_runs += nsets * (nnparas+2)
                                
        print("Total number of model runs required: ",model_runs)


