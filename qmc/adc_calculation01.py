# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 10:06:17 2022

@author: mytum
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 13:45:52 2021

@author: mytum
"""

# from IPython import get_ipython
# get_ipython().magic('reset -sf') 

import qmc.group_symm02 as group
import qmc.adcc_script01 as myadcc
import qmc.adcqmc_integer_representation01 as adcrep
import pyscf
import numpy as np
from pyscf import symm
#import sys
from contextlib import redirect_stdout
from pyscf import symm
from pyscf.geomopt.geometric_solver import optimize



def calculate(nstate,spin, flag_mol_symmetry, mol_basis, mol_unit, mol_name, mol_geo, integral_dir, adc_methods, flag_write_integral, \
              flag_chkfile=False, chkfile=None):
    
    # ###
    mol = pyscf.gto.Mole()
    mol.verbose = 5
    
    mol.symmetry = flag_mol_symmetry

    mol.basis= mol_basis
    mol.unit= mol_unit


    # nstate=26
    # spin=1
    # adc_methods = ["adc2s","adc2x"]
    # mol_name = "c4h4o"
    # flag_write_integral = False
    
    mol.output = 'logHF_'+str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
    # integral_dir = "/gpfs/bwfor/home/hd/hd_hd/hd_mf262/QMCADC/runs/molecules/c4h4o/basis_631g/integrals/"
    integral_files_flnm_root = integral_dir + str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
    flnmroot = str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
    output_flnm = "logADC_"+str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
    

    
    mol.atom = mol_geo
    
    # mol.atom = '''
    # O
    # C 1 1.362
    # C 1 1.362 2 106.7
    # C 2 1.361 1 110.7 3 0.
    # C 3 1.361 1 110.7 2 0.
    # H 2 1.075 1 115.9 3 180.
    # H 3 1.075 1 115.9 2 180.
    # H 4 1.077 2 126.1 1 180.
    # H 5 1.077 3 126.1 1 180.
    # '''
    
    ##########################################
    
    mol.build()
    
    with open(output_flnm, 'w') as f:
        with redirect_stdout(f):
            HFrun, energies, A,B = myadcc.do_HFrun(mol, flag_chkfile, chkfile)
            #mol_eq = optimize(HFrun, maxsteps=100)
            a,b,c,d = 1,A,A+1,B
            sizetot, n1,n2,n3,n4,n5,n6 = adcrep.get_sizetot(spin, a,b,c,d)
            print(mol_name)
            print(mol.output)
            print(output_flnm)
            print(integral_files_flnm_root)
            print("sizetot=",sizetot)
            print("A=",A,"B=",B, "Nvirtual=", B-A)
            print("spin=",spin,"basis=", mol.basis,"symmetry:",mol.symmetry)
            
            if flag_write_integral:
                eri_4fold = mol.ao2mo(HFrun.mo_coeff)
                myadcc.write_integrals2(eri_4fold, HFrun.mo_coeff,energies, integral_files_flnm_root)
        
    
            print("point group: ", mol.topgroup, " supported: ", mol.groupname)
            orbsym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
            
            mygroup = group.group(mol.topgroup)
            
            mydict={}
            for i,a in enumerate(mol.irrep_id):
                mydict[a]= mol.irrep_name[i]
            
            orbsymLetter = []
            temps=""
            temp_mosym = "mosym { "
            for i, os in enumerate(orbsym):
                temps += (str(i+1)+"->"+(mydict[os])).rjust(10)
                if (i+1)%5==0 :
                    print(temps)
                    temps=""
                    if i+1==A:
                        print("*****")
                else:
                    if i+1==A:
                        print(temps)
                        print("*****")
                        count=1
                        temps=""
                        while (i+1+count)%5 !=0:
                            temps +="".rjust(10)
                            count +=1
                
                temp_mosym += "     "+ mydict[os]
                orbsymLetter.append(mydict[os])
            temp_mosym += "  } "
            print(temps)
            print(temp_mosym)
        
        
            for count, method in enumerate(adc_methods):
                print("\n\n**************** ADC METHOD *********** :",method['method'])
                # if spin==3:
                #     if method=="adc2s":
                #         adcstate = myadcc.adcc.adc2(HFrun, n_triplets=nstate)
                #     elif method=="adc2x":
                #         adcstate = myadcc.adcc.adc2x(HFrun, n_triplets=nstate)
                #     else:
                #         assert 0
                        
                # else:
                #     if method=="adc2s":
                #         adcstate = myadcc.adcc.adc2(HFrun, n_singlets=nstate)
                #     elif method=="adc2x":
                #         adcstate = myadcc.adcc.adc2x(HFrun, n_singlets=nstate)     
                #     else :
                #         assert 0
    
                adcstate = myadcc.adcc.run_adc(HFrun, **method) 
    
                print(adcstate.describe())
                print("*******")
                print(adcstate.excitation_energy)
                print("*******")
                print("state  energy(Ha) energy(eV)")
                evolt = 27.211396
                for istate, energy in enumerate(adcstate.excitation_energy):
                    print(istate, energy, energy*evolt)
                
                print("*******")
                print(adcstate.describe_amplitudes())
    
                Y = myadcc.get_eigenvectors_from_adcc(nstate,spin,sizetot, A,B, adcstate)
                nmax=100
                index, maxValues = myadcc.find_top_vector_elements(Y, nmax)
                info_section = "molecule: "+mol_name + " mol_basis: " + mol_basis + " spin: " + str(spin) + " method: "+ method['method']
                myadcc.generate_trial_wf_with_basis_with_group(mygroup, orbsym, index, maxValues, spin, A,B,"trialwf_"+flnmroot+method['method'], info_section)
                for istate in range(nstate):
                    print("istate: ",istate)
                    for counter, imax in enumerate(index[istate]):
                        print(counter, imax, Y[istate,imax])
                    print("***")
                    
            
                assert myadcc.write_Y_vector_to_binary_file(Y, flnmroot+method['method'])
                #M = myadcc.read_Y_vector_from_binary_file(flnmroot+method)

# # def calculate(nstate,spin, flag_mol_symmetry, mol_basis, mol_unit, mol_name, mol_geo, integral_dir, adc_methods, flag_write_integral, \
# #               flag_chkfile=False, chkfile=None):
# def calculate(args):
#     nstate = args['nstate']
#     spin = args['spin']
#     mol_symmetry = args['mol_symmetry']
#     mol_basis = args['mol_basis']
#     mol_unit = args['mol_unit']
#     mol_name = args['mol_name']
#     mol_geo = args['mol_geo']
#     integral_dir = args['integral_dir']
#     adc_methods = args['adc_methods']
#     flag_write_integral = args['flag_write_integral']
#     if 'flag_chkfile' in args:
#         flag_chkfile = args['flag_chkfile'] 
#     else:
#         flag_chkfile = False
        
#      if 'chkfile' in args:
#          chkfile = args['chkfile'] 
#      else:
#          chkfile = None       
    
#     n_guesses = max(4,2*nstate)
#     if 'n_guesses' in args:
#         n_guesses = args['n_guesses']
    
#     max_iter = 70
#     if 'max_iter' in args:
#         max_iter = args['max_iter']    
        
    
#     max_subspace  = 10*nstate
#     if 'max_subspace' in args:
#         max_subspace = args['max_subspace']   
    
#     # ###
#     mol = pyscf.gto.Mole()
#     mol.verbose = 5
    
#     mol.symmetry = flag_mol_symmetry

#     mol.basis= mol_basis
#     mol.unit= mol_unit


#     # nstate=26
#     # spin=1
#     # adc_methods = ["adc2s","adc2x"]
#     # mol_name = "c4h4o"
#     # flag_write_integral = False
    
#     mol.output = 'logHF_'+str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
#     # integral_dir = "/gpfs/bwfor/home/hd/hd_hd/hd_mf262/QMCADC/runs/molecules/c4h4o/basis_631g/integrals/"
#     integral_files_flnm_root = integral_dir + str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
#     flnmroot = str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
#     output_flnm = "logADC_"+str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
    

    
#     mol.atom = mol_geo
    
#     # mol.atom = '''
#     # O
#     # C 1 1.362
#     # C 1 1.362 2 106.7
#     # C 2 1.361 1 110.7 3 0.
#     # C 3 1.361 1 110.7 2 0.
#     # H 2 1.075 1 115.9 3 180.
#     # H 3 1.075 1 115.9 2 180.
#     # H 4 1.077 2 126.1 1 180.
#     # H 5 1.077 3 126.1 1 180.
#     # '''
    
#     ##########################################
    
#     mol.build()
    
#     with open(output_flnm, 'w') as f:
#         with redirect_stdout(f):
#             HFrun, energies, A,B = myadcc.do_HFrun(mol, flag_chkfile, chkfile)
#             #mol_eq = optimize(HFrun, maxsteps=100)
#             a,b,c,d = 1,A,A+1,B
#             sizetot, n1,n2,n3,n4,n5,n6 = adcrep.get_sizetot(spin, a,b,c,d)
#             print(mol_name)
#             print(mol.output)
#             print(output_flnm)
#             print(integral_files_flnm_root)
#             print("sizetot=",sizetot)
#             print("A=",A,"B=",B, "Nvirtual=", B-A)
#             print("spin=",spin,"basis=", mol.basis,"symmetry:",mol.symmetry)
            
#             if flag_write_integral:
#                 eri_4fold = mol.ao2mo(HFrun.mo_coeff)
#                 myadcc.write_integrals2(eri_4fold, HFrun.mo_coeff,energies, integral_files_flnm_root)
        
    
#             print("point group: ", mol.topgroup, " supported: ", mol.groupname)
#             orbsym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
            
#             mygroup = group.group(mol.topgroup)
            
#             mydict={}
#             for i,a in enumerate(mol.irrep_id):
#                 mydict[a]= mol.irrep_name[i]
            
#             orbsymLetter = []
#             temps=""
#             temp_mosym = "mosym { "
#             for i, os in enumerate(orbsym):
#                 temps += (str(i+1)+"->"+(mydict[os])).rjust(10)
#                 if (i+1)%5==0 :
#                     print(temps)
#                     temps=""
#                     if i+1==A:
#                         print("*****")
#                 else:
#                     if i+1==A:
#                         print(temps)
#                         print("*****")
#                         count=1
#                         temps=""
#                         while (i+1+count)%5 !=0:
#                             temps +="".rjust(10)
#                             count +=1
                
#                 temp_mosym += "     "+ mydict[os]
#                 orbsymLetter.append(mydict[os])
#             temp_mosym += "  } "
#             print(temps)
#             print(temp_mosym)
        
        
#             for count, method in enumerate(adc_methods):
#                 print("\n\n**************** ADC METHOD *********** :",method)
#                 if spin==3:
#                     if method=="adc2s":
#                         adcstate = myadcc.adcc.adc2(HFrun, n_triplets=nstate)
#                     elif method=="adc2x":
#                         adcstate = myadcc.adcc.adc2x(HFrun, n_triplets=nstate)
#                     else:
#                         assert 0
                        
#                 else:
#                     if method=="adc2s":
#                         adcstate = myadcc.adcc.adc2(HFrun, n_singlets=nstate)
#                     elif method=="adc2x":
#                         adcstate = myadcc.adcc.adc2x(HFrun, n_singlets=nstate)     
#                     else :
#                         assert 0
    
    
#                 print(adcstate.describe())
#                 print("*******")
#                 print(adcstate.excitation_energy)
#                 print("*******")
#                 print("state  energy(Ha) energy(eV)")
#                 evolt = 27.211396
#                 for istate, energy in enumerate(adcstate.excitation_energy):
#                     print(istate, energy, energy*evolt)
                
#                 print("*******")
#                 print(adcstate.describe_amplitudes())
    
#                 Y = myadcc.get_eigenvectors_from_adcc(nstate,spin,sizetot, A,B, adcstate)
#                 nmax=100
#                 index, maxValues = myadcc.find_top_vector_elements(Y, nmax)
#                 info_section = "molecule: "+mol_name + " mol_basis: " + mol_basis + " spin: " + str(spin) + " method: "+ method
#                 myadcc.generate_trial_wf_with_basis_with_group(mygroup, orbsym, index, maxValues, spin, A,B,"trialwf_"+flnmroot+method, info_section)
#                 for istate in range(nstate):
#                     print("istate: ",istate)
#                     for counter, imax in enumerate(index[istate]):
#                         print(counter, imax, Y[istate,imax])
#                     print("***")
                    
            
#                 assert myadcc.write_Y_vector_to_binary_file(Y, flnmroot+method)
#                 #M = myadcc.read_Y_vector_from_binary_file(flnmroot+method)

    
    
    
    
