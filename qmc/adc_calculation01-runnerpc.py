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

#import .group_symm02 as group
from . import group_symm02 as group
#import .adcc_script01 as myadcc
from . import adcc_script01 as myadcc
#import .adcqmc_integer_representation01 as adcrep
from . import adcqmc_integer_representation01 as adcrep
import pyscf
import numpy as np
from pyscf import symm
import sys
#from contextlib import redirect_stdout
from pyscf import symm
from pyscf.geomopt.geometric_solver import optimize



def calculate(nstate,spin, flag_mol_symmetry, mol_basis, mol_unit, mol_name, mol_geo, \
              integral_dir, adc_methods, flag_write_integral, nmax= 100,\
                  IwriteOpt=2, flag_chkfile=False, chkfile=None):

    assert spin==1
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
    flnmroot = str(mol_name)+"_"+str(mol.basis)+"multiplicity."+str(spin)
    output_flnm = "logADC_"+str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
    Yvectorbinfile_rootname = integral_dir + str(mol_name)+"_"+str(mol.basis)+"_multiplicity."+str(spin)

    
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
    
    print("now HF run begins.")
    sys.stdout.flush()
    HFrun, energies, A,B = myadcc.    (mol, flag_chkfile, chkfile)
    print("now HF run ends.")
    sys.stdout.flush()
    #twobody_integral_file_name, energy_file_name = integral_files_flnm_root+"_4s_int.bin" , integral_files_flnm_root +"_ENERGYbin"
    twobody_integral_file_name, energy_file_name = "None", "None"
    if flag_write_integral:
        print("now ERI transformation to MO integrals.")
        sys.stdout.flush()
        eri_4fold = mol.ao2mo(HFrun.mo_coeff)
        print("ERI transformation is done. IwriteOpt=",IwriteOpt)
        sys.stdout.flush()
        if IwriteOpt==2:
            twobody_integral_file_name, energy_file_name = myadcc.write_integrals2(eri_4fold, HFrun.mo_coeff,energies, integral_files_flnm_root)
        elif IwriteOpt==3:
            twobody_integral_file_name, energy_file_name = myadcc.write_integrals3(eri_4fold, HFrun.mo_coeff,energies, integral_files_flnm_root)     
        elif IwriteOpt==4:
            orbsym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
            #mygroup = group.group(mol.topgroup)
            mygroup = group.group(mol.groupname)
            orbsymm, Nmu, mu_orbs = adcrep.generate_Nmu(mygroup, orbsym)
            integral_flnm_human, integral_flnm_bin = myadcc.write_integrals4(eri_4fold, HFrun.mo_coeff, energies, integral_files_flnm_root, mygroup, Nmu, mu_orbs, orbsym)
        else:
            print("write option ",IwriteOpt, " is not recognised.")
            sys.stdout.flush()
            assert 0
        print("Writing MO integrals ends.")
        sys.stdout.flush()

    print("Now orbital information printing:")
####
    with open(output_flnm+"_"+"QMCADCWFinformation", 'w') as f:
        a,b,c,d = 1,A,A+1,B
        sizetot, n1,n2,n3,n4,n5,n6 = adcrep.get_sizetot(spin, a,b,c,d,'adc2')
        
        print(mol_name, file=f)
        print(mol.output, file=f)
        print(output_flnm, file=f)
        print(integral_files_flnm_root, file=f)
        print(twobody_integral_file_name, file=f)
        print(energy_file_name, file=f)
        print("sizetot=",sizetot," numbers: n1:",n1," n2:",n2," n3:",n3," n4:",n4," n5:",n5," n6:",n6, file=f)
        print("A=",A,"B=",B, "Nvirtual=", B-A, file=f)
        print("spin=",spin,"basis=", mol.basis,"symmetry:",mol.symmetry, file=f)
        
        print("a,b,c,d:",a,b,c,d," A,B=",A,B, file=f)

    
        print("point group: ", mol.topgroup, " supported: ", mol.groupname, file=f)
        orbsym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
        #mygroup = group.group(mol.topgroup)
        mygroup = group.group(mol.groupname)
            
        symmetry_id_numbers = adcrep.get_numbers_with_gp(mygroup, orbsym, spin, A,B)
        
        
        print("numbers with group symmetry", file=f)
        for item in symmetry_id_numbers:
            irrep_id = item[0]
            irrep_name = item[1]
            n1,n2,n3,n4,n5,n6 = tuple(item[2])
    
            sizetot = n1 + n2 + n3 + n4 +n5 + n6
            print("irrep:", irrep_name,irrep_id,"sizetot=",\
                  sizetot," numbers: n1:",n1," n2:",n2," n3:",n3," n4:",n4," n5:",n5," n6:",n6, file=f)
            
    
        
        mydict={}
        for i,q in enumerate(mol.irrep_id):
            mydict[q]= mol.irrep_name[i]
        
        # orbsymLmoletter = []
    
        temp_mosym = "mosym { "
        temp_mosym_id = "mosym_id { "
        for i, os in enumerate(orbsym):                
            temp_mosym += "  "+ mydict[os]
            temp_mosym_id += "  "+str(os)
            if (i+1) % 30==0 and i < len(orbsym)-1:
                temp_mosym += "\n"
                temp_mosym_id += "\n"                    
            # orbsymLetter.append(mydict[os])
        temp_mosym += "  } "
        temp_mosym_id += " } "
        
        temps="\n ---orb symmetries---\n"
        for i, os in enumerate(orbsym[0:A]):
            temps += (str(i+1)+"->"+(mydict[os])).rjust(10)
            if (i+1)%5==0 or i == A-1 :
                temps +="\n"
        temps +="\n*****\n"
        temps +="".rjust(10)*((5 - (A%5))%5)
        for i, os in enumerate(orbsym[A:B]):
            temps += (str(A+i+1)+"->"+(mydict[os])).rjust(10)
            if (i+1)%5==0 or i == len(orbsym[A:B])-1 :
                temps +="\n"
        temps +="*****\n"
    
        print(temps, file=f)
        print(temp_mosym, file=f)
        print(temp_mosym_id, file=f)
        f.flush()
###





    
    sys.stdout.flush()
        
    print(adc_methods)
    adcstate_vec = []
    for count, method in enumerate(adc_methods):
        with open(output_flnm+"_"+method['method'], 'w') as f:
            print(method)
            print("\n\n**************** ADC METHOD *********** :",method['method'], file=f)
            print(method,file=f)
            a,b,c,d = 1,A,A+1,B
            sizetot, n1,n2,n3,n4,n5,n6 = adcrep.get_sizetot(spin, a,b,c,d,method['method'])
            
            print(mol_name, file=f)
            print(mol.output, file=f)
            print(output_flnm, file=f)
            print(integral_files_flnm_root, file=f)
            print(twobody_integral_file_name,file=f)
            print(energy_file_name,file=f)
            print("sizetot=",sizetot," numbers: n1:",n1," n2:",n2," n3:",n3," n4:",n4," n5:",n5," n6:",n6, file=f)
            print("A=",A,"B=",B, "Nvirtual=", B-A, file=f)
            print("spin=",spin,"basis=", mol.basis,"symmetry:",mol.symmetry, file=f)
            
            print("a,b,c,d:",a,b,c,d," A,B=",A,B)
            print("method: ",method['method'],"A:",A,"B:",B,"sizetot:",sizetot,\
                  " numbers: n1:",n1," n2:",n2," n3:",n3," n4:",n4," n5:",n5," n6:",n6)
        
            print("point group: ", mol.topgroup, " supported: ", mol.groupname, file=f)
            orbsym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
            #mygroup = group.group(mol.topgroup)
            mygroup = group.group(mol.groupname)
                
            symmetry_id_numbers = adcrep.get_numbers_with_gp(mygroup, orbsym, spin, A,B)
            
            
            print("numbers with group symmetry", file=f)
            for item in symmetry_id_numbers:
                irrep_id = item[0]
                irrep_name = item[1]
                n1,n2,n3,n4,n5,n6 = tuple(item[2])
                if method['method']=='adc1':
                    sizetot = n1
                else:
                    sizetot = n1 + n2 + n3 + n4 +n5 + n6
                    print("irrep:", irrep_name,irrep_id,"sizetot=",sizetot," numbers: n1:",n1," n2:",n2," n3:",n3," n4:",n4," n5:",n5," n6:",n6, file=f)
                    

            
            mydict={}
            for i,q in enumerate(mol.irrep_id):
                mydict[q]= mol.irrep_name[i]
            
            # orbsymLmoletter = []

            temp_mosym = "mosym { "
            temp_mosym_id = "mosym_id { "
            for i, os in enumerate(orbsym):                
                temp_mosym += "  "+ mydict[os]
                temp_mosym_id += "  "+str(os)
                if (i+1) % 30==0 and i < len(orbsym)-1:
                    temp_mosym += "\n"
                    temp_mosym_id += "\n"                    
                # orbsymLetter.append(mydict[os])
            temp_mosym += "  } "
            temp_mosym_id += " } "
            
            temps="\n ---orb symmetries---\n"
            for i, os in enumerate(orbsym[0:A]):
                temps += (str(i+1)+"->"+(mydict[os])).rjust(10)
                if (i+1)%5==0 or i == A-1 :
                    temps +="\n"
            temps +="\n*****\n"
            temps +="".rjust(10)*((5 - (A%5))%5)
            for i, os in enumerate(orbsym[A:B]):
                temps += (str(A+i+1)+"->"+(mydict[os])).rjust(10)
                if (i+1)%5==0 or i == len(orbsym[A:B])-1 :
                    temps +="\n"
            temps +="*****\n"

            print(temps, file=f)
            print(temp_mosym, file=f)
            print(temp_mosym_id, file=f)
        
        
            adcstate = myadcc.adcc.run_adc(HFrun, **method)
            adcstate_vec.append(adcstate)
            

            print(adcstate.describe(), file=f)
            print("*******", file=f)
            print(adcstate.excitation_energy, file=f)
            print("*******", file=f)
            print("state  energy(Ha) energy(eV)", file=f)
            evolt = 27.211396
            for istate, energy in enumerate(adcstate.excitation_energy):
                print(istate, energy, energy*evolt, file=f)
            
            print("*******", file=f)
            #print(adcstate.describe_amplitudes(), file=f)

            Y = myadcc.get_eigenvectors_from_adcc(nstate,spin,sizetot, A,B, adcstate)
            nmax= min(nmax, Y.shape[1])
            print("Y.shape:",Y.shape)
            print("Y.max and index:", np.argmax(np.abs(Y)))
            index, maxValues = myadcc.find_top_vector_elements(Y, nmax)
            print("index and maxValues", file=f)
            print(index, file=f)
            print(maxValues, file=f)
            info_section = "molecule: "+mol_name + " mol_basis: " + mol_basis + " spin: " + str(spin) + " method: "+ method['method']
            wf_irrep_dict = myadcc.generate_trial_wf_with_basis_with_group(mygroup, orbsym, index, maxValues, spin, A,B,"trialwf_"+flnmroot+method['method'], info_section)
            for istate in range(nstate):
                print("istate: ",istate, file=f)
                for counter, imax in enumerate(index[istate]):
                    print(counter, imax, Y[istate,imax], file=f)
                print("***", file=f)
                
        
            assert myadcc.write_Y_vector_to_binary_file(Y, flnmroot+method['method'])
            #M = myadcc.read_Y_vector_from_binary_file(flnmroot+method)

    assert len(adcstate_vec) == len(adc_methods)
    return adcstate_vec


def make_guess_from_adc1(HFrun, mygroup,nstate, target_irrep_id, orbsym, A, B):
    
    adcstate1 = myadcc.adcc.adc1(HFrun, n_singlets=nstate)
    adcstate2 = myadcc.adcc.adc2(HFrun, n_singlets=nstate, max_iter=1)
    
    guess_vec =[]
    for istate in range(nstate):
        orbs = adcstate1.excitation_vector[istate].ph.select_n_absmax(1)[0][0]
        o  =  ( orbs[0] % A ) +1
        v = (orbs[1] % (B-A)) +1
        if mygroup.get_direct_product_with_id( orbsym[o-1],orbsym[v-1]) == target_irrep_id :  
            print("istate ",istate,"meets the irrep.")
            adcstate2.excitation_vector[istate].ph = adcstate1.excitation_vector[istate].ph.copy()
            adcstate2.excitation_vector[istate].pphh = adcstate2.excitation_vector[istate].pphh.zeros_like()
            guess_vec.append(adcstate2.excitation_vector[istate])
    
    print("make guess: len(guess): ", len(guess_vec), " out of  from ",nstate, "states.")
    return guess_vec    
    

def make_guesses2(stateADC, mygroup, spin, orbsymm, target_irrep_id, A, B):
    
    assert spin==1
    
    nstate = len(stateADC.excitation_vector)
    for istate in range(nstate):
        stateADC.excitation_vector[istate].ph = stateADC.excitation_vector[istate].ph.empty_like()

    if stateADC.method.level > 1:
        for istate in range(nstate):
            stateADC.excitation_vector[istate].pphh = stateADC.excitation_vector[istate].pphh.empty_like()
    #first ph block
    for v in range(A+1,B+1):
        for o in range(1,A+1):
            irrep_o = orbsymm[o-1];
            irrep_v = orbsymm[v-1];
            if mygroup.get_direct_product_with_id(irrep_o, irrep_v) == target_irrep_id:
                for alpha in [True, False]:
                    vorb = adcrep.convert_to_adc_orb(v, it_is_occupied_flag=False, it_is_alpha_spin_flag=alpha, A=A,B=B)
                    oorb = adcrep.convert_to_adc_orb(o, it_is_occupied_flag=True, it_is_alpha_spin_flag=alpha, A=A,B=B)
                    for istate in range(nstate):
                        stateADC.excitation_vector[istate].ph[oorb,vorb] = 1.


    if stateADC.method.level > 1:                
        #now pphh block
        for j in range(A+1,B):  
            for k in range(j+1, B+1):
                for L in range(1,A):
                    for m in range(L+1,A+1):
                        if mygroup.get_direct_product_with_id4(orbsymm[j-1], orbsymm[k-1], orbsymm[L-1],orbsymm[m-1]) == target_irrep_id:
                            
                            #basis2 T1
                            vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 2.
                                
                            #basis2 T2
                            vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 2.                            
                                
                                
                             #basis2 T3
                            vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                              
                                
                                
                             #basis2 T4
                            vorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                                
                               
                            #basis2 T5
                            vorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = -1.                                    
                               
                               
                            #basis2 T6
                            vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = -1.
                                
   
                            #####################  
                            #basis3 T1
                            vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                              
                              
                            #basis3 T2
                            vorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                                   
                              
                              
                            #basis3 T3
                            vorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                                
                              
                              
                              
                            #basis3 T4
                            vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                            vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                            oorb1 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                            oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                            for istate in range(nstate):
                                stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                                
                               
                                
                                
        for j in range(A+1,B+1):  
            for L in range(1,A):
                for m in range(L+1,A+1):
                    if mygroup.get_direct_product_with_id( orbsymm[L-1],orbsymm[m-1]) == target_irrep_id :                            
                                
                        #basis4 T1
                        vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                        vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                        oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                        oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                        for istate in range(nstate):
                            stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                            
                                
                        #basis4 T2
                        vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                        vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                        oorb1 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                        oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                        for istate in range(nstate):
                            stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.                                  
                                
                                
        for j in range(A+1,B):  
            for k in range(j+1, B+1):
                for L in range(1,A+1):
                    if mygroup.get_direct_product_with_id(orbsymm[j-1], orbsymm[k-1]) == target_irrep_id:                          
                                
                        #basis5 T1
                        vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                        vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                        oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                        oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                        for istate in range(nstate):
                            stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.   
                            
                        #basis5 T2
                        vorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                        vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                        oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                        oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                        for istate in range(nstate):
                            stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.      
                            
                            
                            
        for j in range(A+1,B+1):  
            for k in range(1,A+1):
                if  target_irrep_id == 0:  
                    
                    #basis6 T1
                    vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                    vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                    oorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                    oorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)      
                    for istate in range(nstate):
                        stateADC.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] = 1.   
                            
                            
    return stateADC
                            
                            