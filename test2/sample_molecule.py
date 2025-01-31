# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 14:31:44 2023

@author: mytum
"""

from qadc import write_integral
import pyscf
import adcc


spin=1
flag_mol_symmetry = True
mol_basis = "6-31g"   #mol.basis = 'aug-cc-pvdz'   #mol.basis='def2tzvp'
mol_unit='Angstrom'  #mol.unit='Bohr'
mol_name = "h2o"
mol_geo = '''8  0  0.     0
              1  0  -0.757 0.587
              1  0  0.757  0.587'''

flag_chkfile = False
chkfile = ""
integral_dir =""

#######
mol = pyscf.gto.Mole()
mol.verbose = 5

mol.symmetry = flag_mol_symmetry

mol.basis= mol_basis
mol.unit= mol_unit


mol.output = 'logHF_'+str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
integral_files_flnm_root = integral_dir + str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
flnmroot = str(mol_name)+"_"+str(mol.basis)+"multiplicity."+str(spin)
output_flnm = "logADC_"+str(mol_name)+str(mol.basis)+"multiplicity."+str(spin)
Yvectorbinfile_rootname = integral_dir + str(mol_name)+"_"+str(mol.basis)+"_multiplicity."+str(spin)


mol.atom = mol_geo
mol.build()

print(mol.output)    
HFrun = mol.RHF()
HFrun.conv_tol = 1e-13
HFrun.max_cycle=200
if flag_chkfile:
    HFrun.chkfile = chkfile
    HFrun.init_guess = 'chk'
HFrun.kernel()

nstate=30
nmax = 10
flag_Yvector_wbinaryfile = False
method1 = {'method':'adc2','n_singlets':nstate, 'n_guesses':None, 'max_subspace':None, 'max_iter':100}
adcstate = adcc.run_adc(HFrun, **method1)
    
write_integral.write_integrals4(HFrun,integral_files_flnm_root)
write_integral.write_mol_sym(HFrun,output_flnm+"_mosym")
write_integral.write_adc_state(adcstate, HFrun, nstate, nmax, output_flnm+"_ADC", mol_name, flag_Yvector_wbinaryfile)
