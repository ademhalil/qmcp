# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 14:31:44 2023

@author: mytum
"""

from qmc import adc_calculation01 as calculation

nstate=30
spin=1
flag_mol_symmetry = True
mol_basis = "6-31g"   #mol.basis = 'aug-cc-pvdz'   #mol.basis='def2tzvp'
mol_unit='Angstrom'  #mol.unit='Bohr'
mol_name = "h2o"
mol_geo = '''8  0  0.     0
              1  0  -0.757 0.587
              1  0  0.757  0.587'''
integral_dir =""
# method1 = {'method':'adc2','n_singlets':nstate, 'n_guesses':None, 'max_subspace':None, 'max_iter':None}
# method2 =  {'method':'adc2x','n_singlets':nstate, 'n_guesses':None, 'max_subspace':None, 'max_iter':None}
method1 = {'method':'adc2','n_singlets':nstate, 'n_guesses':None, 'max_subspace':None, 'max_iter':100}
#method2 =  {'method':'adc2x','n_singlets':nstate, 'n_guesses':None, 'max_subspace':None, 'max_iter':100}

adc_methods = [method1] #["adc2s","adc2x"]
flag_write_integral = False 
flag_chkfile = True
chkfile = "chkfile_"+mol_name

target_irrep_id = 0

adcstate_vec  = calculation.calculate(nstate,spin, flag_mol_symmetry, mol_basis, mol_unit, mol_name, mol_geo, integral_dir, adc_methods, flag_write_integral, flag_chkfile, chkfile)




