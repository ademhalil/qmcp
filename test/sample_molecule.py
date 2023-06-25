# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 14:31:44 2023

@author: mytum
"""

from qmc import adc_calculation01 as calculation

nstate=24
spin=1
flag_mol_symmetry = True
mol_basis = "6-31g"   #mol.basis = 'aug-cc-pvdz'   #mol.basis='def2tzvp'
mol_unit='Angstrom'  #mol.unit='Bohr'
mol_name = "h2o"
mol_geo = '''8  0  0.     0
              1  0  -0.757 0.587
              1  0  0.757  0.587'''
integral_dir =""
adc_methods = [] #["adc2s","adc2x"]
flag_write_integral = False 
flag_chkfile = True
chkfile = "chkfile_"+mol_name
calculation.calculate(nstate,spin, flag_mol_symmetry, mol_basis, mol_unit, mol_name, mol_geo, integral_dir, adc_methods, flag_write_integral,flag_chkfile, chkfile)