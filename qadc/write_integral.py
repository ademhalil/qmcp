# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:36:03 2024

@author: mytum
"""

from . import group_symm02 as group
from . import adcrep
from . import  myadcc
from pyscf import gto, scf, symm
import adcc
import numpy as np
import pyscf
from numpy.linalg import eigh
import struct 
import math



def func2(line):
    if line[-1]=='\n':
        a=line[:-1].split(" ")
    else:
        a=line.split(" ")
    while a.count('\n')>0:
        a.remove('\n')
    while a.count('')>0:
        a.remove('')
    
    for i in range(len(a)):
        a[i]=a[i].strip('\t')
    return a


	##I is the order of \sum_b \sum_a=0^b, i.e   a<=b
def get_ab(I):
    amin = 0
    bmax = 0.5 * (-1 + math.sqrt(1 - 8 * (amin - I)))
    b = int(bmax)
    a = I - b * (b + 1) // 2

    assert a <= b
    return a, b

def generate_Nmu(mygroup, orbsymm_NUMPY):
    
    orbsymm = [a for a in orbsymm_NUMPY]
    Nmu = [0 for a in range(mygroup.size)]
    mu_orbs= [[] for a in range(mygroup.size)]
    
    for i,s in enumerate(orbsymm):
        Nmu[s] +=1
        mu_orbs[s].append(i)
        
    print("orbsymm:",orbsymm)
    print("Nmu:", Nmu)
    print("mu_orbs:", mu_orbs)
    
    return orbsymm, Nmu, mu_orbs

#def write_integrals4(eri_4fold, orb, energies, integral_flnm, gp, Nmu, mu_orbs, mosym):
def write_integrals4(HFrun,integral_flnm, flag_use_sym=True):   
    mol = HFrun.mol
    orb = HFrun.mo_coeff
    eri_4fold = mol.ao2mo(HFrun.mo_coeff)
    energies = HFrun.mo_energy 
    norb = orb.shape[0]
    assert mol.nelec[0]==mol.nelec[1]
    nocc = mol.nelec[0] 
    A,B = nocc, norb
    
    if flag_use_sym:
        mosym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
        gp = group.group(mol.groupname)
    else:
        mosym = [0 for i in range(B)]
        gp = group.group("C1")
        
    orbsymm, Nmu, mu_orbs = generate_Nmu(gp, mosym)

    
    #####
    integral_test_flnm = integral_flnm +"_test"
    integral_bin_flnm = integral_flnm +"_bin"
    integral_human_flnm = integral_flnm +"_human"
    integral_file =  open(integral_human_flnm,'w')
    integral_file_bin =  open(integral_bin_flnm,'wb')
    integral_file_test = open(integral_test_flnm,'w')
    
    print("integral { group "+ gp.groupName +"  binaryfile "+ integral_bin_flnm +" norb "+str(orb.shape[1]), file=integral_file, end=" ")
    print("irrep_mo_map_id { ", file=integral_file, end=" ")
    for s in mosym:
        print(s, file=integral_file, end=" ")
    print(" } ", file=integral_file, end=" ")
    print(" energies { ", file=integral_file, end=" ")
    for s in energies:
        print(s, file=integral_file, end=" ")
    print(" } ", file=integral_file, end=" ")
    print(" } ", file=integral_file, end=" \n")
    
    
    s = gp.size
    gamma1max = s*(s+1)//2
    Ksize = gamma1max*(gamma1max +1)//2
    
    Ksize_total =0;
    for Ixyzt in range(Ksize):
        Ixy,Izt = get_ab(Ixyzt)
        x,y = get_ab(Ixy)
        z,t = get_ab(Izt)
        
        if gp.does_ABCD_to_lead_to_E(x, y, z, t, 0):
            Jzt_max = (Nmu[t] * (Nmu[t] + 1) // 2) if z == t else (Nmu[t] * Nmu[z])
            Jxy_max = (Nmu[y] * (Nmu[y] + 1) // 2) if x == y else (Nmu[x] * Nmu[y])
            Jxyzt_max = (Jzt_max * (Jzt_max + 1)//2) if Ixy == Izt else (Jzt_max * Jxy_max)
    
            Ksize_total += Jxyzt_max

    #K = np.ndarray((Ksize_total,), dtype='d')
    
    count =0
    for Ixyzt in range(Ksize):
        Ixy,Izt = get_ab(Ixyzt)
        x,y = get_ab(Ixy)
        z,t = get_ab(Izt)
        
        if gp.does_ABCD_to_lead_to_E(x, y, z, t, 0) and min(Nmu[x], Nmu[y], Nmu[z], Nmu[t])>0:  
            Jzt_max = (Nmu[t] * (Nmu[t] + 1) // 2) if z == t else (Nmu[t] * Nmu[z])
            Jxy_max = (Nmu[y] * (Nmu[y] + 1) // 2) if x == y else (Nmu[x] * Nmu[y])
            Jxyzt_max = (Jzt_max * (Jzt_max + 1)//2) if Ixy == Izt else (Jzt_max * Jxy_max)
            
            for Jabcd_xyzt in range(Jxyzt_max):
                if Ixy==Izt:
                    Jab_xy,Jcd_zt = get_ab(Jabcd_xyzt)
                else:
                    Jab_xy = Jabcd_xyzt % Jxy_max
                    Jcd_zt = Jabcd_xyzt // Jxy_max
                
                if x==y:
                    a,b = get_ab(Jab_xy)
                else:
                    a =  Jab_xy % Nmu[x]
                    b = Jab_xy // Nmu[x]
                
                if z==t:
                    c,d = get_ab(Jcd_zt)
                else:
                    c =  Jcd_zt % Nmu[z]
                    d = Jcd_zt // Nmu[z]         
                    
                
                
                # print("Ixyzt/Ksize:",Ixyzt,"/",Ksize,"Jabcd_xyzt/Jxyzt_max:",Jabcd_xyzt, "/",Jxyzt_max, "Jab_xy:",Jab_xy,"Jcd_zt:",Jcd_zt)
                # print("Jxy_max:",Jxy_max,"Jzt_max:",Jzt_max,"Jxyzt_max:",Jxyzt_max)
                # print("Ixy,Izt:", Ixy, Izt)
                # print("a,b,c,d:",a,b,c,d)
                # print("x,y,z,t:",x,y,z,t)
                # print("mu_orbs:", mu_orbs)
                ax = mu_orbs[x][a]
                by = mu_orbs[y][b]
                cz = mu_orbs[z][c]
                dt = mu_orbs[t][d]
                
                i,j = min(ax,by), max(ax,by)
                k,L = min(cz,dt), max(cz, dt)
                
                ij = j*(j+1)//2 + i
                kL = L*(L+1)//2 + k
                
                
                q,w = min(ij,kL), max(ij,kL)
                #K[count] = eri_4fold[q,w]
                integral_file_bin.write(struct.pack('d',eri_4fold[q,w]))
                
                integral_file_test.write(f"xyzt {x} {y} {z} {z}: abcd {a} {b} {c} {d} :\
                                         axbyczdt {ax} {by} {cz} {dt} : eri {eri_4fold[q,w]} \n")
                count +=1
                
    assert count == Ksize_total
    
    integral_file.close()
    integral_file_bin.close()
    integral_file_test.close()
    
    return integral_human_flnm, integral_bin_flnm

def write_mol_sym(HFrun,output_flnm):
    
    mol = HFrun.mol
    orb = HFrun.mo_coeff
    # eri_4fold = mol.ao2mo(HFrun.mo_coeff)
    energies = HFrun.mo_energy 
    orbsym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
    gp = group.group(mol.groupname)
    orbsymm, Nmu, mu_orbs = adcrep.generate_Nmu(gp, orbsym)
    norb = orb.shape[0]
    assert mol.nelec[0]==mol.nelec[1]
    nocc = mol.nelec[0] 
    A,B = nocc, norb        
    spin = mol.spin + 1
    symmetry_id_numbers = adcrep.get_numbers_with_gp(gp, orbsym, spin, A,B)    
    
    with open(output_flnm, 'w') as f:
        a,b,c,d = 1,A,A+1,B
        sizetot, n1,n2,n3,n4,n5,n6 = adcrep.get_sizetot(spin, a,b,c,d,'adc2')
        

        print("sizetot=",sizetot," numbers: n1:",n1," n2:",n2," n3:",n3," n4:",n4," n5:",n5," n6:",n6, file=f)
        print("A=",A,"B=",B, "Nvirtual=", B-A, file=f)
        print("spin=",spin,"basis=", mol.basis,"symmetry:",mol.symmetry, file=f)
        
        print("a,b,c,d:",a,b,c,d," A,B=",A,B, file=f)

    
        print("point group: ", mol.topgroup, " supported: ", mol.groupname, file=f)

        
        
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





def write_adc_state(adcstate, HFrun, nstate, nmax, output_flnm, mol_name, flag_Yvector_wbinaryfile):
    
    
    mol = HFrun.mol
    orb = HFrun.mo_coeff
    # eri_4fold = mol.ao2mo(HFrun.mo_coeff)
    # energies = HFrun.mo_energy 
    orbsym = symm.label_orb_symm(mol,mol.irrep_id, mol.symm_orb, HFrun.mo_coeff)
    gp = group.group(mol.groupname)
    orbsymm, Nmu, mu_orbs = adcrep.generate_Nmu(gp, orbsym)
    norb = orb.shape[0]
    assert mol.nelec[0]==mol.nelec[1]
    nocc = mol.nelec[0] 
    A,B = nocc, norb        
    spin = mol.spin + 1
    # symmetry_id_numbers = adcrep.get_numbers_with_gp(gp, orbsym, spin, A,B)  
        
    with open(output_flnm, 'w') as f:
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

        Y = myadcc.get_eigenvectors_from_adcc(nstate,spin, A,B, adcstate)
        nmax= min(nmax, Y.shape[1])
        print("Y.shape:",Y.shape)
        print("Y.max and index:", np.argmax(np.abs(Y)))
        index, maxValues = myadcc.find_top_vector_elements(Y, nmax)
        print("index and maxValues", file=f)
        print(index, file=f)
        print(maxValues, file=f)
        # info_section = "molecule: "+mol_name + " mol_basis: " + mol.basis + " spin: " + str(spin) + " method: "+ adcstate.method.name
        # wf_irrep_dict = myadcc.generate_trial_wf_with_basis_with_group(gp, orbsym, index, maxValues, spin, A,B,"trialwf_"+output_flnm+adcstate.method.name, info_section)
        for istate in range(nstate):
            print("istate: ",istate, file=f)
            for counter, imax in enumerate(index[istate]):
                print(counter, imax, Y[istate,imax], file=f)
            print("***", file=f)
            
        if flag_Yvector_wbinaryfile:
            assert myadcc.write_Y_vector_to_binary_file(Y, output_flnm+adcstate.method.name)
        #M = myadcc.read_Y_vector_from_binary_file(flnmroot+method)


