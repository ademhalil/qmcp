
# Created on Tue Dec 29 17:11:57 2020

# @author: mytum

#import .group_symm02 as group
from . import group_symm02 as group
from . import adcqmc_integer_representation01 as adcrep
from pyscf import gto, scf
import adcc
import numpy as np
import pyscf
from numpy.linalg import eigh
import struct 
import math
# import qmc.adcqmc_integer_representation01 as adcrep


	##I is the order of \sum_b \sum_a=0^b, i.e   a<=b
def get_ab(I):
    amin = 0
    bmax = 0.5 * (-1 + math.sqrt(1 - 8 * (amin - I)))
    b = int(bmax)
    a = I - b * (b + 1) // 2

    assert a <= b
    return a, b


def write_integrals4(eri_4fold, orb, energies, integral_flnm, gp, Nmu, mu_orbs, mosym):
    
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
                    
                
                
                print("Ixyzt/Ksize:",Ixyzt,"/",Ksize,"Jabcd_xyzt/Jxyzt_max:",Jabcd_xyzt, "/",Jxyzt_max, "Jab_xy:",Jab_xy,"Jcd_zt:",Jcd_zt)
                print("Jxy_max:",Jxy_max,"Jzt_max:",Jzt_max,"Jxyzt_max:",Jxyzt_max)
                print("Ixy,Izt:", Ixy, Izt)
                print("a,b,c,d:",a,b,c,d)
                print("x,y,z,t:",x,y,z,t)
                print("mu_orbs:", mu_orbs)
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

                
    

def write_integrals3(eri_4fold, orb,energies, flnmROOT): 
    
    twobody_integral_file_name = flnmROOT+"_4s_int.bin"
    energy_file_name = flnmROOT+"_ENERGYbin"
    humanEnergyfilnm = flnmROOT+"_humENERGIES"
    N= orb.shape[1]
    twobody_integral_file = open(flnmROOT+"_4s_int.bin",'wb')
    file_twobody_human = open(flnmROOT+"_4s_int.txt",'w')
    file_twobody_human.close()
    file_twobody_human = open(flnmROOT+"_4s_int.txt",'w')
    B = orb.shape[1]
    N1 = B*(B+1)//2
    
    # depo=[]
    # for alpha in range(N1):
    #     for beta in range(alpha,N1):
    #         depo.append(eri_4fold[alpha,beta])
            
    
    flag=False
    a=0
    count=0
    for i in range(orb.shape[1]):
        for j in range(0,i+1):
            
            b=0
            for k in range(orb.shape[1]):
                for L in range(0,k+1):
                    #myint[i,j,k,L]=eri_4fold[a,b]
                    #if abs(eri_4fold[a,b])>10**-10:
                    if True:
                        if i>=j and k>=L:
                            if i<k:
                                flag=True
                            elif i==k:
                                if j<=L:
                                    flag=True
                            else:
                                flag=False
                    
                    if flag:
                        twobody_integral_file.write(struct.pack('d',eri_4fold[a,b]))
                        print(i,j,k,L,eri_4fold[a,b], file=file_twobody_human)
                        # assert eri_4fold[a,b] == depo[count]
                        #print(count, eri_4fold[a,b], depo[count])
                        count +=1
                    flag=False
                   	
                    b +=1
            a +=1
    
    
    
    twobody_integral_file.close()
    file_twobody_human.close()



    energiesBIN = [struct.pack('d',a) for a in energies]


    print(len(energiesBIN), )

    energy_file = open(flnmROOT+"_ENERGYbin",'wb')


    [energy_file.write(a) for a in energiesBIN]
    energy_file.close()
    
    gEn=open(humanEnergyfilnm,'w')
    [gEn.write(str(a)+"\n") for a in energies]    
    gEn.close()
    
    return twobody_integral_file_name, energy_file_name

def write_integrals2(eri_4fold, orb,energies, flnmROOT): 
    humanEnergyfilnm = flnmROOT+"_humENERGIES"
    N= orb.shape[1]
    
    twobody_integral_file_name = flnmROOT+"_4s_int.bin"
    twobody_integral_file = open(twobody_integral_file_name,'wb')
    file_twobody_human = open(flnmROOT+"_4s_int.txt",'w')
    file_twobody_human.close()
    file_twobody_human = open(flnmROOT+"_4s_int.txt",'w')
    B = orb.shape[1]
    N1 = B*(B+1)//2
    
    depo=[]
    for alpha in range(N1):
        for beta in range(alpha,N1):
            depo.append(eri_4fold[alpha,beta])
            
    
    flag=False
    a=0
    count=0
    for i in range(orb.shape[1]):
        for j in range(0,i+1):
            
            b=0
            for k in range(orb.shape[1]):
                for L in range(0,k+1):
                    #myint[i,j,k,L]=eri_4fold[a,b]
                    #if abs(eri_4fold[a,b])>10**-10:
                    if True:
                        if i>=j and k>=L:
                            if i<k:
                                flag=True
                            elif i==k:
                                if j<=L:
                                    flag=True
                            else:
                                flag=False
                    
                    if flag:
                        twobody_integral_file.write(struct.pack('d',eri_4fold[a,b]))
                        print(i+1,j+1,k+1,L+1,eri_4fold[a,b], file=file_twobody_human)
                        assert eri_4fold[a,b] == depo[count]
                        #print(count, eri_4fold[a,b], depo[count])
                        count +=1
                    flag=False
                   	
                    b +=1
            a +=1
    
    
    
    twobody_integral_file.close()
    file_twobody_human.close()



    energiesBIN = [struct.pack('d',a) for a in energies]


    print(len(energiesBIN), )

    energy_file_name = flnmROOT+"_ENERGYbin"
    energy_file = open(energy_file_name,'wb')


    [energy_file.write(a) for a in energiesBIN]
    energy_file.close()
    
    gEn=open(humanEnergyfilnm,'w')
    [gEn.write(str(a)+"\n") for a in energies]    
    gEn.close()
    
    return twobody_integral_file_name, energy_file_name


def write_integtals(eri_4fold, orb,energies, flnmROOT): 
    human4indexflnm = flnmROOT+"_hum4index"
    humanEnergyfilnm = flnmROOT+"_humENERGIES"
    N= orb.shape[1]

    f4=open(human4indexflnm,'w')
    #myint=np.zeros((N,N,N,N))
    depo2a=[]
    depo2b=[]
    flag=False
    a=0
    count=0
    for i in range(orb.shape[1]):
        for j in range(0,i+1):
            
            b=0
            for k in range(orb.shape[1]):
                for L in range(0,k+1):
                    #myint[i,j,k,L]=eri_4fold[a,b]
                    #if abs(eri_4fold[a,b])>10**-10:
                    if True:
                        if i>=j and k>=L:
                            if i<k:
                                flag=True
                            elif i==k:
                                if j<=L:
                                    flag=True
                            else:
                                flag=False
                    
                    if flag:
                        for x in [i+1,j+1,k+1,L+1, count]:
                            depo2a.append(x)
                        depo2b.append(eri_4fold[a,b])
                        count +=1
                        print("{0:3d} {1:3d} {2:3d} {3:3d} {4:3d} {5:3d} {6:5.12f}".\
                        format(i+1,j+1,k+1,L+1,a,b, eri_4fold[a,b]), file=f4)
                    flag=False
                   	
                    b +=1
            a +=1
    
    f4.close()
    
    
    import struct 
    twobodya = [struct.pack('i',a) for a in depo2a]
    twobodyb = [struct.pack('d',a) for a in depo2b]
    print(len(twobodya), len(twobodyb))
    
    twobodya_file = open(flnmROOT+"_ab",'wb')
    twobodyb_file = open(flnmROOT+"_bb",'wb')
    [twobodya_file.write(a) for a in twobodya]
    [twobodyb_file.write(a) for a in twobodyb]
    twobodya_file.close()
    twobodyb_file.close()


    energiesBIN = [struct.pack('d',a) for a in energies]


    print(len(energiesBIN), )

    energy_file = open(flnmROOT+"_ENERGYbin",'wb')


    [energy_file.write(a) for a in energiesBIN]
    energy_file.close()
    
    gEn=open(humanEnergyfilnm,'w')
    [gEn.write(str(a)+"\n") for a in energies]    
    gEn.close()






def do_HFrun(mol, flag_chkfile=False, chkfile=None):
    print(mol.output)    
    HFrun = mol.RHF()
    HFrun.conv_tol = 1e-13
    HFrun.max_cycle=200
    if flag_chkfile:
        HFrun.chkfile = chkfile
        HFrun.init_guess = 'chk'
    HFrun.kernel()
    orb =HFrun.mo_coeff
    #orb1= orb[:,:]
    energies = HFrun.mo_energy 
    #energies1 = energies[:]
    norb = orb.shape[0]
    nocc = int (np.sum(HFrun.get_occ())/2)    
    print("A/B=",nocc, norb)
    A,B = nocc, norb
    
    return HFrun, energies, A,B

def write_4index_integrals(mol,orb,energies, flnmroot):
    eri_4fold = mol.ao2mo(orb)
    write_integtals2(eri_4fold, orb,energies, flnmroot)

def test_integrals(mol,orb,HFrun):
    hcore_ao = mol.intor('int1e_nuc_sph') + mol.intor('int1e_kin_sph')
    hcore_mo = np.matmul(orb.T, np.matmul(hcore_ao, orb))
    #hcore_mo = HFrun.get_hcore()

    mo_ints = pyscf.ao2mo.kernel(mol,orb, aosym=1)


    EHF =0.
    for i in range(nocc):
        EHF +=hcore_mo[i,i]*2
    for i in range(nocc):
        for j in range(nocc):
            EHF += mo_ints[i*norb+i,j*norb+j]*2 - mo_ints[i*norb+j,j*norb+i]
    EHF += HFrun.energy_nuc()
    
    print("EHF from code=",HFrun.e_tot,"  EHF from test=",EHF)


def get_eigenvectors_from_adcc(nstate,spin, sizetot, A,B, state2x):   
    
    assert len(state2x.excitation_vector)==nstate
    a,b,c,d = 1,A,A+1,B
    Y = np.ndarray((nstate,sizetot))

    for order in range(1,sizetot+1):
        #order = i+1
        basis, orbs = adcrep.find_basis_orbs_from_order(order, spin,a,b,c,d)
        if spin==1:
            if basis == 1:
                j, k = orbs[0], orbs[1]
                vorb = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                oorb = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                #print(vorb,oorb,"<---",j,k)
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].ph[oorb,vorb]*np.sqrt(2)
                    
            elif basis ==2:
                j,k,L,m = orbs[0], orbs[1], orbs[2], orbs[3]
                vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                #print(vorb1, vorb2, oorb1, oorb2,"<---",j,k,L,m)
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]*np.sqrt(3)*2         
            elif basis ==3:
                for istate in range(nstate):
                    j,k,L,m = orbs[0], orbs[1], orbs[2], orbs[3]
                    vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                    vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                    oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                    oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                    #print(vorb1, vorb2, oorb1, oorb2,"<---",j,'(',k,')',L,'(',m,')') 
                    temp1 =state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]
                    
                    vorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                    vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                    oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                    oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                    #print(vorb1, vorb2, oorb1, oorb2,"<---",j,'(',k,')',L,'(',m,')') 
                    temp2 =state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]     
                    Y[istate,order-1] = -(temp1+temp2)*2
                    
            elif basis ==4:
                j, k, L = orbs[0], orbs[1], orbs[2]
                vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                oorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                #print(vorb1, vorb2, oorb1, oorb2,"<---",j,'(',j,')',k,'(',L,')') 
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]*(-np.sqrt(2.))*2      

            elif basis ==5:
                j, k, L = orbs[0], orbs[1], orbs[2]
                vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                #print(vorb1, vorb2, oorb1, oorb2,"<---",j,'(',k,')',L,'(',L,')')             
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]*(-np.sqrt(2.))*2 
            elif basis ==6:
                j, k = orbs[0], orbs[1]
                vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                oorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                #print(vorb1, vorb2, oorb1, oorb2,"<---",j,'(',j,')',k,'(',k,')')  
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]*(-1.)*2 
            else:
                assert 0        
                
        elif spin==3:
            if basis == 1:
                j, k = orbs[0], orbs[1]
                vorb = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)   
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].ph[oorb,vorb]*np.sqrt(2)
            elif basis ==2:
                j,k,L,m = orbs[0], orbs[1], orbs[2], orbs[3]
                vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                vorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb1 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb2 = adcrep.convert_to_adc_orb(m, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]*np.sqrt(2)*(-1.)       
            elif basis ==3 or basis==4:
                for istate in range(nstate):
                    i,j,k,L= orbs[0], orbs[1], orbs[2], orbs[3]
                    vorb1 = adcrep.convert_to_adc_orb(i, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                    vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                    oorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                    oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                    temp1 = state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]
                    
                    #j,k,L,m = orbs[0], orbs[1], orbs[2], orbs[3]
                    vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                    vorb2 = adcrep.convert_to_adc_orb(i, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                    oorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                    oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                    temp2 = state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2] 
                    
                    if basis==3:
                        Y[istate,order-1]= (-temp1 + temp2)
                    else: #basis==4
                        Y[istate,order-1]= (temp1 + temp2) 
                        
            elif basis ==5:
                j, k, L = orbs[0], orbs[1], orbs[2]
                vorb1 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                oorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb2 = adcrep.convert_to_adc_orb(L, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]*np.sqrt(2)       
                
            elif basis ==6:
                i,j,k = orbs[0], orbs[1], orbs[2]
                vorb1 = adcrep.convert_to_adc_orb(i, it_is_occupied_flag=False, it_is_alpha_spin_flag=True, A=A,B=B)
                vorb2 = adcrep.convert_to_adc_orb(j, it_is_occupied_flag=False, it_is_alpha_spin_flag=False, A=A,B=B)
                oorb1 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=True, A=A,B=B)
                oorb2 = adcrep.convert_to_adc_orb(k, it_is_occupied_flag=True, it_is_alpha_spin_flag=False, A=A,B=B)
                
                for istate in range(nstate):
                    Y[istate,order-1]=state2x.excitation_vector[istate].pphh[oorb1,oorb2,vorb1,vorb2]*np.sqrt(2)*(-1) 
                
            else :
                assert 0
        else:
            assert 0
    print("Y matrix size:",Y.shape)
    return Y
    
def write_Y_vector_to_binary_file(Y, flnm):
    nstate, L = Y.shape[0], Y.shape[1]
    with open(flnm+".Yvector.bin",'wb') as f:
        f.write(struct.pack('ii',nstate,L))
        for istate in range(nstate):
            for j in range(L):
                f.write(struct.pack('d',Y[istate,j]))
    return 1

def read_Y_vector_from_binary_file(flnm):
    with open(flnm+".Yvector.bin",'rb') as f:
        data = f.read(8)
        (nstate,L) = struct.unpack('@ii',data)
        Y = np.ndarray((nstate,L))
        #assert( len(data)==8+ nstate*L*8)
        for istate in range(nstate):
            for j in range(L):
                data = f.read(8)
                (Y[istate,j],) = struct.unpack('d', data)
        #Y = np.array(struct.unpack('d'*nstate*L,data[8:])).rehsape((nstate,L))
    return Y
        
        

def find_top_vector_elements(Y, nmax):
    nstate, L  = Y.shape[0], Y.shape[1]
    nmax = min(nmax,L)
    index=np.ones((nstate,nmax), dtype=np.int32)*-1
    maxValues = np.ndarray((nstate,nmax))
    for istate in range(nstate):
        for j in range(nmax):
            maxval=0.
            for i in range(L):
                if not (i in index[istate]):
                    if np.abs(Y[istate,i])>=maxval:
                        maxval = np.abs(Y[istate,i])
                        imax = i
            index[istate,j]=imax
            maxValues[istate,j] = Y[istate,imax]
    print("find_top_vector_elements maxValues and index shape:", maxValues.shape, index.shape)
    return index, maxValues
                
def generate_trial_wf(index, maxValues, flnmwf="trialwf"):
    nstate, nconfig  = index.shape[0], index.shape[1]
    with open(flnmwf,'w') as f:
        for istate in range(nstate):
            orderstring = "order { "
            coefstring = "coef { "
            for j in range(nconfig):
                orderstring += str(index[istate,j]+1)+"  "
                coefstring += str(maxValues[istate,j]/np.sign(maxValues[istate,0]))+" "
            orderstring +=" }  "
            coefstring +=" } "
            print(orderstring, file=f)
            print(coefstring,file=f)
            
def generate_trial_wf_with_basis(index, maxValues, spin, A,B, flnmwf="trialwf"):
    nstate, nconfig  = index.shape[0], index.shape[1]
    a,b,c,d = 1,A,A+1,B
    
    with open(flnmwf,'w') as f:
        for istate in range(nstate):
            orderstring = "basis { "
            coefstring = "coef { "
            for j in range(nconfig):
                order = index[istate,j]+1
                basis, orbs = adcrep.find_basis_orbs_from_order(order, spin,a,b,c,d)
                #orderstring += str(index[istate,j]+1)+"  "
                orderstring += str(basis)+" "
                for k in orbs:
                    orderstring +=str(k)+" "
                coefstring += str(maxValues[istate,j]/np.sign(maxValues[istate,0]))+" "
            orderstring +=" }  "
            coefstring +=" } "
            print(orderstring, file=f)
            print(coefstring,file=f)


    with open(flnmwf,'a') as f:
        print("\n\n\n*******", file=f)
        for istate in range(nstate):
            print("state: ",istate, file=f)
            print("i".rjust(4)+"coef".rjust(14)+"basis".rjust(7)+"orbs".rjust(16), file=f)
            for j in range(nconfig):
                order = index[istate,j]+1
                coef = maxValues[istate,j]/np.sign(maxValues[istate,0])
                basis, orbs = adcrep.find_basis_orbs_from_order(order, spin,a,b,c,d)
                temps = str(j).rjust(4)+str(round(coef,7)).rjust(14)+str(basis).rjust(7)
                for k in orbs:
                    temps +=str(k).rjust(4)
                print(temps, file=f)
                #coefstring += str(maxValues[istate,j]/np.sign(maxValues[istate,0]))+" "  
                
def generate_trial_wf_with_basis_with_group(mygroup, orbsym, index, maxValues, spin, A,B, flnmwf="trialwf", info_section=""):
    nstate, nconfig  = index.shape[0], index.shape[1]
    a,b,c,d = 1,A,A+1,B
    
    with open(flnmwf,'w') as f:
        print(info_section, file=f)
        print("#nstate: ",nstate, "nconfig: ",nconfig, file=f)
        print("#******\n\n")
        for istate in range(nstate):
            print("#state ",istate, file=f)
            orderstring = "basis { "
            coefstring = "coef { "
            for j in range(nconfig):
                order = index[istate,j]+1
                basis, orbs = adcrep.find_basis_orbs_from_order(order, spin,a,b,c,d)
                #orderstring += str(index[istate,j]+1)+"  "
                orderstring += str(basis)+" "
                for k in orbs:
                    orderstring +=str(k)+" "
                coefstring += str(maxValues[istate,j]/np.sign(maxValues[istate,0]))+" "
            orderstring +=" }  "
            coefstring +=" } "
            print(orderstring, file=f)
            print(coefstring,file=f)


    with open(flnmwf,'a') as f:
        print("\n\n\n*******", file=f)
        mydict ={}
        for istate in range(nstate):
            j = 0
            order = index[istate,j]+1
            irrep = mygroup.find_irrep_name_orb(orbsym, order, spin, A,B)
            if irrep in mydict:
                mydict[irrep].append(istate)
            else:
                mydict[irrep] = [istate]
        
        print("irreps states", file=f)
        for irrep in mydict:
            temps = "irrep: " + str(irrep) +" --> "
            for istate in mydict[irrep]:
                temps += str(istate)+"  "
            print(temps, file=f)
        temps =""
        
            
        print("\n\n\n*******", file=f)
        for istate in range(nstate):
            print("state: ",istate, file=f)
            print("i".rjust(4)+"irrep".rjust(6)+"coef".rjust(14)+"basis".rjust(7)+"orbs".rjust(16), file=f)
            for j in range(nconfig):
                order = index[istate,j]+1
                irrep = mygroup.find_irrep_name_orb(orbsym, order, spin, A,B)
                coef = maxValues[istate,j]/np.sign(maxValues[istate,0])
                basis, orbs = adcrep.find_basis_orbs_from_order(order, spin,a,b,c,d)
                temps = str(j).rjust(4)+ irrep.rjust(6) +str(round(coef,7)).rjust(14)+str(basis).rjust(7)
                for k in orbs:
                    temps +=str(k).rjust(4)
                print(temps, file=f)
                #coefstring += str(maxValues[istate,j]/np.sign(maxValues[istate,0]))+" "


    return mydict
# aTest =0.

# for i in range(norb):
#     for j in range(norb):
#         aTest += mo_ints[i*norb+i,j*norb+j]*2 - mo_ints[i*norb+j,j*norb+i]
# print("norb:",norb)
# print("(iijj)*2 -(ijji) from 1:norb+1,  aTest: ", aTest)


# Run an ADC(3) calculation, solving for 3 singlets
#state2 = adcc.adc2(HFrun, n_singlets=3)
#state2x = adcc.adc2x(HFrun, n_singlets=3)

#A, B = nocc, norb
#return A,B, HFrun, mo_ints, energies
#return A,B, HFrun, energies, mol


#icase=1
#print("\n\n****** icase=",icase,"  ********\n")
#A,B, HFrun, mo_ints, energies = get_scf(icase)
#A,B, HFrun, energies, mol = get_scf(icase)  
# state1singlet = adcc.adc1(HFrun, n_singlets=30)
# state1triplet = adcc.adc1(HFrun, n_triplets=30)

# state2singlet = adcc.adc2(HFrun, n_singlets=40)
# state2triplet = adcc.adc2(HFrun, n_triplets=40)

# state2xsinglet = adcc.adc2x(HFrun, n_singlets=40)
# state2xtriplet = adcc.adc2x(HFrun, n_triplets=40)

# state3singlet = adcc.adc3(HFrun, n_singlets=40)
# state3triplet = adcc.adc3(HFrun, n_triplets=40)