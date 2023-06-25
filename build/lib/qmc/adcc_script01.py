
# Created on Tue Dec 29 17:11:57 2020

# @author: mytum

import group_symm02 as group
from pyscf import gto, scf
import adcc
import numpy as np
import pyscf
from numpy.linalg import eigh
import struct 
import adcqmc_integer_representation01 as adcrep


def write_integrals3(eri_4fold, orb,energies, flnmROOT): 
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
                        print(i+1,j+1,k+1,L+1,eri_4fold[a,b], file=file_twobody_human)
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

def write_integrals2(eri_4fold, orb,energies, flnmROOT): 
    humanEnergyfilnm = flnmROOT+"_humENERGIES"
    N= orb.shape[1]
    twobody_integral_file = open(flnmROOT+"_4s_int.bin",'wb')
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

    energy_file = open(flnmROOT+"_ENERGYbin",'wb')


    [energy_file.write(a) for a in energiesBIN]
    energy_file.close()
    
    gEn=open(humanEnergyfilnm,'w')
    [gEn.write(str(a)+"\n") for a in energies]    
    gEn.close()


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
                
            else:
                assert 0
                
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