import numpy as np
from . import group_symm02 as group

def generate_config_set(spin,a,b,c,d):
    occ=range(a,b+1)
    vir = range(c,d+1)
    
    basis_set=[]
    orbs_set=[]
    if spin==1:
        for k in occ:
            for j in vir:
                basis=1
                orbs=[j,k]
                basis_set.append(basis)
                orbs_set.append(orbs)

        for L in occ:
            for k in occ:
                for j in vir:
                    for i in vir:
                        if i<j and k<L:
                            basis=2
                            orbs=[i,j,k,L]
                            basis_set.append(basis)
                            orbs_set.append(orbs)
                            
        for L in occ:
            for k in occ:
                for j in vir:
                    for i in vir:
                        if i<j and k<L:
                            basis=3
                            orbs=[i,j,k,L]
                            basis_set.append(basis)
                            orbs_set.append(orbs)    
                            
        for i in vir:
            for L in occ:
                for k in occ:
                    if k<L:
                        basis=4
                        orbs=[i,k,L]
                        basis_set.append(basis)
                        orbs_set.append(orbs)  
        
        for k in occ:
            for j in vir:
                for i in vir:
                    if i<j:
                        basis=5
                        orbs=[i,j,k]
                        basis_set.append(basis)
                        orbs_set.append(orbs)     
                        
        for i in vir:
            for k in occ:
                basis=6
                orbs=[i,k]
                basis_set.append(basis)
                orbs_set.append(orbs) 

    elif spin==3:
        for k in occ:
            for j in vir:
                basis=1
                orbs=[j,k]
                basis_set.append(basis)
                orbs_set.append(orbs)

        for L in occ:
            for k in occ:
                for j in vir:
                    for i in vir:
                        if i<j and k<L:
                            basis=2
                            orbs=[i,j,k,L]
                            basis_set.append(basis)
                            orbs_set.append(orbs)
                            
        for L in occ:
            for k in occ:
                for j in vir:
                    for i in vir:
                        if i<j and k<L:
                            basis=3
                            orbs=[i,j,k,L]
                            basis_set.append(basis)
                            orbs_set.append(orbs)
                            
        for L in occ:
            for k in occ:
                for j in vir:
                    for i in vir:
                        if i<j and k<L:
                            basis=4
                            orbs=[i,j,k,L]
                            basis_set.append(basis)
                            orbs_set.append(orbs)
                            

        for i in vir:
            for L in occ:
                for k in occ:
                    if k<L:
                        basis=5
                        orbs=[i,k,L]
                        basis_set.append(basis)
                        orbs_set.append(orbs)  
        
        for k in occ:
            for j in vir:
                for i in vir:
                    if i<j:
                        basis=6
                        orbs=[i,j,k]
                        basis_set.append(basis)
                        orbs_set.append(orbs)     
                        
    return basis_set, orbs_set
    
def get_numbers(spin,a,b,c,d):
    if spin==1:
        n1= (b-a+1)*(d-c+1)
        n2 = (d-c)*(d-c+1)//2*(b-a)*(b-a+1)//2
        n3=n2
        n4=(b-a+1)*(b-a)//2*(d-c+1)
        n5=(d-c)*(d-c+1)//2*(b-a+1)
        n6 = (b-a+1)*(d-c+1)
    elif spin==3:
        n1= (b-a+1)*(d-c+1)
        n2 = (d-c)*(d-c+1)//2*(b-a)*(b-a+1)//2
        n3=n2
        n4=n2
        n5=(b-a+1)*(b-a)//2*(d-c+1)
        n6 = (d-c)*(d-c+1)//2*(b-a+1)     
    else:
        assert 0

    return n1,n2,n3,n4,n5,n6

def get_blocks_with_gp(mygroup, target_irrep_id):
    s1=[]
    s2=[]
    s4=[]
    for x in mygroup.irrep_id_set:
        for y in mygroup.irrep_id_set: 
            if mygroup.get_direct_product_with_id(x, y)== target_irrep_id:
                s1.append([x,y])
            if x<=y:
                s4.append([x,y])
    
    for pair in s1:
        x,y = tuple(pair)
        
        for i in mygroup.irrep_id_set:
            for j in mygroup.irrep_id_set:
                if j >= i and mygroup.get_direct_product_with_id(i, j)==x:
                    for k in mygroup.irrep_id_set:
                        for L in mygroup.irrep_id_set:
                            if L>=k and mygroup.get_direct_product_with_id(k, L)==y:
                                s2.append([i,j,k,L])
                                
    return s1,s2,s4
                    

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
            
def get_numbers_with_gp(mygroup, orbsymm_NUMPY, spin, A,B):
    
    orbsymm = [a for a in orbsymm_NUMPY]
    symmetry_id_numbers = []
    for target_irrep_id in mygroup.irrep_id_set:
        s1,s2,s4 = get_blocks_with_gp(mygroup, target_irrep_id)
        
        m1 =0
        for pair in s1:
            x,y = tuple(pair)
            nv = orbsymm[A:B].count(x)
            no = orbsymm[0:A].count(y)
            m1 += nv*no

        m1occ =0
        for pair in s4:
            x,y = tuple(pair)
            nv = B-A
            if x < y:
                no = orbsymm[0:A].count(x)*orbsymm[0:A].count(y)
            elif x==y:
                no = orbsymm[0:A].count(x)*(orbsymm[0:A].count(x)-1)//2
            else:
                assert 0
            m1occ += nv*no        
        
        m1vir =0
        for pair in s4:
            x,y = tuple(pair)
            no = A
            if x < y:
                nv = orbsymm[A:B].count(x)*orbsymm[A:B].count(y)
            elif x==y:
                nv = orbsymm[A:B].count(x)*(orbsymm[A:B].count(x)-1)//2
            else:
                assert 0
            m1vir += nv*no      
            
        m2 =0
        for pair in s2:
            x,y = tuple(pair)[:2]
            
            if x < y:
                no = orbsymm[0:A].count(x)*orbsymm[0:A].count(y)
            elif x==y:
                no = orbsymm[0:A].count(x)*(orbsymm[0:A].count(x)-1)//2
            else:
                assert 0
                
            x,y = tuple(pair)[2:4]
            
            if x < y:
                nv = orbsymm[A:B].count(x)*orbsymm[A:B].count(y)
            elif x==y:
                nv = orbsymm[A:B].count(x)*(orbsymm[A:B].count(x)-1)//2
            else:
                assert 0
            m2 += nv*no   
            
        m6 = (A)*(B-A) if target_irrep_id==0 else 0
        
        if spin==1:
            n1 = m1
            n2 = m2
            n3 = m2
            n4 = m1occ
            n5 = m1vir
            n6 = m6
        elif spin==3:
            n1 = m1
            n2 = m2
            n3 = m2
            n4 = m2
            n5 = m1occ
            n6 = m1vir
        else :
            assert 0
        
        symmetry_id_numbers.append([target_irrep_id, mygroup.irrep_set[target_irrep_id], [n1,n2,n3,n4,n5,n6]] )
    return symmetry_id_numbers
        
def get_sizetot(spin,a,b,c,d, adcmethod):
    #a,b,c,d = window[0],window[1], window[2],window[3]    
    n1,n2,n3,n4,n5,n6 = get_numbers(spin,a,b,c,d)
    if adcmethod=='adc1':
        sizetot = n1
    elif adcmethod=='adc2' or adcmethod=='adc2x' or adcmethod=='adc3':
        sizetot = n1 + n2 + n3 + n4 + n5 + n6
    else:
        print(adcmethod, "not defined")
        assert 0
    return sizetot, n1, n2, n3, n4, n5, n6

def find_jk_from_n1(N,a,b,c,d):
    if N%(d-c+1)==0:
        j=d
        k=N//(d-c+1) +a -1
    else:
        j = N%(d-c+1) +c -1
        k = N//(d-c+1) +a 
        
    return j,k

def find_ijkL_from_n2(N,a,b,c,d):
    if N%((d-c+1)*(d-c)//2)==0:
        j=d
        i=d-1
        NkL = N//((d-c+1)*(d-c)//2)
        Lmax = a + 0.5*(1+ np.sqrt(1+8*(NkL-1)))
        L = int(Lmax)
        k = NkL +a-1- (L-a)*(L-a-1)//2
        assert k<L
    else:
        Nij = N%((d-c+1)*(d-c)//2)
        jmax = c + 0.5*(1+np.sqrt(1+8*(Nij-1)))
        j = int(jmax)
        i = Nij +c-1 - (j-c)*(j-c-1)//2
        NkL = N//((d-c+1)*(d-c)//2)
        Lmax = a + 0.5*(1+ np.sqrt(1+8*NkL))
        L = int(Lmax)
        k = NkL +a - (L-a)*(L-a-1)//2
        assert k<L
    
    return i,j,k,L

def find_ikL_from_n4(N,a,b,c,d):
    
    if N%((b-a+1)*(b-a)//2)==0:
        L=b
        k=b-1
        i=N//((b-a+1)*(b-a)//2) +c -1
    else:
        NkL = N%((b-a+1)*(b-a)//2)
        Lmax = a + 0.5*(1+ np.sqrt(1+8*(NkL-1)))
        L = int(Lmax)
        k = NkL +a -1 - (L-a)*(L-a-1)//2
        i = N//((b-a+1)*(b-a)//2) + c
        
    return i,k,L

def find_ijk_from_n5(N,a,b,c,d):
    
    if N%((d-c+1)*(d-c)//2)==0:
        j=d
        i=d-1
        k= a-1 + N//((d-c+1)*(d-c)//2)
    else:
        Nij = N%((d-c+1)*(d-c)//2)
        jmax = c + 0.5*(1+np.sqrt(1+8*(Nij-1)))
        j = int(jmax)
        i = Nij +c-1 - (j-c)*(j-c-1)//2
        
        assert i<j
        k = N//((d-c+1)*(d-c)//2) + a
        
    return i,j,k

def find_ik_from_n6(N,a,b,c,d):
    if N%(b-a+1)==0:
        k=b
        i=N//(b-a+1) + c -1
    else:
        i = N//(b-a+1) + c 
        k = N%(b-a+1) + a -1
    
    return i,k


def get_n1_from_jk(j,k,a,b,c,d):
    n1jk = (k-a)*(d-c+1) + (j-c+1)
    return n1jk

def get_n2_fromijkL(i,j,k,L,a,b,c,d):
    n2ijkL = (L-a)*(L-a-1)//2*(d-c+1)*(d-c)//2 + (k-a)*(d-c+1)*(d-c)//2+ \
        (j-c)*(j-c-1)//2+(i-c+1)
    return n2ijkL

def get_n4_from_ikL(i,k,L,a,b,c,d):
    n4ikL = (i-c)*(b-a+1)*(b-a)//2 + (L-a)*(L-a-1)//2 + k-a +1
    return n4ikL

def get_n5_from_ijk(i,j,k,a,b,c,d):
    n5ijk = (k-a)*(d-c+1)*(d-c)//2 + (j-c)*(j-c-1)//2 + (i-c+1)
    return n5ijk

def get_n6_from_ik(i,k,a,b,c,d):
    n6ik = (i-c)*(b-a+1) + (k-a+1) 
    return n6ik


def find_order_from_basis_orbs(basis, orbs, spin,a,b,c,d):
    
    
    if spin==1:
        if basis==1:
            j,k = orbs
            order = get_n1_from_jk(j,k,a,b,c,d)
        elif basis==2:
            i,j,k,L = orbs
            order = n1 + get_n2_fromijkL(i,j,k,L,a,b,c,d)
        elif basis==3:
            i,j,k,L = orbs
            order = n1 + n2 + get_n2_fromijkL(i,j,k,L,a,b,c,d)            
        elif basis==4:
            i,k,L = orbs
            order = n1 + n2 + n3 +  get_n4_from_ikL(i,k,L,a,b,c,d)
        elif basis==5:
            i,j,k = orbs
            order = n1 + n2 + n3 + n4 + get_n5_from_ijk(i,j,k,a,b,c,d)            
        elif basis==6:
            i,k = orbs
            order = n1 + n2 + n3 + n4 + n5+ get_n6_from_ik(i,k,a,b,c,d)
        else:
            assert 0
            
    elif spin==3:
        if basis==1:
            j,k = orbs
            order = get_n1_from_jk(j,k,a,b,c,d)
        elif basis==2:
            i,j,k,L = orbs
            order = n1 + get_n2_fromijkL(i,j,k,L,a,b,c,d)            
        elif basis==3:
            i,j,k,L = orbs
            order = n1 + n2+ get_n2_fromijkL(i,j,k,L,a,b,c,d)           
        elif basis==4:
            i,j,k,L = orbs
            order = n1 + n2+ n3+ get_n2_fromijkL(i,j,k,L,a,b,c,d)           
        elif basis==5:
            i,k,L = orbs
            order = n1 + n2 + n3 + n4+ get_n4_from_ikL(i,k,L,a,b,c,d)        
        elif basis==6:
            i,j,k = orbs
            order = n1 + n2 + n3 + n4 + n5+ get_n5_from_ijk(i,j,k,a,b,c,d)   
    else:
        assert 0
    
    return order

def find_basis_orbs_from_order(order, spin,a,b,c,d):
    
    n1,n2,n3,n4,n5,n6 = get_numbers(spin,a,b,c,d)
    
    if spin==1:    
        if order <= n1:
            basis=1
            j,k = find_jk_from_n1(order,a,b,c,d)
            orbs=[j,k]
        elif order <= n1 + n2:
            basis = 2
            i,j,k,L = find_ijkL_from_n2(order-n1,a,b,c,d)
            orbs=[i,j,k,L]
        elif order <= n1 + n2+n3:
            basis = 3
            i,j,k,L = find_ijkL_from_n2(order-n1-n2,a,b,c,d)     
            orbs=[i,j,k,L]
        elif order <= n1 + n2 + n3 + n4:
            basis = 4
            i,k,L = find_ikL_from_n4(order-n1-n2-n3,a,b,c,d)  
            orbs=[i,k,L]
        elif order <= n1 + n2 + n3 + n4 + n5:
            basis = 5
            i,j,k = find_ijk_from_n5(order- n1 - n2 - n3 - n4,a,b,c,d)
            orbs=[i,j,k]
        elif order <= n1 + n2 + n3 + n4 + n5 + n6:
            basis = 6
            i,k  = find_ik_from_n6(order- n1 - n2 - n3 - n4-n5,a,b,c,d)
            orbs=[i,k]
        else:
            assert 0
    
    elif spin==3:
        if order <= n1:
            basis=1
            j,k = find_jk_from_n1(order,a,b,c,d)  
            orbs=[j,k]
        elif order <= n1 + n2:
            basis = 2
            i,j,k,L = find_ijkL_from_n2(order-n1,a,b,c,d)  
            orbs=[i,j,k,L]
        elif order <= n1 + n2+n3:
            basis = 3
            i,j,k,L = find_ijkL_from_n2(order-n1-n2,a,b,c,d) 
            orbs=[i,j,k,L]
        elif order <= n1 + n2+n3 + n4:
            basis = 4
            i,j,k,L = find_ijkL_from_n2(order-n1-n2-n3,a,b,c,d)
            orbs=[i,j,k,L]            
        elif order <= n1 + n2 + n3 + n4 +n5:
            basis = 5
            i,k,L = find_ikL_from_n4(order-n1-n2-n3-n4,a,b,c,d)  
            orbs=[i,k,L]
        elif order <= n1 + n2 + n3 + n4 + n5 +n6:
            basis = 6
            i,j,k = find_ijk_from_n5(order- n1 - n2 - n3 - n4 - n5,a,b,c,d)
            orbs=[i,j,k]
    
    return basis, orbs
  
def test_func(basis_set, orbs_set, spin,a,b,c,d):
    n1,n2,n3,n4,n5,n6 = get_numbers(spin,a,b,c,d)
    for  q in range(len(basis_set)):
        order = q+1
        order_ = find_order_from_basis_orbs(basis_set[q], orbs_set[q], spin,a,b,c,d)
        if not order == order_:
            print("test failed in finding order. spin=", spin)
            print("order:", order,"found order:", order_ ,"basis:", basis_set[q], "orbs:", orbs_set[q])
            assert 0
        basis_, orbs_ = find_basis_orbs_from_order(order, spin,a,b,c,d)
        if not (basis_set[q]==basis_ and orbs_set[q]==orbs_):
            print("test failed in finding basis and orbs. spin=",spin)
            print("order:", order, "basis:", basis_set[q], "orbs:", orbs_set[q])  
            print("basis_:", basis_, "orbs_:", orbs_)  
            assert 0
   
def convert_to_adc_orb(j, it_is_occupied_flag, it_is_alpha_spin_flag,A,B, frozen_core=0, frozen_vir=0):
    #note a,b,c,d = 1, A,A+1, B
    
    if it_is_occupied_flag:
        assert  1<=j and j<=A
        if it_is_alpha_spin_flag:
            return j-1-frozen_core
        else:
            return j+A-1-frozen_core*2
    else:
        assert A+1<=j and  j<=B
        if it_is_alpha_spin_flag:
            return j-A-1
        else:
            return j-A-1+(B-A)- frozen_vir
    
        
# def define_globals(a_,b_,c_,d_):       
#     a,b,c,d = a_,b_,c_,d_
#     # print("window :",a,b,c,d)
#     # occ=range(a,b+1)
#     # vir = range(c,d+1)
#     # print("occ:",list(occ))
#     # print("vir:",list(vir))
    
# def get_globals():
#     return a,b,c,d,occ, vir,n1,n2,n3,n4,n5,n6

if 0:
    for i in range(100):
        a=np.random.randint(1,50)
        b= a+ np.random.randint(1,50)
        c=b+1
        d=c+ np.random.randint(1,50)
        
        print("window :",a,b,c,d)
        # occ=range(a,b+1)
        # vir = range(c,d+1)
        
        spin=1
        basis_set, orbs_set = generate_config_set(spin,a,b,c,d)
        assert len(basis_set)==len(orbs_set)
        sizetot = len(basis_set)     
        #sizetot_, n1,n2,n3,n4,n5,n6 = get_sizetot(spin,[a,b,c,d])
        #assert sizetot ==sizetot_
        n1,n2,n3,n4,n5,n6 = get_numbers(spin,a,b,c,d)
        assert sizetot == n1 + n2 + n3 + n4 + n5 + n6   
        print("spin=",spin,"sizetot=", sizetot)
        test_func(basis_set, orbs_set, spin,a,b,c,d)
           
        spin=3
        basis_set, orbs_set = generate_config_set(spin,a,b,c,d)
        assert len(basis_set)==len(orbs_set)
        sizetot = len(basis_set)     
        n1,n2,n3,n4,n5,n6 = get_numbers(spin,a,b,c,d)
        assert sizetot == n1 + n2 + n3 + n4 + n5 + n6   
        print("spin=",spin,"sizetot=", sizetot)
        test_func(basis_set, orbs_set, spin,a,b,c,d)        
        