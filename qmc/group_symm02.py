# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 01:17:03 2022

@author: mytum
"""
import numpy as np
import qmc.adcqmc_integer_representation01 as adcrep

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

class group:
    def __init__(self,groups):
        if groups=='c2' or groups=='C2':
            self.groupName = "C2"
            self.size = 2
            self.irrep_set = ['A','B']
            self.irrep_ir_set = [0,1]
            self.dict={0:'A',1:'B'}
            self.reverse_dict = {'A':0,'B':1}
            self.direct_product_table = np.array([['A','B'],['B','A']])
        elif groups=='cs' or groups=='CS' or groups=="Cs" or groups=="cS":
            self.groupName = "CS"
            self.size = 2
            self.irrep_set = ["A'","A\""]
            self.irrep_ir_set = [0,1]
            self.dict={0:"A'",1:"A\""}
            self.reverse_dict = {"A'":0,"A\"":1}
            self.direct_product_table = np.array([["A'","A\""],["A\"","A'"]])
        elif groups=="D2H" or groups=="d2h" or groups=="D2h" or groups=="d2H":
            self.groupName = "D2H"
            self.size = 8
            self.irrep_set = ["Ag","B1g", "B2g","B3g","Au","B1u","B2u","B3u"]
            self.irrep_ir_set = range(8)
            self.dict={}
            self.reverse_dict={}
            for i in range(self.size):
                self.dict[i]=self.irrep_set[i]
                self.reverse_dict[self.irrep_set[i]]=i
            #self.dict={0:"A'",1:"A\""}
            #self.reverse_dict = {"A'":0,"A\"":1}
            self.direct_product_table = np.array([
                func2("Ag	    B1g	    B2g	    B3g	    Au	    B1u	    B2u	    B3u"),
                func2("B1g	    Ag	    B3g	    B2g	    B1u	    Au	    B3u	    B2u"),
                func2("B2g	    B3g	    Ag	    B1g	    B2u	    B3u	    Au	    B1u"),
                func2("B3g	    B2g	    B1g	    Ag	    B3u	    B2u	    B1u	    Au"),
                func2("Au	    B1u	    B2u	    B3u	    Ag	    B1g	    B2g	    B3g"),
                func2("B1u	    Au	    B3u	    B2u	    B1g	    Ag	    B3g	    B2g"),
                func2("B2u	    B3u	    Au	    B1u	    B2g	    B3g	    Ag	    B1g"),
                func2("B3u	    B2u	    B1u	    Au	    B3g	    B2g	    B1g	    Ag")              
                ])            
        elif groups=="C2V" or groups=="C2v" or groups=="c2v" or groups=="c2V":
            self.groupName = "C2V"
            self.size = 4
            self.irrep_set = ["A1","A2","B1","B2"]
            self.irrep_ir_set = range(self.size)
            self.dict={}
            self.reverse_dict={}
            for i in range(self.size):
                self.dict[i]=self.irrep_set[i]
                self.reverse_dict[self.irrep_set[i]]=i
            #self.dict={0:"A'",1:"A\""}
            #self.reverse_dict = {"A'":0,"A\"":1}
            self.direct_product_table = np.array([
                func2("A1	    A2	    B1	    B2"),
                func2("A2	    A1	    B2	    B1"),
                func2("B1	    B2	    A1	    A2"),
                func2("B2	    B1	    A2	    A1")
                ]) 
        else:
            assert 0
    def get_direct_product_with_id(self,p,q):
        assert p<self.size and q<self.size
        product_name = self.direct_product_table[p,q]
        product_id = self.reverse_dict[product_name]
        return product_id

    def get_direct_product_with_id4(self,p,q,r,s):
        assert p<self.size and q<self.size and r < self.size and s < self.size
        product_name1 = self.direct_product_table[p,q]
        product_id1 = self.reverse_dict[product_name1]
        product_name2 = self.direct_product_table[r,s]
        product_id2 = self.reverse_dict[product_name2]
        final_name = self.direct_product_table[product_id1,product_id2]
        final_id = self.reverse_dict[final_name]
        return final_id

    def get_direct_product_with_name(self,p,q):
        pID, qID = self.reverse_dict[p], self.reverse_dict[q]
        product_name = self.direct_product_table[pID,qID]
        return product_name
      
    
    def find_irrep_id_orb(self, orbsym, order, spin, A,B):
        a,b,c,d = 1,A,A+1,B
        basis, orbs = adcrep.find_basis_orbs_from_order(order, spin,a,b,c,d)
        if spin==1:
            if basis==1:
                orb1, orb2 = orbs[0], orbs[1]
                irrep1, irrep2 = orbsym[orb1-1], orbsym[orb2-1]
                return self.get_direct_product_with_id(irrep1, irrep2)
            elif basis==2 or basis==3:
                orb1, orb2, orb3, orb4 = orbs[0], orbs[1], orbs[2], orbs[3]
                irrep1, irrep2, irrep3, irrep4 = orbsym[orb1-1], orbsym[orb2-1], orbsym[orb3-1],orbsym[orb4-1]
                return self.get_direct_product_with_id4(irrep1, irrep2, irrep3, irrep4)
            elif basis==4:
                orb1, orb2, orb3 = orbs[0], orbs[1], orbs[2]
                irrep1, irrep2, irrep3 = orbsym[orb1-1], orbsym[orb2-1], orbsym[orb3-1]
                return self.get_direct_product_with_id(irrep2, irrep3)        
            elif basis==5:
                orb1, orb2, orb3 = orbs[0], orbs[1], orbs[2]
                irrep1, irrep2, irrep3 = orbsym[orb1-1], orbsym[orb2-1], orbsym[orb3-1]
                return self.get_direct_product_with_id(irrep1, irrep2)     
            else :
                return self.get_direct_product_with_id(0, 0) 
        else:
            if basis==1:
                orb1, orb2 = orbs[0], orbs[1]
                irrep1, irrep2 = orbsym[orb1-1], orbsym[orb2-1]
                return self.get_direct_product_with_id(irrep1, irrep2)
            elif basis==2 or basis==3 or basis==4:
                orb1, orb2, orb3, orb4 = orbs[0], orbs[1], orbs[2], orbs[3]
                irrep1, irrep2, irrep3, irrep4 = orbsym[orb1-1], orbsym[orb2-1], orbsym[orb3-1],orbsym[orb4-1]
                return self.get_direct_product_with_id4(irrep1, irrep2, irrep3, irrep4)
            elif basis==5:
                orb1, orb2, orb3 = orbs[0], orbs[1], orbs[2]
                irrep1, irrep2, irrep3 = orbsym[orb1-1], orbsym[orb2-1], orbsym[orb3-1]
                return self.get_direct_product_with_id(irrep2, irrep3)        
            elif basis==6:
                orb1, orb2, orb3 = orbs[0], orbs[1], orbs[2]
                irrep1, irrep2, irrep3 = orbsym[orb1-1], orbsym[orb2-1], orbsym[orb3-1]
                return self.get_direct_product_with_id(irrep1, irrep2)     
            else :
                return self.get_direct_product_with_id(0, 0)            
    def find_irrep_name_orb(self, orbsym, order, spin, A,B):
         irrep_id = self.find_irrep_id_orb(orbsym, order, spin, A,B)
         return self.dict[irrep_id]
        
# g = group('c2')    
# print(g.get_direct_product_with_id(0,1))
# print(g.get_direct_product_with_name('B','A'))