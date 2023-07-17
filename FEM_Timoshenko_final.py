#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:04:45 2023
@author: abdiel

### FEA based on Timoshenko beam theory###
Element stiffness matrix was compared to other solution and is ok!
"""
import numpy
import sympy as sym
#from sympy import MatrixSymbol, Matrix
#import matplotlib.pyplot as plt
#from array import *

class FEM_2D_element:
    
    def __init__(self,n_elem,*args):
       self.L_value = args[0]
       self.E_value = args[1]
       self.I_value = args[2]
       self.G_value = args[3] 
       self.A_value = args[4]*5/6 # 5/6 is the shear reducxtion factor in timoshenko beams
       self.n_elem = n_elem
       self.dof = n_elem*2 + 2
        
       self.L = sym.Symbol('L')        
       self.E = sym.Symbol('E')        
       self.I = sym.Symbol('I')        
       self.G = sym.Symbol('G')          
       self.A = sym.Symbol('A') 
       
    def Compute_Shape_Functions(self):
        
       self.N = sym.zeros(2, 4)
      
       self.psi = sym.Symbol('psi')


       # shape functions
       self.N1 = 1./2.*(1. - self.psi)
       self.N2 = 1./2.*(1. + self.psi)

       
       self.N[0,0] = self.N1; self.N[0,2] = self.N2; 
       self.N[1,1] = self.N1; self.N[1,3] = self.N2; 

       # derivatives w.r.t. psi and eta
       
       self.N1_psi = sym.diff(self.N1, self.psi)
       self.N2_psi = sym.diff(self.N2, self.psi) 
       return 

    def Compute_J(self):
        
       det_J = self.L/2
       return det_J

    def Compute_B_matrix(self,det_J):
       #self.B = numpy.zeros((3,6))
  
       # array does not work with symbolic math
       B = sym.zeros(2, 4)#
        
       # J^(-1) = |x,psi y,psi|
       #          |x,eta y,eta|   
       
       self.N1_x = -1/self.L
       self.N2_x = 1/self.L     


       
       # loading strain-disp matrix
       #phi = cross section rotation
       #w = vertical displacement
       #gamma = phi + dw/dw = shear angle
       
       #k (curvature) = d_phi/dx 
       B[0,1] = self.N1_x; B[0,3] = self.N2_x; 
       B[1,0] = self.N1_x; B[1,1] = self.N1;  
       B[1,2] = self.N2_x; B[1,3] = self.N2;  


       return B
            
 #       return p,fb    


    def Compute_C(self):

      self.C = sym.zeros(2,2)      
      self.C[0,0] = self.E*self.I;             
      self.C[1,1] = self.G*self.A;  
    

       
    def Calculate_ke(self):
                
       aux_ke  = sym.zeros(4,2)
       
       det_J = self.Compute_J()
       B = self.Compute_B_matrix(det_J)
       aux_ke = B.T*self.C
       ke = aux_ke*B*self.L/2

       ke_int = sym.integrate(ke, (self.psi, -1, 1))
       ke_int = sym.simplify(ke_int)
       #evaluate the stiffness integral for the different values of area, etc.
       ke_int = ke_int.evalf(subs={self.A:self.A_value,self.L:self.L_value,self.E:self.E_value,self.I:self.I_value,self.G:self.G_value}) 
     #  print(ke_int)
       return ke_int
    



    def Compute_Kg(self):
        self.K = numpy.zeros((self.dof,self.dof)) 
        ke = self.Calculate_ke()
        Id = numpy.zeros((self.n_elem*4))
        
        for k in range(self.n_elem):
             for i in range(4):           
                 Id[i+4*k] = i+2*k
# convert to list         
        Id = Id.tolist() # list   
        Id = [int(x) for x in Id]
        for k in range(self.n_elem):

            for i in range(4):
                for j in range(4):
                    self.K[Id[i + 4*k],Id[j + 4*k]] = ke[i,j] + self.K[Id[i + 4*k],Id[j + 4*k]]

        return   self.K
 

    def Compute_Kred(self,K):
       self.restrictions  =  2 
       self.free_dof = self.dof-self.restrictions
       self.dof_block= numpy.zeros(self.free_dof)        
       self.Kred = numpy.zeros((self.free_dof,self.free_dof))    
       
       # Id matrix containing restrictions
       for i in range(self.free_dof):
          self.dof_block[i] = i+2
# convert to list         
       self.dof_block = self.dof_block.tolist() # list   
       self.dof_block = [int(x) for x in self.dof_block]
        
        
       for i in range(self.free_dof):
           for j in range(self.free_dof):
                       
                   self.Kred[i,j] = K[self.dof_block[i],self.dof_block[j]]
                   
      # print(self.Kred)
       return self.Kred   
         

    def Compute_Qred(self):
      
       self.Qred = numpy.zeros(self.dof-self.restrictions)    
       self.Qred[self.n_elem*2 - 2] = -1000
      # for i in range(self.dof-self.restrictions):   
                   
       #            self.Qred[i] = f[self.dof_block[i]]        
       return self.Qred       

    def Solve_System(self):
        self.u = numpy.dot(numpy.linalg.inv(self.Kred),self.Qred)
        return self.u

   # def Update_positions(self):                  

    def Solve(self):
        self.Compute_Shape_Functions()
        self.Compute_C()
        K = self.Compute_Kg()
        self.Compute_Kred(K)         
        self.Compute_Qred()
       # print("K_red",self.Kred)
        u = self.Solve_System()
        print("u_FEM",u)
        print("u_tip-point-load", ["3.35e-4 m","1.007e-3 rads"])
       
        #print(C)       
  
    
def main():
 number_elements = 1000
 Length = 0.5/number_elements
 E_modulus = 69e9
 I = 1.8e-6
 v = 0.33
 G_modulus = E_modulus/(2*(1+v))
 Cross_area = 0.1*0.06
 args = (Length,E_modulus,I,G_modulus,Cross_area)   
 
 Beam = FEM_2D_element(number_elements,*args) 
 Beam.Solve() 

 #pf, fbf = pressure.Compute_fb()   
 #plt.plot(pf, fbf)
 
 
main() 