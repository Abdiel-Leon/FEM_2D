#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:04:45 2023
@author: abdiel

### Code to produce FEM analyses in 2D for very simple geometries###
Benchmark adopted from lecture notes of Prof. Suvrano
https://homepages.rpi.edu/~des/IFEA2007Fall.html

"""
import numpy
import sympy as sym
from sympy import MatrixSymbol, Matrix
import matplotlib.pyplot as plt

class FEM_2D_element_quad:
    
    def __init__(self,n_elem,*args):
        self.L = args[0]
        self.H = args[1]
        self.q = args[2]
        self.E = args[3]
        self.v = args[4]
        self.n_elem = n_elem
        self.dof = 8       
        self.t = 0.5 #  for plane-strain problems, thickness asssumed as 1 m.
        self.q_vector = sym.zeros(2, 1) # traction vector in Pa       
        self.q_vector[1] = self.q  
        
        


    def Compute_C(self):

     #Material Matrix (plane stress)
      self.C = numpy.zeros((3,3))      
      self.C[0,0] = 1.;             self.C[1,1] = 1.;   self.C[2,2] = (1.-self.v)/2.; 
      self.C[0,1] = self.v; self.C[1,0] = self.v; 
      self.C = self.E/(1-self.v*self.v)*self.C
      return 0   

    
    def Compute_k1(self):
        B = numpy.zeros((3, 6))        
        x1=self.L;  x2=self.L; x3=0.0;
        y1=0;  y2=self.H; y3=0.0;  
        
        B[0,0] = y2-y3; B[0,2] = y3-y1; B[0,4] = y1-y2;
        B[1,1] = x3-x2; B[1,3] = x1-x3; B[1,5] = x2-x1;                   
        B[2,0] = x3-x2; B[2,1] = y2-y3; B[2,2] = x1-x3;\
        B[2,3] = y3-y1; B[2,4] = x2-x1; B[2,5] = y1-y2; 
        
        J = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
        A = J/2
        
        B = 1/(2*A)*B
        aux_ke_1 = numpy.matmul(numpy.transpose(B),self.C)
        ke_1 = self.t*A*numpy.matmul(aux_ke_1,B)
        return ke_1   

    def Compute_k2(self):
        B = numpy.zeros((3, 6))    
        self.N = numpy.zeros((2, 6))        
        
        x1=0;       x2=0;   x3=self.L;
        y1=self.H;  y2=0.0; y3=self.H;  
        
        
        
        B[0,0] = y2-y3; B[0,2] = y3-y1; B[0,4] = y1-y2;
        B[1,1] = x3-x2; B[1,3] = x1-x3; B[1,5] = x2-x1;                   
        B[2,0] = x3-x2; B[2,1] = y2-y3; B[2,2] = x1-x3;\
        B[2,3] = y3-y1; B[2,4] = x2-x1; B[2,5] = y1-y2; 
        
        J = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
        A = J/2
        
        B = 1/(2*A)*B
        
        aux_ke_2 = numpy.matmul(numpy.transpose(B),self.C)
        ke_2 = self.t*A*numpy.matmul(aux_ke_2,B)
        return ke_2      


    
    def Define_Local_Load_Vector(self):
        # acting on the edge of element 2
       self.N = sym.zeros(2, 6)
       fpoint = numpy.zeros((6,1))
       fpoint[5] = -1000 #point load of 1000 pounds

       self.psi = sym.Symbol('psi')
       self.eta = sym.Symbol('eta')

       # shape functions
       self.N1 = 1-self.psi-self.eta
       self.N2 = self.psi
       self.N3 = self.eta
      
       self.N[0,0] = self.N1;        self.N[0,2] = self.N2;        self.N[0,4] = self.N3;
       self.N[1,1] = self.N1;        self.N[1,3] = self.N2;        self.N[1,5] = self.N3;

       f=self.t*self.N.T*self.q_vector
       #evaluation of the integral at psi = 0 (top edge of element)
      # psi1 =0
       f = f.subs([(self.psi, 0)])
       f = sym.integrate(f,self.eta)
       f = f.subs([(self.eta, 1)]) - f.subs([(self.eta, 0)])
       # deta_dx -> from parametric to natural coordinates
       x = self.N1*0 + self.N2*0+ self.N3*3 
       dx_deta = sym.diff(x, self.eta)

       f = f*dx_deta + fpoint
       # add point load of 1000 pounds
       print(f)
       
       return f
    

    def Compute_Kg(self):
        self.K = numpy.zeros((self.dof,self.dof)) 
        Id = numpy.zeros((8))
        Id = [0,1,2,3,6,7,4,5,6,7,2,3]

        for k in range(self.n_elem):
            
            if (k==0):
                ke = self.Compute_k1()
                #print("k1 = ", ke)
            if (k==1):
                ke = self.Compute_k2()  
                #print("k2 = ", ke)                             
                
                
            for i in range(6):
                for j in range(6):
                    self.K[Id[i + 6*k],Id[j + 6*k]] = ke[i,j] + self.K[Id[i + 6*k],Id[j + 6*k]]    
          
        return   self.K
    
   


    def check_symmetric(self, tol=1e-8):
        
        return numpy.all(numpy.abs(self.K-self.K.T) < tol)

 
    def Define_Global_Load_Vector(self,f): 
        self.Q = numpy.zeros(self.dof) 

        Id_load = numpy.zeros((6))
        Id_load = [4,5,6,7,2,3]
                        
                                
        for i in range(6):
                    self.Q[Id_load[i]] = f[i]    
        # point load
        return  
    
    

    def Compute_Kred(self):
       self.restrictions  =  5   
       self.dof_block= numpy.zeros(self.dof-self.restrictions)        
       self.Kred = numpy.zeros((self.dof-self.restrictions,self.dof-self.restrictions))    
       
       # Id matrix with free dof
       self.dof_block = [0,2,3]
       
       for i in range(self.dof-self.restrictions):
           for j in range(self.dof-self.restrictions):
                  # print(self.dof_block[i],self.dof_block[j])    
                   self.Kred[i,j] = self.K[self.dof_block[i],self.dof_block[j]]
                   
      # print(self.K)
       return self.Kred    

    def Compute_Qred(self):
      
       self.Qred = numpy.zeros(self.dof-self.restrictions)    
       
       for i in range(self.dof-self.restrictions):   
                  # print(self.dof_block[i], self.Qred)
                   self.Qred[i] = self.Q[self.dof_block[i]]        
       return self.Qred       

    def Solve_System(self,Kred,Qred):
        Kinv = numpy.linalg.inv(Kred)
        self.u = numpy.dot(Kinv,Qred)
        return self.u

   # def Update_positions(self):                  

    def Solve(self):
        self.Compute_C()
        self.Compute_Kg()        
        symmetry = self.check_symmetric()
        if(symmetry==False):
          print("IS Not Symmetric, please check!!")
        
        Kred = self.Compute_Kred()
        print("Kred = ", Kred)
        f = self.Define_Local_Load_Vector()
        self.Define_Global_Load_Vector(f)      
        Qred= self.Compute_Qred()
        print("Qred = ", Qred)
        # Calculate displacements
        u = self.Solve_System(Kred,Qred) 
        print("u = ", u)
        print("u_ref = ",[0.2337e-4,0.1069e-4,-0.9084e-4])
        
       
        #print(C)       
                
     


    
    #def Plot_displacements(self):
        # plots original and deformed bar lengths
       # plt.ylim(0.5,1.1)
        #y = [1] * (self.n + 1)
        #yupd = [1.05] * (self.n + 1)        
        #plt.plot(self.x0,y,label="Original bar",marker=8)
        #plt.plot(self.xupd,yupd, label="Stressed bar",marker=8)
        #plt.legend(loc="lower right")    
    
def main():
 number_elements = 2
 Length = 3
 Height = 2
 q_load = -300 #distribuited load of 300 psi
 E_modulus = 30e6
 v_poisson = 0.25
 args = (Length,Height,q_load,E_modulus,v_poisson)   
 
 Beam = FEM_2D_element_quad(number_elements,*args) 
 Beam.Solve() 

 #pf, fbf = pressure.Compute_fb()   
 #plt.plot(pf, fbf)
 
 
main() 