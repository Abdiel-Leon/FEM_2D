#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:04:45 2023
@author: abdiel

### Code to produce FEM analyses in 2D for very simple geometries###
"""
import numpy
import sympy as sym
from sympy import MatrixSymbol, Matrix
import matplotlib.pyplot as plt

class FEM_2D_element:
    
    def __init__(self,n_elem,*args):
        self.L = args[0]
        self.H = args[1]
        self.q = args[2]
        self.E = args[3]
        self.v = args[4]
        self.n_elem = n_elem
        self.dof = self.n_elem*8        
        self.t = 1. #  for plane-strain problems, thickness asssumed as 1 m.
        self.q_vector = sym.zeros(2, 1) # traction vector in Pa       
        self.q_vector[1] = self.q  
        
        
        
    
    def Compute_Geometry(self):
        x = numpy.zeros(int(self.n_elem + 1))
        y = numpy.zeros(int(self.n_elem + 1))
        x = numpy.linspace(0, self.L, num=numpy.math.ceil(self.n_elem/2+1))   
        y = numpy.linspace(0, self.H, num=numpy.math.ceil(self.n_elem/2+1))    
       # xv, yv = numpy.meshgrid(x, y)
        
        self.Elem = numpy.zeros((self.n_elem,8))
       
        if(self.n_elem==1):
                    print('!NEED TO BE DEBUGGED!')
                    print('!NEED TO BE DEBUGGED!')            
                    print('DEbug line 41: y coordinates should be next to their corresponding x, e.g., x1,y1, and not after the last x')
                    for i in range(2):
                        self.Elem[0,i] = x[i]
                        self.Elem[0,i+2] = x[1-i] 
                   
                        # y coordinates
                        self.Elem[0,i+4] = y[0]                   
                        self.Elem[0,i+6] = y[1]                    
                  # 
        
        if(self.n_elem>1):
       #For now, computed only for 4 elements 2 x 2 grid 
       # For finer grids (e.g, 4x4) I still have to improve this part
       # passing throug the elemnts
            n_row = 2#2 # number of divisions along y axis
            n_hor = 2#2 # number of divisions along x axis
            path = self.n_elem
            path = numpy.math.ceil(path/n_row) 
            for row in range(n_row):
                c = 0 # counter
                for k in range(path): 
                    for i in range(2):
                        self.Elem[k+n_hor*row,i] = x[i+c]
                        self.Elem[k+n_hor*row,i+2] = x[1-i+c]                    
                        # y coordinates
                        self.Elem[k+n_hor*row,i+4] = y[row]                   
                        self.Elem[k+n_hor*row,i+6] = y[row+1]                    
                    c = c+1  
                  # 
                   
        print(self.Elem)
        return 0

    
    def Compute_Shape_Functions(self):
        
       self.N = sym.zeros(2, 8)
      
       self.psi = sym.Symbol('psi')
       self.eta = sym.Symbol('eta')

       # shape functions
       self.N1 = 1./4.*(1. - self.psi)*(1. - self.eta)
       self.N2 = 1./4.*(1. + self.psi)*(1. - self.eta)
       self.N3 = 1./4.*(1. + self.psi)*(1. + self.eta)
       self.N4 = 1./4.*(1. - self.psi)*(1. + self.eta)
       
       self.N[0,0] = self.N1; self.N[0,2] = self.N2; self.N[0,4] = self.N3; self.N[0,6] = self.N4;
       self.N[1,1] = self.N1; self.N[1,3] = self.N2; self.N[1,5] = self.N3; self.N[1,7] = self.N4;

       # derivatives w.r.t. psi and eta
       
       self.N1_psi = sym.diff(self.N1, self.psi)
       self.N1_eta = sym.diff(self.N1, self.eta)
       self.N2_psi = sym.diff(self.N2, self.psi)
       self.N2_eta = sym.diff(self.N2, self.eta)
       self.N3_psi = sym.diff(self.N3, self.psi)
       self.N3_eta = sym.diff(self.N3, self.eta)       
       self.N4_psi = sym.diff(self.N4, self.psi)
       self.N4_eta = sym.diff(self.N4, self.eta)   
       return 

    def Compute_J(self,i):
        
       # Interpolation geometry of element 
       J = numpy.zeros((2,2))
       J_inv = numpy.zeros((2,2))

       
       x = self.Elem[i,0]*self.N1 + self.Elem[i,1]*self.N2 + self.Elem[i,2]*self.N3 + self.Elem[i,3]*self.N4
       y = self.Elem[i,4]*self.N1 + self.Elem[i,5]*self.N2 + self.Elem[i,6]*self.N3 + self.Elem[i,7]*self.N4
       
       # Entries of the Jacobian
       x_psi = sym.diff(x, self.psi) 
       
       x_eta = sym.diff(x, self.eta) 
       
       y_psi = sym.diff(y, self.psi) 
       
       y_eta = sym.diff(y, self.eta)     
       
       det_J = x_psi*y_eta - (y_psi*x_eta)
       
       # J = |x,psi y,psi|
       #     |x,eta y,eta|
       J[0,0] = x_psi; J[0,1] = y_psi
       J[1,0] = x_eta; J[1,1] = y_eta
       
       J_inv = numpy.linalg.inv(J)
       
        # J^(-1) = |psi,x eta,x|
       #           |psi,y eta,y|      

      # print(J_inv)
       return J_inv, det_J

    def Compute_B_matrix(self,J_inv):
       #self.B = numpy.zeros((3,6))
  
       # array does not work with symbolic math
       B = sym.zeros(3, 8)#Matrix([[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]]);
        
       # J^(-1) = |x,psi y,psi|
       #          |x,eta y,eta|   
       
       self.N1_x = self.N1_psi * J_inv[0,0] + self.N1_eta * J_inv[1,0]
       self.N2_x = self.N2_psi * J_inv[0,0] + self.N2_eta * J_inv[1,0]     
       self.N3_x = self.N3_psi * J_inv[0,0] + self.N3_eta * J_inv[1,0]
       self.N4_x = self.N4_psi * J_inv[0,0] + self.N4_eta * J_inv[1,0]        
       
       self.N1_y = self.N1_psi * J_inv[0,1] + self.N1_eta * J_inv[1,1]
       self.N2_y = self.N2_psi * J_inv[0,1] + self.N2_eta * J_inv[1,1]     
       self.N3_y = self.N3_psi * J_inv[0,1] + self.N3_eta * J_inv[1,1]
       self.N4_y = self.N4_psi * J_inv[0,1] + self.N4_eta * J_inv[1,1]    
       
       #print(self.N1_x)
       
       # loading strain-disp matrix
       B[0,0] = self.N1_x; B[0,2] = self.N2_x; B[0,4] = self.N3_x; B[0,6] = self.N4_x;
       B[1,1] = self.N1_y; B[1,3] = self.N2_y; B[1,5] = self.N3_y; B[1,7] = self.N4_y; 
       B[2,0] = self.N1_y; B[2,2] = self.N2_y; B[2,4] = self.N3_y; B[2,6] = self.N4_y;       
       B[2,1] = self.N1_x; B[2,3] = self.N2_x; B[2,5] = self.N3_x; B[2,7] = self.N4_x; 
       


       return B
            
 #       return p,fb    


    def Compute_C(self):

     #Material Matrix
      self.C = numpy.zeros((3,3))      
      self.C[0,0] = 1.-self.v;             self.C[1,1] = 1.-self.v;   self.C[2,2] = (1.-2.*self.v)/2.; 
      self.C[0,1] = self.v; self.C[1,0] = self.v; 
      self.C = self.E/((1.+self.v)*(1.-2.*self.v))*self.C
      return 0      
    
 
   # def Clear_Matrices(self):
    #   J_inv = numpy.zeros((2,2))  
     #  B =  numpy.zeros((3,8))  
      # det_J = 0
       #ke = numpy.zeros((8,8))
       #aux_ke = numpy.zeros((8,3))  
       #return J_inv,B,det_J,ke, aux_ke 
       
    def Calculate_ke(self,i):
       aux_ke  = sym.zeros(8,3)
      # I have to be aware that I'm mixing here sympy with numpy operations and it's really not consistent, so
      # some operations allowed with numpy funcitons are not allowed with sympy and viceversa
      # J_inv,B,det_J,ke, aux_ke  = self.Clear_Matrices()
       
       J_inv, det_J = self.Compute_J(i)
       B = self.Compute_B_matrix(J_inv)
       aux_ke = B.T*self.C
       ke = aux_ke*B*det_J

       
       #Integration Gauss quadrature (2 x 2 for quads)
       psi_1 = -1./numpy.sqrt(3); psi_2 = 1./numpy.sqrt(3); eta_1 = -1./numpy.sqrt(3); eta_2 = 1./numpy.sqrt(3);
       w1 = 1.; w2 = 1.;
       
       ke = ke.evalf(subs={self.psi:psi_1, self.eta:eta_1})*w1*w1 + ke.evalf(subs={self.psi:psi_1, self.eta:eta_2})*w1*w2 +\
       ke.evalf(subs={self.psi:psi_2, self.eta:eta_1})*w2*w1 + ke.evalf(subs={self.psi:psi_2, self.eta:eta_2})*w2*w2
       
       # ke was verified and is symm!

       #print("ke ==",ke[0,0])
       return ke
    

    
    def Define_Local_Load_Vector(self,i):
        
       psi_1 = -1./numpy.sqrt(3); psi_2 = 1./numpy.sqrt(3); 
       w1 = 1.; w2 = 1.;       
       J_inv, det_J = self.Compute_J(i)  
       aux_N = sym.zeros(1,8)             
       f = sym.zeros(1,8)       
       aux_N = (self.N.T)*self.q_vector*det_J              
      # f = aux_N.evalf(subs={self.eta:1,self.psi:psi_1})*w1 + aux_N.evalf(subs={self.eta:1,self.psi:psi_2})*w2
       f[5]=-1
       return f
   
    

    def Compute_Kg(self):
        self.K = numpy.zeros((self.dof,self.dof)) 

        for i in range(self.n_elem):
            ke = self.Calculate_ke(i)
            self.K = ke

        return   self.K
 
    def Define_Global_Load_Vector(self): 
       self.Q = numpy.zeros(self.dof) 
       for i in range(self.n_elem):
            q_local = self.Define_Local_Load_Vector(i)
            self.Q = q_local     
       #print("load   ", self.Q)
       return    

    def Compute_Kred(self):
       self.restrictions  =  4   
       self.dof_block= numpy.zeros(self.dof-self.restrictions)        
       self.Kred = numpy.zeros((self.dof-self.restrictions,self.dof-self.restrictions))    
       
       # Id matrix containing restrictions
       self.dof_block = [2,3,4,5]
       
       for i in range(self.dof-self.restrictions):
           for j in range(self.dof-self.restrictions):
                       
                   self.Kred[i,j] = self.K[self.dof_block[i],self.dof_block[j]]
                   
       
       return self.Kred    

    def Compute_Qred(self):
      
       self.Qred = numpy.zeros(self.dof-self.restrictions)    
       
       for i in range(self.dof-self.restrictions):   
                   
                   self.Qred[i] = self.Q[self.dof_block[i]]        
       return self.Qred       

    def Solve_System(self):
        self.u = numpy.dot(numpy.linalg.inv(self.Kred),self.Qred)
        return self.u

   # def Update_positions(self):                  

    def Solve(self):
        self.Compute_Geometry()
        self.Compute_Shape_Functions()
        self.Compute_C()
        self.Compute_Kg()
        self.Compute_Kred()  
        self.Define_Global_Load_Vector()
        self.Compute_Qred()
        u = self.Solve_System()
        print(u)
       
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
 number_elements = 1
 Length = 3.0
 Height = 0.5
 q_load = 1
 E_modulus = 30000
 v_poisson = 0.0
 args = (Length,Height,q_load,E_modulus,v_poisson)   
 
 Beam = FEM_2D_element(number_elements,*args) 
 Beam.Solve() 

 #pf, fbf = pressure.Compute_fb()   
 #plt.plot(pf, fbf)
 
 
main() 