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
       # self.q_vector = sym.zeros(2, 1) # traction vector in Pa       
      #  self.q_vector[1] = self.q 
      
        self.a = 3
        self.b= 0.5
        self.E = 30000
        self.h = 1
        self.v = 0
        self.c = self.b/self.a
        self.ka = 4*self.c*(1-self.v)
        self.kb = (2/self.c)*(1-2*self.v)
        self.kc = (4/self.c)*(1-self.v)
        self.kd = (2*self.c)*(1-2*self.v)      
        
        
    
    
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

    def Compute_J(self):
        
       # Interpolation geometry of element 
       J = numpy.zeros((2,2))
       J_inv = numpy.zeros((2,2))

       i =0
       # nodal coordinates
       x1 = 3; x2 = 0; x3 =0; x4 = 3;
       y1 = 0.5; y2 = 0.5; y3 =0.; y4 = 0.;
       
       x = x1*self.N1 + x2*self.N2 + x3*self.N3 + x4*self.N4
       y = y1*self.N1 + y2*self.N2 + y3*self.N3 + y4*self.N4       
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
       B = sym.zeros(3, 8)#
        
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

     #Material Matrix (plane strain)
      self.C = numpy.zeros((3,3))      
      self.C[0,0] = 1.-self.v;             self.C[1,1] = 1.-self.v;   self.C[2,2] = (1.-2.*self.v)/2.; 
      self.C[0,1] = self.v; self.C[1,0] = self.v; 
      self.C = self.E/((1.+self.v)*(1.-2.*self.v))*self.C
      return 0      
    
 

       
    def Calculate_ke(self):
       aux_ke  = sym.zeros(8,3)
      # I have to be aware that I'm mixing here sympy with numpy operations and it's really not consistent, so
      # some operations allowed with numpy funcitons are not allowed with sympy and viceversa
      # J_inv,B,det_J,ke, aux_ke  = self.Clear_Matrices()
       
       J_inv, det_J = self.Compute_J()
       B = self.Compute_B_matrix(J_inv)
       aux_ke = B.T*self.C
       ke = self.t*aux_ke*B*det_J

       
       #Integration Gauss quadrature (2 x 2 for quads)
       psi_1 = -1./numpy.sqrt(3); psi_2 = 1./numpy.sqrt(3); eta_1 = -1./numpy.sqrt(3); eta_2 = 1./numpy.sqrt(3);
       w1 = 1.; w2 = 1.;
       
       ke = ke.evalf(subs={self.psi:psi_1, self.eta:eta_1})*w1*w1 + ke.evalf(subs={self.psi:psi_1, self.eta:eta_2})*w1*w2 +\
       ke.evalf(subs={self.psi:psi_2, self.eta:eta_1})*w2*w1 + ke.evalf(subs={self.psi:psi_2, self.eta:eta_2})*w2*w2
       
       # ke was verified and is symm!

      # print(self.C)
       return ke
    

    
    def Define_Local_Load_Vector(self):
        
       #psi_1 = -1./numpy.sqrt(3); psi_2 = 1./numpy.sqrt(3); 
       #w1 = 1.; w2 = 1.;       
       #J_inv, det_J = self.Compute_J()  
       #aux_N = sym.zeros(1,8)             
       f = sym.zeros(1,8)       
      # aux_N = (self.N.T)*self.q_vector*det_J              
      # f = aux_N.evalf(subs={self.eta:1,self.psi:psi_1})*w1 + aux_N.evalf(subs={self.eta:1,self.psi:psi_2})*w2
       f[5]=-1
       return f
   
 
    def Compute_KG(self):
        
        K = numpy.zeros((8,8)) 
        
        k00 =self.ka + self.kb; k01 = 1.5; k02 = -self.ka + self.kb/2; k03 = 6*self.v- 1.5; k04 = -self.ka/2 -self.kb/2; k05 = -1.5; k06 = self.ka/2 - self.kb; k07 = 1.5 - 6*self.v;
        k10 =1.5; k11 = self.kc + self.kd; k12 = 1.5-6*self.v; k13 = self.kc/2-self.kd; k14 = -1.5; k15 = -self.kc/2-self.kd/2; k16 = 6*self.v - 1.5; k17 = -self.kc + self.kd/2;
        k20 = -self.ka + self.kb/2; k21 = 1.5-6*self.v; k22 = self.ka + self.kb; k23 = -1.5; k24 = self.ka/2 - self.kb; k25 = 6*self.v - 3/2; k26 = -self.ka/2 - self.kb/2; k27 =1.5;
        k30 = 6*self.v - 3/2; k31 = self.kc/2 - self.kc; k32 = -1.5; k33 = self.kc + self.kd; k34 = 3/2 -6*self.v; k35 =-self.kc + self.kd/2; k36 = 1.5; k37 = -self.kc/2 - self.kd/2;
        k40 = -self.ka/2 -self.kb/2; k41 = -1.5; k42 = self.ka/2 -self.kb; k43 = 3/2 - 6*self.v; k44 = self.ka + self.kb; k45 =1.5; k46 = -self.ka + self.kb/2; k47 = 6*self.v - 1.5;
        k50 = -1.5; k51 = -self.kc/2 -self.kd/2; k52 = 6*self.v - 3/2; k53 = -self.kc + self.kd/2; k54 = 1.5; k55 = self.kc + self.kd; k56 =3/2 - 6*self.v; k57 =self.kc/2 - self.kd;
        k60 = self.ka/2 - self.kb; k61 =6*self.v - 1.5; k62 = -self.ka/2 -self.kb/2; k63 =1.5; k64 =-self.ka + self.kb/2;k65 = 3/2 - 6*self.v; k66 =self.ka + self.kb; k67 =-1.5;
        k70 =1.5-6*self.v; k71 =-self.kc + self.kd/2; k72 = 1.5; k73 = -self.kc/2 -self.kd/2; k74 = 6*self.v -1.5; k75 = self.kc/2 - self.kd; k76 = -1.5; k77 = self.kc + self.kd;
        
        K[0,0] = k00; K[0,1] = k01; K[0,2] = k02; K[0,3] = k03; K[0,4] = k04; K[0,5] = k05;  K[0,6] = k06; K[0,7] = k07; 
        K[1,0] = k10; K[1,1] = k11; K[1,2] = k12; K[1,3] = k13; K[1,4] = k14; K[1,5] = k15;  K[1,6] = k16; K[1,7] = k17; 
        K[2,0] = k20; K[2,1] = k21; K[2,2] = k22; K[2,3] = k23; K[2,4] = k24; K[2,5] = k25;  K[2,6] = k26; K[2,7] = k27; 
        K[3,0] = k30; K[3,1] = k31; K[3,2] = k32; K[3,3] = k33; K[3,4] = k34; K[3,5] = k35;  K[3,6] = k36; K[3,7] = k37; 
        K[4,0] = k40; K[4,1] = k41; K[4,2] = k42; K[4,3] = k43; K[4,4] = k44; K[4,5] = k45;  K[4,6] = k46; K[4,7] = k47;         
        K[5,0] = k50; K[5,1] = k51; K[5,2] = k52; K[5,3] = k53; K[5,4] = k54; K[5,5] = k55;  K[5,6] = k56; K[5,7] = k57;
        K[6,0] = k60; K[6,1] = k61; K[6,2] = k62; K[6,3] = k63; K[6,4] = k64; K[6,5] = k65;  K[5,6] = k66; K[6,7] = k67; 
        K[7,0] = k70; K[7,1] = k71; K[7,2] = k72; K[7,3] = k73; K[7,4] = k74; K[7,5] = k75;  K[7,6] = k76; K[7,7] = k77;        
        K = self.E*self.h/(12*(1+self.v)*(1-2*self.v))*K
        return K

    

    def Compute_Kg(self):
        self.K = numpy.zeros((self.dof,self.dof)) 


        ke = self.Calculate_ke()
        self.K = ke

        return   self.K
 

    def Compute_Kred(self,K):
       self.restrictions  =  4   
       self.dof_block= numpy.zeros(self.dof-self.restrictions)        
       self.Kred = numpy.zeros((self.dof-self.restrictions,self.dof-self.restrictions))    
       
       # Id matrix containing restrictions
       self.dof_block = [4,5,2,3]
       
       for i in range(self.dof-self.restrictions):
           for j in range(self.dof-self.restrictions):
                       
                   self.Kred[i,j] = K[self.dof_block[i],self.dof_block[j]]
                   
      # print(self.Kred)
       return self.Kred    

    def Compute_Qred(self,f):
      
       self.Qred = numpy.zeros(self.dof-self.restrictions)    
       self.Qred[3] = -1
      # for i in range(self.dof-self.restrictions):   
                   
       #            self.Qred[i] = f[self.dof_block[i]]        
       return self.Qred       

    def Solve_System(self):
        self.u = numpy.dot(numpy.linalg.inv(self.Kred),self.Qred)
        return self.u

   # def Update_positions(self):                  

    def Solve(self):
        #self.Compute_Shape_Functions()
        #self.Compute_C()
       # K = self.Compute_Kg()
#        print("fem", Kg)
        K = self.Compute_KG()
        self.Compute_Kred(K)  
        f = self.Define_Local_Load_Vector()
        self.Compute_Qred(f)
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
 Length = 3
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