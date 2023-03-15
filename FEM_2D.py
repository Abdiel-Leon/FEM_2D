#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:04:45 2023
@author: abdiel

### Code to produce FEM analyses in 2D for very simple geometries###
"""
import numpy
import sympy as sym
import matplotlib.pyplot as plt

class FEM_2D_element:
    
    def __init__(self,n_elem,*args):
        self.L = args[0]
        self.H = args[1]
        self.q = args[2]
        self.E = args[3]
        self.v = args[4]
        self.n_elem = n_elem
        
    
    def Compute_Geometry(self):
        x = numpy.zeros(int(self.n_elem + 1))
        y = numpy.zeros(int(self.n_elem + 1))
        x = numpy.linspace(0, self.L, num=numpy.math.ceil(self.n_elem/2+1))   
        y = numpy.linspace(0, self.H, num=numpy.math.ceil(self.n_elem/2+1))    
       # xv, yv = numpy.meshgrid(x, y)
        
        self.Elem = numpy.zeros((self.n_elem,8))
       
        if(self.n_elem==1):
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
        return 

    
    def Compute_Shape_Functions(self):
       psi = sym.Symbol('psi')
       eta = sym.Symbol('eta')

       # shape functions
       self.N1 = 1./4.*(1. - psi)*(1. - eta)
       self.N2 = 1./4.*(1. + psi)*(1. - eta)
       self.N3 = 1./4.*(1. + psi)*(1. + eta)
       self.N4 = 1./4.*(1. - psi)*(1. + eta)

       # derivatives w.r.t. psi and eta
       
       self.N1_psi = sym.diff(self.N1, psi)
       self.N1_eta = sym.diff(self.N1, eta)
       self.N2_psi = sym.diff(self.N2, psi)
       self.N2_eta = sym.diff(self.N2, eta)
       self.N3_psi = sym.diff(self.N3, psi)
       self.N3_eta = sym.diff(self.N3, eta)       
       self.N4_psi = sym.diff(self.N4, psi)
       self.N4_eta = sym.diff(self.N4, eta)  
      # print(self.N1_psi)        
      # print(self.N1_eta)        

       #N1.evalf(subs={psi:0.3, eta:0}) # evaluation of shape function at selected points
       return 

    def Compute_J(self,i):
        
       # Interpolation geometry of element 
       self.x = self.Elem[i,0]*self.N1 + self.Elem[i,1]*self.N2 + self.Elem[i,2]*self.N3 + self.Elem[i,3]*self.N4
       self.y = self.Elem[i,4]*self.N1 + self.Elem[i,5]*self.N2 + self.Elem[i,6]*self.N3 + self.Elem[i,7]*self.N4
       
       # Entries of the Jacobian
       self.x_psi = self.Elem[i,0]*self.N1_psi + self.Elem[i,1]*self.N2_psi + \
       self.Elem[i,2]*self.N3_psi + self.Elem[i,3]*self.N4_psi
       
       self.x_eta = self.Elem[i,0]*self.N1_eta + self.Elem[i,1]*self.N2_eta + \
       self.Elem[i,2]*self.N3_eta + self.Elem[i,3]*self.N4_eta
       
       self.y_psi = self.Elem[i,4]*self.N1_psi + self.Elem[i,5]*self.N2_psi + \
       self.Elem[i,6]*self.N3_psi + self.Elem[i,7]*self.N4_psi
       
       self.y_eta = self.Elem[i,4]*self.N1_eta + self.Elem[i,5]*self.N2_eta + \
       self.Elem[i,6]*self.N3_eta + self.Elem[i,7]*self.N4_eta     
       
       self.det_J = self.x_psi*self.y_eta - (self.y_psi*self.x_eta)
       
       return

            
 #       return p,fb    


  # def Compute_C(self):

            
     #   return p,fb       
    

 #  def Calculate_ke(self):

    
  # def Calculate_K(self):
     #   self.Calculate_ke()
     #   self.K = numpy.zeros((self.n + 1,self.n + 1))        

    
  # def Define_Load_Vector(self):
        # global load vector
     #   self.F = numpy.zeros(self.n) 
     #   self.F[self.n-1] = 5000        
      #  return self.F
    
   # def Calculate_Kred(self):
      
                
   # def Solve_System(self):
     #   self.u = numpy.dot(numpy.linalg.inv(self.Kred),self.F)
       # return self.u

   # def Update_positions(self):
#
    def Solve_FEM(self):
        self.Compute_Geometry()
        self.Compute_Shape_Functions()       
        
        for i in range(self.n_elem):
                self.Compute_J(i)
                #print("elem = ", i,"x = ", x)
    
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
 v_poisson = 0.3
 args = (Length,Height,q_load,E_modulus,v_poisson)   
 
 Beam = FEM_2D_element(number_elements,*args) 
 Beam.Solve_FEM() 

 #pf, fbf = pressure.Compute_fb()   
 #plt.plot(pf, fbf)
 
 
main() 