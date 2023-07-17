# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:04:45 2023
@author: abdiel

### Code to produce FEM analyses in 2D for very simple geometries###
"""
import numpy
import sympy as sym
#import matplotlib.pyplot as plt

class FEM_2D_element_quad:

    def __init__(self, n_elem, *args):
        self.L = args[0]
        self.H = args[1]
        self.q = args[2]
        self.E = args[3]
        self.v = args[4]
        self.n_elem = n_elem
        self.Le =(self.L)/(self.n_elem) # number of horizontal divisions (2 triangoles make 1 quad)
        self.dof = (self.n_elem + 1)*4
        self.t = 0.1 #  for plane-strain problems, thickness asssumed as 1 m.
        #self.q_vector = sym.zeros(2, 1) # traction vector in Pa
       # self.q_vector[1] = self.q
            

    def Compute_C(self):

     #Material Matrix (plane strain)
      #self.C = numpy.zeros((3,3))
      #self.C[0,0] = 1.-self.v;             self.C[1,1] = 1.-self.v;   self.C[2,2] = (1.-2.*self.v)/2.;
      #self.C[0,1] = self.v; self.C[1,0] = self.v;
      #self.C = self.E/((1.+self.v)*(1.-2.*self.v))*self.C

     #Material Matrix (plane stress)
      self.C = numpy.zeros((3,3))      
      self.C[0,0] = 1.;             self.C[1,1] = 1.;   self.C[2,2] = (1.-self.v)/2.; 
      self.C[0,1] = self.v; self.C[1,0] = self.v; 
      self.C = self.E/(1-self.v*self.v)*self.C      
  #    print(self.C)
      return 

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

    def Compute_J(self,x1,y1,x2,y2,x3,y3,x4,y4):
        
       # Interpolation geometry of element 
       J = numpy.zeros((2,2))
       J_inv = numpy.zeros((2,2))

       
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


       return J_inv, det_J

    def Compute_B_matrix(self,J_inv):
       #self.B = numpy.zeros((3,6))
  
       # array does not work with symbolic math
       B = sym.zeros(3, 8)
        
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
       
       
       # loading strain-disp matrix
       B[0,0] = self.N1_x; B[0,2] = self.N2_x; B[0,4] = self.N3_x; B[0,6] = self.N4_x;
       B[1,1] = self.N1_y; B[1,3] = self.N2_y; B[1,5] = self.N3_y; B[1,7] = self.N4_y; 
       B[2,0] = self.N1_y; B[2,2] = self.N2_y; B[2,4] = self.N3_y; B[2,6] = self.N4_y;       
       B[2,1] = self.N1_x; B[2,3] = self.N2_x; B[2,5] = self.N3_x; B[2,7] = self.N4_x; 
       


       return B


    def Compute_ke(self,x1,y1,x2,y2,x3,y3,x4,y4):
       self.Compute_Shape_Functions()
       J_inv, det_J = self.Compute_J(x1,y1,x2,y2,x3,y3,x4,y4)
       B = self.Compute_B_matrix(J_inv)
        
       aux_ke = B.T*self.C
       ke = self.t*aux_ke*B*det_J

       
       #Integration Gauss quadrature (2 x 2 for quads)
       psi_1 = -1./numpy.sqrt(3); psi_2 = 1./numpy.sqrt(3); eta_1 = -1./numpy.sqrt(3); eta_2 = 1./numpy.sqrt(3);
       w1 = 1.; w2 = 1.;
  # too slow!     
      # ke = ke.evalf(subs={self.psi:psi_1, self.eta:eta_1})*w1*w1 + ke.evalf(subs={self.psi:psi_1, self.eta:eta_2})*w1*w2 +\
       #ke.evalf(subs={self.psi:psi_2, self.eta:eta_1})*w2*w1 + ke.evalf(subs={self.psi:psi_2, self.eta:eta_2})*w2*w2
      
       ke_num = sym.lambdify([self.psi,self.eta],ke,"numpy")
       ke_num = ke_num(psi_1,eta_1)*w1*w1  + ke_num(psi_1,eta_2)*w1*w2 + ke_num(psi_2,eta_1)*w2*w1 +  ke_num(psi_2,eta_2)*w2*w2  
       return ke_num 


    def Compute_Kg(self):
        self.K = numpy.zeros((self.dof,self.dof))
        Id = numpy.zeros((self.n_elem*8)) # 8 dof per element    
        elem = 0 
        #loop over elemetns
        for e in range(self.n_elem):
            elem = elem + 1
            # geometry of the elements
            x1 = self.Le*elem-self.Le; x2 = self.Le*elem; x3 = self.Le*elem; x4 = self.Le*elem-self.Le;
            y1 = 0;                    y2 = 0;            y3 = self.H;       y4 = self.H;  
            
            ke = self.Compute_ke(x1,y1,x2,y2,x3,y3,x4,y4) 


            for i in range(4):
               index = 8*e + i
               Id[index] = 2*e +i 
               Id[index + 4] = (self.dof - 4) + i -2*e                   

#
            for i in range(8):
                for j in range(8):
                    self.K[int(Id[i + 8*e]),int(Id[j + 8*e])] = ke[i,j] + self.K[int(Id[i + 8*e]),int(Id[j + 8*e])]
        return   self.K




    def check_symmetric(self, tol=1e-4):

        return numpy.all(numpy.abs(self.K-self.K.T) < tol)


    def Define_Global_Load_Vector(self):
       self.Q = numpy.zeros(self.dof)
       dof_load = self.n_elem*2+3
       self.Q[dof_load] = self.q

       #for i in range(self.n_elem):
        #    q_local = self.Define_Local_Load_Vector(i)
         #   self.Q = q_local
       #print("load   ", self.Q)
       return


    def Compute_Kred(self):
       self.restrictions  =  4 # number of restrictions
       self.dof_fix= numpy.zeros(self.restrictions)
       self.dof_full= numpy.zeros(self.dof)
       self.dof_block= numpy.zeros(self.dof-self.restrictions)      
       self.Kred = numpy.zeros((self.dof-self.restrictions,self.dof-self.restrictions))
       
       #List of dof
       for i in range(self.dof):
           self.dof_full[i] = i
           
       #fixed dof (number of the fixed dof)
       self.dof_fix = [0,1,(self.n_elem)*4+2,(self.n_elem)*4+3]
       # free dof using list comprehension
       self.dof_block =  [x for x in self.dof_full if x not in self.dof_fix]
       for i in range(self.dof-self.restrictions):
           for j in range(self.dof-self.restrictions):
                   self.Kred[i,j] = self.K[int(self.dof_block[i]),int(self.dof_block[j])]
                  
       return self.Kred

    def Compute_Qred(self):       
       self.Qred = numpy.zeros(self.dof-self.restrictions)

       for i in range(self.dof-self.restrictions):
                   self.Qred[i] =self.Q[int(self.dof_block[i])]
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
        self.Define_Global_Load_Vector()
        Qred= self.Compute_Qred()    
        u = self.Solve_System(Kred,Qred)
        print("u_tip-point-load", ["3.35e-4 m"])   
        u_sci = "{:e}".format(u[self.n_elem*2 + 3])
        print("u =", u_sci)

        return #u





    #def Plot_displacements(self):
        # plots original and deformed bar lengths
       # plt.ylim(0.5,1.1)
        #y = [1] * (self.n + 1)
        #yupd = [1.05] * (self.n + 1)
        #plt.plot(self.x0,y,label="Original bar",marker=8)
        #plt.plot(self.xupd,yupd, label="Stressed bar",marker=8)
        #plt.legend(loc="lower right")

def main():
 number_elements, Length, Height, q_load = 50,0.5,0.06,-1000
 E_modulus,v_poisson = 69e9, 0.0
 args = (Length,Height,q_load,E_modulus,v_poisson)

 Beam = FEM_2D_element_quad(number_elements,*args)
 u = Beam.Solve()
 #with open('u_4ele.txt', 'w') as f:
  #  f.write(str(u))
   # f.close

 #pf, fbf = pressure.Compute_fb()
 #plt.plot(pf, fbf)


main()