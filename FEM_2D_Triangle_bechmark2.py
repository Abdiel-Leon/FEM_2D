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
Benchmark adopted from Finite Element Method SME 3033 notes by Dr. Nazri Kamsah 
"""
import numpy
#from operator import itemgetter
#import matplotlib.pyplot as plt

class FEM_2D_element_quad:

    def __init__(self, n_elem, *args):
        self.L = args[0]
        self.H = args[1]
        self.q = args[2]
        self.E = args[3]
        self.v = args[4]
        self.n_elem = n_elem
        self.Le =(self.L)/(self.n_elem/2) # number of horizontal divisions (2 triangoles make 1 quad)
        self.dof = self.n_elem*2 +4
        self.t = 0.5 #  for plane-strain problems, thickness asssumed as 1 m.
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

      return 

    def Compute_ke(self,x1,y1,x2,y2,x3,y3):
        B = numpy.zeros((3, 6))            

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


    def Compute_Kg(self):
        self.K = numpy.zeros((self.dof,self.dof))
        Id = numpy.zeros((self.n_elem*6)) # nelem*6 (6 dof per elem)
        index =0     
        even = 0 # even elements n elem = 2,4,6...
        odd = 0
        for e in range(self.n_elem):              
            # even elements
            if (e%2==0):
                even = even + 1
                # geometry of the elements
                x3 = self.Le*even; x1 = self.Le*even-self.Le; x2 = self.Le*even-self.Le;
                y3 = self.H;            y1 = self.H;       y2 = 0;                
                ke = self.Compute_ke(x1,y1,x2,y2,x3,y3)             
                # Id_matrix
                for i in range(6):
                       index = 6*e + i
                       Id[index] = 2*e + i  

            else:
                # odd elements
                odd = odd + 1
                # geometry
                x1 = self.Le*odd; x2 = self.Le*odd; x3 = self.Le*odd-self.Le;
                y1 = 0;            y2 = self.H;       y3 = 0;
                ke = self.Compute_ke(x1,y1,x2,y2,x3,y3)
                
                #Id_matrix              
                for i in range(3):
                       index = 6*e + 2*i
                       Id[index] = 2*e + 4 - 2*i  
                       Id[index+1] = 2*e + 5 - 2*i                        

            for i in range(6):
                for j in range(6):
                    self.K[int(Id[i + 6*e]),int(Id[j + 6*e])] = ke[i,j] + self.K[int(Id[i + 6*e]),int(Id[j + 6*e])]
        return   self.K




    def check_symmetric(self, tol=1e-6):

        return numpy.all(numpy.abs(self.K-self.K.T) < tol)


    def Define_Global_Load_Vector(self):
       self.Q = numpy.zeros(self.dof)
       dof_load = 5
       self.Q[dof_load] = self.q

       return


    def Compute_Kred(self):
       self.restrictions  =  5 # number of restrictions
       self.dof_fix= numpy.zeros(self.restrictions)
       self.dof_full= numpy.zeros(self.dof)
       self.dof_block= numpy.zeros(self.dof-self.restrictions)

       
       self.Kred = numpy.zeros((self.dof-self.restrictions,self.dof-self.restrictions))
       
       #List of dof
       for i in range(self.dof):
           self.dof_full[i] = i
           
       #fixed dof (indicate the number of the fixed dof)
       self.dof_fix = [0,1,2,3,7]
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
       
        print("u =", u)

        return u





    #def Plot_displacements(self):
        # plots original and deformed bar lengths
       # plt.ylim(0.5,1.1)
        #y = [1] * (self.n + 1)
        #yupd = [1.05] * (self.n + 1)
        #plt.plot(self.x0,y,label="Original bar",marker=8)
        #plt.plot(self.xupd,yupd, label="Stressed bar",marker=8)
        #plt.legend(loc="lower right")

def main():
 number_elements, Length, Height, q_load = 2,3,2,-1000
 E_modulus,v_poisson = 30e6, 0.25
 args = (Length,Height,q_load,E_modulus,v_poisson)

 Beam = FEM_2D_element_quad(number_elements,*args)
 u = Beam.Solve()
 #with open('u_4ele.txt', 'w') as f:
  #  f.write(str(u))
   # f.close

 #pf, fbf = pressure.Compute_fb()
 #plt.plot(pf, fbf)


main()