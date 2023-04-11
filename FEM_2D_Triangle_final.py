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

class FEM_2D_element_triangle:

    def __init__(self, n_elem, *args):
        self.L = args[0]
        self.H = args[1]
        self.q = args[2]
        self.E = args[3]
        self.v = args[4]
        self.n_elem = n_elem
        self.Le =(self.L)/(self.n_elem/2)
        self.dof = self.n_elem*2 +4
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
        ke_1 = A*numpy.matmul(aux_ke_1,B)
        return ke_1


    def Compute_Kg(self):
        self.K = numpy.zeros((self.dof,self.dof))
        Id = numpy.zeros((120))
        #Id = [2,3,24,25,0,1,22,23,0,1,24,25,4,5,26,27,2,3,24,25,2,3,26,27,\
         #     6,7,28,29,4,5,26,27,4,5,28,29,8,9,30,31,6,7,28,29,6,7,30,31,\
          #   10,11,32,33,8,9,30,31,8,9,32,33,12,13,34,35,10,11,32,33,10,11,34,35,\
           #  14,15,36,37,12,13,34,35,12,13,36,37,16,17,38,39,14,15,36,37,14,15,38,39,\
            #18,19,40,41,16,17,38,39,16,17,40,41,20,21,42,43,18,19,40,41,18,19,42,43 ]

        Id = [2,3,4,5,0,1,6,7,0,1,4,5,8,9,10,11,2,3,4,5,2,3,10,11,12,13,14,15,8,9,10,11,8,9,14,15,\
              16,17,18,19,12,13,14,15,12,13,18,19,20,21,22,23,16,17,18,19,16,17,22,23,24,25,26,27,20,21,\
              22,23,20,21,26,27,28,29,30,31,24,25,26,27,24,25,30,31,32,33,34,35,28,29,30,31,28,29,34,35,\
              36,37,38,39,32,33,34,35,32,33,38,39,40,41,42,43,36,37,38,39,36,37,42,43]

        even = 0
        odd = 0
        for e in range(self.n_elem):
            if (e==0 or e==2 or e==4 or e==6 or e==8 or e==10 or e==12 or e==14\
                or e==16 or e==18):
                even = even + 1
                x1 = self.Le*even; x2 = self.Le*even; x3 = self.Le*even-self.Le;
                y1 = 0;            y2 = self.H;       y3 = 0;
                ke = self.Compute_ke(x1,y1,x2,y2,x3,y3)
                print(e,ke)

            if (e==1 or e==3 or e==5 or e==7 or e==9 or e==11 or e==13 or e==15\
                or e==17 or e==19):
                odd = odd + 1
                x3 = self.Le*odd; x1 = self.Le*odd-self.Le; x2 = self.Le*odd-self.Le;
                y3 = self.H;            y1 = self.H;       y2 = 0;
                ke = self.Compute_ke(x1,y1,x2,y2,x3,y3)
               # print(e,ke)

            for i in range(6):
                for j in range(6):
                    #print(e,Id[i + 6*e],Id[j + 6*e])
                    self.K[Id[i + 6*e],Id[j + 6*e]] = self.t*ke[i,j] + self.K[Id[i + 6*e],Id[j + 6*e]]

        return   self.K




    def check_symmetric(self, tol=1e-6):

        return numpy.all(numpy.abs(self.K-self.K.T) < tol)


    def Define_Global_Load_Vector(self):
       self.Q = numpy.zeros(self.dof)
       self.Q[self.dof-1] = self.q

       #for i in range(self.n_elem):
        #    q_local = self.Define_Local_Load_Vector(i)
         #   self.Q = q_local
       #print("load   ", self.Q)
       return

    def Compute_Kred(self):
       self.restrictions  =  4
       self.dof_block= numpy.zeros(self.dof-self.restrictions)
       self.Kred = numpy.zeros((self.dof-self.restrictions,self.dof-self.restrictions))

       # Id matrix with free dof
       self.dof_block = [2,3,4,5,8,9,10,11,12,13,14,15,16,17,18,19,20,\
                         21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43]

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
                   self.Qred[i] =self.Q[self.dof_block[i]]
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
        print(u)

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
 number_elements = 20
 Length = 0.5
 Height = 0.06
 q_load = -1000
 E_modulus = 69e9
 v_poisson = 0.33
 args = (Length,Height,q_load,E_modulus,v_poisson)

 Beam = FEM_2D_element_triangle(number_elements,*args)
 u = Beam.Solve()
 #with open('u_4ele.txt', 'w') as f:
  #  f.write(str(u))
   # f.close

 #pf, fbf = pressure.Compute_fb()
 #plt.plot(pf, fbf)


main()