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

class K_analythical:
    
    def __init__(self):
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

    def Compute_Ke(self):
        
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
 
    
def main():
    Stiffness = K_analythical()
    K= Stiffness.Compute_Ke()
    print("analythical",K)

main()    