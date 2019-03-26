# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 21:51:59 2019

@author: Emanuel
"""
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 09:54:54 2018
@author: Emanuel
"""

import numpy
import fp
import mf

l = 0;
Sp = 0;

e = .62; #porosity unitless
r = .22e-6; #mean pore radius meters.
t  = ((2-e)**2)/e; # unitless 
delta = 126e-6 ;#Using PTFE - 1. Thickness in meters.

#Initial bulk fluid conditions
Tbf = 60+273.15;
Tbp = 20+273.15;
k_ptfe = .19

Dh = 3e-3
vf = 2.13
vp = 3.5
L_both = 1e-1
n,tmf,tmp = mf.Membrane(Tbf,Tbp,e,r,t,delta,Sp,k_ptfe,Dh,L_both,vf,vp)
print("The mass flux  ", n)
print("The feed side temperature", tmf);
print("The permeate side temperature",tmp);

