# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 22:44:07 2019

MEMBRANE FUNCTIONS


@author: Emanuel
"""

import fp
import numpy

Kb = 1.38064852e-23;
Mwater = 18.01528e-3; #kg/mol
R = 8.31;
P = 10**5;

def get_memconstants(e,r,t):
        B0 = e*(r**2)/(8*t);
        K0 = 2*e*r/(3*t);
        K1 = e/t;
        return B0,K0,K1;


def get_knudsen(T_avg,p_avg,r):
#Returns the knudsen number 
        mfp = Kb*T_avg/(1.41421*p_avg*(2.641e-10)**2);
        Kn  = mfp/(2*r);
        return Kn;
    
def Isotrain(e,ks,T):
    kg = 2.72e-3+7.77e-5*T; #Thermal conductivity of water vapor. 
    km = (1-e)*ks+e*kg;
    return km;

def Isostress(e,ks,T):
    kg = 2.72e-3+7.77e-5*T; #Thermal conductivity of water vapor. 
    km = (e/kg + (1-e)/ks)**-1 ;
    return km;    


def get_Bm(e,r,t,Kn,T_m,vmean,delta):
    B0,K0,K1 = get_memconstants(e,r,t)
    if (Kn >1):        
        Bm = K0*vmean/(R*T_m); 
    elif (Kn < .01):
        Bm = ( e/(t*delta) ) * ((1.895e-5 *T_m**(2.071))/(101.3e3)) *(Mwater*R/T_m);
    else:
        Bm = (  (3*t*delta*vmean)/(2*e*r) + \
              (t*delta/e)*( (101.3e3)/ (1.895e-5 *T_m**(2.071)) ) *(R*T_m/Mwater))**-1;  
    return Bm;


def Membrane(Tbf,Tbp,e,r,t,delta,Sp,ks,Dh,L,vf,vp) :
    Tmf = Tbf ;
    Tmp = Tbp ;
    N = 0 ;
    error  = 1 ; 
    while error > 1e-5:
        tempf = Tmf;
        tempp = Tmp;
        Tmean_feed = (Tbf + Tmf)/2;
        Tmean_perm =(Tbp + Tmp)/2;
#Obtain _veat transfer coeffecient for feed side    
        Re_f = fp.get_Re(Dh,vf,Sp,Tmean_feed); 
        Pr_bf  = fp.get_Pr(Sp,Tbf);
        Pr_mf  = fp.get_Pr(Sp,Tmf);
        
        Nu_f = fp.get_Nu(Re_f,Pr_bf,Pr_mf,Dh,L)
        k_f = fp.get_k(Sp,Tmean_feed);
        hf = Nu_f*k_f/Dh;
#Obtain heat transfer coeffecient for permeate side           
#For the permeate assume pure water
        Re_p = fp.get_Re(Dh,vp,0,Tmean_perm); 
        Pr_bp  = fp.get_Pr(0,Tbp);
        Pr_mp  = fp.get_Pr(0,Tmp);
        
        Nu_p = fp.get_Nu(Re_p,Pr_bp,Pr_mp,Dh,L)
        k_p = fp.get_k(0,Tmp);
        hp = Nu_p*k_p/Dh;
#Determination of the membrane heat transfer coeffecient.
#Isostress,Isotrain, or Maxwell corellations can be used.
#Isostress used as default
        T_m = (Tmf+Tmp)/2; 
        km = Isostress(e,ks,T_m);
        hm = km/delta # W/m2k 
#Determination of the membrane interface temperature. 
        H_v = .018*fp.get_heat_vaporization(Sp,Tmf)
      # print(H_v)
        Tmf = (Tbf*hf +hm*(Tbp + Tbf*(hf/hp)) - N*H_v)/(hf*(1+hm/hp) +hm)
        Tmp = (Tbp*hp +hm*(Tbf + Tbp*(hp/hf)) + N*H_v)/(hp*(1+hm/hf) +hm)

#Calculation of flux N
        
#Antoinne equation to calculate interfacial membrane temperatures.
        pmf = fp.get_pressure(Sp,Tmf);
        pmp = fp.get_pressure(0,Tmp);
        pavg = (pmf+pmp)/2;
        dP = pmf-pmp;
#Dusty Gas model

#Determination of knudsen number to determine regime
        Kn = get_knudsen(T_m,pavg,r);  
        vmean = numpy.sqrt(8*T_m*R/(3.141592*Mwater)) #m/s
#Determination of the membrane mass trasnfer coefficient
        Bm = get_Bm(e,r,t,Kn,T_m,vmean,delta);
#Determination of the mass flux
        N_temp = N
        N_molar = Bm*(dP)
        N = N_molar*Mwater;

        errorTmf = numpy.abs(tempf - Tmf)
        errorTmp = numpy.abs(tempp - Tmp)
        errorN = numpy.abs(N_temp - N)
        error = numpy.amax([errorTmf,errorTmp,errorN])
      # print(Kn,T_m,vmean)
    return N,Tmf,Tmp
