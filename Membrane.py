# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 09:54:54 2018
@author: Emanuel
"""

import numpy

def get_cp(Sp,T):
#Returns the specific heat of water in joules/Kg*K calculated at salinitiy
#Sp and temperature T in Kelvin.
#Sp is the salinity in grams of dissolved NaCl per Kg of water. 
    
    A = 5.328 - 9.76e-2*Sp + 4.04e-4*Sp**2;
    B =  -6.913e-3 + 7.351e-4*Sp - 3.15e-6*Sp**2;
    C =  9.6e-6 - 1.927e-6*Sp + 8.23e-9*Sp**2;
    D =  2.5e-9 + 1.666e-9*Sp - 7.125e-12*Sp**2;
    cp_bf = 1000*(A +B*T + C*T**2 + D*T**3);
    return cp_bf;


def get_k(Sp,T):    
#Returns the thermal conductivity of water in Watts/m*K for a given
#salinity Sp, grams of dissolved Nacl per kg of water, and T
    logk = numpy.log10(240 +.0002*Sp) +  \
    .434*(2.3 - (343.5+.037*Sp)/(T))*(1 - (T)/(647 + .03*Sp))**(.33333);
    k = 10**(logk) / 1000;
    return k;


def get_mu_pure(T): 
#https://syeilendrapramuditya.wordpress.com/2011/08/20/water-thermodynamic-properties/
#Returns the viscoity of pure water in Pa*s or kg/m*s
    mu_p =  9.67e-2 - 8.207e-4*(T) + 2.344e-6*(T)**2 - 2.244e-9*(T)**3;
    return mu_p;

def get_mu(Sp,T):
#Returns the density of seawater in kg/m*s
    Sp = Sp/1000; #Sp used in this correlation is in kg/kg
    mu_p = get_mu_pure(T); #Get the viscosity of pure water
    I  = (19.915*Sp)/(1-1.00487*Sp) #The ionic strength
    log10mus = .0428*I + .00123*I**2 + .000131*I**3 + \
    (-.03724*I + .01859*I**2 -  .00271*I**3)*numpy.log10(1e3*mu_p)
    mu = mu_p*10**( log10mus );
    return mu;


def get_rho(Sp,T):
#Returns the density of seawater in kg/m^3
#T is given in Kelvin, and Sp is g/kg.
    T=T-273;
    B = (2*Sp -150)/150;
    G1 = .5;
    G2 = B;
    G3 = 2*B**2 -1;
    A1 = 4.032*G1 + .115*G2 + 3.26e-4*G3;
    A2 = -.108*G1 + 1.571e-3*G2 - 4.23e-4*G3;
    A3 = -.012*G1 + 1.74e-3*G2 - 9e-6*G3;
    A4 = 6.92e-4*G1 - 8.7e-5*G2 - 5.3e-5*G3;
    
    A = (2*T - 200)/160;
    F1 = .5;
    F2 = A;
    F3 = 2*A**2 -1;
    F4 = 4*A**3 -3*A;
    rho = 1000*(A1*F1 + A2*F2 + A3*F3 + A4*F4);
    return rho;

def get_pressure(Sp,T):
#Determines the pressure in Pa based on the Antoinne equation
#Corrects for the salinity of water.
        pmf = numpy.exp(23.5377 - 4016.3632/(T- 38.6339))
        #Adjust feed temperature according to water activity. 
        m  =  Sp/58.4;
        aw = 1 - .03112*m-.001482*m**2;
        pmf = pmf*aw;
        return pmf;

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

def Maxwell(e,ks,T):
    kg = 2.72e-3+7.77e-5*T; #Thermal conductivity of water vapor. 
    B = (ks-kg)/(ks+2*kg);
    km = ( kg*(1 +2*B(1-e)) )/(1-B(1-e));
    return km; 

def get_Re(Dh,v,Sp,T):
    mu = get_mu(Sp,T);
    rho = get_rho(Sp,T);
    Re = v*Dh*rho/mu; 
    return Re;

def get_Pr(Sp,T):
    cp = get_cp(Sp,T);
    mu = get_mu(Sp,T);
    k = get_k(Sp,T);
    Pr = mu*cp/k;
    return Pr;

def get_Nu(Re,Pr_b,Pr_m,D,L):
    if (Re > 10000):
        Nu = .023*(1+ (6*D/L))*(Re**.8)*(Pr_b**.33333) 
    else:
        Nu  = 4.36 + ( (.036)*Re*Pr_b*(D/L) ) / (1 +( .0011* (Re*Pr_b*(D/L))**.8) )
    return Nu;

def get_Bm(e,r,t,Kn,T_m,vmean):
    
    if (Kn >1):        
        Bm = (2*e*r/R*T_m*3*t*delta) *vmean; 
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
#Obtain heat transfer coeffecient for feed side    
        Re_f = get_Re(Dh,vf,Sp,Tmean_feed); 
        Pr_bf  = get_Pr(Sp,Tbf);
        Pr_mf  = get_Pr(Sp,Tmf);
        
        Nu_f = get_Nu(Re_f,Pr_bf,Pr_mf,Dh,L)
        k_f = get_k(Sp,Tmean_feed);
        hf = Nu_f*k_f/Dh;
#Obtain heat transfer coeffecient for permeate side           
#For the permeate assume pure water
        Re_p = get_Re(Dh,vp,0,Tmean_perm); 
        Pr_bp  = get_Pr(0,Tbp);
        Pr_mp  = get_Pr(0,Tmp);
        
        Nu_p = get_Nu(Re_p,Pr_bp,Pr_mp,Dh,L)
        k_p = get_k(0,Tmp);
        hp = Nu_p*k_p/Dh;
#Determination of the membrane heat transfer coeffecient.
#Isostress,Isotrain, or Maxwell corellations can be used.
#Isostress used as default
        T_m = (Tmf+Tmp)/2; 
        km = Isostress(e,ks,T_m);
        hm = km/delta # W/m2k 
#Determination of the membrane interface temperature.        
        Tmf = (Tbf*hf +hm*(Tbp + Tbf*(hf/hp)) - N*H_v)/(hf*(1+hm/hp) +hm)
        Tmp = (Tbp*hp +hm*(Tbf + Tbp*(hp/hf)) + N*H_v)/(hp*(1+hm/hf) +hm)

#Calculation of flux N
        
#Antoinne equation to calculate interfacial membrane temperatures.
        pmf = get_pressure(Sp,Tmf);
        pmp = get_pressure(0,Tmp);
        pavg = (pmf+pmp)/2;
        dP = pmf-pmp;
#Dusty Gas model

#Determination of knudsen number to determine regime
        Kn = get_knudsen(T_m,pavg,r);  
        vmean = numpy.sqrt(8*T_m*R/(3.141592*Mwater)) #m/s
#Determination of the membrane mass trasnfer coefficient
        Bm = get_Bm(e,r,t,Kn,T_m,vmean);
#Determination of the mass flux
        N_temp = N
        N_molar = Bm*(dP)
        N = N_molar*Mwater;

        errorTmf = numpy.abs(tempf - Tmf)
        errorTmp = numpy.abs(tempp - Tmp)
        errorN = numpy.abs(N_temp - N)
        error = numpy.amax([errorTmf,errorTmp,errorN])
        print(Kn,T_m,vmean)
    return N,Tmf,Tmp

H_v = 2.257 # j/kg of water at 44 celsius
P = 10**5;
Kb = 1.38064852e-23;
Mwater = 18.01528e-3; #kg/mol
R = 8.31 ;
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
n,tmf,tmp = Membrane(Tbf,Tbp,e,r,t,delta,Sp,k_ptfe,Dh,L_both,vf,vp)
print(n,tmf,tmp)
