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

H_v = 43172 # j/mol of water at 44 celsius
P = 10**5;
Kb = 1.38064852e-23;
Mwater = 18.01528e-3; #kg/mol
R = 8.31 ;
l = 0;

def Membrane(Tbf,Tbp,e,r,t,delta,Sp,ks) :
    Tmf = Tbf ;
    Tmp = Tbp ;
    N = 0 ;
    error  = 1 ; 
    while error > 1e-5:
        tempf = Tmf;
        tempp = Tmp;
        
        Tmean_feed = (Tbf + Tmf)/2;
        Tmean_perm =(Tbp + Tmp)/2;
        
        cp_bf = get_cp(Sp,Tbf);
        cp_mf = get_cp(Sp,Tmf);
        
        k_bf = get_k(Sp,Tbf);
        k_mf = get_k(Sp,Tmf);
        k_f = get_k(Sp,Tmf);
        
        mu_bf = get_mu(Sp,Tbf);
        mu_mf = get_mu(Sp,Tmf);
        mu_f = get_mu(Sp,Tmean_feed);
        
        rho_f = get_rho(Sp,Tmean_feed);
        
        Dh = 9.09e-3  
        vf = 2.87;
        Re_f = vf*Dh*rho_f/mu_f; 
        
        Pr_bf  = mu_bf*cp_bf/k_bf
        Pr_mf = mu_mf*cp_mf/k_mf
        
        Nu_f = .097*(Re_f**.73)*(Pr_bf**.13)*(Pr_bf/Pr_mf)**.25
        #Nu_f  = 4.36 + (.036)*Re_bf*Pr_bf*(df/L)/(1 +( .0011*Re_bf*Pr_bf*(df/L))**.8)
        #Nu_f = .023*(1+ (6*D/L))*(Re_bf**.8)*(Pr_bf**.33333) 
        hf = Nu_f*k_f/Dh;
##For the permeate assume pure water
        cp_bp = get_cp(Sp,Tbp);
        cp_mp = get_cp(Sp,Tmp);
        
        k_bp = get_k(Sp,Tbp);
        k_mp = get_k(Sp,Tmp);
        k_p = get_k(Sp,Tmp);
        
        mu_bp = get_mu(Sp,Tbp);
        mu_mp = get_mu(Sp,Tmp);
        mu_p = get_mu(Sp,Tmean_perm);
        
        rho_p = get_rho(Sp,Tmean_perm);
        
        vp = 4.67
        Re_p = vp*Dh*rho_p/mu_p 
        Pr_bp  = mu_bp*cp_bp/k_bp
        Pr_mp  = mu_mp*cp_mp/k_mp
        
        Nu_p = .097*(Re_p**.73)*(Pr_bp**.13)*(Pr_bp/Pr_mp)**.25
        #Nu_p  = 4.36 + (.036*Re_bp*Pr_bp*(dp/L))/(1 +( .0011*Re_bp*Pr_bp*(dp/L))**.8)
        #Nu_p = .023*(1+ (6*D/L))*(Re_bp**.8)*(Pr_bp**.33333)
        hp = Nu_p*k_p/Dh;
#Determination of the membrane heat transfer coeffecient.
        T_m = (Tmf+Tmp)/2; 
        km = Isostress(e,ks,T_m);
        hm = km/delta # W/m2k 
#Determination of the membrane interface temperature.        
        Tmf = (Tbf*hf +hm*(Tbp + Tbf*(hf/hp)) - N*H_v)/(hf*(1+hm/hp) +hm)
        Tmp = (Tbp*hp +hm*(Tbf + Tbp*(hp/hf)) + N*H_v)/(hp*(1+hm/hf) +hm)

#########################################################################
#Now we need to calculate the flux N

#Assume pure water and use the antoinne equation.
        pmf = get_pressure(Sp,Tmf);
        pmp = get_pressure(0,Tmp);
        pavg = (pmf+pmp)/2;
        paf = P+pmf;
        pap = P+pmp;
        pa_avg =(pap+pap)/2;
#Dusty Gas model
#Constants we need. depend on the membrane structure but this is a good guess
        B0,K0,K1 = get_memconstants(e,r,t);
        Kn = get_knudsen(T_m,pavg,r);
        
        vmean = numpy.sqrt(8*T_m*R/(3.141592*Mwater)) #m/s
        Dij = 27.425e-6#m/s
        Dijm = K1*pa_avg*Dij
        Dijk = K0*vmean
        paln = (paf-pap)/numpy.log(paf/pap)
        
        N_temp = N
        if Kn > 1:
            N =  K0*vmean/(R*T_m) *  (pmf-pmp)/(delta);
            mu_mem = get_mu(0,T_m);
            N = N + (1/(R*T_m*delta))*(Dijk + (pavg*B0/mu_mem)*(pmp-pmf))

        if Kn <.01:
            N =  Dijm/(R*T_m*delta) * (pmf-pmp)/(paln)
            mu_mem = get_mu(0,T_m);
            N = N + (1/(R*T_m*delta))*(Dijk + (pavg*B0/mu_mem)*(pmp-pmf))
        if Kn <1.0 and Kn>.01:
            N = Dijm/(R*T_m*delta) * numpy.log((Dijk*pap +Dijm)/(Dijk*pap +Dijm))
            mu_mem = get_mu(0,T_m);
            N = N + (1/(R*T_m*delta))*(Dijk + (pavg*B0/mu_mem)*(pmp-pmf))
      
        errorTmf = numpy.abs(tempf - Tmf)
        errorTmp = numpy.abs(tempp - Tmp)
        errorN = numpy.abs(N_temp - N)
        error = numpy.amax([errorTmf,errorTmp,errorN])
        print(errorTmf , errorTmf ,errorN );
        return N

e = .62; #porosity unitless
r = .22e-6/2; #mean pore radius meters.
t  = ((2-e)**2)/e; # unitless 
delta = 126e-6 ;#Using PTFE - 1. Thickness in meters.
#Initial bulk fluid conditions
Tbf = 60+273.15;
Tbp = 20+273.15;
Sp = 0;
k_ptfe = .25
n = Membrane(Tbf,Tbp,e,r,t,delta,Sp,k_ptfe )
