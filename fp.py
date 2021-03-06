# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 23:33:18 2019

FLUID PROPERTIES functions


@author: Emanuel

Correlations for thermophysical properties of seawater. 

Sharqawy, Mostafa H., John H. Lienhard, and Syed M. Zubair. 
"Thermophysical properties of seawater: a review of existing correlations and data." 
Desalination and water Treatment 16.1-3 (2010): 354-380.

For all functions 
Sp is to be given in units of grams of dissolved NaCl per kg of water
T in Kelvin
"""

import numpy
import matplotlib

def get_cp(Sp,T):
#Returns the specific heat of water in joules/Kg*K calculated at salinitiy
#Validity: csw in (kJ/kg K); 273.15 < T68 < 453.15 K; 0 < SP < 180 g/kg
#Accuracy: ±0.28 %

    A = 5.328 - 9.76e-2*Sp + 4.04e-4*Sp**2;
    B =  -6.913e-3 + 7.351e-4*Sp - 3.15e-6*Sp**2;
    C =  9.6e-6 - 1.927e-6*Sp + 8.23e-9*Sp**2;
    D =  2.5e-9 + 1.666e-9*Sp - 7.125e-12*Sp**2;
    cp = 1000*(A +B*T + C*T**2 + D*T**3);
    return cp;


def get_k(Sp,T):    
#Returns the thermal conductivity of water in Watts/m*K for a given
#Validity: ksw in (W/m K); 0 < t68 < 180 C; 0 < SP < 160 g/kg
#Accuracy: ±3 %
    logk = numpy.log10(240 +.0002*Sp) +  \
    .434*(2.3 - (343.5+.037*Sp)/(T))*(1 - (T)/(647 + .03*Sp))**(1/3);
    k = 10**(logk) / 1000;
    return k;


def get_viscosity_pure(T): 
#Returns the viscoity of pure water in Pa*s or kg/m*s
#Validity: μw in (kg/m.s); 273 ≤ t ≤ 453 K;
#Accuracy: ±0.05 % (best fit to IAPWS 2008 [73] data)
    T = T-273;
    mu_p = 4.2844e-5 + (.157*(T+64.993)**(2)  - 91.296)**(-1)
    return mu_p;

def get_viscosity(Sp,T):
#Returns the dynamic viscoisty of seawater in kg/m*s
#I is the ionic strength I = (19.915*Sp)/(1-1.00487*Sp) 
#Validity: 293 < t < 453 K; 15 < Sp < 130 kg/kg
#Accuracy: ±0.4%
    if(Sp != 0):
        Sp = Sp/1000; #Sp used in this correlation is in kg/kg
        mu_p = get_viscosity_pure(T); #Get the viscosity of pure water
        I  = (19.915*Sp)/(1-1.00487*Sp) #The ionic strength
        log10mus = .0428*I + .00123*I**2 + .000131*I**3 + \
        (-.03724*I + .01859*I**2 -  .00271*I**3)*numpy.log10(1e3*mu_p)
        mu = mu_p*10**( log10mus );
    else:
        mu = get_viscosity_pure(T);
    return mu;

def get_density(Sp,T):
#Returns the density of seawater in kg/m^3
#T is given in Kelvin, and Sp is g/kg.
#Validity: ρ sw in (kg/m3); 293 < t < 453 K; 10 < SP < 160 g/kg
#Accuracy: ±0.1 %
    if(Sp != 0):
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
    else:
        rho = get_density_pure(T);
    return rho;

def get_density_pure(T):
#Returns the density of water in kg/m^3
#Validity: ρw in (kg/m3); 273 ≤ t ≤ 453 K
#Accuracy: ±0.01 % (best fi t to IAPWS 1995 [24] data)
    T = T-273;
    a1 = 9.999e2;
    a2 = 2.034e-2;
    a3 = -6.162e-3;
    a4 = 2.261e-5;
    a5 = -4.657e-8;
    p = a1 + a2*T +  a3*T**2 +a4*T**3 + a5*T**4
    return p;

def get_pressure_a(Sp,T):
#Determines the pressure in Pa based on the Antoinne equation
#Corrects for the salinity of water.
#A.J.J. Fontana (Ed.), Water Activity in Foods: Fundamentals and Applications,
#Wiley, 2007 (Appendix B).
        pmf = numpy.exp(23.5377 - 4016.3632/(T- 38.6339))
        #Adjust feed temperature according to water activity. 
        m  =  Sp/58.4;
        aw = 1 - .03112*m-.001482*m**2;
        pmf = pmf*aw;
        return pmf;

def get_pressure(Sp,T):
#Determines the pressure in Pa
#Validity: pv,sw in (atm); 273 < T48 < 313 K; 0 < SP < 40 g/kg;
#Accuracy: ±0.015 %    
    ln_pv_sw = 24.4543 - 67.4509*(100/T) - 4.8489*numpy.log(T/100) - 5.44e-4*Sp;
    pv_sw = numpy.exp(ln_pv_sw);
    pv_sw = pv_sw * 101325;
    return pv_sw;

def get_kinematic_viscosity(Sp,T):
#Returns the kinematic viscosity in m/s^2
    mu = get_viscosity(Sp,T);
    rho = get_density(Sp,T);
    kin_visc  = mu/rho;
    return kin_visc;

def get_heat_vaporization(Sp,T):
#Returns the latent Heat of Vaporization in Joules/kg
#validity: hfg,w in (J/kg); 0 ≤ t ≤ 200 oC
#Accuracy: ±0.01 % (best fi t to IAPWS 1995 [24] data)
    T = T-273
    hfw = 2.501e6 -2.369e3*T + 2.678e-1*T**2 -8.103e-3*T**3 -2.079e-5*T**4
    sal = (1-Sp/1000)
    hsw = hfw*sal
    return hsw

def get_specific_enthalpy(Sp,T):
#Validity: hsw and hw in (J/kg K); 283 ≤ t ≤ 393 K; 0 ≤ S ≤ 0.12 kg/kg;
#Accuracy: ±0.5 % from IAPWS 2008 
    S = Sp/1000;
    hw = get_enthalpy_pure(T)
    T = T-273;
    a1 = -2.348e4
    a2 = 3.152e5
    a3 = 2.803e6
    a4 = -1.446e7
    a5 = 7.826e3
    a6 = -4.417e1
    a7 = 2.139e-1
    a8 = -1.991e4
    a9 = 2.778e4
    a10 = 9.728e1
    hsw = hw - S*(a1 + a2*S + a3*S**2 + a4*S**3 + a5*T + a6*T**2 + a7*T**3 + a8*S*T + a9*S**2*T + a10*S*T**2)
    return hsw

def get_enthalpy_pure(T):
#Validity: hw in (J/kg); 278 ≤ t ≤ 473 K
#Accuracy: ±0.02% (best fi t to IAPWS 1995 data)    
    T = T-273
    hw = 141.335 + 4202.07*T -.535*T**2 + 8.193e-5*T**3
    return hw

def get_Re(Dh,v,Sp,T):
    kv = get_kinematic_viscosity(Sp,T);
    Re = v*Dh/kv; 
    return Re;

def get_Pr(Sp,T):
    cp = get_cp(Sp,T);
    mu = get_viscosity(Sp,T);
    k = get_k(Sp,T);
    Pr = mu*cp/k;
    return Pr;

def get_Nu(Re,Pr_b,Pr_m,D,L):
    if (Re > 10000):
        Nu = .023*(1+ (6*D/L))*(Re**.8)*(Pr_b**.33333) 
    else:
        Nu  = 4.36 + ( (.036)*Re*Pr_b*(D/L) ) / (1 +( .0011* (Re*Pr_b*(D/L))**.8) )
    return Nu;

#For Graphing
##############################
#Ts = numpy.linspace(278,473,1000);
#G = get_heat_vaporization(120,Ts);
#A = get_heat_vaporization(100,Ts);
#B = get_heat_vaporization(80,Ts);
#C = get_heat_vaporization(60,Ts);
#D = get_heat_vaporization(40,Ts);
#E = get_heat_vaporization(20,Ts);
#F = get_heat_vaporization(0,Ts);
#G = G/1000
#B = B/1000
#D = D/1000
#F = F/1000
#matplotlib.pyplot.plot(Ts,G,Ts,B,Ts,D,Ts,F)
#matplotlib.pyplot.xlabel("Temperature (K)")
#matplotlib.pyplot.ylabel("Latent Heat of Vaporization (kJ/kg) ")
