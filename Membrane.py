# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 20:33:59 2018

@author: Emanuel
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 20:36:21 2018

@author: Emanuel
"""

import numpy

Tbf = 50+273.15
Tbp = 20+273.15

#Initial guess
Tmf = Tbf
Tmp = Tbp 
N = 0 #Start with this assumption then continue until steady state. 

#Membrane properties
e = .62 #porosity unitless
r = .22e-6/2 #mean pore radius meters.
t  = ((2-e)**2)/e # unitless
delta = 126e-6 #Using PTFE - 1. Thickness in meters.
#We also need to know some values of the bulk fluid in order to calculate the prandtl
# and reynolds numbers for the calculation of the nusselt number.

##For the Feed
H_v = 43172 # j/mol of water at 44 celsius
#values taken for saltwater at 30 celsius
muf = .436e-3 #dynamic viscosity in Ns/m2
rhof = 1022 # desity 
kinvis = muf/rhof # The kinematic viscoisty
cp_f = 4042.9 #Specific heat capacity in J/kgK
k_f = .66 #W/mK
d = .4 # m. the characteristic size. I actually dont know.  
v = 3.33e-5 # 2l/m in m3/s

#isotrain for the heat transfer coeff. is what they used in the literature
#to derive their coeff so thats what ill assume
km = .041
hm = km/delta # W/m2k
#Lets calculate the reynolds number, prandtl, and nusselt for the feed.
L = 1E-1
error = 1
l = 0
#Caculating the interfactial temps
for n in range(2) :
    while error >1e-5:
        l = l+1
        Tmfi = Tmf
        Tmpi = Tmp
        
        ##For the Feed
        H_v = 43172# Saltwater of 30 g/kg
        #values taken for saltwater at 30 celsius
        rhof = 1022 # density 
        kinvis = muf/rhof # The kinematic viscoisty
        
        #Caclulation of cp_bf
        Sp = 0
        A = 5.328 - 9.76e-2*Sp + 4.04e-4*Sp**2
        B =  -6.913e-3 + 7.351e-4*Sp - 3.15e-6*Sp**2
        C =  9.6e-6 - 1.927e-6*Sp + 8.23e-9*Sp**2
        D =  2.5e-9 + 1.666e-9*Sp - 7.125e-12*Sp**2
        
        cp_bf = 1000*(A +B*Tbf + C*Tbf**2 + D*Tbf**3)
        #Caclulaton k_bf
        logk_bf = numpy.log10(240 +.0002*Sp) +.434*(2.3 - (343.5+.037*Sp)/(Tbf))*(1 - (Tbf)/(647 + .03*Sp))**.333
 
        k_bf = 10**(logk_bf) / 1000
        ##Caclulation of the viscosity
        mu_pwb =  9.67e-2 -8.207e-4*(Tbf) +2.344e-6*(Tbf**2)-2.244e-9*(Tbf**3)
        Tbf   = Tbf - 273.15 #This correlation uses celcsius so lets convert back for a second
        Sp = 0.0
        I  = 19.915*Sp/(1-1.00487*Sp)
        log10mus = .0428*I + .00123*I**2 + .000131*I**3 + (-.03724*I + .01859*I**2 -  .00271*I**3)*numpy.log10(1e3*mu_pwb)
        mu_bf = mu_pwb*10**( log10mus )
        Tbf   = Tbf + 273.15

        ##For the wall.
        
        
        #Caclulation of cp_mf
        Sp = 0 #g/kg
        A = 5.328 - 9.76e-2*Sp + 4.04e-4*Sp**2
        B =  -6.913e-3 + 7.351e-4*Sp - 3.15e-6*Sp**2
        C =  9.6e-6 - 1.927e-6*Sp + 8.23e-9*Sp**2
        D =  2.5e-9 + 1.666e-9*Sp - 7.125e-12*Sp**2
        
        cp_mf = 1000*(A +B*Tmf + C*Tmf**2 + D*Tmf**3)
        
        #Caclulation k_mf
        logk_mf = numpy.log10(240 +.0002*Sp) + .434*(2.3 - (343.5+.037*Sp)/(Tmf))*(1 - (Tmf)/(647 + .03*Sp))**.333
 
        k_mf = 10**(logk_mf) / 1000
        
        #Caclulation of the viscosity at the wall
        mu_pwb_wall =  9.67e-2 -8.207e-4*(Tmf) +2.344e-6*(Tmf**2)-2.244e-9*(Tmf**3)
        Tmf   = Tmf - 273.15 #This correlation uses celcsius so lets convert back for a second
        Sp = .00 #kg/kg normal sp divided by 1000
        I  = 19.915*Sp/(1-1.00487*Sp)
        log10muswall = .0428*I + .00123*I**2 + .000131*I**3 + (-.03724*I + .01859*I**2 -  .00271*I**3)*numpy.log10(1e3*mu_pwb_wall)
        mu_mf = mu_pwb*10**( log10muswall )
        Tmf   = Tmf + 273.15
        
        df = 9.09e-3  
        vf = 3.14 

        Re_bf = vf*df*rhof/muf 
     
        Pr_bf  = mu_bf*cp_bf/k_bf
        Pr_mf = mu_mf*cp_mf/k_mf
       # Nu_f = .097*(Re_bf**.73)*(Pr_bf**.13)*(Pr_bf/Pr_mf)**.25
       # Nu_f  = 4.36 + (.036)*Re_bf*Pr_bf*(df/L)/(1 +( .0011*Re_bf*Pr_bf*(df/L))**.8)
        Nu_f = .023*(1+ (6*D/L))*(Re_bf**.8)*(Pr_bf**.33333)
 
        hf = Nu_f*k_f/df


##For the permeate assume pu water
        #For the bulk
        rho_bp = y = 765.33+1.8142*Tbp-.0035*Tbp**2
        
        Sp = 0
        A = 5.328 - 9.76e-2*Sp + 4.04e-4*Sp**2
        B =  -6.913e-3 + 7.351e-4*Sp - 3.15e-6*Sp**2
        C =  9.6e-6 - 1.927e-6*Sp + 8.23e-9*Sp**2
        D =  2.5e-9 + 1.666e-9*Sp - 7.125e-12*Sp**2
        
        cp_bp = 1000*(A +B*Tbp + C*Tbp**2 + D*Tbp**3)
        
        mu_bp =  9.67e-2 -8.207e-4*Tbp +2.344e-6*Tbp**2-2.244e-9*Tbp**3
        k_bp = -.5752+6.397e-3*Tbp-8.151e-6*Tbp**2
        
        #For the wall
       
        Sp = 0
        A = 5.328 - 9.76e-2*Sp + 4.04e-4*Sp**2
        B =  -6.913e-3 + 7.351e-4*Sp - 3.15e-6*Sp**2
        C =  9.6e-6 - 1.927e-6*Sp + 8.23e-9*Sp**2
        D =  2.5e-9 + 1.666e-9*Sp - 7.125e-12*Sp**2
        
        cp_mp = 1000*(A +B*Tmp + C*Tmp**2 + D*Tmp**3)
        
        mu_mp =  9.67e-2 -8.207e-4*Tmp +2.344e-6*Tmp**2-2.244e-9*Tmp**3
        k_mp = -.5752+6.397e-3*Tmp-8.151e-6*Tmp**2
        
        dp = df
        vp = 4.67
        Re_bp = vp*dp*rho_bp/mu_bp 
        Pr_bp  = mu_bp*cp_bp/k_bp
        Pr_mp  = mu_mp*cp_mp/k_mp
        #Nu_p = .097*(Re_p**.73)*(Pr_bp**.13)*(Pr_bp/Pr_mp)**.25
        #Nu_p  = 4.36 + (.036*Re_bp*Pr_bp*(dp/L))/(1 +( .0011*Re_bp*Pr_bp*(dp/L))**.8)
        Nu_p = .023*(1+ (6*D/L))*(Re_bp**.8)*(Pr_bp**.33333)

        hp = Nu_p*k_bp/dp

        Tmf = (Tbf*hf +hm*(Tbp + Tbf*(hf/hp)) - N*H_v)/(hf*(1+hm/hp) +hm)
        Tmp = (Tbp*hp +hm*(Tbf + Tbp*(hp/hf)) + N*H_v)/(hp*(1+hm/hf) +hm)

#########################################################################
#Now we need to calculate the flux N

#First we need to link the pressure to the temp. P(T).
#Assume pure water and use the antoinne equation.

        pmf =numpy.exp(23.5377 - 4016.3632/(Tmf- 38.6339))
        #Adjust feed temperature according to water activity. 
        m  =  0
        aw = 1 - .03112*m-.001482*m**2
        pmf = pmf*aw
        pmp = numpy.exp(23.5377 - 4016.3632/(Tmp- 38.6339))
        pavg = (pmf+pmp)/2

#Lets use the Dusty gas models

#Constants we need. depend on the membrane structure but this is a good guess
        B0 = e*(r**2)/(8*t)
        K0 = 2*e*r/(3*t)
        K1 = e/t

        T = (Tmf+Tmp)/2 
        P = 10**5
        Kb = 1.38064852e-23
        mfp = Kb*T/(1.41421*pavg*(2.641e-10)**2)
        Kn  = mfp/(2*r)
        Mwater = 18.01528e-3 #kg/mol
        R = 8.31
        Dij = 27.425e-6#m/s
        paf = P-pmf
        pap = P-pmp

        vmean = numpy.sqrt(8*T*R/(3.141592*Mwater)) #m/s

        Dijm = K1*P*Dij
        Dijk = K0*vmean
        paln = (paf-pap)/numpy.log(paf/pap)
    
        Ni = N
        if Kn > 1:
            N =  K0*vmean/(R*T) *  (pmf-pmp)/(delta)
        if Kn <.01:
            N =  Dijm/(R*T*delta) * (pmf-pmp)/(paln)
        if Kn <1.0 and Kn>.01:
            N = Dijm/(R*T*delta) * numpy.log((Dijk*pap +Dijm)/(Dijk*pap +Dijm))
        errorTmf = numpy.abs(Tmf - Tmfi)
        errorTmp = numpy.abs(Tmp - Tmpi)
        errorN = numpy.abs(Ni - N)
        error = numpy.amax([errorTmf,errorTmp,errorN])
        

tau = (Tmf-Tmp)/(Tbf-Tbp)
print(tau)
