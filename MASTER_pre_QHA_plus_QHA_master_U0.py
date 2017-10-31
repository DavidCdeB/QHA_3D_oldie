#
# QHA -Master- program. David Carrasco de Busturia, 13 October 2017 
# Please read the documentation and istructions on: https://github.com/DavidCdeB/QHA
# This program is under the GNU General Public License v3.0. 
# davidcarrascobustur@gmail.com
# d.carrasco-de-busturia@imperial.ac.uk

import sympy as sym
import numpy as np
import sys


from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt

import time
start = time.time()


################### CONSTANTS: they are 
# going to be called constantly by functions, so let's 
# assign them as global variables:

global KB, h, speed_of_light, U0, V0, B0, B0_prime, N_k, conv_fac_nu, conv_fac_denu

# Value of the constants:

# KB = boltmann cte, KB = 1.38064852(79)x10-23 J/K
KB = 1.38064852E-23

# h = plank constant, h = 6.626070040(81)x10-34 J s
h = 6.626070040E-34

# c = speed of light, c = 2.99792458E8 m/s
speed_of_light = 2.99792458E+8

# conversion factor: conv_fac; when the fit is performed in freqs(cm-1) vs volume(A^3):

# h     *  (c*V**2 + d*V + f) * conv_fac   = [J]
# [Js]     [       cm^-1    ]    [1E+2 * c]
     
# h     *  (2*c*V + d*V    ) * conv_fac_denu  = [J]
# [Js]     [ cm^-1 * V^-1  ]   [1E+2 * c * 1E+30]
 

conv_fac_nu =  1E+2 * speed_of_light  
conv_fac_denu =  1E+2 * speed_of_light * 1E+30 

# We import U0, V0, B0 and B0_prime, N_k and n_F_u variables from the pre_QHA script:
import pre_QHA_plus_QHA_master_U0 

U0 =   pre_QHA_plus_QHA_master_U0.U0
V0 =   pre_QHA_plus_QHA_master_U0.V0
B0 =   pre_QHA_plus_QHA_master_U0.B0
B0_prime =  pre_QHA_plus_QHA_master_U0.B0_prime
N_k =  pre_QHA_plus_QHA_master_U0.N_k
n_F_u = pre_QHA_plus_QHA_master_U0.n_F_u

print 'U0 =  ', U0
print 'V0 =  ', V0
print 'B0 =  ', B0
print 'B0_prime =  ', B0_prime
print 'N_k =  ', N_k
print 'n_F_u =  ', n_F_u



##### For each "i" frequency, there is a c_i, d_i and f_i :

Cs, Ds, Fs, mode = np.loadtxt('./done_modes_sorted.dat', skiprows = 1).T


#### Here I set the values of "T" and "V" for which I would like
# to calculate G.
Ts = np.linspace(10.0, 2000.0, 100)
print Ts

# The same Vs as we computed on DFT:
Vs = np.array([ 116.573346 , 118.139505 , 119.713653 , 121.297131 , 122.894958 , 124.510211, \
  124.512598,  125.884132 , 127.054446 , 127.265886 , 128.656314,  130.054927,\
  131.463313 , 132.880152 , 134.309582 , 135.750582 , 137.200672])

########### Functions:

# Eq. 9:
def U(V):
       return  U0+ (2.293710449E+17)*(1E-21)*( (9.0/16.0)*(V0*B0) * (  (((V0/V)**(2.0/3.0)-1.0)**3.0)*B0_prime  + ((V0/V)**(2.0/3.0)-1)**2  * (6.0-4.0*(V0/V)**(2.0/3.0))  ))

# Eq. 2:
def P(V):
    f0=(3.0/2.0)*B0
    f1=((V0/V)**(7.0/3.0))-((V0/V)**(5.0/3.0))
    f2=((V0/V)**(2.0/3.0))-1
    pressure= f0*f1*(1+(3.0/4.0)*(B0_prime-4)*f2)
    return pressure 

def nu(V, c, d, f):
   return (c*V**2 + d*V + f) * conv_fac_nu

def denu(V, c, d):
   return (2*c*V + d) * conv_fac_denu

def P_sq(V, T, c, d, f):
   return (0.5 + (1.0 / (( np.exp(h * nu(V, c, d, f)  / (KB*T) )) - 1.0 ))) * h * denu(V, c, d) 

def ET_sq(V, T, c, d, f):
   return  h * nu(V, c, d, f) / (( np.exp( h * nu(V, c, d, f)  /(KB*T))) - 1.0)  

def S_sq(V, T, c, d, f):
   return  np.log (1.0 - np.exp((-h * nu(V, c, d, f))/(KB*T))  )

def ZPE_sq(V, c, d, f):
   return  h * nu(V, c, d, f)  


Gs = []
for V in Vs:
        aux = []
        for T in Ts:
                ET_sq_sum = 0.0
                P_sq_sum = 0.0
                S_sq_sum = 0.0
                ZPE_sq_sum = 0.0
                for c, d, f in zip(Cs, Ds, Fs):        
                        ET_sq_sum += ET_sq(V, T, c, d, f)
        	        P_sq_sum += P_sq(V, T, c, d, f)
                	S_sq_sum += S_sq(V, T, c, d, f)
                	ZPE_sq_sum += ZPE_sq(V, c, d, f)
        	aux.append( U(V) \
        		     +  (1.0/N_k) * 0.5 * ZPE_sq_sum * (1/(4.3597482E-18))\
                             + (1.0/N_k) * ET_sq_sum * (1/(4.3597482E-18))\
                             + (P(V) - (1.0/N_k) * P_sq_sum  * 1E-9 )* V * 1E+9 * 1E-30 * (1/(4.3597482E-18))\
                             - T * (1.0/N_k) * ( (-KB * S_sq_sum * (1/(4.3597482E-18))) + ((1.0/T) * ET_sq_sum * (1/(4.3597482E-18))) ))

        Gs.append(aux)

# converting Vs to Ps:  
Ps = []
for V in Vs:
        aux = []
        for T in Ts:
                P_sq_sum = 0.0
                for c,d,f in zip(Cs, Ds, Fs):
                        P_sq_sum += P_sq(V, T, c, d, f)
                aux.append(P(V) - (1.0/N_k) * P_sq_sum * 1E-9)
        Ps.append(aux)

print 'Gs = ', Gs
print 'Ts = ', Ts
print 'Vs = ', Vs
print 'Ps = ', Ps
print ' np.shape(Ts) = ', np.shape(Ts)
print ' np.shape(Vs) = ', np.shape(Vs)
print ' np.shape(Gs) = ' , np.shape(Gs)
print ' np.shape(Ps) = ', np.shape(Ps)
len_Ts = len(Ts)
len_Vs = len(Vs)

#n_F_u = 2.0
print 'type(Vs) = ', type(Vs) 
Vs = Vs / n_F_u

print 'type(Gs) = ', type(Gs) 
print 'type(Ps) = ', type(Ps) 


# Transform Gs and Ps from list -> numpy array:
Gs = np.asarray(Gs)
Ps = np.asarray(Ps)
print 'type(Gs) = ', type(Gs) 
print 'type(Ps) = ', type(Ps)
# Now that Gs is numpy array, I can divide it by n_F_u: 
Gs = Gs / n_F_u
print 'Gs = ', Gs

# Flatten both Gs and Ps:
Gs = Gs.flatten()
Ps = Ps.flatten()
print 'Gs = ', Gs
print 'type(Gs) = ', type(Gs) 
print 'Ps = ', Ps
print 'type(Ps) = ', type(Ps) 


Vs = np.repeat(Vs, len_Ts)
Ts = np.tile(Ts,len_Vs)

print 'Gs = ', Gs
print 'Ts = ', Ts
print 'Vs = ', Vs
print 'Ps = ', Ps
print ' np.shape(Ts) = ', np.shape(Ts)
print ' np.shape(Vs) = ', np.shape(Vs)
print ' np.shape(Gs) = ' , np.shape(Gs)
print ' np.shape(Ps) = ', np.shape(Ps)


output_array = np.vstack((Vs, Ps, Ts, Gs)).T 
np.savetxt('Vs_Ps_Gs.dat', output_array, header="Vs  Ps  Ts  Gs", fmt="%0.13f")


##### Plotting:

# Load data:
y_data, z_data, x_data  = np.loadtxt('./solid_1__xyz_sorted_as_P_wise.dat').T

y_data_2, z_data_2, x_data_2  = np.loadtxt('./solid_1__xyz_sorted_as_P_wise.dat').T


####### Calcite I scattered:
# In a new figure, each surface separately:
# set "fig" and "ax" varaibles
fig = plt.figure()
ax = fig.gca(projection='3d')

# Plot the initial scattered points
ax.scatter(x_data, y_data, z_data, color='r', marker='o') # 'ro') #color='r', marker='o')
ax.scatter(Ts, Ps, Gs, color='b', marker='o') # 'ro') #color='r', marker='o')


ax.set_xlabel('\nT (K)')
ax.set_ylabel('P (GPa)')
ax.set_zlabel('\nGibbs free energy / F.unit (a.u.)', linespacing=3)
ax.set_title('\n\nCalcite I', linespacing=3)
#plt.xticks(x, labels, rotation='vertical')
xlabels=[0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]
ax.set_xticklabels(xlabels,rotation=90,
                  verticalalignment='baseline',#)#,
                  horizontalalignment='left')

ylabels = [-4, -2, 0, 2, 4, 6, 8, 10]
ax.set_yticklabels(ylabels,rotation=0,
                  verticalalignment='baseline')#,
#                 horizontalalignment='left')

fig.savefig("Calcite_I_scattered_DFT_data_plus_GTP.pdf",  bbox_inches='tight', pad_inches=0.3)#, tight_layout() )#, bbox_inches=bbox)

end = time.time()
print end-start , "seconds//"

plt.show()
#end = time.time()

