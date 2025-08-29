# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 18:26:34 2025

@author: akira
"""
import numpy as np
import matplotlib.pyplot as plt

# Link asteroid parameters with the CS parameter
d_nea = 10 # assume a value for the diameter of the asteroid
d_sigma = 0 # assume a value for one sigma (uncertainty on the asteroid size)

def CS_mass_estimation(d_nea, d_sigma):
    # Material volume/surface calculation -----------------------
    
    # --- Cylindrical portion -----
    R_cyl = (d_nea + 3*d_sigma)/2 * 1.5 # m
    R_ext = 0.5/2 # external radius of the exoskelton tubes (m) 
    e = 0.005 # thickness of the exoskelton (m) 
    h_cyl = R_cyl /7.5 * 8.5 # height of the cylinder (m)
    S_cyl_exoskelton = (R_ext*2*np.pi)*(6*h_cyl + R_cyl*2*np.pi)
    #print('Surface cyl exosketon: ' + str(S_cyl_exoskelton) + ' m2')
    
    R_airbag = R_cyl /7.5 * 5.3/2 #external radius of the airbag (m)
    #print('airbag radius ' + str(R_airbag) + ' m')
    e = 0.005 # thickness of the exoskelton (m) 
    S_airbag = R_airbag*np.pi*4*6
    #print('Surface airbag: ' + str(S_airbag) + ' m2')
    
    # ----Conical portion ---------
    h_cone = R_cyl /7.5 * 1.5 #because total height = 10 m 
    S_conic_exoskelton = np.pi*R_cyl*(np.sqrt(h_cone**2+R_cyl**2)+R_cyl)
    #print('Surface conic exosketon: ' + str(S_conic_exoskelton) + ' m2')
    
    e_arms = 0.01 # thickness of each plate of the arms (m)
    h_arms = 0.1
    l_arms = 0.5
    V_arms = 6*(h_arms*l_arms/2*R_cyl-(h_arms-e_arms)*(l_arms-e_arms)/2*(R_cyl-e_arms))
    
    # Mass calculation ------------
    
    density_kevlar = 0.049 * 0.001*10000 #kg/m2
    mass_exoskelton = density_kevlar*(S_airbag+S_conic_exoskelton+S_cyl_exoskelton)
    #print('mass of exoskeltons: ' + str(mass_exoskelton) + ' kg')
    
    density_CarbonFiber = 1.79*0.001*1000000 #kg/m3
    mass_arms = V_arms*density_CarbonFiber
    #print('mass of arms: ' + str(mass_arms) + ' kg')
    
    # Actuator ---------------
    # ------- torque applied on the actuator of the arms
    mu_sun = 1.327*10**20
    distance_SunSC = 149.6*10**9 #took distance earth sun
    T = R_cyl/2 * (mass_exoskelton/6 + mass_arms / 6) * mu_sun / distance_SunSC**2
    #print('torque: ' + str(T) + ' Nm')
    
    # ------- mass
    mass_actuator_arms = 0.9 #kg / unit 
    
    # Total mass -------------------
    mass_total = mass_arms + mass_exoskelton + mass_actuator_arms*6
    #print('total mass of the bag capture system: ' + str(mass_total) + ' kg')
    
    return mass_total
print('total mass of the bag capture system: ' + str(CS_mass_estimation(d_nea, d_sigma)) + ' kg')



# Plot -----------
# Generate some data
d_nea_int = np.linspace(5, 15, 100)
d_sigma_int = np.linspace(0, 10, 100)

mass_total_int = np.zeros(100)
mass_total_inf = np.zeros(100)
mass_total_sup = np.zeros(100)

for i in range(len(d_nea_int)):
    #print(d_nea_int[i])
    #print(mass_total_int[i])
    mass_total_int[i] = CS_mass_estimation(d_nea_int[i], 0)
    mass_total_inf[i] = CS_mass_estimation(d_nea_int[i]-1, 0)
    mass_total_sup[i] = CS_mass_estimation(d_nea_int[i]+1, 0)
    #mass_total_int[i] = CS_mass_estimation(10, d_sigma_int[i])

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(d_nea_int, mass_total_int, label='Capture system mass', color='blue', linestyle='-', linewidth=2)
#plt.plot(d_nea_int, mass_total_inf, label='Capture system min mass', color='blue', linestyle='--', linewidth=2)
#plt.plot(d_nea_int, mass_total_sup, label='Capture system max mass', color='blue', linestyle='-.', linewidth=2)

# Add titles and labels
plt.title('Capture system mass estimation with respect to the diameter of the asteroid')
plt.xlabel('diameter of the asteroid + uncertainty')
plt.ylabel('Mass of the capture system')
plt.xlim(5,15)

# Add a legend
plt.legend()

# Show the plot
plt.grid(True)
plt.show()


