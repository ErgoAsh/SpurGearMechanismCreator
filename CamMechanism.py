import seaborn as sns
import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt

dTheta = 1 # Resolution of the x axis (degrees)

""" Initial data """
h = 23 # Stroke of the follower
Ro = 102.67 # Base circle calculated in SolidWorks
start_angle = 90 + 5.08

alpha_samples = [64.97, 15.03, 130, 150]
alpha_angles = [alpha_samples[0], 
                alpha_samples[0] + alpha_samples[1],
                alpha_samples[0] + alpha_samples[1] + alpha_samples[2],
                alpha_samples[0] + alpha_samples[1] + alpha_samples[2] + alpha_samples[3]]

samples = [65, 15, 130, 150] # i-th Theta angles
rise_angle = np.radians(samples[0])
return_angle = np.radians(samples[2])

theta_1 = np.arange(0, samples[0], dTheta)
theta_2 = np.arange(0, samples[1], dTheta)
theta_3 = np.arange(0, samples[2], dTheta)
theta_4 = np.arange(0, samples[3], dTheta)
theta = np.linspace(0, 360, len(theta_1) + len(theta_2) + len(theta_3) + len(theta_4))

theta_1_rad = np.linspace(0, np.radians(samples[0]), len(theta_1))
theta_2_rad = np.linspace(0, np.radians(samples[1]), len(theta_2))
theta_3_rad = np.linspace(0, np.radians(samples[2]), len(theta_3))
theta_4_rad = np.linspace(0, np.radians(samples[3]), len(theta_4))
theta_rad = np.linspace(0, 2 * math.pi, len(theta_1) + len(theta_2) + len(theta_3) + len(theta_4))

""" Acceleration """
acc_rise = (2 * h * math.pi / rise_angle**2) \
              * np.sin(2 * math.pi * theta_1_rad / rise_angle)

acc_dwell_top = np.repeat(0, len(theta_2))

acc_return = np.flipud(-(2 * h * math.pi / return_angle**2) \
                         * np.sin(2 * math.pi * theta_3_rad / return_angle))

acc_dwell_bottom = np.repeat(0, len(theta_4))

acceleration = np.concatenate((acc_rise, acc_dwell_top, acc_return, acc_dwell_bottom))

""" Velocity """ 
vel_rise = (h / rise_angle) * (1 - np.cos(2 * math.pi * theta_1_rad / rise_angle))

vel_dwell_top = np.repeat(0, len(theta_2))

vel_return = np.flipud(-(h / return_angle) \
                  * (1 - np.cos(2 * math.pi * theta_3_rad / return_angle)))

vel_dwell_bottom = np.repeat(0, len(theta_4))

velocity = np.concatenate((vel_rise, vel_dwell_top, vel_return, vel_dwell_bottom))

""" Displacement """
h_rise = h * ((theta_1_rad / rise_angle) - (1 / (2 * math.pi)) \
                * np.sin(2 * math.pi * theta_1_rad / rise_angle))

h_dwell_top = np.repeat(h, len(theta_2))

h_return = np.flipud((h / math.pi) * ((math.pi * theta_3_rad / return_angle) - (1 / 2) \
                * np.sin(2 * math.pi * theta_3_rad / return_angle)))

h_dwell_bottom = np.repeat(0, len(theta_4))

displacement = np.concatenate((h_rise, h_dwell_top, h_return, h_dwell_bottom))

""" Visualization """
motion_data = pd.DataFrame({
   'theta': theta, 
   'disp': displacement,
   'vel': velocity,
   'acc': acceleration
});

sns.set(rc = {'text.usetex': True})
fig, axs = plt.subplots(nrows=3, figsize=(11.7, 8.27)) # A4 size plots
plt.subplots_adjust(hspace = 0.5)

sns.lineplot(data = motion_data, x = 'theta', y = 'acc', ax = axs[0])
sns.set_style("whitegrid")
axs[0].set(xlabel = r'$\theta$ [$deg$]', ylabel = r'Acceleration [$\frac{mm}{rad^2}$]')

sns.lineplot(data = motion_data, x = 'theta', y = 'vel', ax = axs[1])
sns.set_style("whitegrid")
axs[1].set(xlabel = r'$\theta$ [$deg$]', ylabel = r'Velocity [$\frac{mm}{rad}$]')

sns.lineplot(data = motion_data, x = 'theta', y = 'disp', ax = axs[2])
sns.set_style("whitegrid")
axs[2].set(xlabel = r'$\theta$ [$deg$]', ylabel = r'Displacement [$mm$]')

for i in range(0, 3):
   axs[i].set_xlim(0, 360)
   axs[i].set_xticks(np.arange(0, 360.1, 30))

""" Data export """
# X axis: displacement, Y axis: velocity analogue
vel_diagram_rise = pd.DataFrame({'x': -vel_rise,
                                 'y': np.linspace(0, 23, len(theta_1)),
                                 'z': np.repeat(0, len(theta_1))})

vel_diagram_return = pd.DataFrame({'x': -vel_return,
                                   'y': np.linspace(0, 23, len(theta_3)),
                                   'z': np.repeat(0, len(theta_3))})

#Export 2 velocity splines;  %2.5f - 2 integers digits, 5 decimal digits
np.savetxt('vel_data_1.txt', vel_diagram_rise.values, delimiter = ' ', fmt = '%2.5f')
np.savetxt('vel_data_2.txt', vel_diagram_return.values, delimiter = ' ', fmt = '%2.5f')

# 0.01 offset is added to differentiate 2 equal points
alpha = np.concatenate([np.linspace(0, alpha_angles[0], len(theta_1)),
                        np.linspace(alpha_angles[0] + 0.01, alpha_angles[1], len(theta_2)),
                        np.linspace(alpha_angles[1] + 0.01, alpha_angles[2], len(theta_3)),
                        np.linspace(alpha_angles[2] + 0.01, alpha_angles[3], len(theta_4))], axis = 0)
                        
# Translate polar coordinates to x y coordinates,
displacement_cartesian = pd.DataFrame({
   'x': (motion_data['disp'] + Ro) * np.cos(np.radians(start_angle - alpha)),
   'y': (motion_data['disp'] + Ro) * np.sin(np.radians(start_angle - alpha)),
   'z': np.repeat(0, len(motion_data['theta']))
})

np.savetxt('cam_profile.txt', displacement_cartesian.values, delimiter = ' ', fmt = '%2.5f')
plt.figure(figsize=(25,20))

np.savetxt('data.txt', motion_data.values, delimiter = ' ', fmt = '%2.5f')

plot = sns.scatterplot(data = displacement_cartesian, x = 'x', y = 'y', size = 0.5)
plt.xlim(-125, 125)
plt.ylim(-125, 125)
sns.set_style("darkgrid")
plt.rcParams['figure.figsize'] = (1, 1)
#axs[3].set_xticks(np.arange(-30, 30, 1))
