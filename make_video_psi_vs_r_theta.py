import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.animation as manimation
import time
import sys
from mpl_toolkits.mplot3d import Axes3D

start_time = time.time()

plt.rcParams['animation.ffmpeg_path']='C:/Users/d-w-h/Downloads/ffmpeg-20200818-1c7e55d-win64-static/ffmpeg-20200818-1c7e55d-win64-static/bin/ffmpeg.exe'
writer=manimation.FFMpegWriter(bitrate=20000, fps=15)

fig = plt.figure()
ax = fig

file_timesteps = 'number_of_timesteps.txt'
Nt = np.genfromtxt(file_timesteps, unpack=True)

file_limits = 'limits.txt'
max_real, min_real, max_im, min_im, R, max_pd = np.genfromtxt(file_limits, unpack=True)

file_dims = 'dims.txt'
n_r, n_theta, n_phi = np.genfromtxt(file_dims, unpack=True)

n_r = int(n_r)
n_theta = int(n_theta)
n_phi = int(n_phi)

X_3d = np.zeros((n_r, n_theta, n_phi))
Y_3d = np.zeros((n_r, n_theta, n_phi))
Y_3d_theta = np.zeros((n_r, n_theta, n_phi))
Z_3d = np.zeros((n_r, n_theta, n_phi))
Z_3d_im = np.zeros((n_r, n_theta, n_phi))

X_2d = np.zeros((n_r, n_theta))
Y_2d = np.zeros((n_r, n_theta))
Z_2d = np.zeros((n_r, n_theta))

x_arr = np.zeros((n_r))
z_arr = np.zeros((n_r))
z_arr_im = np.zeros((n_r))

max_val = 0
min_val = 0
if(max_real > max_im):
    max_val = max_real
else:
    max_val = max_im

if(min_real < min_im):
    min_val = min_real
else:
    min_val = min_im

theta_index = int(n_theta/2)
phi_index = int(n_phi/2)
    
def animate(i):
    my_file = 'psi_vs_t_' + str(i) + '.txt'
    print(i)
    fig.clear()
    r_p, theta_p, phi_p, psi_real, psi_im = np.genfromtxt(my_file, unpack=True)
    it = 0
    for I in range(0, n_r):
        for j in range(0, n_theta):
            for k in range(0, n_phi):
                X_3d[I][j][k] = r_p[it]
                Y_3d[I][j][k] = phi_p[it]
                Y_3d_theta[I][j][k] = theta_p[it]
                Z_3d[I][j][k] = psi_real[it]
                Z_3d_im[I][j][k] = psi_im[it]
                it = it + 1

    for I in range(0, n_r):
        for j in range(0, n_theta):
            X_2d[I][j] = X_3d[I][j][phi_index]
            Y_2d[I][j] = Y_3d_theta[I][j][phi_index]
            Z_2d[I][j] = Z_3d[I][j][phi_index]
                
    cont = plt.contourf(X_2d, Y_2d, Z_2d, linspace(min_val, max_val, 30))
    return cont

size_t = int(Nt)
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('SE_solution_psi_vs_r_theta.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
