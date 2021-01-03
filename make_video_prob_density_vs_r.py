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
Z_3d = np.zeros((n_r, n_theta, n_phi))

x_arr = np.zeros((n_r))
z_arr = np.zeros((n_r))

theta_index = int(n_theta/2)
phi_index = int(n_phi/2)

def animate(i):
    my_file = 'pd_vs_t_' + str(i) + '.txt'
    print(i)
    fig.clear()
    r_p, theta_p, phi_p, pd_real, pd_im = np.genfromtxt(my_file, unpack=True)
    ax = plt.axes(xlim=(0.0, R), ylim=(0.0, max_pd))
    it = 0
    for I in range(0, n_r):
        for j in range(0, n_theta):
            for k in range(0, n_phi):
                X_3d[I][k] = r_p[it]
                Y_3d[I][k] = phi_p[it]
                Z_3d[I][k] = pd_real[it]
                it = it + 1

    for I in range(0, n_r):
        x_arr[I] = X_3d[I][theta_index][phi_index]
        z_arr[I] = Z_3d[I][theta_index][phi_index]
            
    cont = plt.plot(x_arr, z_arr)
    return cont

size_t = int(Nt)
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('prob_density_solution_vs_r.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
