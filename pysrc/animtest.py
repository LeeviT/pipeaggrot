import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

print(plt.style.available)
plt.style.use('seaborn')


def readfile(filename):
    with open(filename, 'r') as data:
        t = []
        visc = []
        for line in data:
            p = line.split()
            t.append(float(p[0]))
            visc.append(float(p[1]))
    return t, visc


filelist = []

t0, visc0 = readfile("visctot.dat")

fig = plt.figure(figsize=(8, 4))
base = fig.add_subplot(1, 2, 2)
base.plot(t0, visc0)
base.set_xlabel('time (s)')
base.set_ylabel('viscosity')
velo = fig.add_subplot(1, 2, 1)
velo.set_xlim(0, 1)
velo.set_ylim(0, 1)
velo.set_ylabel("flow velocity")
velo.set_xlabel("radius")
axes = fig.add_subplot(1, 2, 2)
axes.set_xlim(1e-4, 0.1)
axes.set_ylim(min(visc0)-0.02, max(visc0)+0.02)
# axes.semilogx()
plt.subplots_adjust(top=0.97, bottom=0.12, right=0.97, left=0.07, hspace=0.05, wspace=0.2)

point, = axes.plot([t0[0]], [visc0[0]], 'o', markersize=5)

filelist = []
r0, vx0 = readfile("r_vx0.dat")
vxnorm = vx0[0]

for j in range(0, 1000):
    filelist.append("r_vx%s.dat" % j)

curve, = velo.plot([readfile(filelist[0])[0]], [readfile(filelist[0])[1]], '-')


def readfile(filename):
    with open(filename, 'r') as data:
        x = []
        y = []
        for line in data:
            p = line.split()
            x.append(float(p[0]))
            y.append(float(p[1]))
    return x, y


def ani(coords, filelist, i):
    point.set_data([coords[0]], [coords[1]])
    curve.set_data([readfile(filelist[i])[0]], [readfile(filelist[i])[1]])
    return point


def frames1():
    for i in range(0, 1000):
        r, vx = readfile(filelist[i])
        velo.plot(r, np.array(vx)/vxnorm)
        yield t0[i], visc0[i]


def frames2():
    for i in range(100, 300):
        r, vx = readfile(filelist[i])
        velo.plot(r, np.array(vx)/vxnorm)
        yield t0[i], visc0[i]


def frames3():
    for i in range(300, 1000):
        r, vx = readfile(filelist[i])
        velo.plot(r, np.array(vx)/vxnorm)
        yield t0[i], visc0[i]


ani1 = FuncAnimation(fig, ani, frames=frames1(), interval=35, repeat=False, save_count=10)
ani1.save("anim1.mp4", dpi=300)
# ani2 = FuncAnimation(fig, ani, frames=frames2(), interval=35, repeat=False, save_count=300)
# ani2.save("anim2.mp4", dpi=300)
# ani3 = FuncAnimation(fig, ani, frames=frames3(), interval=35, repeat=False, save_count=600)
# ani3.save("anim3.mp4", dpi=300)

plt.show()
