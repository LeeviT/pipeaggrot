import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

print(plt.style.available)
plt.style.use('bmh')


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
base.set_ylabel('viscosity (mPa*s)')
velo = fig.add_subplot(1, 2, 1)
axes = fig.add_subplot(1, 2, 2)
axes.set_xlim(1e-4, 1)
axes.set_ylim(min(visc0)-0.02, max(visc0)+0.02)
axes.semilogx()
plt.subplots_adjust(top=0.97, bottom=0.05, right=0.97, left=0.04, hspace=0.05, wspace=0.2)

point, = axes.plot([t0[0]], [visc0[0]], 'o', markersize=5)

filelist = []
r0, vx0 = readfile("r_vx0.dat")
vxnorm = vx0[0]

for i in range(0, 30):
    filelist.append("r_vx%s.dat" % i)

def readfile(filename):
    with open(filename, 'r') as data:
        t = []
        visc = []
        for line in data:
            p = line.split()
            t.append(float(p[0]))
            visc.append(float(p[1]))
    return t, visc


def ani(coords):
    point.set_data([coords[0]], [coords[1]])
    return point


def frames1():
    for i in range(0, 100):
        yield t0[i], visc0[i]


def frames2():
    for i in range(100, 300):
        yield t0[i], visc0[i]


def frames3():
    for i in range(300, 1000):
        yield t0[i], visc0[i]


ani1 = FuncAnimation(fig, ani, frames=frames1(), interval=100, repeat=False, save_count=100)
ani1.save("anim1.mp4", dpi=200)
#ani2 = FuncAnimation(fig, ani, frames=frames2(), interval=50, repeat=True, save_count=300)
#ani2.save("anim2.mp4", dpi=200)
#ani3 = FuncAnimation(fig, ani, frames=frames3(), interval=25, repeat=True, save_count=600)
#ani3.save("anim3.mp4", dpi=200)

plt.show()
