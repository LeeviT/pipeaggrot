import matplotlib.pyplot as plt
import numpy as np

print(plt.style.available)
plt.style.use('bmh')
plt.rcParams['figure.dpi'] = 300


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
# r0, vx0 = readfile("r_vx0.dat")
# vxnorm = vx0[0]

# for i in range(0, 30):
#    filelist.append("r_vx%s.dat" % i)

# for fname in filelist:
#    r, vx = readfile(fname)
#    plt.plot(r, np.array(vx)/vxnorm, '-')

plt.show()

plt.xlim(1e-3, 1)
# plt.ylim(1, 1.3)
t0, visc0 = readfile("visctot1.dat")
plt.semilogx()
plt.plot(t0, visc0, '-')
plt.title('viscosity as a funtion of t in pipe y-midpoint')
plt.xlabel('time (s)')
plt.ylabel('viscosity')
plt.show()
