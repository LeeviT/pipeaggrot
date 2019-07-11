import matplotlib.pyplot as plt
import numpy as np

print(plt.style.available)
plt.style.use('bmh')
plt.rcParams['figure.dpi'] = 200


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
# r0, vx0 = readfile("r_vx10.dat")
# vxnorm = vx0[0]

for i in range(0, 9):
    filelist.append("r_vx%s.dat" % i)

# for fname in filelist:
#    r, vx = readfile(fname)
#    plt.plot(r, np.array(vx)/vxnorm, '-', linewidth=1.3)

# plt.plot(r0, np.array(vx0)/vxnorm, '--k', linewidth=1.3)

# plt.xlabel('radius')
# plt.ylabel('velocity')
# plt.xlim(0.0, 1.0)
# plt.show()

plt.xlim(1e-2, 2)
# plt.ylim(1, 1.3)
t0, visc0 = readfile("visctot0.dat")
t1, visc1 = readfile("visctot4.dat")
# plt.semilogx()
plt.plot(t0, visc0, '-', label='r=0.0')
plt.plot(t1, visc1, '-', label='r=0.8')
plt.title('viscosity as a funtion of t, only orientation included')
# plt.xlabel('time (s)')
# plt.ylabel('viscosity')
plt.legend()
plt.show()
