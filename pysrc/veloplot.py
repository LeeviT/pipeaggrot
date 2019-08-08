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
r0, vx0 = readfile("r_vx1000.dat")
vxnorm = vx0[0]

for i in range(0, 1000):
    filelist.append("r_vx%s.dat" % i)

for fname in filelist:
    r, vx = readfile(fname)
    # plt.plot(r, np.array(vx)/vxnorm, '-', linewidth=1.3)

r, vx = readfile("r_vx1.dat")
plt.plot(r, np.array(vx)/vxnorm, '-', label="t=0.001 s")
r, vx = readfile("r_vx10.dat")
plt.plot(r, np.array(vx)/vxnorm, '-', label="t=0.01 s")
r, vx = readfile("r_vx50.dat")
plt.plot(r, np.array(vx)/vxnorm, '-', label="t=0.05 s")
r, vx = readfile("r_vx250.dat")
plt.plot(r, np.array(vx)/vxnorm, '-', label="t=0.25 s")
r, vx = readfile("r_vx1000.dat")
plt.plot(r, np.array(vx)/vxnorm, '-', label="t=1.0 s")
plt.title("dp = 10 MPa")
plt.xlabel('radius')
plt.ylabel('velocity')
plt.xlim(0.0, 1e-3)
plt.ylim(0, 1.03)
plt.legend()
plt.show()

