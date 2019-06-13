import matplotlib.pyplot as plt

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

# plt.xlim(1e-3, 10)
# plt.ylim(1, 1.3)
plt.semilogx()

# r0, vx0 = readfile("r_vx0.dat")
# vxnorm = vx0[0]

# for i in range(0, 10001):
    # filelist.append("r_vx%s.dat" % i)

# for fname in filelist:
    # r, vx = readfile(fname)
    # plt.plot(r, np.array(vx)/vxnorm, '-')


plt.plot(t0, visc0, '-')
# plt.ylim(0.515, 0.525)
plt.show()
