import matplotlib.pyplot as plt


def readfile(filename):
    with open(filename, 'r') as data:
        t = []
        visc = []
        for line in data:
            p = line.split()
            t.append(float(p[0]))
            visc.append(float(p[1]))

    return t, visc


t0, visc0 = readfile('r_vx0.dat')
t1, visc1 = readfile('r_vx1.dat')
t2, visc2 = readfile('r_vx2.dat')
# t2, visc2 = readfile('viscoraw2.dat')

# print(r)
# print(vx)

plt.plot(t0, visc0, '-g')
plt.plot(t1, visc1, '-b')
plt.plot(t2, visc2, '-r')
# plt.yscale('log')
#plt.xlim(0.0, 1.5)
plt.ylim(0, 0.4)
plt.show()
