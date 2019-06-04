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


t1, visc1 = readfile('viscoraw1.dat')
t2, visc2 = readfile('viscoraw2.dat')

# print(r)
# print(vx)

plt.plot(t1, visc1, '-bo')
plt.plot(t2, visc2, '-r')
# plt.xscale('log')
plt.xlim(0.005, 0.1)
plt.ylim(0.072, 0.0783)
plt.show()
