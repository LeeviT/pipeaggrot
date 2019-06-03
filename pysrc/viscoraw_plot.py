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


r, vx = readfile('r_vx.dat')

print(r)
print(vx)

plt.plot(r, vx, '-bo')
# plt.xscale('log')
plt.xlim(0, 0.12)
plt.ylim(0, 2)
plt.show()
