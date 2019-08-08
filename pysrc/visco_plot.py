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


plt.xlim(1e-3, 1.5)
# plt.ylim(1, 1.3)
t0, visc0 = readfile("visctot0.dat")
t1, visc1 = readfile("visctot9.dat")
plt.semilogx()
plt.plot(t0, visc0, '-', label='r=0.0')
plt.plot(t1, visc1, '-', label='r=0.9')
plt.title('viscosity as a funtion of t, a&o included')
# plt.xlabel('time (s)')
# plt.ylabel('viscosity')
plt.legend()
plt.show()
