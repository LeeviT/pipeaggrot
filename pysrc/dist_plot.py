import matplotlib.pyplot as plt
import numpy as np

plt.style.use('bmh')
plt.rcParams['figure.dpi'] = 200


def readfile(filename):
    with open(filename, 'r') as data:
        t = []
        phi = []
        theta = []
        for line in data:
            p = line.split()
            if not (float(p[0]) in t):
                t.append(float(p[0]))
            phi.append(float(p[1]))
            theta.append(float(p[2]))
    return t, phi, theta


time_steps = 5001
dt = 0.001
ny = 10

t0, phi0, theta0 = readfile("../cppsrc/dist/class_0.dat")
plt.xlim(dt, time_steps*dt)

phi_avg = []
phi_upper = []
phi_lower = []
theta_avg = []
theta_upper = []
theta_lower = []

for i in range(1, time_steps + 1):
    phi_sum = 0
    phi_max = 0
    phi_min = 10
    theta_sum = 0
    theta_max = 0
    theta_min = 10

    for j in range((i - 1) * (ny + 1), (i - 1) * (ny + 1) + ny):
        phi_sum += phi0[j]
        theta_sum += theta0[j]
        if phi0[j] > phi_max:
            phi_max = phi0[j]
        if theta0[j] > theta_max:
            theta_max = theta0[j]
        if phi0[j] < phi_min:
            phi_min = phi0[j]
        if theta0[j] < theta_min:
            theta_min = theta0[j]

    phi_avg.append(phi_sum / ny)
    phi_upper.append(phi_max)
    phi_lower.append(phi_min)
    theta_avg.append(theta_sum / ny)
    theta_upper.append(theta_max)
    theta_lower.append(theta_min)

print(phi_upper)
plt.xscale('log')
plt.xlim(1e-3, 1)
plt.plot(dt*np.array(t0), phi_avg)
plt.fill_between(dt*np.array(t0), phi_upper, phi_lower, facecolor='0.85', label="min-max-interval")
# plt.plot(dt*np.array(t0), theta_avg)
# plt.fill_between(dt*np.array(t0), theta_upper, theta_lower, facecolor='0.85', label="min-max-interval")
# plt.plot(t1, phi1, label="also aggr. incl.")
plt.xlabel("time (s)")
# plt.ylabel("phi")
# plt.title("average orientation (phi) as a function of timesteps")
plt.title("10 MPa")
plt.legend()
plt.show()
