import numpy as np
from matplotlib import pyplot as plt

target = [39.387247, -8.286478] # lat/lon

pos = [39.381157, -8.278333] # lat/lon
z = 1000 # altitude
vh = 5 # horizontal vel
vv = -3 # vertical vel
t = 0 
dt = 0.1
head = 33 * np.pi / 180 # heading, deg
rho = 1.225 # assume sea level for now
m = 3 # kg
g = 9.81 # m/s
bank = 0

# max limits
bankmax = 20 * np.pi / 180 
ratemax = 5 * np.pi / 180

# wind
wn = 3
we = -2

# pos in m
xn = (pos[0] - target[0]) * 111320
ye = (pos[1] - target[1]) * 111320 * np.cos(np.radians(target[0]))
# target is at origin
tn = 0
te = 0

# parafoil params
CL = 0.8
CD = 0.2
A = 1.5 # m^2

# control params
Kp = 1
Ki = 0
integ = 0
bias_known = 10 * np.pi / 180
bias_unknown = 10 * np.pi / 180

# plotting lists
lats = []
lons = []
zs = []
vs = []
ts = []
xns = []
yes = []
banks = []

# main loop - no control
while z > 0 and t < 1000:

    t = t + dt
    bias_unknown = t / 120 * (np.pi / 180) # constantly inceasing bias, 0.5deg/min
    dN = (vh * np.cos(head) + wn) * dt
    dE = (vh * np.sin(head) + we) * dt
    pos[0] += dN / 111320
    pos[1] += dE / (111320 * np.cos(np.radians(pos[0])))
    xn += dN
    ye += dE
    z = z + vv * dt

    # forces
    L = CL * A * 0.5 * rho * (vv ** 2 + vh ** 2)
    D = CD * A * 0.5 * rho * (vv ** 2 + vh ** 2)
    # resolve with angles
    gam = np.arctan2(vv, vh)
    Fv = L * np.cos(gam) * np.cos(bank + bias_unknown) - D * np.sin(gam)
    Fh = -L * np.sin(gam) * np.cos(bank + bias_unknown)- D * np.cos(gam)

    # explicit euler
    vv += (Fv / m - g) * dt
    vh += Fh / m * dt
    head += L * np.sin(bank + bias_unknown) / (m * vh) * dt

    # add to lists
    lats.append(pos[0])
    lons.append(pos[1])
    zs.append(z)
    vs.append(np.sqrt(vv ** 2 + vh ** 2))
    ts.append(t)
    xns.append(xn)
    yes.append(ye)
    banks.append(bank)

    # control
    rn = tn - xn # can get these from GPS lat/lon via the conversions above
    re = te - ye
    htarg = np.arctan2(re, rn)
    dist = np.hypot(rn, re) + 1e-6 # hypotenuse for magnitude
    # head we can either get from GPS or from differentiating
    err = (htarg - head + np.pi) % (2 * np.pi) - np.pi # wrapping
    
    com_bank = Kp * err + Ki * integ + bias_known
    if abs(com_bank) > bankmax:
        integ += err * dt
    com_bank = np.clip(com_bank, -bankmax, bankmax)
    # rate limiter
    bank += np.clip(com_bank - bank, -ratemax * dt, ratemax * dt)
    # print(f"pos{pos} z{z} vv{vv} vh{vh} t{t}s")

# matplotlib
fig, axs = plt.subplots(2, 2)
axs[0,0].plot(lons, lats)
axs[0,0].set_xlabel("lon")
axs[0,0].set_ylabel("lat")
axs[0,1].plot(ts, zs)
axs[0,1].set_xlabel("t")
axs[0,1].set_ylabel("alt")
axs[1,0].plot(ts, np.multiply(banks, 180/np.pi))
axs[1,0].set_xlabel("t")
axs[1,0].set_ylabel("bank")

axs[1,1].remove()
ax3d = fig.add_subplot(2,2,4,projection='3d')
ax3d.plot(xns, yes, zs)
ax3d.set_xlabel("lon")
ax3d.set_ylabel("lat")
ax3d.set_zlabel("alt")
plt.show()