from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import time
from spot_disk_crossing import *
import sys

spots = []
for i in range(5):
	theta = -np.pi + 2*np.pi*np.random.random()
	phi = -np.pi/2 + np.pi *np.random.random()
	r = 0.2 + 0.1*np.random.random()
	temp = 6400 - 800 * np.random.random() 
	spots.append(spot(theta,phi,r,temp))


#host_star = star(period = 0.5,temp=6000,limb_coeffs=[0.2,0.1],theta = -np.pi/2.0,phi = np.pi/2.0)
#disk = disk(r1 = 30,r2 = 45,inclination = 89.5,opacity = 1)
times = np.linspace(0,1,500)


model = spd_model(res=300,bandpass=[400,900])
model.host_star = star(period = 0.5,temp=6000,limb_coeffs=[0.2,0.1],theta = 0,phi = 0)
model.spots = spots
model.disk = disk(r1 = 30,r2 = 45,inclination = 89.5,opacity = 1)

t0 = time.time()
model.setup()

lightcurve = model.lightcurve(times)
print(time.time()-t0)

plt.plot(times,lightcurve / np.median(lightcurve),label = '300')

model.show_all()
plt.show()