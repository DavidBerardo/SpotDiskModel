from spot_disk_crossing import *

#make some random spots
spot1 = spot(theta = -np.pi + 2*np.pi*np.random.random(), phi = -np.pi/2.0 + np.pi*np.random.random(), r = 0.2, temp = 6100)
spot2 = spot(theta = -np.pi + 2*np.pi*np.random.random(), phi = -np.pi/2.0 + np.pi*np.random.random(), r = 0.2, temp = 5800)
spot3 = spot(theta = -np.pi + 2*np.pi*np.random.random(), phi = -np.pi/2.0 + np.pi*np.random.random(), r = 0.2, temp = 5800)

#compile the spots
spots = [spot1,spot2,spot3]
#make the star object
host_star = star(period = 0.5,temp=6000,limb_coeffs=[0.2,0.1],theta = 2*np.pi*np.random.random(),phi = np.pi*np.random.random())
#make the dist object
disk = disk(r1 = 30,r2 = 45,inclination = 89.5,opacity = 1)

#make the model
model = spd_model(res=400,bandpass=[400,900])

#assign attributes
model.spots = spots
model.host_star = host_star
model.disk = disk

#setup model
model.setup()

times = np.linspace(0,0.5,300)
lightcurve = model.lightcurve(times = times)

plt.plot(times,lightcurve)

model.show_all()
plt.show()
