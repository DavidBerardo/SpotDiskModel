Coordinate system:
	+x points towards observer
	+y points to the right
	+z points up

Main script to import is spot_disk_crossing.py

Spot:
	params:
		theta: 
			latitude, from -pi to pi
		phi: 
			longitude, from -pi/2 to pi/2
		r: 
			projected radius if you were looking at the spot face on, in units of r_star, from 0 to 1
		temp: 
			temperature of spot in kelvin

Star:
	parameters:
		period: 
			rotation period of star. Can be in any units, but time for fluxes must be in same unit.
		theta,phi: 
			angles that determine where the stellar rotation pole points. 0,0 points in [0,0,1] in xyz
		temp:
			temperature of star in kelvin
		limb_coeffs:
			two coefficients for a quadratic limb darkening model
Disk:
	parameters:
		inclination:
			inclination of disk plane relative to observer in degrees. Follows exoplanet orbit conventions
		opacity:
			opaicty of disk. 0 = completely transparent, 1 = completely opaque
	
Model:
	res:
		resolution of theta elements. Phi elements is set to res / 2. Default is 300, which seems to be good to about 50 ppm
	bp:
		limits of bandpass, in nanometers. Given as a 1d array with two entries



The smallest time resolution is orbital period / (theta resolution)

Under current flood-fill way of drawing spots, the area of the spots varies by about 2% (might not even be this bad, should check further)

