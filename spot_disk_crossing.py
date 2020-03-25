from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import time
from scipy.spatial.transform import Rotation as Rot
from scipy.integrate import quad

#TODO 
# - put in realistic temperatures/fluxes
# - disk line
# - deal with time steps

class spot:
	""" 
    Class for holding a spot object 
      
    Attributes: 
        theta (float): latitude of stellar spot, varies from -pi to pi
        phi (float): longitude of spot, varies from -pi/2 to pi/2
        r (float): projected radius of spot when looking face on, units of r_star, goes from 0 to 1
		temp: temperature of the spot in Kelvin, gets converted to flux
    """

	def __init__(self,theta,phi,r,temp):
		self.theta = theta
		self.phi = phi
		self.r = r
		#calculate proper flux based on temperature
		self.temp = temp
		#parameters used for drawing the spot
		self.g = np.sqrt(1 - self.r**2)
		self.a = np.cos(self.theta) * np.cos(self.phi)
		self.b = np.sin(self.theta) * np.cos(self.phi)
		self.c = np.sin(self.phi)

class star:
	""" 
    Class for holding a star object 
      
    Attributes:
    	period (float): orbital period in days
    	temp (float): surface temperature of star, gets converted to flux
		limb_coeffs (1d array of floats): c1, c2, limb darkening coefficients for a law of the form 1 - c1*(1-u) - c2*(1-u)^2
		theta (float): on angle for orienting the rotation axis
		phi (float): other angle for orienting the rotation axis
    """
	def __init__(self,period,temp,limb_coeffs,theta,phi):
		self.period = period
		self.temp = temp
		self.limb_coeffs = limb_coeffs
		self.theta = theta 
		self.phi = phi

class disk:
	""" 
    Class for holding a disk object
      
    Attributes:
    	r1 (float): inner edge of disk, in stellar radii
    	r2 (float): outer edge of disk, in stellar radii
    	inc (float): inclination of disk, relative to observer
    	opacity (float): opacity of disk. 0 - 1, 1 == completely opaque, 0 == completely transparent
    """
	def __init__(self,r1,r2,inclination,opacity):
		self.r1 = r1
		self.r2 = r2
		self.inc = inclination * np.pi/180.0
		self.opacity = opacity

class spd_model:
	def __init__(self,res,bandpass):

		self.n_phi = int(res/2.0)
		self.n_theta = res

		self.dtheta = 2 * np.pi / self.n_theta
		self.dphi = np.pi / self.n_phi
		self.bp = bandpass

	#everything to get the model ready to calculate fluxes 
	def setup(self):
		#convert temperature to fluxes, relative to the stars flux
		self.calc_specific_fluxes()

		#make underlying surface map
		self.draw_surface()

		#get matrix to rotate star 
		self.calc_rotation_matrix()
		#print(time.time() - t0)

		#get surface normals and rotate them so that rotation axis of star points the right way
		self.calculate_normal_vectors()

		#calculate the areas and projection factors
		self.calc_geometry()

		#calculate the limb darkening map
		self.calc_limb_darkening()
		#print(time.time() - t0)

		#calculate the disk map
		self.calc_disk_map()

		#paints the spots onto the surface of the star
	def draw_surface(self):
		#the background star flux is normalized to 1
		self.surface_flux = np.ones((self.n_phi,self.n_theta))
		def fill(s,dx,dy):
			i_0 = self.theta_index(s.theta)
			j_0 = self.phi_index(s.phi)
			
			j = j_0
			while True:
				phi = self.get_phi(j)
				i = i_0
				while True:
					theta = self.get_theta(i)

					x = np.cos(theta) * np.cos(phi)
					y = np.sin(theta) * np.cos(phi)
					z = np.sin(phi)

					if s.a*x + s.b*y + s.c*z >= s.g and not(abs(i - i_0) == self.n_theta+1):
						self.surface_flux[j][i%self.n_theta] = s.flux
						i += dx
					else:
						break
				j += dy
				if not(0 <= j < self.n_phi) or i == i_0:
					break
			return 

		
		for s in self.spots:
			fill(s,1,1)
			fill(s,-1,1)
			fill(s,1,-1)
			fill(s,-1,-1)

	def calc_specific_fluxes(self):
		#stellar flux
		u1 = 0.01439 / (self.host_star.temp * self.bp[0])
		u2 = 0.01439 / (self.host_star.temp * self.bp[1])

		self.host_star.flux = self.band_flux(self.host_star.temp)
		for s in self.spots:
			s.flux = self.band_flux(s.temp) / self.host_star.flux
		return

	#calculate the flux of blacbody in a given waveband
	def band_flux(self,temp):
		def integrand(x):
			return x**3 / (np.exp(x)-1.0)
		u1 = 0.01439 / temp / (self.bp[0] * 1e-9) 
		u2 = 0.01439 / temp / (self.bp[1] * 1e-9) 
		I = quad(integrand,u2,u1)
		return I[0]

	#rotate everything into the correct frame of reference re the stellar rotation axis
	#theta phi = 0,0, rotation axis is (0,0,1)
	def calc_rotation_matrix(self):
		x = np.cos(self.host_star.theta)*np.sin(self.host_star.phi) 
		y = np.sin(self.host_star.theta)*np.sin(self.host_star.phi) 
		z = np.cos(self.host_star.phi)
		
		if z == 1:
			self.rotation_matrix = Rot.from_rotvec(self.host_star.theta * np.asarray([0,0,1]))
		elif z == -1:
			self.rotation_matrix = Rot.from_rotvec(np.pi * np.asarray([0,1,0]))
		else:
			n_mag = np.sqrt(x**2 + y**2) 
			theta = np.arcsin(n_mag) 
			if z < 0: 
				theta = np.pi - theta 

			self.rotation_matrix = Rot.from_rotvec(theta / n_mag * np.asarray([y,x,0]))
		return

	#go across the sphere and calculate the normal vectors, rotated into correct plane
	def calculate_normal_vectors(self):
		normal_vectors = []
		for i in range(self.n_phi):
			phi = self.get_phi(i)
			cphi = np.cos(phi)
			sphi = np.sin(phi)
			for j in range(self.n_theta):
				theta = self.get_theta(j)
				x,y,z = cphi * np.cos(theta), cphi * np.sin(theta), sphi
				normal_vectors.append(np.asarray([x,y,z]))

		normal_vectors = self.rotation_matrix.apply(normal_vectors)
		self.normal_vectors = normal_vectors.reshape((self.n_phi,self.n_theta,3))
		return

	#calculate the u variable, the angle between normal and observer line of sight
	#used for limb darkening, and also projecting the surface areas
	#also calculate the areas
	def calc_geometry(self):
		self.umap = self.normal_vectors[:,:,0].clip(min = 0)
		
		self.areas = np.zeros((self.n_phi,self.n_theta))
		for j in range(self.n_phi):
			phi = self.get_phi(j)
			self.areas[j][:] = np.sin(phi) * self.dtheta * self.dphi

		return

	def calc_limb_darkening(self):
		self.limb_darkening = 1 - self.host_star.limb_coeffs[0] * (1 - self.umap) - self.host_star.limb_coeffs[1] * (1 - self.umap)**2
		return

	def calc_(self):
		self.projection = np.ones((self.n_phi,self.n_theta))
		return 

	#calculated the parts of the map that are occulted by the disk
	#is always a straight horizontal line from the observer point of view
	def calc_disk_map(self):
		b1 = -self.disk.r1*np.cos(self.disk.inc)
		b2 = -self.disk.r2*np.cos(self.disk.inc)

		z = self.normal_vectors[:,:,2]
		self.disk_map = np.logical_or(z>b1,z<b2).astype(int) + np.logical_and(z<b1,z>b2).astype(int) * (1 - self.disk.opacity)
		return

	#utilitiy functions to convert between angles and indices
	def get_theta(self,index):
		return -np.pi + (index%self.n_theta)*self.dtheta
	def get_phi(self,index):
		return -np.pi/2.0 + (index)*self.dphi
	def phi_index(self,phi):
		return int((phi + np.pi/2) / self.dphi)
	def theta_index(self,theta):
		return int((theta + np.pi) / self.dtheta)


	def flux(self,time):
		#figure out how much to shift the theta map by
		roll_ind = int(np.round(time  * (self.n_theta / self.host_star.period)))
		flux = np.sum(np.roll(self.surface_flux,roll_ind,axis=1) * self.limb_darkening * self.disk_map * self.umap * self.areas)
		return flux

	def lightcurve(self,times):
		return [self.flux(i) for i in times]

	def show_map(self,image = None):
		if image is None:
			plt.imshow(self.surface_flux)
		else:
			plt.imshow(image)
		plt.colorbar()

	def show_all(self):
		fig = plt.figure()
		plt.subplot(2,2,1)
		plt.imshow(self.surface_flux)
		plt.title("spots")
		plt.xlabel('θ')
		plt.ylabel('ϕ')

		plt.subplot(2,2,2)
		plt.imshow(self.limb_darkening)
		plt.title("limb_darkening")
		plt.xlabel('θ')
		plt.ylabel('ϕ')

		plt.subplot(2,2,3)
		plt.imshow(self.disk_map)
		plt.title("disk")
		plt.xlabel('θ')
		plt.ylabel('ϕ')
		plt.colorbar()
		
		plt.subplot(2,2,4)
		plt.imshow(self.surface_flux * self.limb_darkening * self.disk_map * self.umap)
		plt.title("All")
		plt.xlabel('θ')
		plt.ylabel('ϕ')
		plt.colorbar()

