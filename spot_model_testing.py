from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import sph_harm
import time
import sys

n = 500
m = 250

dtheta = 2 * np.pi / n
dphi = np.pi / m

def phi_index(phi):
	return int((phi + np.pi/2) / dphi)
def theta_index(theta):
	return int((theta + np.pi) / dtheta)

def get_theta(index):
	return -np.pi + (index%n)*dtheta
def get_phi(index):
	return -np.pi/2.0 + (index)*dphi

def disk_map():
	a1,b1 = 1,0.1
	a2,b2 = 1,-0.1

	y11 = (-2*a1*b1 + np.sqrt((2*a1*b1)**2 - 4 * (a1 + 1) * (b1**2 -1))) / (2 * (a1 + 1))
	y12 = (-2*a1*b1 - np.sqrt((2*a1*b1)**2 - 4 * (a1 + 1) * (b1**2 -1))) / (2 * (a1 + 1))

	y21 = (-2*a2*b2 + np.sqrt((2*a2*b2)**2 - 4 * (a2 + 1) * (b2**2 -1))) / (2 * (a2 + 1))
	y22 = (-2*a2*b2 - np.sqrt((2*a2*b2)**2 - 4 * (a2 + 1) * (b2**2 -1))) / (2 * (a2 + 1))

	z11 = a1*y11 + b1
	z12 = a1*y12 + b1
	
	z21 = a2*y21 + b2
	z12 = a2*y22 + b2

	phi_11 = np.arcsin(z11)
	phi_12 = np.arcisn(z12)

	#for i in range(phi_index(phi_11),phi_index(phi_12)):
	#	phi_cur = 



	return

def area_map():
	areas = []
	#area only depends on phi
	for i in range(m):
		phi = get_phi(i)
		areas.append(n * [np.cos(phi)*dphi*dtheta])
	return

def draw_spot(flux_grid,theta_0,phi_0,r):
	area = 0
	g = np.sqrt(1 - r**2)
	
	a = np.cos(theta_0) * np.cos(phi_0)
	b = np.sin(theta_0) * np.cos(phi_0)
	c = np.sin(phi_0)

	i_0 = theta_index(theta_0)
	j_0 = phi_index(phi_0)
	print(i_0,j_0)

	
	for j in range(m):
		for i in range(n):
			theta = get_theta(i)
			phi = get_phi(j)

			x = np.cos(theta) * np.cos(phi)
			y = np.sin(theta) * np.cos(phi)
			z = np.sin(phi)

			if a*x + b*y + c*z >= g:
				flux_grid[j][i] = 1
				area += np.cos(phi)*dphi*dtheta
	print(area)
	return
	
	def fill1(flux_grid,i_0,j_0,di,dy):
		i = i_0
		while True:
			theta = get_theta(i)
			j = j_0
			while True:
				phi = get_phi(j)

				x = np.cos(theta) * np.cos(phi)
				y = np.sin(theta) * np.cos(phi)
				z = np.sin(phi)

				if a*x + b*y + c*z >= g and 0 <= j < m:
					flux_grid[j][i] = 1
					j += int(dy)
				else:
					break
			i += int(di)
			if abs(i - i_0) == m or j == j_0:
				break
		return

	def fill2(flux_grid,i_0,j_0,dx,dy):
		j = j_0
		while True:
			phi = get_phi(j)
			i = i_0
			while True:
				theta = get_theta(i)

				x = np.cos(theta) * np.cos(phi)
				y = np.sin(theta) * np.cos(phi)
				z = np.sin(phi)

				if a*x + b*y + c*z >= g and not(abs(i - i_0) == n):
					flux_grid[j][i%n] = 1
					i += int(dx)
				else:
					break
			j += int(dy)
			if not(0 <= j < m) or i == i_0:
				break
		return

	fill2(flux_grid,i_0,j_0,1,1)
	fill2(flux_grid,i_0,j_0,-1,1)
	fill2(flux_grid,i_0,j_0,1,-1)
	fill2(flux_grid,i_0,j_0,-1,-1)
	print(area)

	return

def u_map():
	umap = np.zeros((m,n))
	for i in range(n):
		theta = -np.pi + dtheta * i
		if abs(theta) > np.pi/2:
			continue
		for j in range(m):
			phi = -np.pi/2 + dphi * j
			umap[j][i] = np.cos(theta) * np.cos(phi)
	return umap

flux_grid = np.zeros((m,n))

t0 = time.time()
#for i in range(5):
#draw_spot(flux_grid,0,(-3/4)*np.pi/2.0,0.2)
#draw_spot(flux_grid,0,(3/4)*np.pi/2.0,0.2)
draw_spot(flux_grid,0,(4/4)*np.pi/2.0,0.2)

print(time.time() - t0)

plt.imshow(flux_grid)
plt.colorbar()
plt.show()