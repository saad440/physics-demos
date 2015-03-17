#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  compton-scattering.py version 1.0
#  A demonstration of Compton Scattering using VPython and matplotlib
#  
#  Copyright 2014 Muhammad Saad <muhammad.saad@yandex.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


import visual
import random
import matplotlib.pyplot as plt


# Plank's constant
h = 6.62606957 * 10**(-34)
# Typical X-ray wavelength would be ~10**-9
# But to make the shift more visible, we use Gamma wavelengths
#lmda = 0.1 * 10**-9
lmda = 1 * 10**-12  # meters
c = 3.0 * 10**8  # meters/second
freq_o = (c / lmda)
# Scaling factor for lengths, to make the wavelength viewable
scalefactor_length = 10.0**20  # ratio of displayed length to actual length
freq_s = freq_o / scalefactor_length
# Scaling factor for time, to make movement viewable
scalefactor_time = 2 * 10.0**7  # ratio of elapsed time to actual time
# Time interval per step
dt = 0.05  # seconds
# Maximum times per second to render a frame when processing loops
animation_run_rate = 100

# Initialize the Scene
scene = visual.display(title="Compton Scattering (M. Saad)", x = 0, y = 0, width=800, height=600)
scene.autoscale = False
scene.range = (200,200,200)
scene.fullscreen = False

# The Photon's Frame of Reference
photon_frame = visual.frame()
# Initialize the Photon
photon = visual.curve( x = visual.arange(-200,-180), y=0, display=scene, frame=photon_frame, color=visual.color.blue)
photon.radius = 0.5
photon.wavelength_init = lmda  # meters
photon.freq_init = freq_o  # hertz
photon.energy_init = h * photon.freq_init  # joules
# Make it wavy
photon.y = 7*visual.sin((freq_s*photon.x)/(2*visual.pi))

# The Electron's Frame of Reference
electron_frame = visual.frame()
# Initialize the Electron
electron = visual.sphere(x=0, y=0, display=scene, frame=electron_frame, radius=10)
electron.mass = 9.10938291 * 10**(-31)  # kilograms
electron.energy_rest = electron.mass * (c**2)  # joules

# Define momenta
# Electron is at rest
electron.momentum_init = visual.vector(0,0,0)
# Photon is moving along X axis
photon.momentum_init = visual.vector(h/photon.wavelength_init,0,0)

# Display data for electron and photon
electron.label = visual.label(pos=(30,100,0), xoffset=1, height=15, border=10, font='STIX', box=False, text="Electron")
photon.label = visual.label(pos=(-160,100,0), xoffset=1, height=15, border=10, font='STIX', box=False, text="Photon")

# Plot the Graph of Wavelength Shift as a function of Scattering Angle
plt.xlim(-visual.pi, visual.pi)
plt.ylim( 0, 6 * 10**(-12))
plt.title('Wavelength Shift as a Function of Scattering Angle')
plt.xlabel('Scattering Angle (radians)')
plt.ylabel('Increase in Wavelength (x 10^-12 meters)')
plt.ion()
plt.show()
plt.draw()


def reset():
	# Reset everything to initial conditions
	electron.pos = (0,0,0)
	electron.momentum_fin = visual.vector(0,0,0)
	electron.velocity = electron.momentum_fin / electron.mass
	electron.speed = visual.mag( electron.velocity )
	electron.mass_rel = electron.mass
	electron.energy_rel = electron.mass_rel * (c**2)
	electron_frame.axis = (1,0,0)
	electron.label.text = u'Electron\nv = {0:.5} m/s\n'.format(electron.speed)
	
	photon.wavelength_fin = lmda
	photon.momentum_fin = visual.vector(h/photon.wavelength_fin,0,0)
	photon.freq_fin = ( c / photon.wavelength_fin )
	photon.energy_fin = h * photon.freq_fin
	photon.x = visual.arange(-200,-180)
	photon.y = 7*visual.sin((freq_s*photon.x)/(2*visual.pi))
	photon.color = visual.color.blue
	photon_frame.axis = (1,0,0)
	photon.label.text = u'Photon\nλ = {0:.5} meters\n∡ = {1:.5} π radians'.format(photon.wavelength_init,0.0)


def collide():
	# Move the photon right till it reaches the electron
	while photon.pos[-1,0] < electron.pos.x:
		visual.rate(animation_run_rate)
		photon.pos += [(c*dt)/scalefactor_time,0,0]
		photon.y = 7*visual.sin((freq_s*photon.x)/(2*visual.pi))
		photon.label.text = u'Photon\nλ = {0:.5} meters\n∡ = {1:.5} π radians'.format(photon.wavelength_init,0.0)	
	# Change color after collision to indicate "redshift"
	photon.color = visual.color.red
	# Deflect at random angle
	rand_angle = random.uniform(-visual.pi,visual.pi)
	# Rotate the photon's frame at that angle to avoid complicated calculations for drawing and velocity
	photon_frame.axis = (visual.cos(rand_angle), visual.sin(rand_angle), 0)
	
	# Compton wavelength of electron = h/mc = 2.4263102389 * 10**(-12)
	compton_wavelength_electron = 2.4263102389 * 10**(-12)  # meters
	# Compton Shift Formula
	del_lmda = compton_wavelength_electron * (1-visual.cos(rand_angle))
	
	lmda2 = lmda + del_lmda
	freq2 = ( c / lmda2 )
	freq_s2 = freq2 / ( scalefactor_length )
	
	# After collision
	photon.y = 7*visual.sin((freq_s2*photon.x)/(2*visual.pi))
	photon.wavelength_fin = lmda2
	photon.momentum_fin = visual.vector((h/photon.wavelength_fin)*visual.cos(rand_angle),(h/photon.wavelength_fin)*visual.sin(rand_angle),0)
	photon.freq_fin = ( c / photon.wavelength_fin )
	photon.energy_fin = h * photon.freq_fin
	photon.label.text = u'Photon\nλ = {0:.5} meters\n∡ = {1:.5} π radians'.format(photon.wavelength_fin, rand_angle/visual.pi)
	
	electron.momentum_fin = ( photon.momentum_init + electron.momentum_init ) - photon.momentum_fin
	electron.energy_rel = ( photon.energy_init - photon.energy_fin ) + electron.energy_rest
	electron.mass_rel = electron.energy_rel / (c**2)
	electron.velocity = electron.momentum_fin / electron.mass_rel
	electron.speed = visual.mag( electron.momentum_fin / electron.mass_rel )
	electr_angle = visual.arctan( electron.velocity.y / electron.velocity.x )
	electron.label.text = u'Electron\nv = {0:.5} m/s\n∡ = {1:.5} π radians'.format(electron.speed, electr_angle/visual.pi)
	
	# Move both the electron and scattered photon for a while
	for i in range(300):
		visual.rate(animation_run_rate)
		photon.pos += [(c*dt)/scalefactor_time,0,0]
		photon.y = 7*visual.sin((freq_s2*photon.x)/(2*visual.pi))
		electron.pos += (electron.velocity * dt)/scalefactor_time
	
	# Plot the data point
	plt.plot(rand_angle,del_lmda, 'ro')
	plt.draw()
	

reset()
# Perform a finite number of collisions
for i in range(100):
	collide()
	reset()
