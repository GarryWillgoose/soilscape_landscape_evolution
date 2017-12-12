#  Line load out of plane loading on a thin plate
#  Elastic support

#  This code is designed to be used in conjunction with the explanatory details
#  in Willgoose (2018) "Principles of Soilscape and Landscape Evolution", Cambridge University Press
#
#    Copyright (C) 2017  Garry Willgoose, The University of Newcastle
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#   length units km

# D is flexural rigidity
# K is the mantle bouyancy force


import numpy
import matplotlib.pyplot

# For figure 12.7 set to True
POINT = True

# For figure 12.6 set to False
#POINT = False

# Mauna Lau load spread over about 100km, 50km one side of centre 50 km the other.
LOAD_WIDTH_DISTRIBUTED = 50.e3          #  distributed load spread over 50km either side opf centre
# NARROW LINE LOAD in Y DIRECTION RATHER THAN STRICTLY A POINT LOAD
#   (this value needs to be greater than the resolution of the grid)
LOAD_WIDTH_POINT = 5.e3           #  point load spread over 5km either side of centre

#DX = 5.0e3                      #  grid spacing 5 km
DX = 1.0e3                      #  grid spacing 1 km
DY = 100.0e3                      #  grid spacing 100 km  (Mauna Lau load spread over about 100km)
#DOMAIN = 101
DOMAIN = 501
THICKNESS = 30.e3                  # nominal crust thickness 30km
DENSITY_CRUST = 2700.
DENSITY_MANTLE = 3300.
DENSITY_ICE = 1000.

#  10000 cubic km additonal (to the self weight of the crust) load (approx volume of Mauna Lau(N)
VOLUME = 10000.e9

#  plate flexural rigidity E 50.e9 Pa (granite) typical of most igneous rocks,
#  http://community.dur.ac.uk/~des0www4/cal/dams/geol/mod.htm
#  I 2e12 (30km thick)
E = 50.e9

#  these are the sensitivity test multipliers

# order is thickness,K (mantle density),p (load density)
#PARAMETERS = [[0.01,1,1],[1,1,1],[100,1,1],[1,0.5,1],[1,2,1],[1,0.1,1],[1,0.05,1],[1,0.02,1],[1,0.01,1]]
#   variability in  crustal thickness
PARAMETERS = [[0.333,1,1],[0.6666,1,1],[1,1,1],[2.16666,1,1],[3.3333,1,1]]
#   no  bouyancy
#PARAMETERS = [[0.333,0,1],[0.6666,0,1],[1,0,1],[2.16666,0,1],[3.3333,0,1]]
#   variable  bouyancy
#PARAMETERS = [[1,0.82,1],[1,0.91,1],[1,1,1],[1,1.09,1],[1,1.18,1]]
#   variable  bouyancy (1st number is the correctino for submerged weight (i.e. in the ocean)
#PARAMETERS = [[1,1,0.63],[1,1,0.82],[1,1,1],[1,1,1.09],[1,1,1.18]]
w = numpy.zeros((len(PARAMETERS),DOMAIN), 'float64')
p = numpy.zeros((len(PARAMETERS),DOMAIN), 'float64')
bouyancy = numpy.zeros((len(PARAMETERS),DOMAIN), 'float64')
amatrix = numpy.zeros((DOMAIN,DOMAIN),'float64')
wmink = numpy.zeros((len(PARAMETERS)),'float64')
wmaxk = numpy.zeros((len(PARAMETERS)),'float64')

for k in range(len(PARAMETERS)):
  pload=DENSITY_CRUST * VOLUME * 9.81* PARAMETERS[k][2]
  D = E * (THICKNESS* PARAMETERS[k][0])**3*DY/12.0
  DX4 = 1./DX**4
  DY4 = 1./DY**4
  DX2Y2 = 1./(DX**2*DY**2)
#  mantle bouyancy force
  K = DX*DY*DENSITY_MANTLE*9.81 * PARAMETERS[k][1]
  print 'D,K',D,K

#  w = numpy.zeros((DOMAIN), 'f')
#  p = numpy.zeros((DOMAIN), 'f')
  centre_x=int(DOMAIN*0.5)
  if POINT:
#  sensivitity case with width of 5km either way
    ix = int(round(LOAD_WIDTH_POINT/DX))
  else:
# Mauna Lau load spread over about 100km, 50km one 50 km the other.
    ix = int(round(LOAD_WIDTH_DISTRIBUTED/DX))
  span = ix*2+1
  p[k,:]=0
  for i in range(-ix,ix+1):
#  (assumed uniformly distributed for the example)
    p[k,centre_x+i] = -pload/span
  amatrix[:,:]=0
  for i in range(DOMAIN):
    amatrix[i,i]=6*DX4*D + K
    if i == 0:
      amatrix[i,i+1]=-4*DX4*D
      amatrix[i,i+2]=1*DX4*D
    elif i == 1:
      amatrix[i,i-1]=-4*DX4*D
      amatrix[i,i+1]=-4*DX4*D
      amatrix[i,i+2]=1*DX4*D
    elif i == DOMAIN-2:
      amatrix[i,i-2]=1*DX4*D
      amatrix[i,i-1]=-4*DX4*D
      amatrix[i,i+1]=-4*DX4*D
    elif i == DOMAIN-1:
      amatrix[i,i-2]=1*DX4*D
      amatrix[i,i-1]=-4*DX4*D
    else:
      amatrix[i,i-2]=1*DX4*D
      amatrix[i,i-1]=-4*DX4*D
      amatrix[i,i+1]=-4*DX4*D
      amatrix[i,i+2]=1*DX4*D

  amatrixinv = numpy.linalg.inv(amatrix)
  w[k,:]= numpy.dot(amatrixinv,p[k,:])
  bouyancy[k,:] = K*w[k,:]
  wmink[k] = numpy.min(w[k,:])
  wmaxk[k] = numpy.max(w[k,:])
print 'min and max deflections ',wmink, wmaxk

fig = matplotlib.pyplot.figure()
contents = fig.add_subplot(111)
for k in range(len(PARAMETERS)):
  stuff = contents.plot(w[k],label=str(int(round(PARAMETERS[k][0]*THICKNESS/1000.)))+'km')
contents.set_title('Deflection (m)\n thickness='+
      str(THICKNESS)+' mantle='+str(DENSITY_MANTLE)+' crust='+str(DENSITY_CRUST)+
      '\n(Min,Max) Max=('+str(min(wmaxk))+','+str(max(wmaxk))+')')
contents.legend(loc='lower right')
matplotlib.pyplot.xlabel('Distance (km)')
matplotlib.pyplot.ylabel('Incremental Deflection (m)')

fig2 = matplotlib.pyplot.figure()
contents2 = fig2.add_subplot(111)
for k in range(len(PARAMETERS)):
  stuff2 = contents2.plot(-bouyancy[k],label=str(int(round(PARAMETERS[k][0]*THICKNESS/1000.)))+'km')
contents2.set_title('Bouyancy force (N)')
contents2.legend(loc='upper right')
matplotlib.pyplot.xlabel('Distance (km)')
matplotlib.pyplot.ylabel('Net Bouyancy (N)')

matplotlib.pyplot.show()

