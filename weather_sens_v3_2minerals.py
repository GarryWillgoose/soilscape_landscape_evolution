
#  This is the transport limited weathering simulator with the added weathering product state
#  for 2 minerals simultaneously weathering
#
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


import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import time
import os
import os.path
import multiprocessing

global zaxis,a_state_store,s_state_store,l_state_store,r_rate_store
global last_out,last_outk,s_state_ts,s_statemean_ts
global s_state_store_2,s_state_ts_2,s_statemean_ts_2,r_rate_store_2

#========================================
#========================================
#  CHANGE RUN HERE
#NEWDIR = 'bioturb 0.005 depletion test-parallel1 print'.lower()
NEWDIR = 'k1 depletion print 2minerals 1'.lower()
#========================================
#========================================



def simulation(factor,nominal_index, nominal):
  a_state = numpy.zeros((Z_DIM+2),dtype = 'f')
  a_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  l_state = numpy.zeros((Z_DIM+2),dtype = 'f')
  l_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  r_rate = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state_2 = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state_new_a_2 = numpy.zeros((Z_DIM+2),dtype = 'f')
  r_rate_2 = numpy.zeros((Z_DIM+2),dtype = 'f')
  temp_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state_ts = numpy.zeros((len(z_index),len(TIMES_OUT)), dtype = 'f')
  s_statemean_ts = numpy.zeros((len(TIMES_OUT)), dtype = 'f')
  s_state_ts_2 = numpy.zeros((len(z_index),len(TIMES_OUT)), dtype = 'f')
  s_statemean_ts_2 = numpy.zeros((len(TIMES_OUT)), dtype = 'f')
  last_out = TIMES_OUT[0]
  last_outk = last_out

  if TIME_SERIES_RUN:
    a_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    s_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    s_state_store_2 = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    l_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    r_rate_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    r_rate_store_2 = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  else:
    a_state_store = numpy.zeros((Z_DIM+2),dtype = 'f')
    s_state_store = numpy.zeros((Z_DIM+2),dtype = 'f')
    s_state_store_2 = numpy.zeros((Z_DIM+2),dtype = 'f')
    l_state_store = numpy.zeros((Z_DIM+2),dtype = 'f')
    r_rate_store = numpy.zeros((Z_DIM+2),dtype = 'f')
    r_rate_store_2 = numpy.zeros((Z_DIM+2),dtype = 'f')

  V = V_NOMINAL
  DIFFUSIVITY = DIFFUSIVITY_NOMINAL
  K1 = K1_NOMINAL
  K2 = K2_NOMINAL
  K3 = K3_NOMINAL
  ERODE = ERODE_NOMINAL
  LEQUILIB = LEQUILIB_NOMINAL
  DZ = DZ_NOMINAL
  A_BC1 = A_BC_NOMINAL
#  etarun=False
  if (NEWDIR == 'v' or 'v leachate' in NEWDIR or 'v const mass' in NEWDIR or
        'v depletion' in NEWDIR or 'v ' in NEWDIR):
    V = V_NOMINAL*factor
    if 'v const mass' in NEWDIR:
      A_BC1 = A_BC_NOMINAL/factor
  elif NEWDIR == 'dispersivity':
    DIFFUSIVITY = DIFFUSIVITY_NOMINAL*factor
  elif (NEWDIR == 'k1' or 'k1 leachate' in NEWDIR or 'k1 ' in NEWDIR or
        'k1 fixedco2' in NEWDIR or 'k1 leachate fixedco2' in NEWDIR):
    K1 = K1_NOMINAL*factor
    K1_2 = K1_NOMINAL*factor*4
  elif NEWDIR == 'k2' or 'k2 depletion' in NEWDIR:
    K2 = K2_NOMINAL*factor
  elif NEWDIR == 'k3':
    K3 = K3_NOMINAL*factor
  elif (NEWDIR == 'erode' or 'erode ' in NEWDIR or
        'erode fixedco2' in NEWDIR or 'erode leachate fixedco2' in NEWDIR or
        'erode depletion' in NEWDIR):
    ERODE = ERODE_NOMINAL*factor
  elif NEWDIR == 'lequilib':
    LEQUILIB = LEQUILIB_NOMINAL*factor
  elif 'depth' in NEWDIR:
    DZ = DZ_NOMINAL*factor
  elif NEWDIR == 'eta':
    pass
#    etarun=True
  elif 'erode-velocity depletion' in NEWDIR:
    ERODE = ERODE_NOMINAL*((factor-1)*2+1.0)
    V = V_NOMINAL*factor
  elif 'timeseries' in NEWDIR or 'numerics' in NEWDIR:
    pass
  elif 'lumped_test' in NEWDIR:
    V = V_NOMINAL*factor
    ERODE=0
  elif 'bioturb' in NEWDIR:
    BIOTURB = BIOTURB_NOMINAL*factor
  else:
    print 'ERROR: wrong sensitivity mode selected'
    import sys
    sys.exit()
  DIFFUSIVITY_1 = DIFFUSIVITY/DZ**2
  if bioturb_active:
    bioturb_1 = BIOTURB/DZ**2
  V_1 = V/DZ
  ERODE_1 = ERODE/DZ

  a_state[:] = 0
  a_state_new_a[:] = 0
  s_state[:] = 0
  s_state_new_a[:] = 0
  s_state_2[:] = 0
  s_state_new_a_2[:] = 0
  l_state[:] = 0
  l_state_new_a[:] = 0

  a_state[0] = a_state[1] = A_BC1
  a_state_new_a[0] = a_state_new_a[1] = A_BC1
  s_state[:] = S_BC
  s_state_new_a[:] = S_BC
  s_state_2[:] = S_BC
  s_state_new_a_2[:] = S_BC
  l_state[:] = L_BC
  l_state_new_a[:] = L_BC

  for t in range(1,T_DIM+1):
    a_state[end-1] = 2*a_state[end-2]-a_state[end-3]
    s_state_new_a[0] = s_state[1]
    s_state_new_a_2[0] = s_state_2[1]
#     the follwing statement causes numerical problems for certain parameters
#    l_state[end-1] = 2*l_state[end-2]-l_state[end-3]
    temp_a[1:Z_DIM+1] = (K1 * a_state[1:Z_DIM+1] * s_state[1:Z_DIM+1]*
                          (LEQUILIB-l_state[1:Z_DIM+1])/LEQUILIB)
    temp_a_2[1:Z_DIM+1] = (K1_2 * a_state[1:Z_DIM+1] * s_state_2[1:Z_DIM+1]*
                          (LEQUILIB-l_state[1:Z_DIM+1])/LEQUILIB)
    if etarun:
      for i in range(1,Z_DIM+1):
        temp_a[i] = temp_a[i]*(0.5*Z_DIM/(i+1))**factor
        temp_a_2[i] = temp_a_2[i]*(0.5*Z_DIM/(i+1))**factor
    a_state_new_a[1:Z_DIM+1] = a_state[1:Z_DIM+1] + DT*(
                -V_1 * (a_state[1:Z_DIM+1]-a_state[0:Z_DIM])
                +DIFFUSIVITY_1 *
                      (a_state[2:Z_DIM+2]-2*a_state[1:Z_DIM+1]+a_state[0:Z_DIM])
                -K2 * (temp_a[1:Z_DIM+1]+temp_a_2[1:Z_DIM+1]))
    if bioturb_active:
      s_state_new_a[1:Z_DIM+1] = s_state[1:Z_DIM+1] + DT*(
                ERODE_1 * (s_state[2:Z_DIM+2]-s_state[1:Z_DIM+1])
                +bioturb_1 *
                      (s_state[2:Z_DIM+2]-2*s_state[1:Z_DIM+1]+s_state[0:Z_DIM])
                -temp_a[1:Z_DIM+1])
      s_state_new_a_2[1:Z_DIM+1] = s_state_2[1:Z_DIM+1] + DT*(
                ERODE_1 * (s_state_2[2:Z_DIM+2]-s_state_2[1:Z_DIM+1])
                +bioturb_1 *
                      (s_state_2[2:Z_DIM+2]-2*s_state_2[1:Z_DIM+1]+s_state_2[0:Z_DIM])
                -temp_a_2[1:Z_DIM+1])
    else:
      s_state_new_a[1:Z_DIM+1] = s_state[1:Z_DIM+1] + DT*(
                ERODE_1 * (s_state[2:Z_DIM+2]-s_state[1:Z_DIM+1])
                -temp_a[1:Z_DIM+1])
      s_state_new_a_2[1:Z_DIM+1] = s_state_2[1:Z_DIM+1] + DT*(
                ERODE_1 * (s_state_2[2:Z_DIM+2]-s_state_2[1:Z_DIM+1])
                -temp_a_2[1:Z_DIM+1])

    s_state_new_a.clip(min=0.000001)
    s_state_new_a_2.clip(min=0.000001)
    l_state_new_a[1:Z_DIM] = l_state[1:Z_DIM] + DT*(
                -V_1 * (l_state[1:Z_DIM]-l_state[0:Z_DIM-1])
                +DIFFUSIVITY_1 *
                      (l_state[2:Z_DIM+1]-2*l_state[1:Z_DIM]+l_state[0:Z_DIM-1])
                +K3 * (temp_a[1:Z_DIM]+temp_a_2[1:Z_DIM]))
    l_state_new_a[Z_DIM] = l_state[Z_DIM] + DT*(
                -V_1 * (l_state[Z_DIM]-l_state[Z_DIM-1])
                +DIFFUSIVITY_1 * (-l_state[Z_DIM]+l_state[Z_DIM-1])
                +K3 * (temp_a[Z_DIM]+temp_a_2[Z_DIM]))
    r_rate = temp_a.clip(min=1.001e-3)
    r_rate_2 = temp_a_2.clip(min=1.001e-3)
    invlgth = 1.0/(len(a_state[1:Z_DIM+1]))
    a_diff = a_state[1:Z_DIM+1]-a_state_new_a[1:Z_DIM+1]
    a_sls = numpy.dot(a_diff,a_diff)*invlgth
    s_diff = s_state[1:Z_DIM+1]-s_state_new_a[1:Z_DIM+1]
    s_sls = numpy.dot(s_diff,s_diff)*invlgth
    s_diff_2 = s_state_2[1:Z_DIM+1]-s_state_new_a_2[1:Z_DIM+1]
    s_sls_2 = numpy.dot(s_diff_2,s_diff_2)*invlgth
    l_diff = l_state[1:Z_DIM+1]-l_state_new_a[1:Z_DIM+1]
    l_sls = numpy.dot(l_diff,l_diff)*invlgth
    a_state = a_state_new_a.copy()
    s_state = s_state_new_a.copy()
    s_state_2 = s_state_new_a_2.copy()
    l_state = l_state_new_a.copy()
#    l_state[Z_DIM+1] = l_state[Z_DIM]   #  lowerBC in dispersiivity
#    if a_sls < 1.e-14 and s_sls < 1.e-14 and l_sls < 1.e-14:
    if (a_sls < 1.e-14 and s_sls < 1.e-14 and l_sls < 1.e-14) or t == T_DIM:
      if not TIME_SERIES_RUN:
        a_state_store[:] = a_state[:]
        s_state_store[:] = s_state[:]
        s_state_store_2[:] = s_state_2[:]
        l_state_store[:] = l_state[:]
        r_rate_store[:] = r_rate[:]
        r_rate_store_2[:] = r_rate_2[:]
      time_done=t
      print 'time finished:',t,a_sls,s_sls,s_sls_2,l_sls
      break
    if t in TIMES_OUT:
#      print 'time:',t,a_sls,s_sls,l_sls
      j = TIMES_OUT.index(t)
      if TIME_SERIES_RUN:
        time_done=t
        last_out = j
        last_outk = j
        a_state_store[j,:] = a_state[:]
        s_state_store[j,:] = s_state[:]
        s_state_store_2[j,:] = s_state_2[:]
        l_state_store[j,:] = l_state[:]
        r_rate_store[j,:] = r_rate[:]
        r_rate_store_2[j,:] = r_rate_2[:]
      s_statemean_ts[j] = numpy.sum(s_state[:])/len(s_state[:])
      s_statemean_ts_2[j] = numpy.sum(s_state_2[:])/len(s_state_2[:])
      last_outk = j
      if nominal:
        time_done=t
        last_out = j
        for i in range(len(z_index)):
          s_state_ts[i,j] = min(s_state[z_index[i]],S_BC)
          s_state_ts_2[i,j] = min(s_state_2[z_index[i]],S_BC)
  return(a_state_store, s_state_store, l_state_store, r_rate_store,
          last_out, last_outk, s_state_ts, s_statemean_ts, s_state_store_2,
          s_state_ts_2, s_statemean_ts_2, r_rate_store_2)


def loop(kkk):
  global zaxis,a_state_store,s_state_store,l_state_store,r_rate_store
  global last_out,last_outk,s_state_ts,s_statemean_ts
  global s_state_store_2,s_state_ts_2,s_statemean_ts_2,r_rate_store_2
  
  time_start = time.time()
  print 'SIM=',kkk,'\nstart:',time.ctime()
  nominal = (kkk == NOMINAL_INDEX)
  fresult = simulation(SENSITIVITY[kkk],kkk,nominal)
  if TIME_SERIES_RUN:
    a_state_store = fresult[0]
    s_state_store = fresult[1]
    s_state_store_2 = fresult[8]
    l_state_store = fresult[2]
    r_rate_store = fresult[3]
    r_rate_store_2 = fresult[11]
  else:
    a_state_store[kkk,:] = fresult[0]
    s_state_store[kkk,:] = fresult[1]
    s_state_store_2[kkk,:] = fresult[8]
    l_state_store[kkk,:] = fresult[2]
    r_rate_store[kkk,:] = fresult[3]
    r_rate_store_2[kkk,:] = fresult[11]
  if nominal:
    last_out = fresult[4]
    s_state_ts = fresult[6]
    s_state_ts_2 = fresult[9]
  last_outk[kkk]=fresult[5]
  s_statemean_ts[kkk,:] = fresult[7]
  s_statemean_ts_2[kkk,:] = fresult[10]
  time_end = time.time()
  print 'end: CPU=',time_end-time_start
  return(kkk)


def loop1(kkk):
  global zaxis,a_state_store,s_state_store,l_state_store,r_rate_store
  global last_out,last_outk,s_state_ts,s_statemean_ts
  global s_state_store_2,s_state_ts_2,s_statemean_ts_2,r_rate_store_2
  
  time_start = time.time()
  print 'SIM=',kkk,'\nstart:',time.ctime()
  nominal = (kkk == NOMINAL_INDEX)
  fresult = simulation(SENSITIVITY[kkk],kkk,nominal)
  time_end = time.time()
  print 'end: CPU=',time_end-time_start
  return(fresult)




if 'timeseries' in NEWDIR:
  print 'TIMESERIES RUN: ',NEWDIR
elif 'numerics' in NEWDIR:
  print 'NUMERICS RUN: ', NEWDIR
else:
  print 'SENSITIVITY RUN: '+NEWDIR
print 40*'-'
VERSIONDIR = 'figures'
LOCAL = VERSIONDIR+os.sep+NEWDIR+os.sep
if os.path.exists(LOCAL):
  pass
#  import shutil
#  shutil.rmtree(LOCAL, ignore_errors=True)
else:
  os.mkdir(LOCAL[:-1])
#Z_DIM = 100
#T_DIM = 100
Z_DIM = 50
#T_DIM_S = 2000
T_DIM_S = 4000
NUM_TIMES = 200
if 'numerics' in NEWDIR:
  T_DIM_S = 1000
  NUM_TIMES = 1
elif 'lumped_test' in NEWDIR:
  T_DIM_S = 1000
  NUM_TIMES = 100

T_DIM = T_DIM_S*NUM_TIMES
TIMES_OUT = []
for i in range(NUM_TIMES+1):
  TIMES_OUT.append(i*T_DIM_S)
last_out = TIMES_OUT[0]

if 'timeseries' in NEWDIR or 'numerics' in NEWDIR:
  TIME_SERIES_RUN = True
else:
  TIME_SERIES_RUN = False

if NEWDIR == 'eta':
  etarun=True
else:
  etarun=False

V = .1
V_NOMINAL = V
#DIFFUSIVITY = 0.01
DIFFUSIVITY = 0.000001
DIFFUSIVITY_NOMINAL = DIFFUSIVITY
if 'depletion' in NEWDIR:
  K1 = 0.5
else:
  K1 = 0.1
K1_NOMINAL = K1
K2 = 2.0
#  Pedometrics
#K2 = 1.0
K2_NOMINAL = K2
K3 = 1.0
K3_NOMINAL = K3
ERODE = 0.1
ERODE_NOMINAL = ERODE
BIOTURB = 0.005
BIOTURB_NOMINAL = BIOTURB
bioturb_active = False
if 'bioturb' in NEWDIR:
  bioturb_active = True


LEQUILIB = 1000
if 'leachate' in NEWDIR:
  LEQUILIB = 0.4
LEQUILIB_NOMINAL = LEQUILIB

A_BC = 1.0
A_BC_NOMINAL = A_BC
if 'depletion' in NEWDIR:
  S_BC = 0.25
#  pedometrics
#  S_BC = 0.5
else:
  S_BC = 1.0
#S_BC = 0.25
L_BC = 0.0

DT = 0.001                     #  timestep
DZ = 2.0/(Z_DIM-1)          #  2m deep soil
DZ_NOMINAL = DZ          #  2m deep soil

if 'detail' in NEWDIR:
  SENSITIVITY = [0.8, 0.9, 1.0, 1.1, 1.2]        #  more detail around the nominal value
else:
  SENSITIVITY = [1.0]        #  most  runs
#  SENSITIVITY = [0.5, 0.666, 1.0, 1.5, 2.0]        #  most  runs
NOMINAL_INDEX = 2         # which index of the sensitivity runs to use for the time series plots

if TIME_SERIES_RUN:
  SENSITIVITY = [1.0]       # for time series only want to do it for nominal parameters
  NOMINAL_INDEX = 0         # which index of the sensitivity runs to use for the time series plots
if 'dispersivity' in NEWDIR:
  SENSITIVITY = [1.0, 1000, 10000, 100000, 1000000]        #  dispersivity
  NOMINAL_INDEX = 2         # which index of the sensitivity runs to use for the time series plots
if 'bioturb' in NEWDIR:
#  SENSITIVITY = [50]        #  dispersivity
#  NOMINAL_INDEX = 0         # which index of the sensitivity runs to use for the time series plots
  SENSITIVITY = [0.0, 1.0, 5, 10, 50, 100]        #  dispersivity
  NOMINAL_INDEX = 2         # which index of the sensitivity runs to use for the time series plots
if 'k2 depletion' in NEWDIR:
  SENSITIVITY = [1.0, 1.5, 2, 3, 4]        #  k2 depletion test
  NOMINAL_INDEX = 2         # which index of the sensitivity runs to use for the time series plots
if 'eta' == NEWDIR:
  SENSITIVITY = [0,0.5,1.0, 2]         # for area sensitivity runs
  NOMINAL_INDEX = 2         # which index of the sensitivity runs to use for the time series plots

#  SENSITIVITY = [1.8,2.05, 2.06, 2.1, 2.4]        #  k2 runs
#  NOMINAL_INDEX = 1         # which index of the sensitivity runs to use for the time series plots

LINE_STYLE = ['--', '-', ':', '-.', '-','--', '-', ':', '-.', '-']
LINE_WIDTH = [2.0, 2.0, 2.0, 2.0, 2.0,2.0, 2.0, 2.0, 2.0, 2.0]

#  make states 2 nodes bigger so can apply flux BC
#   [0] is the min BC node,
#   [1] is the first computational node
#   [Z_DIM] is the last computational node
#   [Z_DIM+1] is the max BC node
a_state = numpy.zeros((Z_DIM+2),dtype = 'f')
a_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
s_state = numpy.zeros((Z_DIM+2),dtype = 'f')
s_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
s_state_2 = numpy.zeros((Z_DIM+2),dtype = 'f')
s_state_new_a_2 = numpy.zeros((Z_DIM+2),dtype = 'f')
l_state = numpy.zeros((Z_DIM+2),dtype = 'f')
l_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
r_rate = numpy.zeros((Z_DIM+2),dtype = 'f')
r_rate_2 = numpy.zeros((Z_DIM+2),dtype = 'f')
temp_a = numpy.zeros((Z_DIM+2),dtype = 'f')
temp_a_2 = numpy.zeros((Z_DIM+2),dtype = 'f')

z_index = [1,int(0.25*Z_DIM),int(0.5*Z_DIM),int(0.75*Z_DIM), Z_DIM]
s_state_ts = numpy.zeros((len(z_index),len(TIMES_OUT)), dtype = 'f')
s_statemean_ts = numpy.zeros((len(SENSITIVITY),len(TIMES_OUT)), dtype = 'f')
s_state_ts_2 = numpy.zeros((len(z_index),len(TIMES_OUT)), dtype = 'f')
s_statemean_ts_2 = numpy.zeros((len(SENSITIVITY),len(TIMES_OUT)), dtype = 'f')
last_outk = numpy.zeros((len(SENSITIVITY)),dtype='i')
last_outk[:] = last_out
end = len(l_state)

if TIME_SERIES_RUN:
  a_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  s_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  s_state_store_2 = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  l_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  r_rate_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  r_rate_store_2 = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  NOMINAL_INDEX = 0
  if len(SENSITIVITY) != 1:
    print '#### For a time series run there can only be one sensitivity value',repr(SENSITIVITY)
    import sys
    sys.exit()
  zaxis = numpy.zeros((1,Z_DIM+2),dtype = 'f')
else:
  a_state_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  s_state_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  s_state_store_2 = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  l_state_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  r_rate_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  r_rate_store_2 = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  zaxis = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')

a_state_store[0,:] = a_state[:]
s_state_store[0,:] = s_state[:]
s_state_store_2[0,:] = s_state_2[:]
l_state_store[0,:] = l_state[:]
r_rate_store[0,:] = 0.0001
r_rate_store_2[0,:] = 0.0001
for i in range(len(z_index)):
  s_state_ts[i,0] = S_BC
  s_statemean_ts[:,0] = S_BC
  s_state_ts_2[i,0] = S_BC
  s_statemean_ts_2[:,0] = S_BC

time_done = T_DIM_S*NUM_TIMES

if 'parallel' in NEWDIR:
  for kkk in range(len(SENSITIVITY)):
    for j in range(1,Z_DIM+1):
      zaxis[kkk,j] = -(j-1)*DZ
    zaxis[kkk,0] = 0
    zaxis[kkk,Z_DIM+1] = zaxis[kkk,Z_DIM]
  num_cpus = multiprocessing.cpu_count()
  threads =max(1,int(num_cpus*0.75))
  print 'Number of CPUs available: ',num_cpus
  print 'Number of CPUs used     : ',threads
  pool = multiprocessing.Pool(processes=threads)
  fresult = pool.map(loop1,range(len(SENSITIVITY)))
  for kkk in range(len(SENSITIVITY)):
    nominal = (kkk == NOMINAL_INDEX)
    if TIME_SERIES_RUN:
      a_state_store = fresult[kkk][0]
      s_state_store = fresult[kkk][1]
      s_state_store_2 = fresult[kkk][8]
      l_state_store = fresult[kkk][2]
      r_rate_store = fresult[kkk][3]
      r_rate_store_2 = fresult[kkk][11]
    else:
      a_state_store[kkk,:] = fresult[kkk][0]
      s_state_store[kkk,:] = fresult[kkk][1]
      s_state_store_2[kkk,:] = fresult[kkk][8]
      l_state_store[kkk,:] = fresult[kkk][2]
      r_rate_store[kkk,:] = fresult[kkk][3]
      r_rate_store_2[kkk,:] = fresult[kkk][11]
    if nominal:
      last_out = fresult[kkk][4]
      s_state_ts[:,1:] = fresult[kkk][6][:,1:]
      s_state_ts_2[:,1:] = fresult[kkk][9][:,1:]
    last_outk[kkk]=fresult[kkk][5]
    s_statemean_ts[kkk,1:] = fresult[kkk][7][1:]
    s_statemean_ts_2[kkk,1:] = fresult[kkk][10][1:]
else:
  for kkk in range(len(SENSITIVITY)):
    for j in range(1,Z_DIM+1):
      zaxis[kkk,j] = -(j-1)*DZ
    zaxis[kkk,0] = 0
    zaxis[kkk,Z_DIM+1] = zaxis[kkk,Z_DIM]
    loop(kkk)

#  for kkk in range(len(SENSITIVITY)):
#    for j in range(1,Z_DIM+1):
#      zaxis[kkk,j] = -(j-1)*DZ
#    zaxis[kkk,0] = 0
#    zaxis[kkk,Z_DIM+1] = zaxis[kkk,Z_DIM]
#    time_start = time.time()
#    print 'SIM=',kkk,'\nstart:',time.ctime()
#    nominal = (kkk == NOMINAL_INDEX)
#    fresult = simulation(SENSITIVITY[kkk],kkk,nominal)
#    if TIME_SERIES_RUN:
#      a_state_store = fresult[0]
#      s_state_store = fresult[1]
#      l_state_store = fresult[2]
#      r_rate_store = fresult[3]
#    else:
#      a_state_store[kkk,:] = fresult[0]
#      s_state_store[kkk,:] = fresult[1]
#      l_state_store[kkk,:] = fresult[2]
#      r_rate_store[kkk,:] = fresult[3]
#    if nominal:
#      last_out = fresult[4]
#      s_state_ts = fresult[6]
#    last_outk[kkk]=fresult[5]
#    s_statemean_ts[kkk,:] = fresult[7]
#    time_end = time.time()
#    print 'end: CPU=',time_end-time_start
try:
  if 'numerics' in NEWDIR:
    ff = open('cpu.txt','a')
    ff.write(str(time_end-time_start)+'\n')
    ff.close()
    print 'A=',a_state
    raise SystemExit

  pstring = (NEWDIR+' Sensitivity: V='+str(V_NOMINAL)+' D='+str(DIFFUSIVITY_NOMINAL)+
              ' K1='+str(K1_NOMINAL)+' K2='+str(K2_NOMINAL)+
              ' E='+str(ERODE_NOMINAL)+' K3='+str(K3_NOMINAL)+
              ' LEQUILIB='+str(LEQUILIB_NOMINAL)+' BIOTURB='+str(BIOTURB_NOMINAL))
  if TIME_SERIES_RUN:
    pstring=pstring+'\n Max Time = '+str(time_done)+'@'+str(T_DIM_S)
  if 'print' in NEWDIR:
#  size to make the print job fit Cambridge Books 125mm wide 300dpi
    fig1 = matplotlib.pyplot.figure(figsize=(15,10))
    resolution = 300
    matplotlib.rcParams.update({'font.family': 'sans-serif'})
    matplotlib.rcParams.update({'font.size': 18})
    matplotlib.rcParams.update({'axes.labelsize': 'large'})
    matplotlib.rcParams.update({'legend.fontsize': 'medium'})
    matplotlib.rcParams.update({'legend.labelspacing': 0.2})
  else:
#  screen size
    fig1 = matplotlib.pyplot.figure(figsize=(15,10))
    resolution = 300
  fig1.suptitle(pstring)

  contents1=fig1.add_subplot(231)
  if TIME_SERIES_RUN:
    stuff1 = contents1.plot(a_state_store[0,:], zaxis[0,:],
                              label = str(SENSITIVITY[0]))
    for i in range(1,len(TIMES_OUT)-1):
      stuff1 = contents1.plot(a_state_store[i,:], zaxis[0,:])
  else:
    for i in range(len(SENSITIVITY)):
      stuff1 = contents1.plot(a_state_store[i,:-1], zaxis[i,:-1],
                              label = str(SENSITIVITY[i]),
                              ls = LINE_STYLE[i],
                              lw = LINE_WIDTH[i])
  title1 = '(a) A: Acid Conc.'
  contents1.set_title(title1)
  if not 'clean' in NEWDIR:
    if len(SENSITIVITY) > 1:
      contents1.legend(loc='lower right')
#      contents1.legend(loc='lower center')

  contents2=fig1.add_subplot(232)
  if TIME_SERIES_RUN:
    stuff2 = contents2.plot(s_state_store[0,:], zaxis[0,:],
                              label = str(SENSITIVITY[0]))
    for i in range(1,len(TIMES_OUT)):
      stuff2 = contents2.plot(s_state_store[i,:], zaxis[0,:])
  else:
#    for i in range(len(SENSITIVITY)):
#      stuff2 = contents2.plot(s_state_store[i,:-1], zaxis[i,:-1],
#                            label = str(SENSITIVITY[i]),
#                              ls = LINE_STYLE[i],
#                              lw = LINE_WIDTH[i])
      stuff2 = contents2.plot(s_state_store[0,:-1], zaxis[0,:-1],
                            label = '0.5',
                              ls = LINE_STYLE[0],
                              lw = LINE_WIDTH[0])
      stuff2 = contents2.plot(s_state_store_2[0,:-1], zaxis[0,:-1],
                            label = '2',
                              ls = LINE_STYLE[1],
                              lw = LINE_WIDTH[1])
  title2 = '(b) S: Substrate'
  contents2.set_title(title2)
#  if not 'clean' in NEWDIR:
#    if len(SENSITIVITY) > 1:
  contents2.legend(loc='lower left')
#      contents2.legend(loc='lower center')

  contents5=fig1.add_subplot(233)
  if TIME_SERIES_RUN:
    stuff5 = contents5.plot(l_state_store[0,:], zaxis[0,:],
                              label = str(SENSITIVITY[0]))
    for i in range(1,len(TIMES_OUT)):
      stuff5 = contents5.plot(l_state_store[i,:], zaxis[0,:])
  else:
    for i in range(len(SENSITIVITY)):
      stuff5 = contents5.plot(l_state_store[i,:-1], zaxis[i,:-1],
                            label = str(SENSITIVITY[i]),
                              ls = LINE_STYLE[i],
                              lw = LINE_WIDTH[i])
  title5 = '(c) L: Leachate Conc.'
  contents5.set_title(title5)
#  if not 'clean' in NEWDIR:
#    if len(SENSITIVITY) > 1:
#  #    contents5.legend(loc='upper right')
#  #    contents5.legend(loc='lower left')
#      contents5.legend(loc='lower center')

  contents3=fig1.add_subplot(234)
  matplotlib.pyplot.locator_params(axis = 'x', nbins=4)
  if TIME_SERIES_RUN:
    stuff3 = contents3.plot(r_rate_store[1,:], zaxis[0,:],
                              label = str(SENSITIVITY[0]))
    for i in range(2,len(TIMES_OUT)):
      stuff3 = contents3.plot(r_rate_store[i,:], zaxis[0,:])
  else:
#    for i in range(len(SENSITIVITY)):
    stuff3 = contents3.plot(r_rate_store[0,:-1], zaxis[0,:-1],
                              ls = LINE_STYLE[0],
                              lw = LINE_WIDTH[0])
    stuff3 = contents3.plot(r_rate_store_2[0,:-1], zaxis[0,:-1],
                              ls = LINE_STYLE[1],
                              lw = LINE_WIDTH[1])
  title3 = '(d) W: Weathering'
  contents3.set_title(title3)
#  if not 'clean' in NEWDIR:
#    if len(SENSITIVITY) > 1:
#  #    contents3.legend(loc='upper left')
#  #    contents3.legend(loc='upper right')
#  #    contents3.legend(loc='lower right')
#      contents3.legend(loc='center right')

  contents4=fig1.add_subplot(235)
  if TIME_SERIES_RUN:
    stuff4 = contents4.semilogx(r_rate_store[1,:], zaxis[0,:],
                              label = str(SENSITIVITY[0]))
    for i in range(2,len(TIMES_OUT)):
      stuff4 = contents4.semilogx(r_rate_store[i,:], zaxis[0,:])
  else:
#    for i in range(len(SENSITIVITY)):
    stuff4 = contents4.semilogx(r_rate_store[0,:-1], zaxis[0,:-1],
                              ls = LINE_STYLE[0],
                              lw = LINE_WIDTH[0])
    stuff4 = contents4.semilogx(r_rate_store_2[0,:-1], zaxis[0,:-1],
                              ls = LINE_STYLE[0],
                              lw = LINE_WIDTH[0])
  title4 = '(e) W: log(Weathering)'
  contents4.set_title(title4)
#  if not 'clean' in NEWDIR:
#    if len(SENSITIVITY) > 1:
#  #    contents4.legend(loc='lower right')
#  #    contents4.legend(loc='upper left')
#  #    contents4.legend(loc='upper right')
#      contents4.legend(loc='center right')

#  print '\a\a\a\a'
  #matplotlib.pyplot.savefig(LOCAL+'figure.png')
#  matplotlib.pyplot.savefig(LOCAL+'figure.tiff', dpi=resolution)
  matplotlib.pyplot.savefig(LOCAL+'figure.tiff', dpi=resolution)
  matplotlib.pyplot.show()
except:
  if not 'numerics' in NEWDIR:
    import traceback
    print traceback.format_exc()

