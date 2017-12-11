
#  This is the transport limited weathering simulator with the added weathering product state
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
import sys

print ' Chemical Weathering Simulator V4: single mineral version'
print ' --------------------------------------------------------'
print ' Copyright (C) 2017 Garry Willgoose, The University of Newcastle'
print ''
print ' Input Options are (all of these can be used):\n'
print ' parallel: run the code in parallel across 75% of the available cores'
print ' print: output figures in a form suitable for Cambridge University Press'
print ' clean: minimal annotation for Cambridge University Press'
print ' depletion: parameter set that results in full depletion of the acid and substrate'
print ' fixedco2: acid concentration is the same down the profile'
print ' numerics: a simple run to test the numerics and do CPU timing'
print ' timeseries: a single run that outputs a time series of output'
print ' treethrow:  reduces the porosity by 0.01 half way down the profile'
print ' leachate: a run where the concentration of the leachate is limited to 0.4'
print ''
print ' Sensitivity Study options are (only one option can be used):\n'
print ' v: velocity'
print ' k1: parameter k1'
print ' k2: parameter k2'
print ' k3: parameter k3'
print ' erode: erosion rate'
print ' bioturb: bioturbation rate'
print ' depth: depth of the profile\n'

global zaxis,a_state_store,s_state_store,l_state_store,r_rate_store
global last_out,last_outk,s_state_ts,s_statemean_ts

#========================================
#  GLOBAL PARAMETERS AND CONSTANTS
#========================================

K1_NORMAL = 0.1
K1_DEPLETION = 0.5
K2 = 2.0
K3 = 1.0
BIOTURB = 0.005
ERODE = 0.1
DIFFUSIVITY = 0.000001
V = 0.1
LEQUILIB_NORMAL = 1000
LEQUILIB_LTD = 0.4

A_INITIALCONDITION = 1.0                      # this is also the upper BC
S_INITIALCONDITION_NORMAL = 1.0               # this is also the lower BC
S_INITIALCONDITION_DEPLETION = 0.5            # this is also the lower BC
L_INITIALCONDITION = 0.0

PARALLEL = True         # uses 75% of the available cores, if this causes problems set to False
Z_DIM = 50              # number of nodes in the soil profile
TREETHROW_POROSITY_FACTOR = 10
SOIL_DEPTH = 2.0

# The multipliers for the sensitivities are gathered here to make it easy to modify them
# The nominal index is which index of the sensitivities the substrate time series will be calculated for
SENSITIVITY_DETAIL = [0.8, 0.9, 1.0, 1.1, 1.2]
NOMINAL_INDEX_DETAIL = 2

SENSITIVITY_NORMAL = [0.5, 0.666, 1.0, 1.5, 2.0]
NOMINAL_INDEX_NORMAL = 2

SENSITIVITY_DISPERSIVITY = [1.0, 1000, 10000, 100000, 1000000]
NOMINAL_INDEX_DISPERSIVITY = 2

SENSITIVITY_EETA = [0,0.5,1.0, 2]
NOMINAL_INDEX_EETA = 2

SENSITIVITY_K2_DEPLETION = [1.0, 1.5, 2, 3, 4]
NOMINAL_INDEX_K2_DEPLETION = 2

SENSITIVITY_BIOTURBATION = [0.0, 1.0, 5, 10, 50, 100]
NOMINAL_INDEX_BIOTURBATION = 2

#========================================
book = {  'fig8.2': 'timeseries',
          'fig8.3': 'timeseries depletion',
          'fig8.4': 'k1',
          'fig8.5': 'k1 depletion',
          'fig8.6': 'k2',
          'fig8.7': 'k3',
          'fig8.8': 'v',
          'fig8.9': 'v const mass',
          'fig8.10': 'erode',
          'fig8.11': 'dispersivity',
          'fig8.12': 'k1 leachate',
          'fig8.13': 'v leachate',
          'fig8.14': 'v depletion',
          'fig8.15': 'v depletion detail',
          'fig8.16': 'erode depletion',
          'fig8.17': 'erode depletion detail',
          'fig8.19': 'depth',
          'fig8.20': 'depth depletion',
          'fig8.21': 'eeta',
          'fig8.22': 'erode fixedco2',
          'fig8.23': 'erode fixedco2 leachate',
          'fig8.24': 'bioturb depletion',
          'fig8.25': 'bioturb',
          'fig11.1': 'treethrow depletion',
        }
# for figure 8.26 see the two minerals script
#========================================
#  SET RUN HERE
#========================================

if len(sys.argv) >= 2:
#  setting parameters by book figure number

  if not sys.argv[1].lower() in book:
    print 'The command line argument must be a figure number from the book'
    print sys.argv[1].lower(),' not in ',book.keys()
    quit()
  else:
    NEWDIR = book[sys.argv[1].lower()]
    if PARALLEL:
      NEWDIR = NEWDIR+' parallel'
else:
# setting parameters by explicit command arguments

  #NEWDIR = 'timeseries print'.lower()
  #NEWDIR = 'timeseries depletion print'.lower()
  #NEWDIR = 'bioturb 0.005 depletion test-parallel1 print'.lower()
  #NEWDIR = 'erode parallel print'.lower()
  #NEWDIR = 'erode depletion parallel print'.lower()
  #NEWDIR = 'erode depletion detail parallel print'.lower()
  NEWDIR = 'k1 parallel print'.lower()
  #NEWDIR = 'k1 depletion parallel print'.lower()
  #NEWDIR = 'k1 leachate parallel print'.lower()
  #NEWDIR = 'k2 parallel print'.lower()
  #NEWDIR = 'k3 parallel print'.lower()
  #NEWDIR = 'v parallel print'.lower()
  #NEWDIR = 'v const mass parallel print'.lower()
  #NEWDIR = 'v leachate parallel print'.lower()
  #NEWDIR = 'v depletion parallel print'.lower()
  #NEWDIR = 'v depletion detail parallel print'.lower()
  #NEWDIR = 'dispersivity parallel print'.lower()
  #NEWDIR = 'depth parallel print'.lower()
  #NEWDIR = 'depth depletion parallel print'.lower()
  #NEWDIR = 'eeta parallel print'.lower()
  #NEWDIR = 'erode fixedco2 parallel print'.lower()
  #NEWDIR = 'erode leachate fixedco2 parallel print'.lower()
  #NEWDIR = 'bioturb parallel print'.lower()
  #NEWDIR = 'bioturb depletion parallel print'.lower()
  #NEWDIR = 'treethrow depletion parallel print'.lower()

#========================================
#   DONE GLOBALS
#========================================



def simulation(factor,nominal_index, nominal):
  global S_BC
  a_state = numpy.zeros((Z_DIM+2),dtype = 'f')
  a_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  l_state = numpy.zeros((Z_DIM+2),dtype = 'f')
  l_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  r_rate = numpy.zeros((Z_DIM+2),dtype = 'f')
  temp_a = numpy.zeros((Z_DIM+2),dtype = 'f')
  s_state_ts = numpy.zeros((len(z_index),len(TIMES_OUT)), dtype = 'f')
  s_statemean_ts = numpy.zeros((len(TIMES_OUT)), dtype = 'f')
  vel_state = numpy.zeros((Z_DIM+2),dtype = 'f')
  diff_state = numpy.zeros((Z_DIM+2),dtype = 'f')

  last_out = TIMES_OUT[0]
  last_outk = last_out
  midz = int(round(Z_DIM*0.5))

  if TIME_SERIES_RUN:
    a_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    s_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    l_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
    r_rate_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  else:
    a_state_store = numpy.zeros((Z_DIM+2),dtype = 'f')
    s_state_store = numpy.zeros((Z_DIM+2),dtype = 'f')
    l_state_store = numpy.zeros((Z_DIM+2),dtype = 'f')
    r_rate_store = numpy.zeros((Z_DIM+2),dtype = 'f')

  V = V_NOMINAL
  DIFFUSIVITY = DIFFUSIVITY_NOMINAL
  K1 = K1_NOMINAL
  K2 = K2_NOMINAL
  K3 = K3_NOMINAL
  ERODE = ERODE_NOMINAL
  LEQUILIB = LEQUILIB_NOMINAL
  DZ = DZ_NOMINAL
  A_BC1 = A_BC_NOMINAL
  if (NEWDIR == 'v' or 'v leachate' in NEWDIR or 'v const mass' in NEWDIR or
        'v depletion' in NEWDIR or 'v ' in NEWDIR):
    V = V_NOMINAL*factor
    if 'v const mass' in NEWDIR:
      A_BC1 = A_BC_NOMINAL/factor
  elif 'dispersivity' in NEWDIR:
    DIFFUSIVITY = DIFFUSIVITY_NOMINAL*factor
  elif (NEWDIR == 'k1' or 'k1 leachate' in NEWDIR or 'k1 ' in NEWDIR or
        'k1 fixedco2' in NEWDIR or 'k1 leachate fixedco2' in NEWDIR):
    K1 = K1_NOMINAL*factor
  elif 'k2' in NEWDIR or 'k2 depletion' in NEWDIR:
    K2 = K2_NOMINAL*factor
  elif 'k3' in NEWDIR:
    K3 = K3_NOMINAL*factor
  elif (NEWDIR == 'erode' or 'erode ' in NEWDIR or
        'erode fixedco2' in NEWDIR or 'erode leachate fixedco2' in NEWDIR or
        'erode depletion' in NEWDIR):
    ERODE = ERODE_NOMINAL*factor
  elif NEWDIR == 'lequilib':
    LEQUILIB = LEQUILIB_NOMINAL*factor
  elif 'depth' in NEWDIR:
    DZ = DZ_NOMINAL*factor
  elif 'eeta ' in NEWDIR:
    pass
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
  elif 'treethrow' in NEWDIR:
    ERODE = ERODE_NOMINAL*factor
  else:
    print 'ERROR: wrong sensitivity mode selected'
    import sys
    sys.exit()
  if bioturb_active:
    bioturb_1 = BIOTURB/DZ**2
  if 'fixedco2' in NEWDIR:
    fixedco2 = True
  else:
    fixedco2 = False
  if 'treethrow' in NEWDIR:
    vel_state[:midz]= V/DZ
    vel_state[midz:]= TREETHROW_POROSITY_FACTOR*V/DZ
    diff_state[:midz] = DIFFUSIVITY/DZ**2
    diff_state[midz:] = TREETHROW_POROSITY_FACTOR*DIFFUSIVITY/DZ**2
  else:
    vel_state[:] = V/DZ
    diff_state[:] = DIFFUSIVITY/DZ**2
  if 'timeseries' in NEWDIR and 'depletion' in NEWDIR:
    S_BC = S_BC*0.5
  ERODE_1 = ERODE/DZ

  a_state[:] = 0
  a_state_new_a[:] = 0
  s_state[:] = 0
  s_state_new_a[:] = 0
  l_state[:] = 0
  l_state_new_a[:] = 0

  a_state[0] = a_state[1] = A_BC1
  a_state_new_a[0] = a_state_new_a[1] = A_BC1
  s_state[:] = S_BC
  s_state_new_a[:] = S_BC
  s_statemean_ts[0] = S_BC
  s_state_ts[0] = S_BC
  l_state[:] = L_BC
  l_state_new_a[:] = L_BC
  if fixedco2:
    a_state[:] = A_BC1
    a_state_new_a[:] = A_BC1

  for t in range(1,T_DIM+1):
    a_state[end-1] = 2*a_state[end-2]-a_state[end-3]
    s_state_new_a[0] = s_state[1]
    if fixedco2:
      a_state_new_a[1:Z_DIM+1] = A_BC1
      temp_a[1:Z_DIM+1] = (K1 * a_state[1:Z_DIM+1] * s_state[1:Z_DIM+1]*
                            (LEQUILIB-l_state[1:Z_DIM+1])/LEQUILIB)
    else:
      temp_a[1:Z_DIM+1] = (K1 * a_state[1:Z_DIM+1] * s_state[1:Z_DIM+1]*
                            (LEQUILIB-l_state[1:Z_DIM+1])/LEQUILIB)
      if etarun:
        for i in range(1,Z_DIM+1):
          temp_a[i] = temp_a[i]*(0.5*Z_DIM/(i+1))**factor
      a_state_new_a[1:Z_DIM+1] = a_state[1:Z_DIM+1] + DT*(
                  -vel_state[0:Z_DIM] * (a_state[1:Z_DIM+1]-a_state[0:Z_DIM])
                  +diff_state[0:Z_DIM] *
                        (a_state[2:Z_DIM+2]-2*a_state[1:Z_DIM+1]+a_state[0:Z_DIM])
                  -K2 * temp_a[1:Z_DIM+1])
    if bioturb_active:
      s_state_new_a[1:Z_DIM+1] = s_state[1:Z_DIM+1] + DT*(
                ERODE_1 * (s_state[2:Z_DIM+2]-s_state[1:Z_DIM+1])
                +bioturb_1 *
                      (s_state[2:Z_DIM+2]-2*s_state[1:Z_DIM+1]+s_state[0:Z_DIM])
                -temp_a[1:Z_DIM+1])
    else:
      s_state_new_a[1:Z_DIM+1] = s_state[1:Z_DIM+1] + DT*(
                ERODE_1 * (s_state[2:Z_DIM+2]-s_state[1:Z_DIM+1])
                -temp_a[1:Z_DIM+1])

    s_state_new_a.clip(min=0.000001)
    l_state_new_a[1:Z_DIM] = l_state[1:Z_DIM] + DT*(
                -vel_state[1:Z_DIM] * (l_state[1:Z_DIM]-l_state[0:Z_DIM-1])
                +diff_state[1:Z_DIM] *
                      (l_state[2:Z_DIM+1]-2*l_state[1:Z_DIM]+l_state[0:Z_DIM-1])
                +K3 * temp_a[1:Z_DIM])
    l_state_new_a[Z_DIM] = l_state[Z_DIM] + DT*(
                -vel_state[Z_DIM] * (l_state[Z_DIM]-l_state[Z_DIM-1])
                +diff_state[Z_DIM] * (-l_state[Z_DIM]+l_state[Z_DIM-1])
                +K3 * temp_a[Z_DIM])
    r_rate = temp_a.clip(min=1.001e-3)
    invlgth = 1.0/(len(a_state[1:Z_DIM+1]))
    a_diff = a_state[1:Z_DIM+1]-a_state_new_a[1:Z_DIM+1]
    a_sls = numpy.dot(a_diff,a_diff)*invlgth
    s_diff = s_state[1:Z_DIM+1]-s_state_new_a[1:Z_DIM+1]
    s_sls = numpy.dot(s_diff,s_diff)*invlgth
    l_diff = l_state[1:Z_DIM+1]-l_state_new_a[1:Z_DIM+1]
    l_sls = numpy.dot(l_diff,l_diff)*invlgth
    a_state = a_state_new_a.copy()
    s_state = s_state_new_a.copy()
    l_state = l_state_new_a.copy()
#    if a_sls < 1.e-14 and s_sls < 1.e-14 and l_sls < 1.e-14:
    if (a_sls < 1.e-14 and s_sls < 1.e-14 and l_sls < 1.e-14) or t == T_DIM:
      if not TIME_SERIES_RUN:
        a_state_store[:] = a_state[:]
        s_state_store[:] = s_state[:]
        l_state_store[:] = l_state[:]
        r_rate_store[:] = r_rate[:]
      time_done=t
      print 'time finished:',t,a_sls,s_sls,l_sls
      break
    if t in TIMES_OUT:
      j = TIMES_OUT.index(t)
      if TIME_SERIES_RUN:
        time_done=t
        last_out = j
        last_outk = j
        a_state_store[j,:] = a_state[:]
        s_state_store[j,:] = s_state[:]
        l_state_store[j,:] = l_state[:]
        r_rate_store[j,:] = r_rate[:]
      s_statemean_ts[j] = numpy.sum(s_state[:])/len(s_state[:])
      last_outk = j
      if nominal:
        time_done=t
        last_out = j
        for i in range(len(z_index)):
          s_state_ts[i,j] = min(s_state[z_index[i]],S_BC)
  return(a_state_store, s_state_store, l_state_store, r_rate_store,
          last_out, last_outk, s_state_ts, s_statemean_ts)


def loop(kkk):
  global zaxis,a_state_store,s_state_store,l_state_store,r_rate_store
  global last_out,last_outk,s_state_ts,s_statemean_ts
  time_start = time.time()
  print 'SIM=',kkk,'\nstart:',time.ctime()
  nominal = (kkk == NOMINAL_INDEX)
  fresult = simulation(SENSITIVITY[kkk],kkk,nominal)
  if TIME_SERIES_RUN:
    a_state_store = fresult[0]
    s_state_store = fresult[1]
    l_state_store = fresult[2]
    r_rate_store = fresult[3]
  else:
    a_state_store[kkk,:] = fresult[0]
    s_state_store[kkk,:] = fresult[1]
    l_state_store[kkk,:] = fresult[2]
    r_rate_store[kkk,:] = fresult[3]
  if nominal:
    last_out = fresult[4]
    s_state_ts = fresult[6]
  last_outk[kkk]=fresult[5]
  s_statemean_ts[kkk,:] = fresult[7]
  time_end = time.time()
  print 'end: CPU=',time_end-time_start
  return(kkk)


def loop1(kkk):
  global zaxis,a_state_store,s_state_store,l_state_store,r_rate_store
  global last_out,last_outk,s_state_ts,s_statemean_ts
  time_start = time.time()
  print 'SIM=',kkk,'\nstart:',time.ctime()
  nominal = (kkk == NOMINAL_INDEX)
  fresult = simulation(SENSITIVITY[kkk],kkk,nominal)
  time_end = time.time()
  print 'end: CPU=',time_end-time_start
  return(fresult)


# ============================================
# ============================================
#    SETTING RUN PARAMETERS
# ============================================
# ============================================


if 'timeseries' in NEWDIR:
  print 'TIMESERIES RUN: ',NEWDIR
elif 'numerics' in NEWDIR:
  print 'NUMERICS RUN: ', NEWDIR
else:
  print 'SENSITIVITY RUN: '+NEWDIR
print 40*'-'

# the code automatically saves of copy of the graphics output to a directory called 'figures'
VERSIONDIR = 'figures'
LOCAL = VERSIONDIR+os.sep+NEWDIR+os.sep
if os.path.exists(LOCAL):
  pass
else:
  os.mkdir(LOCAL[:-1])

#  run times
T_DIM_S = 4000
NUM_TIMES = 200
if 'numerics' in NEWDIR:
  T_DIM_S = 1000
  NUM_TIMES = 1
elif 'lumped_test' in NEWDIR:
  T_DIM_S = 1000
  NUM_TIMES = 100

# ------------------------
# setting parameter values
# ------------------------

T_DIM = T_DIM_S*NUM_TIMES
TIMES_OUT = []
for i in range(NUM_TIMES+1):
  TIMES_OUT.append(i*T_DIM_S)
TIMES_OUT[0] = 1
last_out = TIMES_OUT[0]

if 'timeseries' in NEWDIR or 'numerics' in NEWDIR:
  TIME_SERIES_RUN = True
else:
  TIME_SERIES_RUN = False

if 'eeta' in NEWDIR:
  etarun=True
else:
  etarun=False

V_NOMINAL = V
DIFFUSIVITY_NOMINAL = DIFFUSIVITY
if 'depletion' in NEWDIR:
  K1 = K1_DEPLETION
else:
  K1 = K1_NORMAL
K1_NOMINAL = K1
K2_NOMINAL = K2
K3_NOMINAL = K3
if 'deposition' in NEWDIR:
  ERODE = -ERODE
ERODE_NOMINAL = ERODE
BIOTURB_NOMINAL = BIOTURB
bioturb_active = False
if 'bioturb' in NEWDIR:
  bioturb_active = True


LEQUILIB = LEQUILIB_NORMAL
if 'leachate' in NEWDIR:
  LEQUILIB = LEQUILIB_LTD
LEQUILIB_NOMINAL = LEQUILIB

A_BC = A_INITIALCONDITION
A_BC_NOMINAL = A_BC
if 'depletion' in NEWDIR:
  S_BC = S_INITIALCONDITION_DEPLETION
else:
  S_BC = S_INITIALCONDITION_NORMAL
L_BC = L_INITIALCONDITION

if 'dispersivity' in NEWDIR:
  DT = 0.0001                     #  timestep
else:
  DT = 0.001                     #  timestep

DZ = SOIL_DEPTH/(Z_DIM-1)
DZ_NOMINAL = DZ

if 'detail' in NEWDIR:
  SENSITIVITY = SENSITIVITY_DETAIL        #  more detail around the nominal value
  NOMINAL_INDEX = NOMINAL_INDEX_DETAIL
else:
  SENSITIVITY = SENSITIVITY_NORMAL        #  most  runs
  NOMINAL_INDEX = NOMINAL_INDEX_NORMAL      # which index of the sensitivity runs to use for the time series plots

if TIME_SERIES_RUN:
  SENSITIVITY = [1.0]       # for time series only want to do it for nominal parameters
  NOMINAL_INDEX = 0         # which index of the sensitivity runs to use for the time series plots

if 'dispersivity' in NEWDIR:
  SENSITIVITY = SENSITIVITY_DISPERSIVITY        #  dispersivity
  NOMINAL_INDEX = NOMINAL_INDEX_DISPERSIVITY      # which index of the sensitivity runs to use for the time series plots

if 'bioturb' in NEWDIR:
  SENSITIVITY = SENSITIVITY_BIOTURBATION        #  BIOTURBATION dispersivity
  NOMINAL_INDEX = NOMINAL_INDEX_BIOTURBATION      # which index of the sensitivity runs to use for the time series plots

if 'k2 depletion' in NEWDIR:
  SENSITIVITY = SENSITIVITY_K2_DEPLETION        #  k2 depletion test
  NOMINAL_INDEX = NOMINAL_INDEX_K2_DEPLETION         # which index of the sensitivity runs to use for the time series plots

if 'eeta' in NEWDIR:
  SENSITIVITY = SENSITIVITY_EETA         # for area sensitivity runs
  NOMINAL_INDEX = NOMINAL_INDEX_EETA         # which index of the sensitivity runs to use for the time series plots

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
l_state = numpy.zeros((Z_DIM+2),dtype = 'f')
l_state_new_a = numpy.zeros((Z_DIM+2),dtype = 'f')
r_rate = numpy.zeros((Z_DIM+2),dtype = 'f')
temp_a = numpy.zeros((Z_DIM+2),dtype = 'f')

z_index = [1,int(0.25*Z_DIM),int(0.5*Z_DIM),int(0.75*Z_DIM), Z_DIM]
s_state_ts = numpy.zeros((len(z_index),len(TIMES_OUT)), dtype = 'f')
s_statemean_ts = numpy.zeros((len(SENSITIVITY),len(TIMES_OUT)), dtype = 'f')
last_outk = numpy.zeros((len(SENSITIVITY)),dtype='i')
last_outk[:] = last_out
end = len(l_state)

if TIME_SERIES_RUN:
  a_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  s_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  l_state_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  r_rate_store = numpy.zeros((len(TIMES_OUT),Z_DIM+2),dtype = 'f')
  NOMINAL_INDEX = 0
  if len(SENSITIVITY) != 1:
    print '#### For a time series run there can only be one sensitivity value',repr(SENSITIVITY)
    import sys
    sys.exit()
  zaxis = numpy.zeros((1,Z_DIM+2),dtype = 'f')
else:
  a_state_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  s_state_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  l_state_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  r_rate_store = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')
  zaxis = numpy.zeros((len(SENSITIVITY),Z_DIM+2),dtype = 'f')

a_state_store[0,:] = a_state[:]
s_state_store[0,:] = s_state[:]
l_state_store[0,:] = l_state[:]
r_rate_store[0,:] = 0.0001
for i in range(len(z_index)):
  s_state_ts[i,0] = S_BC
s_statemean_ts[:,0] = S_BC

time_done = T_DIM_S*NUM_TIMES

# ============================================
# ============================================
#    CALCULATIONS
# ============================================
# ============================================


if 'parallel' in NEWDIR:
#  doing the calculation in parallel mode

  for kkk in range(len(SENSITIVITY)):
    if 'depth' in NEWDIR:
      DZ11 = DZ*SENSITIVITY[kkk]
    else:
      DZ11 = DZ
    for j in range(1,Z_DIM+1):
      zaxis[kkk,j] = -(j-1)*DZ11
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
      l_state_store = fresult[kkk][2]
      r_rate_store = fresult[kkk][3]
    else:
      a_state_store[kkk,:] = fresult[kkk][0]
      s_state_store[kkk,:] = fresult[kkk][1]
      l_state_store[kkk,:] = fresult[kkk][2]
      r_rate_store[kkk,:] = fresult[kkk][3]
    if nominal:
      last_out = fresult[kkk][4]
      s_state_ts[:,1:] = fresult[kkk][6][:,1:]
    last_outk[kkk]=fresult[kkk][5]
    s_statemean_ts[kkk,1:] = fresult[kkk][7][1:]
else:
#  doing the calculation in serial mode

  for kkk in range(len(SENSITIVITY)):
    if 'depth' in NEWDIR:
      DZ11 = DZ*SENSITIVITY[kkk]
    else:
      DZ11 = DZ
    for j in range(1,Z_DIM+1):
      zaxis[kkk,j] = -(j-1)*DZ11
    zaxis[kkk,0] = 0
    zaxis[kkk,Z_DIM+1] = zaxis[kkk,Z_DIM]
    loop(kkk)

try:
  if 'numerics' in NEWDIR:
    ff = open('cpu.txt','a')
    ff.write(str(time_end-time_start)+'\n')
    ff.close()
    print 'A=',a_state
    raise SystemExit
  
# ============================================
# ============================================
#    GRAPHICS OUTPUT
# ============================================
# ============================================

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
  contents1.set_ylabel('depth (m)')
  matplotlib.pyplot.locator_params(axis = 'x', nbins=4)
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
    for i in range(len(SENSITIVITY)):
      stuff2 = contents2.plot(s_state_store[i,:-1], zaxis[i,:-1],
                            label = str(SENSITIVITY[i]),
                              ls = LINE_STYLE[i],
                              lw = LINE_WIDTH[i])
  title2 = '(b) S: Substrate'
  contents2.set_title(title2)
  contents2.yaxis.set_ticklabels([])
  matplotlib.pyplot.locator_params(axis = 'x', nbins=4)
#  if not 'clean' in NEWDIR:
#    if len(SENSITIVITY) > 1:
#  #    contents2.legend(loc='lower left')
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
  contents5.yaxis.set_ticklabels([])
  matplotlib.pyplot.locator_params(axis = 'x', nbins=4)
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
    for i in range(len(SENSITIVITY)):
      stuff3 = contents3.plot(r_rate_store[i,:-1], zaxis[i,:-1],
                            label = str(SENSITIVITY[i]),
                              ls = LINE_STYLE[i],
                              lw = LINE_WIDTH[i])
  title3 = '(d) W: Weathering'
  contents3.set_title(title3)
  contents3.set_ylabel('depth (m)')
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
    for i in range(len(SENSITIVITY)):
      stuff4 = contents4.semilogx(r_rate_store[i,:-1], zaxis[i,:-1],
                                label = str(SENSITIVITY[i]),
                              ls = LINE_STYLE[i],
                              lw = LINE_WIDTH[i])
  title4 = '(e) W: log(Weathering)'
  contents4.set_title(title4)
  contents4.yaxis.set_ticklabels([])
#  if not 'clean' in NEWDIR:
#    if len(SENSITIVITY) > 1:
#  #    contents4.legend(loc='lower right')
#  #    contents4.legend(loc='upper left')
#  #    contents4.legend(loc='upper right')
#      contents4.legend(loc='center right')

  contents6=fig1.add_subplot(236)
#  contents6.set_ylabel('Substrate')
#  plot6,axis6 = matplotlib.pyplot.subplots()
#  print contents6,fig1
#  major_locator = matplotlib.ticker.MultipleLocator(10000)
  matplotlib.pyplot.locator_params(axis = 'x', nbins=4)
  xaxis = []
  for i in range(s_state_ts.shape[1]):
    xaxis.append(TIMES_OUT[i])
  for i in range(len(z_index)):
    if 'lumped_test' in NEWDIR:
      temp = s_state_ts[i,:last_out].clip(min=1.001e-2)
      stuff6 = contents6.semilogy(xaxis[:last_out], temp,
                            label = 'i='+str(z_index[i]),
                              ls = LINE_STYLE[i],
                              lw = LINE_WIDTH[i])
    else:
      stuff6 = contents6.plot(xaxis[:last_out], s_state_ts[i,:last_out],
                            label = 'i='+str(z_index[i]),
                              ls = LINE_STYLE[i],
                              lw = LINE_WIDTH[i])
  import copy
  if len(SENSITIVITY) > 1:
    k=0
    xxaxis = copy.copy(xaxis)
    if last_out > last_outk[k]:
      xxaxis[last_outk[k]+1]=xxaxis[last_out]
      s_statemean_ts[k,last_outk[k]+1] = s_statemean_ts[k,last_outk[k]]
      last_outk[k]=last_outk[k]+1
    else:
      last_outk[k]=last_out
    if 'lumped_test' in NEWDIR:
      temp = s_statemean_ts[k,:last_outk[k]].clip(min=1.001e-2)
      stuff6 = contents6.semilogy(xxaxis[:last_outk[k]], temp,
                              label = 'ave '+str(SENSITIVITY[k]),
                                ls = LINE_STYLE[0],
                                lw = 2*LINE_WIDTH[0], color='k')
    else:
      stuff6 = contents6.plot(xxaxis[:last_outk[k]], s_statemean_ts[k,:last_outk[k]],
                              label = 'ave '+str(SENSITIVITY[k]),
                                ls = LINE_STYLE[0],
                                lw = 2*LINE_WIDTH[0], color='k')
  k=NOMINAL_INDEX
  xxaxis = copy.copy(xaxis)
  if last_out > last_outk[k]:
    xxaxis[last_outk[k]+1]=xxaxis[last_out]
    s_statemean_ts[k,last_outk[k]+1] = s_statemean_ts[k,last_outk[k]]
    last_outk[k]=last_outk[k]+1
  else:
    last_outk[k]=last_out
  if 'lumped_test' in NEWDIR:
    temp = s_statemean_ts[k,:last_outk[k]].clip(min=1.001e-2)
    stuff6 = contents6.semilogy(xxaxis[:last_outk[k]], temp,
                            label = 'ave '+str(SENSITIVITY[k]),
                              ls = LINE_STYLE[1],
                              lw = 2*LINE_WIDTH[1], color='k')
  else:
    stuff6 = contents6.plot(xxaxis[:last_outk[k]], s_statemean_ts[k,:last_outk[k]],
                            label = 'ave '+str(SENSITIVITY[k]),
                              ls = LINE_STYLE[1],
                              lw = 2*LINE_WIDTH[1], color='k')

  if len(SENSITIVITY) > 1:
    k=len(SENSITIVITY)-1
    xxaxis = copy.copy(xaxis)
    if last_out > last_outk[k]:
      xxaxis[last_outk[k]+1]=xxaxis[last_out]
      s_statemean_ts[k,last_outk[k]+1] = s_statemean_ts[k,last_outk[k]]
      last_outk[k]=last_outk[k]+1
    else:
      last_outk[k]=last_out
    if 'lumped_test' in NEWDIR:
      temp = s_statemean_ts[k,:last_outk[k]].clip(min=1.001e-2)
      stuff6 = contents6.semilogy(xxaxis[:last_outk[k]], temp,
                              label = 'ave '+str(SENSITIVITY[k]),
                                ls = LINE_STYLE[3],
                                lw = 2*LINE_WIDTH[3], color='k')
    else:
      stuff6 = contents6.plot(xxaxis[:last_outk[k]], s_statemean_ts[k,:last_outk[k]],
                              label = 'ave '+str(SENSITIVITY[k]),
                                ls = LINE_STYLE[3],
                                lw = 2*LINE_WIDTH[3], color='k')
#    contents6.axis.set_major_locator(major_locator)
    matplotlib.pyplot.locator_params(axis = 'x', nbins=3)

  title6 = '(f) Substrate Time Series'
  contents6.set_title(title6)
  if not 'clean' in NEWDIR:
#    contents6.legend(loc='lower right')
    contents6.legend(loc='upper right')
#    contents6.legend(loc='lower left')

  print '\a\a\a\a'
  #matplotlib.pyplot.savefig(LOCAL+'figure.png')
#  matplotlib.pyplot.savefig(LOCAL+'figure.tiff', dpi=resolution)
  matplotlib.pyplot.savefig(LOCAL+'figure.tiff', dpi=resolution)
  matplotlib.pyplot.show()
except:
  if not 'numerics' in NEWDIR:
    import traceback
    print traceback.format_exc()

