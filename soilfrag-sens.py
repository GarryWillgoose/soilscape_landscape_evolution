#   The analytic matrix fragmentation model and its statistics change with time
#   copyright 2015 Garry Willgoose, The University of Newcastle, Australia
#                   garry.willgoose@newcastle.edu.au
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

try:
  import numpy
except:
  print 'ERROR: The extension library numpy is not installed '
  quit()
try:
  import matplotlib
  import matplotlib.pyplot
except:
  print 'ERROR: The extension library matplotlib is not installed '
  quit()
import sys
import copy
import os
import os.path

version = '1.04'

#  the rationale for the daughter fragmentation shorthand is each number is
#  % of the mass in the next smaller size fraction. The 1st number after the
#  'p' is the parent particle size fraction and '-' seperates the size
#   fractions as they get progressively smaller

daughter_list = ['p2-stoch','p2-0-100','p2-0-50-50','p2-0-50-0-50',
                'p2-0-0-0-100','p2-0-50-25-12.5-12.5',
                'px-98-1-2*0-1','px-98-1-3*0-1','px-98-1-n*0-1','px-98-1',
                'px-98-1-x','px-98-1-spread1','px-98-1-spread2']
daughter_help = [
        '(1) 2 fragments; 2 fragments uniformly randomly distributed between '+
            'the original volume and zero (random fracture of a cube)',
        '(2) 2 fragments; 2 x 1/2 volume',
        '(3) 3 fragments; one x 1/2 the volume, two x 1/4 the volume',
        '(4) 5 fragments; one x 1/2 the size, four x 1/8 volume', 
        '(5) 8 fragments; 8 x 1/8 the volume (the Finke model)',
        '(6) 5 fragments; 1 x 1/2 volume, 1 x 1/4 volume, 1 x 1/8 volume, '+
            '2 x 1/16 volume',
        '(7) X fragments; particles spall, lots of small particles and large'+
            'particles decline slightly in size (model 1)',
        '(8) X fragments; particles spall, lots of small particles and large'+
            'particles decline slightly in size (model 2)',
        '(9) X fragments; particles spall, all small particles to smallest '+
            'grading and large particles decline slightly in size '+
            '(model 3)',
        '(10) X fragments; particles spall, all small particles are dissolved '+
            'and large particles decline slightly in size '+
            '(model 4)',
        '(11) X fragments; as for model 4 but the finest fraction chemically  '+
            'weathers (model 5)',
        '(12) X fragments; as for model 4 but the finest fraction spread over  '+
            'the fine fractions (model 6)',
        '(13) X fragments; as for model 6 but with chemical dissolution (model 7)'
                  ]
daughter = 'p2-0-100'
lgth = len(daughter_list)
matplotlib.rcParams['figure.max_open_warning'] = 2*lgth+5
#   the diameter that we want the PSD for
#largest_dia = 100.0      # the starting diameter for the Coulthard grading runs
#largest_dia = 19.0      # the starting diameter for many of the Ranger expts
#largest_dia = 3.5       # the average diameter for Ranger
largest_dia = 16.0       # close to Ranger and close to phi grading boundary
dia_obj = 0.1
#num_timesteps = 3000
num_timesteps = 1000
debug = False
cum_obj = {}
calibrate_data = False
calibrate_auto = False

OBJ_TEXT = ['','mean','d50','','','','','','','','',
            'd10','d20','d30','d40','d50','d60','d70','d80','d90','d95','d975']
obj = 2                 # default calibration objective function is to use d50

w_rate = 0.03
cw_rate = 0.01
cw_power = 1.0
rate_range = [0.5,0.6667,1.0,1.5,2.0]
spread = 20
nominal_rate_index = 2
ts_ave = {}
ts_stone = {}
cal_filename = ''
cal_d = {}
cal_mass = {}
cal_time = {}
cal_timesort = []
cal_sort = {}
cal_stats = {'d50':{}, 'ave':{}}
cal_stats_file = 'cal_stats'
cumpsdt_save = {}
time_output = num_timesteps
fig_size = (20,10)
fig_size0 = (15,10)
fig_size_cambridge = (15,10)
fig_cambridge_def = {'font.family': 'sans-serif', 'font.size': 18,
                  'axes.labelsize': 'large','legend.fontsize': 'medium',
                  'legend.labelspacing': 0.2}
fig_resolution_cambridge = 300
cambridge = False
fig_save = False
fig_name = ['figure','.tiff']
fig_number = 1

init_psd_filename = ''
init_psd_d = {}
init_psd_mass = {}
init_psd_time = {}


#  ============================================================================
#  ========================== END GLOBALS =====================================
#  ============================================================================



def cal_read(filelist,time=False):
#  global cal_d, cal_mass, cal_time
  cal_d = {}
  cal_mass = {}
  cal_time = {}
  for filename in filelist:
    if filename.strip() != '':
      cal_filename = filename.replace('\n','')
#      cal_filename = filename[:-1]
      if cal_filename[1] =='-':
        raise ValueError('Invalid filename'+cal_filename)
      input = open(cal_filename,'rU')
      try:
        cal_d[cal_filename] = []
        cal_mass[cal_filename] = []
        cal_data = input.readlines()
        for i in range(len(cal_data)):
          cal_data[i] = cal_data[i].replace('\t',' ')
          cal_data[i] = cal_data[i].replace(',',' ')
          cal_data[i] = cal_data[i].replace('\n','')
#     read time for the file     ... line 4
        try:
          cal_time[cal_filename] = int(cal_data[3].split()[1])
        except:
          if not time:
            cal_time[cal_filename] = 0
          else:
           print('ERROR: No time specified in file: '+str(cal_filename))
           return ()
#  read columns that the data is in   ... line 5
        if cal_data[4].strip() == '':
          col_d = 0
          col_mass = 1
        else:
          stuff = cal_data[4].split()
          col_d = int(stuff[0])-1
          col_mass = int(stuff[1])-1
#  read the data
        for line in cal_data[5:]:        #  skip the first 5 lines of header
          if line.strip() != '':
            try:
              ll = line.strip().split(' ')
              cal_d[cal_filename].append(float(ll[col_d]))
              cal_mass[cal_filename].append(float(ll[col_mass]))
            except:
              pass
  # normalisation to [0,1] no matter the units of input
        for i in range(len(cal_mass[cal_filename])):
          cal_mass[cal_filename][i] = cal_mass[cal_filename][i]/cal_mass[cal_filename][-1]
        input.close()
      except:
        import traceback
        traceback.print_exc()
        input.close()
        print 'ERROR: In read of the calibration data from: ',cal_filename
  return(cal_d, cal_mass, cal_time)

def cambridge_fig_defs():
  matplotlib.rcParams.update(fig_cambridge_def)
  return()

def output_stats(stats):
  stats_list = stats.keys()
  for stat in stats_list:
    outfile = open(cal_stats_file+'_'+str(stat)+'.txt','w')
    outfile.write(' STATS from weather-d\n\n')
    frag_models = stats[stat].keys()
    frag_models.sort()
    times = stats[stat][frag_models[0]].keys()
    times.sort()
    outfile.write(' frag_model ')
    for time in times:
      outfile.write(str(time)+' ')
    outfile.write('\n')
    for model in frag_models:
      outfile.write(' '+model+' ')
      for time in times:
        outfile.write(str(stats[stat][model][time])+' ')
      outfile.write('\n')
    outfile.close()

def init_psd(dia, psd, num_psd, daughter, init_psd_filename, init_psd_d, init_psd_mass):
  """
  initialise the psd
  """
  import numpy
  cum_psd = True
  
  if init_psd_d == {}:
    if daughter == 'p2-0-0-0-100':
  # ugly hack for Finke model to remove oscillations
      psd[num_psd-1] = 0.333
      psd[num_psd-2] = 0.334
      psd[num_psd-3] = 0.333
    else:
  #  all the material starts as the coarsest grading
      psd[num_psd-1] = 1.0
    psd1=psd
  else:
    psd[:]=1.0
    for i in range(num_psd):
      d = dia[i]
      if d < init_psd_d[init_psd_filename][0]:
        if d < 0.0001:
          prop=0
        else:
          interp = d/init_psd_d[init_psd_filename][0]
          prop = interp*(init_psd_mass[init_psd_filename][0])
      elif d >= init_psd_d[init_psd_filename][len(init_psd_d[init_psd_filename])-1]:
        prop=1.0
      else:
        for j in range(len(init_psd_d[init_psd_filename])):
          jj = j
          if d <= init_psd_d[init_psd_filename][j]:
            break
        interp = ((d-init_psd_d[init_psd_filename][jj-1])/
                  (init_psd_d[init_psd_filename][jj]-init_psd_d[init_psd_filename][jj-1]))
        prop = (init_psd_mass[init_psd_filename][jj-1]+
                interp*(init_psd_mass[init_psd_filename][jj]-init_psd_mass[init_psd_filename][jj-1]))
      psd[i] = prop
      psd1=numpy.zeros((psd.shape[0]),'f')
      if cum_psd:
        psd1[0]=psd[0]
        for i in range(1,psd.shape[0]):
          psd1[i]=psd[i]-psd[i-1]
  return(psd1)

def psd_stats(dia,cumsum,num_psd):
  for i in range(num_psd):
    if cumsum[i] > 0.1:
      if i > 0:
        d10 = (dia[i-1]+(0.1-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d10 = max(0.0,(0.1-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.2:
      if i > 0:
        d20 = (dia[i-1]+(0.2-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d20 = max(0.0,(0.2-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.3:
      if i > 0:
        d30 = (dia[i-1]+(0.3-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d30 = max(0.0,(0.3-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.4:
      if i > 0:
        d40 = (dia[i-1]+(0.4-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d40 = max(0.0,(0.4-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.5:
      if i > 0:
        d50 = (dia[i-1]+(0.5-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d50 = max(0.0,(0.5-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.6:
      if i > 0:
        d60 = (dia[i-1]+(0.6-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d60 = max(0.0,(0.6-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.7:
      if i > 0:
        d70 = (dia[i-1]+(0.7-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d70 = max(0.0,(0.7-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.8:
      if i > 0:
        d80 = (dia[i-1]+(0.8-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d80 = max(0.0,(0.8-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.9:
      if i > 0:
        d90 = (dia[i-1]+(0.9-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d90 = max(0.0,(0.9-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.95:
      if i > 0:
        d95 = (dia[i-1]+(0.95-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d95 = max(0.0,(0.95-cumsum[i])/cumsum[i]*dia[i])
      break
  for i in range(num_psd):
    if cumsum[i] > 0.975:
      if i > 0:
        d975 = (dia[i-1]+(0.975-cumsum[i-1])/
                (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
      else:
        d975 = max(0.0,(0.975-cumsum[i])/cumsum[i]*dia[i])
      break
  return(d10,d20,d30,d40,d50,d60,d70,d80,d90,d95,d975)


#  ============== START PROCESSING COMMAND LINE ARGUMENTS =====================

if '-help' in sys.argv:
  print 
  print(60*'-'+'\n'+
        ' Physical fragmentation particle size distribution model '+version+'\n'+
        ' 2015 Garry Willgoose, Uni of Newcastle (Australia)\n'+
        60*'-'+'\n'+
        ' Fragmentation model options:')
  for i in range(len(daughter_list)):
    print daughter_list[i],':: ',daughter_help[i]
  print('\n'+
   'Command line options:\n'+
   ' -help            : Output some brief help on using the model\n'+
   ' -cdaughter <arg> : Input of a custom fragmentation model (not yet implemented)\n'+
   ' -daughter <arg>  : Select one of the built-in fragmentation models (overrides the default='+repr(daughter)+')\n'+
   ' -all_daughters   : Select all fragmentation models (overrides the default='+repr(daughter)+')\n'+
   '                  : (overrides the default='+repr(OBJ_TEXT[obj])+')\n'+
   ' -dia_max <arg>   : Set the maximum diameter used in the model (overrides the default='+str(largest_dia)+')\n'+
   ' -time <arg>      : Set the number of timesteps for the simulation (overrides the default='+str(num_timesteps)+')\n'+
   ' -w_rate <arg>    : Set the nominal fragmentationweathering rate (overrides the default='+str(w_rate)+')\n'+
   ' -cw_rate <arg>   : Set the nominal chemical weathering rate (overrides the default='+str(cw_rate)+')\n'+
   ' -spread <arg>    : Set the spread of sines in spalling (overrides the default='+str(spread)+')\n'+
   ' -init_psd <file name> : The grading distribution to be used as the initial condition (default=100% largest size)\n'+
   ' -cal_dia <arg>   : Set the diameter for which full PSD stats are output (overrides the default='+str(dia_obj)+')\n'+
   ' -cal_obj <arg>   : Set what statistic to calibrate to in -cal_dia. Options \'mean\' and \'d50\'\n'+
   '                    and \'d10\' to \'d90\' in increments of 10\n'+
   ' -cal_psd <filename> : Input the file in <filename> as data to calibrate the model to.\n'+
   '                         If -cal_auto is OFF then just display this data at time given by -cal_dia\n'+
   ' -cal_list <filename> : This is a file that includes a list of file names to compare against,\n'+
   '                         one file name per line. The data files are in the format used by -cal_psd.\n'+
   '                         The fourth line of the data file must be \n'+
   '                         time <time value>\n'+
   '                         NOTE: this option is different to -cal_auto because -cal_psd outputs the grading at a given\n'+
   '                         time, matching the time in the files rather than the best fit, as done by -cal_auto.\n'+
   '                         It is assumed that the data files are listed in order from early time to late time\n'+
   ' -cal_auto        : Automatically calibrate the model to data input '+
          'in -cal_psd and using the objective in -cal_obj.\n'+
   '                    The value of -cal_dia is ignored and the value '+
          'derived from the calibration data is used \n'+
   '                    (the default is cal_auto ON).\n'+
   ' -debug           : Turn on debug mode (default is debugging OFF)\n'+
   ' -fig_leg_fsize <points> : The font size for the legend font '
   ' -fig_leg_lspace <arg> : The line spacing for the legend, line spacing=(arg*fsize)'
   ' -fig_fsize <points> : The font size for the figure'
   ' -fig_save        : Save figures in a high resolution tiff file\n'+
   ' -fig_cambridge   : Set figure defaults for 1 column format for Cambridge Press in -fig output\n'+
   60*'-'+'\n'+
   ' The rationale for the daughter fragmentation shorthand is each number is the % of the mass \n'+
   ' in the next smaller size fraction. The 1st number after the "p" is the parent particle size \n'+
   ' fraction "-" seperates the size fractions as they get progressively smaller\n'+
   '\n'
   ' The output filenames for the diameter objective have two strings embedded in their\n'+
   ' filename. The first is \'_model####\' where #### is the fragmentation model used for that\n'+
   ' file. The second is \'_dia####\' where #### is the objective diameter that the file contains.\n'+
   ' The code overwrites any pre-existing files.\n'+
   '\n'+
   60*'-'+'\n')

else:
  if '-cdaughter' in sys.argv:
    print 'Custom fragmentation model not yet implemented'
    index = sys.argv.index('-cdaughter')
    print 'Model input = ',sys.argv[index+1]
    print 'Using default model = ',daughter
  elif '-daughter' in sys.argv:
    index = sys.argv.index('-daughter')
    daughter1=sys.argv[index+1]
    if daughter1 in daughter_list:
      print 'Valid Model requested = ',daughter1
      daughter=daughter1
    else:
      print 'ERROR: Requested model does not exist = ',daughter1
      print 'Using default model = ',daughter

  if '-cal_dia' in sys.argv:
    index = sys.argv.index('-cal_dia')
    try:
      dia_obj = float(sys.argv[index+1])
    except:
      print 'ERROR: Requested final diameter invalid = ',sys.argv[index:]
      print 'Using default final diameter = ',dia_obj

  if '-cal_obj' in sys.argv:
    index = sys.argv.index('-cal_obj')
    try:
      obj = OBJ_TEXT.index(sys.argv[index+1].lower())
    except:
      print 'ERROR: Invalid option on -cal_obj:',sys.argv[index+1]
      print 'Option must be one of ',OBJ_TEXT

  if '-dia_max' in sys.argv:
    index = sys.argv.index('-dia_max')
    try:
      largest_dia = float(sys.argv[index+1])
    except:
      print 'ERROR: Requested final diameter invalid = ',sys.argv[index:]
      print 'Using default final diameter = ',largest_dia

  if '-time' in sys.argv:
    index = sys.argv.index('-time')
    temp = sys.argv[index+1]
    try:
      num_timesteps=int(round(float(temp)))
      time_output=num_timesteps
    except:
      print 'ERROR: Requested number of timesteps invalid = ',sys.argv[index:]
      print 'Using default number of timesteps = ',num_timesteps

  if '-w_rate' in sys.argv:
    index = sys.argv.index('-w_rate')
    temp = sys.argv[index+1]
    try:
      w_rate=float(temp)
    except:
      print 'ERROR: Requested weathering rate invalid = ',sys.argv[index:]
      print 'Using default weathering rate = ',w_rate

  if '-cw_rate' in sys.argv:
    index = sys.argv.index('-cw_rate')
    temp = sys.argv[index+1]
    try:
      cw_rate=float(temp)
    except:
      print 'ERROR: Requested weathering rate invalid = ',sys.argv[index:]
      print 'Using default weathering rate = ',w_rate

  if '-spread' in sys.argv:
    index = sys.argv.index('-spread')
    temp = sys.argv[index+1]
    try:
      spread=int(temp)
    except:
      print 'ERROR: Requested spread invalid = ',sys.argv[index:]
      print 'Using default spread = ',spread

  if '-debug' in sys.argv:
    debug=True

  if '-fig_cambridge' in sys.argv:
    cambridge = True
    fig_size = fig_size_cambridge
    fig_size0 = fig_size_cambridge

  if '-fig_save' in sys.argv:
    fig_save = True

  if '-fig_leg_lspace' in sys.argv:
    index = sys.argv.index('-fig_leg_lspace')
    temp = sys.argv[index+1]
    matplotlib.rcParams.update({'legend.labelspacing':temp})

  if '-fig_leg_fsize' in sys.argv:
    index = sys.argv.index('-fig_leg_fsize')
    temp = sys.argv[index+1]
    matplotlib.rcParams.update({'legend.fontsize':temp})

  if '-fig_fsize' in sys.argv:
    index = sys.argv.index('-fig_fsize')
    temp = sys.argv[index+1]
    matplotlib.rcParams.update({'font.size':temp})

  if '-init_psd' in sys.argv:
    index = sys.argv.index('-init_psd')
    init_psd_filename = sys.argv[index+1].strip()
    temp1 = cal_read([init_psd_filename,],time=False)
    init_psd_d = temp1[0]
    init_psd_mass = temp1[1]
    init_psd_time = temp1[2]

  if '-cal_psd' in sys.argv:
    index = sys.argv.index('-cal_psd')
    cal_filename = sys.argv[index+1]
    temp = cal_read([cal_filename,])
    cal_d = temp[0]
    cal_mass = temp[1]
    cal_time = temp[2]
    cal_sort[time_output] = cal_filename
    cal_timesort = [cal_time[cal_filename],]
    calibrate_data = True

  if '-cal_list' in sys.argv:
    index = sys.argv.index('-cal_list')
    cal_filename = sys.argv[index+1]
    try:
      if cal_filename[1] =='-':
        raise ValueError('Invalid filename'+cal_filename)
      input = open(cal_filename,'rU')
      cal_filelist = input.readlines()
      temp = cal_read(cal_filelist)
      cal_d = temp[0]
      cal_mass = temp[1]
      cal_time = temp[2]
      input.close()
      cal_sort = {}
      for file in cal_time.keys():
        cal_sort[cal_time[file]] = file
      cal_timesort=cal_sort.keys()
      cal_timesort.sort()
      calibrate_data = True
      calibrate_auto = False
    except:
      import traceback
      traceback.print_exc()
      input.close()
      print 'ERROR: In read of the calibration data from: ',cal_filename
#    print 'after cal_list',cal_d,cal_mass,cal_time,cal_sort,cal_timesort

  if '-cal_auto' in sys.argv:
# this option MUST come after dealing with -cal_psd, -cal_obj, and
# -cal_dia because these options are either needed or may overwrite the
# results of the following code. Incompatible with -cal_list.
#    cal_file = cal_d.keys()[0]
    if len(cal_d) > 1:
      print 'ERROR: Only calibration data file is allowed in -cal_auto'
      print '       Do not use -cal_list in combination with -cal_auto'
    else:
      if  calibrate_data:
        cal_file = cal_d.keys()[0]
        if obj == 1:
  # this code assumes the data is normalised between [0,1] in -cal_psd
          sum1 = 0
          for i in range(1,len(cal_d[cal_file])):
            sum1 = sum1 + cal_d[cal_file][i]*(cal_mass[cal_file][i]-cal_mass[cal_file][i-1])
          dia_obj = sum1
          calibrate_auto = True
        elif obj == 2:
          for i in range(1,len(cal_mass[cal_file])):
            if cal_mass[cal_file][i] >= 0.5:
  #            dia_obj = cal_d[i]
  #            break
              if i > 0:
                dia_obj = (cal_d[cal_file][i-1]+(0.5-cal_mass[cal_file][i-1])/
                            (cal_mass[cal_file][i]-cal_mass[cal_file][i-1])*
                            (cal_d[cal_file][i]-cal_d[cal_file][i-1]))
              else:
                dia_obj = max(0.0,(0.5-cal_mass[cal_file][i])/cal_mass[cal_file][i]*cal_d[cal_file][i])
              break
          calibrate_auto = True
        elif obj >= 11 and obj <=21:
          threshold = 0.1*(obj-10)
          for i in range(1,len(cal_mass[cal_file])):
            if cal_mass[cal_file][i] >= threshold:
              if i > 0:
                dia_obj = (cal_d[cal_file][i-1]+(threshold-cal_mass[cal_file][i-1])/
                            (cal_mass[cal_file][i]-cal_mass[cal_file][i-1])*
                            (cal_d[cal_file][i]-cal_d[cal_file][i-1]))
              else:
                dia_obj = max(0.0,(threshold-cal_mass[cal_file][i])/cal_mass[cal_file][i]*cal_d[cal_file][i])
              break
          calibrate_auto = True
        else:
          print 'ERROR: Invalid objective function in -cal_auto processing'
      else:
  #  this option calibrates to the median/mean diameter in dia_obj
        print ('Calibrating to diameter:'+str(dia_obj))
        calibrate_auto = True


#  ============== END PROCESSING COMMAND LINE ARGUMENTS =====================

#  the rationale for the shorthand is each number is % of the mass in the next
#  smaller size fraction. The 1st number number is the parent particle size
#  fraction and '-' seperates the size fractions as they get progressively
#   smaller

  print 'Physical fragmentation grading model ',version
  print 60*'-'
  print ' 2016 Garry Willgoose, Uni of Newcastle (Australia)'
  print 60*'-'
  print ' Fragmentation model options:'
  print daughter_list
  num_psd = 60
  num_saves = 50
  LINE_STYLE = ['-', '-', '-', '-.', '-']
  #LINE_STYLE = ['--', '-', ':', '-.', '-']
  num_sens = len(rate_range)
  scaling = 1.0/0.7937
  smallest_dia = largest_dia/scaling**(num_psd-1)
  savesf = []
  for i  in range(num_saves+1):
    savesf.append(1.0/num_saves*i)
  saves = []
  for s in savesf:
    saves.append(int(s*num_timesteps))

  cum = numpy.zeros((num_psd),dtype='f')
  psd = numpy.zeros((num_psd),dtype='f')
  psd_save = numpy.zeros((num_psd,len(saves)),dtype='f')
  cumpsd_save = numpy.zeros((num_psd,len(saves)),dtype='f')
  ave_save = numpy.zeros((num_sens,len(saves)),dtype='f')
  sa_save = numpy.zeros((num_sens,len(saves)),dtype='i')
  stone_save = numpy.zeros((num_sens,len(saves)),dtype='f')

  median_save = numpy.zeros((len(saves)),dtype='f')
  d10_save = numpy.zeros((len(saves)),dtype='f')
  d90_save = numpy.zeros((len(saves)),dtype='f')
  diac = numpy.zeros((num_psd),dtype='f')
  W = numpy.zeros((num_psd,num_psd),dtype = 'f')
  CW = numpy.zeros((num_psd,num_psd),dtype = 'f')

# if -all_daughters option is selected then run the code for all gradings, if
  if '-all_daughters' in sys.argv:
    daughter_list1 = daughter_list
  else:
    if isinstance(daughter, (list,tuple)):
      daughter_list1 = daughter
    elif isinstance(daughter, str):
      if daughter.strip() == '':
        daughter_list1 = daughter_list
      else:
        daughter_list1 = [daughter.strip()]
    else:
      print 'ERROR: Specified fragmentation model is not a string or list:',daughter
      print '   --- Continuing using a null daughter string'
      daughter_list1 = daughter_list
  if calibrate_data:
    for t in cal_sort:
      cumpsdt_save[t]={}
  else:
    cumpsdt_save[time_output]={}
    cal_timesort=[time_output]
  print 'Selected Fragmentation Model =',daughter_list1
  print 'Selected Final Diameter =',dia_obj
  print 'Selected Objective Function =',OBJ_TEXT[obj]
  print 'Selected Weathering Rate =',w_rate
  for daughter in daughter_list1:
    psd[:] = 0
    psd_save[:,:] = 0
    cumpsd_save[:,:] = 0
    ave_save[:,:] = 0
    sa_save[:,:] = 0
    stone_save[:,:] = 0
    cal_stats['d50'][daughter] ={}
    cal_stats['ave'][daughter] ={}

    median_save[:] = 0
    d10_save[:] = 0
    d90_save[:] = 0
    diac[:] = 0
    W[:,:] = 0
    CW[:,:] = 0

    file_output = False
    try:
      filename_rate = ('weather_psd_model'+str(daughter).strip()+
                        '_dia-'+str(dia_obj).strip()[:8]+'.txt')
    except:
      filename_rate = ('weather_psd_model'+str(daughter).strip()+
                        '_dia-'+str(dia_obj).strip()+'.txt')
    filename_rate = filename_rate.replace('*','')
    for k in range(num_sens):
      W[:,:] = 0
      CW[:,:] = 0
      diac[:] = 0
      psd[:] = 0
      rate = rate_range[k]*w_rate
      rate_cw = rate_range[k]*cw_rate
      
      for i in range(num_psd):
        diac[i] = i
        W[i,i] = 1
        CW[i,i] = 1

      dia = []
      sa_dia = []
      for c in diac:
        diameter = smallest_dia*scaling**c
        dia.append(diameter)
        sa_dia.append(6.0/diameter)
        if diameter <= 2.0:
          stone_index = c
#      print 'stone index: ',stone_index

      if daughter.lower() == 'p2-0-100':
#  one size of daughter particles
        for i in range(1,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = rate
        W[0,0] = 1.0             #  the smallest fraction captures all the fines
      elif daughter == 'p2-0-50-50':
        for i in range(2,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 0.5*rate
          W[i-2,i] = 0.5*rate
        W[0,1] = rate
        W[1,1] = W[1,1]-rate
        W[0,0] = 1.0             #  the smallest fraction captures all the fines
      elif daughter == 'p2-0-50-0-50':       # spalling
        for i in range(3,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 0.5*rate
          W[i-3,i] = 0.5*rate
        W[0,1] = rate
        W[0,2] = 0.5*rate
        W[1,2] = 0.5*rate
        W[1,1] = W[1,1]-rate
        W[2,2] = W[2,2]-rate
        W[0,0] = 1.0             #  the smallest fraction captures all the fines
      elif daughter == 'p2-0-0-0-100':       # Finke model
        for i in range(3,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-3,i] = rate
        W[0,0] = 1.0             #  the smallest fraction captures all the fines
        W[0,1] = rate
        W[1,1] = W[1,1]-rate
        W[0,2] = rate
        W[1,2] = 0
        W[2,2] = W[2,2]-rate
      elif daughter == 'p2-0-50-25-12.5-12.5':       # scaling model
        for i in range(4,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 0.5*rate
          W[i-2,i] = 0.25*rate
          W[i-3,i] = 0.125*rate
          W[i-4,i] = 0.125*rate
        W[0,0] = 1.0             #  the smallest fraction captures all the fines
        W[0,1] = rate
        W[1,1] = W[1,1]-rate
        W[0,2] = 0.5*rate
        W[1,2] = 0.5*rate
        W[2,2] = W[2,2]-rate
        W[0,3] = 0.25*rate
        W[1,3] = 0.25*rate
        W[2,3] = 0.5*rate
        W[3,3] = W[3,3]-rate
      elif daughter == 'p2-stoch':
#         resolution 5 fractions (same resolution as scaling model) ...
#         can be extended if need be
#
# the proportion of the mass that is transferred out of parent grading
        transfer = 0.69262
        fractions=[-transfer, 0.5*transfer, 0.25*transfer, 0.125*transfer, 0.125*transfer]
        for i in range(len(fractions)-1,num_psd):
          for kk in range(len(fractions)):
            W[i-kk,i] = W[i-kk,i]+fractions[kk]*rate
        for i in range(len(fractions)):
          for kk in range(1,i+1):
            W[kk,i] = W[kk,i]+fractions[i-kk]*rate
          W[0,i] = W[0,i]+sum(fractions[i:])*rate
        W[0,0]=1.0        # captures all the smallest particles
      elif daughter == 'px-98-1-2*0-1':    # spalling model 1
        for i in range(4,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
          W[i-2,i] = 0
          W[i-3,i] = 0
          W[i-4,i] = rate/2.62
        W[0,0] = 1.0             #  the smallest fraction captures all the fines
        W[0,1] = rate
        W[1,1] = W[1,1]-rate
        W[0,2] = rate/2.62
        W[1,2] = 1.62*rate/2.62
        W[2,2] = W[2,2]-rate
        W[0,3] = rate/2.62
        W[1,3] = 0
        W[2,3] = 1.62*rate/2.62
        W[3,3] = W[3,3]-rate
      elif daughter == 'px-98-1-3*0-1':    # spalling model 2
        for i in range(5,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
          W[i-2,i] = 0
          W[i-3,i] = 0
          W[i-4,i] = 0
          W[i-5,i] = rate/2.62
        W[0,0] = 1.0             #  the smallest fraction captures all the fines
        W[0,1] = rate
        W[1,1] = W[1,1]-rate
        W[0,2] = rate/2.62
        W[1,2] = 1.62*rate/2.62
        W[2,2] = W[2,2]-rate
        W[0,3] = rate/2.62
        W[1,3] = 0
        W[2,3] = 1.62*rate/2.62
        W[3,3] = W[3,3]-rate
        W[0,4] = rate/2.62
        W[1,4] = 0
        W[2,4] = 0
        W[3,4] = 1.62*rate/2.62
        W[4,4] = W[4,4]-rate
      elif daughter == 'px-98-1-n*0-1':    # spalling model 3
        for i in range(2,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
        W[0,2:] = 1.0/1.62*rate        #  the smallest fraction captures all the fines
        W[1,1] = W[1,1]-rate
        W[0,1] = rate
      elif daughter == 'px-98-1':    # spalling model 4
        for i in range(2,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
        W[0,2:] = 1.0/1.62*rate       # the smallest fraction captures all the 
                                      # fines and chemically weathers them
        W[1,1] = W[1,1]-rate
        W[0,1] = rate
        W[0,0] = 0
      elif daughter == 'px-98-1-x':    # spalling model 5
        alpha =0.98
        for i in range(2,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
        W[0,2:] = alpha/1.62*rate       # the smallest fraction captures all 
                                        # the fines and chemically weathers them
        W[1,1] = W[1,1]-rate
        W[0,1] = alpha*rate
        W[0,0] = alpha
      elif daughter == 'px-98-1-spread1':    # spalling model 6
                                            # uniformly spread over the next 'spread' smaller fractions
        for i in range(spread+1,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
          for j in range(2,2+spread):
            W[i-j,i] = rate/(2.62*spread)
        for i in range(2,spread+1):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
          for j in range(2,i):
            W[i-j,i] = rate/(2.62*(i-1))
        W[1,1] = W[1,1]-rate
        W[0,1] = rate
      elif daughter == 'px-98-1-spread2':    # spalling model 7 with chemical weathering
                                            # uniformly spread over the next 'spread' smaller fractions
        for i in range(spread+1,num_psd):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
          for j in range(2,2+spread):
            W[i-j,i] = rate/(2.62*spread)
        for i in range(2,spread+1):
          W[i,i] = W[i,i]-rate
          W[i-1,i] = 1.62*rate/2.62
          for j in range(2,i):
            W[i-j,i] = rate/(2.62*(i-1))
        W[1,1] = W[1,1]-rate
        W[0,1] = rate
        d_norm = 0.001
        CW_norm = d_norm**cw_power      # normalised by the rate at 1 micron
        print 'CW_norm',CW_norm
        for kk in range(num_psd):
          if d_norm <= dia[kk]:
            CW[kk,kk] = CW[kk,kk]-rate_cw*CW_norm/dia[kk]**cw_power
          else:
            CW[kk,kk] = CW[kk,kk]-rate_cw
        W = numpy.dot(W,CW)
      else:
        print 'Invalid Daughter model Option'
        break
      if debug:
        print 'W bottom 6'
        print W[:6,:6]
        print 'W bottom 6-12'
        print W[:12,6:12]
        print 'W top'
        print W[-6:,-6:]
        print 'CW bottom 6'
        print CW[:6,:6]
        print 'CW bottom 6-12'
        print CW[:12,6:12]
        print 'CW top'
        print CW[-6:,-6:]

      psd = init_psd(dia, psd, num_psd, daughter, init_psd_filename, init_psd_d, init_psd_mass)

      ave = numpy.sum(numpy.dot(psd,dia))
      if debug:
        print 'average diameter START = ',k,ave
      save_index = 1
      cumsum = numpy.cumsum(psd)
      psd_save[:,0]=psd[:]
      cumpsd_save[:,0]=cumsum[:]

      for t in range(1,num_timesteps+1):
#
#  START calcs
# ==============================================================================
        psdt = numpy.dot(W,psd)
        ave = numpy.sum(numpy.dot(psdt,dia))
#        sa = numpy.sum(numpy.dot(psdt,sa_dia))
        cumsum = numpy.cumsum(psdt)
        cumsum[:] = cumsum[:]/cumsum[cumsum.shape[0]-1]
        stone = 100*(1.0-cumsum[stone_index])
        if stone < 0:
          print '#### stone',stone,stone_index,cumsum[stone_index-5:]
        stats = psd_stats(dia,cumsum,num_psd)
        d50 = stats[4]
        if obj >=11 and obj <= 21:
          ddd= stats[obj-11]
#        for i in range(num_psd):
#          if cumsum[i] > 0.5:
#            if i > 0:
#              d50 = (dia[i-1]+(0.5-cumsum[i-1])/
#                      (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
#            else:
#              d50 = max(0.0,(0.5-cumsum[i])/cumsum[i]*dia[i])
#            break
#
#  END calcs
# ==============================================================================
#
#  START output to file calibration results when its for the nominal value of
#   the weathering rate
# ==============================================================================
        if calibrate_auto:
          if k == nominal_rate_index:
            converged = False
            if obj == 1:
              converged = (dia_obj >= ave and not file_output)
            elif obj == 2:
              converged = (dia_obj >= d50 and not file_output)
            elif obj >= 11 and obj <= 21:
              converged = (dia_obj >= ddd and not file_output)
            if converged:
              if cal_filename.strip() != '':
                cal_dir = os.path.splitext(cal_filename)
              else:
                cal_dir = ('weather_psd','')
              try:
                cum_obj[daughter] = copy.copy(cumsum)
                cal_timesort = [t,]
                cal_sort = {t: cal_filename}
                cumpsdt_save[t] = {}
                cumpsdt_save[t][daughter]=copy.copy(cumsum)
                if not os.path.exists(cal_dir[0]):
                  os.mkdir(cal_dir[0])
                else:
                  if not os.path.isdir(cal_dir[0]):
                    print 'ERROR: a file already exists with a name required: ',cal_dir[0]
                    raise ValueError()
                filename = cal_dir[0]+os.sep+filename_rate
                fout = open(filename,'w')
                fout.write('WEATHERING_PSD '+str(version)[:3]+'\n\n')
                fout.write(' Weathering Rate '+str(w_rate)+'\n')
                fout.write(' Time '+str(t)+'\n')
                fout.write(' diameter cdf-'+daughter.strip()+' psd-'+
                            daughter.strip()+'\n')
                for iii in range(num_psd):
                  fout.write(' '+str(dia[iii])+' '+str(cumsum[iii])+' '+
                                  str(psdt[iii])+'\n')
                fout.close()
                file_output = True
                print '#### Objective grading output at ',t,' in ',filename
              except:
                import traceback
                traceback.print_exc()
                file_output = True
                print 'ERROR: problem outputting the output file: ',cal_dir[0]+os.sep+filename_rate
#
#  END output to file calibration results when its for the nominal value of the
#   weathering rate
# ==============================================================================
#
#  START output to file data matches for specified times if calibration data is
#   ON but auto calibration is OFF.
# ==============================================================================
        elif calibrate_data:
          if k == nominal_rate_index:
            if t in cal_sort:
              try:
                cum_obj[daughter] = copy.copy(cumsum)
                cumpsdt_save[t][daughter]=copy.copy(cumsum)
                cal_stats['d50'][daughter][t] = d50
                cal_stats['ave'][daughter][t] = ave
                if cal_filename.strip() != '':
                  cal_dir = os.path.splitext(cal_filename)
                else:
                  cal_dir = ('weather_psd','')
                if not os.path.exists(cal_dir[0]):
                  os.mkdir(cal_dir[0])
                else:
                  if not os.path.isdir(cal_dir[0]):
                    print 'ERROR: a file already exists with a name required for a directory: ',cal_dir[0]
                    raise ValueError()
                filename_split = os.path.splitext(filename_rate)
                filename = (cal_dir[0]+os.sep+filename_split[0]+'_t'+str(t)+
                            filename_split[1])
                fout = open(filename,'w')
                fout.write('WEATHERING_PSD '+str(version)[:3]+'\n\n')
                fout.write(' Weathering Rate '+str(w_rate)+'\n')
                fout.write(' Time '+str(t)+'\n')
                fout.write(' diameter cdf-'+daughter.strip()+' psd-'+
                              daughter.strip()+'\n')
                for iii in range(num_psd):
                  fout.write(' '+str(dia[iii])+' '+str(cumsum[iii])+' '+
                              str(psdt[iii])+'\n')
                fout.close()
                file_output = True
                print '#### Objective grading output in ',filename
              except:
                import traceback
                traceback.print_exc()
                file_output = True
                print 'ERROR: problem outputting the output file: ',filename
#
#  START output to file data matches for specified times if calibration data is
#   ON but auto calibration is OFF.
# ==============================================================================
        else:
          if k == nominal_rate_index:
            if t+1 == time_output:
              cum_obj[daughter] = copy.copy(cumsum)
              cumpsdt_save[time_output][daughter]=copy.copy(cumsum)
#
#  END output to file data matches for specified times
# ==============================================================================
#
#  START save results for time series graphics output
# ==============================================================================
        if t in saves:
          if k == nominal_rate_index:
            psd_save[:,save_index]=psdt
            cumpsd_save[:,save_index]=cumsum[:]
          median_save[save_index] = d50
          done = False
          for i in range(num_psd):
            if cumsum[i] > 0.1:
              if i > 0:
                d10 = (dia[i-1]+(0.1-cumsum[i-1])/
                        (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
              else:
                d10 = max(0.0,(0.1-cumsum[i])/cumsum[i]*dia[i])
              break
          d10_save[save_index] = d10
          for i in range(num_psd):
            if cumsum[i] > 0.9:
              if i > 0:
                d90 = (dia[i-1]+(0.9-cumsum[i-1])/
                        (cumsum[i]-cumsum[i-1])*(dia[i]-dia[i-1]))
              else:
                d90 = max(0.0,(0.9-cumsum[i])/cumsum[i]*dia[i])
              break
          d90_save[save_index] = d90
          ave_save[k,save_index] = ave
          sa_save[k,save_index] = numpy.sum(numpy.dot(psdt,sa_dia))
          stone_save[k,save_index] = stone
          save_index = save_index + 1
        psd = psdt
        
      if k == nominal_rate_index:
        ts_ave[daughter] =copy.copy(ave_save[k,:])
        ts_stone[daughter] =copy.copy(stone_save[k,:])
      if debug:
        print 'average diameter END = ',t,ave
#
#  END save results for time series graphics output
# ==============================================================================
#
#  END timestep loop
# ==============================================================================
#  OUTPUT PSD, SA and grading stats with time for each fragmentation model
# ==============================================================================

    pstring = ('Particle size weathering\nnum_psd='+str(num_psd)+
              ', daughter='+daughter+', weathering rate='+str(w_rate))
    fig1 = matplotlib.pyplot.figure(figsize=fig_size)
    fig1.suptitle(pstring)
    if cambridge:
      cambridge_fig_defs()

    try:
      contents5=fig1.add_subplot(141)
      for k in range(num_sens):
        stuff4 = contents5.semilogy(saves, ave_save[k,:],
                                label='rate ='+str(rate_range[k]*w_rate)[:5],
                                ls = LINE_STYLE[k], lw=2)
      contents5.set_title('Mean Diameter')
      contents5.legend()
      matplotlib.pyplot.xlabel('Timesteps')
      matplotlib.pyplot.ylabel('Diameter (mm)')
    except:
      import traceback
      traceback.print_exc()

    try:
      contents5=fig1.add_subplot(142)
      for k in range(num_sens):
        stuff4 = contents5.semilogy(saves, stone_save[k,:],
                                label='rate ='+str(rate_range[k]*w_rate)[:5],
                                ls = LINE_STYLE[k], lw=2)
      contents5.set_title('Stoniness')
      contents5.legend()
      matplotlib.pyplot.xlabel('Timesteps')
      matplotlib.pyplot.ylabel('Stoniness')
    except:
      import traceback
      traceback.print_exc()

    try:
      contents6=fig1.add_subplot(143)
      for k in range(num_sens):
        stuff6 = contents6.semilogy(saves, sa_save[k,:],
                                label='rate ='+str(rate_range[k]*w_rate)[:5],
                                ls = LINE_STYLE[k], lw=2)
      contents6.set_title('Specific Area')
      contents6.legend()
      matplotlib.pyplot.xlabel('Timesteps')
      matplotlib.pyplot.ylabel('Specific Area (mm$^3$/mm$^2$)')
      if fig_save:
        temp = fig_name[0]+str(fig_number)+fig_name[1]
        matplotlib.pyplot.savefig(temp, dpi=fig_resolution_cambridge)
        fig_number=fig_number+1
    except:
      import traceback
      traceback.print_exc()

    try:
      contents1=fig1.add_subplot(144)
#      stuff1 = contents1.semilogx(dia[:], psd_save[:,:], lw=2)
      stuff1 = contents1.semilogx(dia[:], cumpsd_save[:], lw=2)
      contents1.set_title('normalised PSD')
      matplotlib.pyplot.xlim([1.e-4,None])
    except:
      import traceback
      traceback.print_exc()

    try:
      fig = matplotlib.pyplot.figure(figsize=fig_size0)
      if cambridge:
        cambridge_fig_defs()
      for t in cumpsdt_save.keys():
        try:
          contents1=fig.add_subplot(111)
          stuff1 = contents1.semilogx(dia[:], cumpsd_save[:], lw=2)
          contents1.set_title('time vs PSD: daughter='+daughter+
                '\nEvery '+str(saves[1])+' timesteps')
          matplotlib.pyplot.xlim([1.e-4,None])
        except:
          import traceback
          traceback.print_exc()
    except:
      import traceback
      traceback.print_exc()

#
#  END different fragmentation models loop
# ==============================================================================
#  OUTPUT PSDs for all the fragmentation models used
# ==============================================================================
  try:
    if calibrate_auto:
#  Results have been saved when results have converged and so are for different
#  times for each of the fragmentation models
#      print cum_obj
#      print cal_sort
#      print calibrate_data
#      print cal_d
#      print cal_mass
      temp = cal_d.keys()
      if len(temp) == 1:
        name = temp[0]
      else:
        name = ''
      if init_psd_filename != '' :
        temp1 = init_psd_d.keys()
        if len(temp1) == 1:
          init_psd_name = temp1[0]
        else:
          init_psd_name = ''
      
      fig = matplotlib.pyplot.figure(figsize=fig_size0)
      if cambridge:
        cambridge_fig_defs()
      contents=fig.add_subplot(111)
      order =daughter_list1
      order.sort()
      for dataset in order:
        if 'stoch' in dataset:
          stuff6 = contents.semilogx(dia, cum_obj[dataset],
#          stuff6 = contents.semilogx(dia, cumpsdt_save[t][dataset],
                                    label=dataset,
                                    ls='--', lw=2)
        else:
          stuff6 = contents.semilogx(dia, cum_obj[dataset],
#          stuff6 = contents.semilogx(dia, cumpsdt_save[t][dataset],
                                    label=dataset,
                                    ls=LINE_STYLE[k], lw=2)
      if calibrate_data and name != '' :
        stuff6 = contents.semilogx(cal_d[name], cal_mass[name],
                                    label='Calibration data',
                                    ls='None', lw=4,
                                    marker='o', mec='k', mfc='k')
      
      if init_psd_filename != '' :
        stuff6 = contents.semilogx(init_psd_d[init_psd_name], init_psd_mass[init_psd_name],
                                    label='Initial data',
                                    ls='None', lw=4,
                                    marker='o', mec='k', mfc='r')
      if obj == 1:
        titlestr = 'Fitted PSD, mean='+str(dia_obj)+' mm, t='+str(t)
      elif obj == 2:
        titlestr = 'Fitted PSD, d50='+str(dia_obj)+' mm, t='+str(t)
      elif obj >= 11 and obj <= 21:
        titlestr = 'Fitted PSD, '+OBJ_TEXT[obj]+'='+str(dia_obj)+' mm, t='+str(t)
      if calibrate_auto and calibrate_data:
        titlestr=titlestr+'\n'+name
      contents.set_title(titlestr)
      contents.legend(loc='best')
      matplotlib.pyplot.ylabel('Proportion of mass less than')
      matplotlib.pyplot.xlabel('Diameter (mm)')
      matplotlib.pyplot.ylim([0.0,1.0])
      matplotlib.pyplot.xlim([1.e-4,None])
      if fig_save:
        temp = fig_name[0]+str(fig_number)+fig_name[1]
        matplotlib.pyplot.savefig(temp, dpi=fig_resolution_cambridge)
        fig_number=fig_number+1
    else:
# This is different from calibrate auto because we know exactly what times ALL
# results output are saved at (in timesort)
      for t in cal_timesort:
        name = ''
        fig = matplotlib.pyplot.figure(figsize=fig_size)
        if cambridge:
          cambridge_fig_defs()
        contents=fig.add_subplot(111)
        order = cumpsdt_save[t].keys()
        order.sort()
        for dataset in order:
          if 'stoch' in dataset:
  #          stuff6 = contents.semilogx(dia, cum_obj[dataset],
            stuff6 = contents.semilogx(dia, cumpsdt_save[t][dataset],
                                      label=dataset,
                                      ls='--', lw=2)
          else:
  #          stuff6 = contents.semilogx(dia, cum_obj[dataset],
            stuff6 = contents.semilogx(dia, cumpsdt_save[t][dataset],
                                      label=dataset,
                                      ls=LINE_STYLE[k], lw=2)
        if calibrate_data:
          name = cal_sort[t]
          stuff6 = contents.semilogx(cal_d[name], cal_mass[name],
                                      label='Calibration data',
                                      ls='None', lw=4,
                                      marker='o', mec='k', mfc='k')
        if obj == 1:
          titlestr = 'Fitted PSD, mean='+str(dia_obj)+' mm, t='+str(t)
        elif obj == 2:
          titlestr = 'Fitted PSD, d50='+str(dia_obj)+' mm, t='+str(t)
        elif obj >= 11 and obj <= 21:
          titlestr = 'Fitted PSD, '+OBJ_TEXT[obj]+'='+str(dia_obj)+' mm, t='+str(t)
        if calibrate_auto and calibrate_data:
          titlestr=titlestr+'\n'+name
        contents.set_title(titlestr)
        contents.legend(loc='best')
        matplotlib.pyplot.ylabel('Proportion of mass less than')
        matplotlib.pyplot.xlabel('Diameter (mm)')
        matplotlib.pyplot.ylim([0.0,1.0])
        matplotlib.pyplot.xlim([1.e-4,None])
        if fig_save:
          temp = fig_name[0]+str(fig_number)+fig_name[1]
          matplotlib.pyplot.savefig(temp, dpi=fig_resolution_cambridge)
          fig_number=fig_number+1
  except:
    import traceback
    traceback.print_exc()
    
#
#  ===================================================
#                       PLOT FINAL OUTPUT
#  ===================================================
#

  try:
    fig = matplotlib.pyplot.figure(figsize=fig_size0)
    if cambridge:
      cambridge_fig_defs()
    contents=fig.add_subplot(111)
    order = cum_obj.keys()
    order.sort()
    for dataset in order:
      if 'stoch' in dataset:
        stuff6 = contents.semilogy(saves, ts_ave[dataset],
                                  label=dataset,
                                  ls='--', lw=2)
      else:
        stuff6 = contents.semilogy(saves, ts_ave[dataset],
                                  label=dataset,
                                  ls=LINE_STYLE[k], lw=2)
    contents.set_title('Average diameter with time, weathering rate='+
                            str(w_rate))
    contents.legend(loc='best')
    matplotlib.pyplot.ylabel('Diameter (mm)')
    matplotlib.pyplot.xlabel('Time')
    matplotlib.pyplot.ylim([1.e-4,None])
    if fig_save:
      temp = fig_name[0]+str(fig_number)+fig_name[1]
      matplotlib.pyplot.savefig(temp, dpi=fig_resolution_cambridge)
      fig_number=fig_number+1
  except:
    import traceback
    traceback.print_exc()

  try:
    fig = matplotlib.pyplot.figure(figsize=fig_size0)
    if cambridge:
      cambridge_fig_defs()
    contents=fig.add_subplot(111)
    order = cum_obj.keys()
    order.sort()
    for dataset in order:
      if 'stoch' in dataset:
        stuff6 = contents.semilogy(saves, ts_stone[dataset],
                                  label=dataset,
                                  ls='--', lw=2)
      else:
        stuff6 = contents.semilogy(saves, ts_stone[dataset],
                                  label=dataset,
                                  ls=LINE_STYLE[k], lw=2)
    contents.set_title('Average stoniness with time, weathering rate='+
                            str(w_rate))
    contents.legend(loc='best')
    matplotlib.pyplot.ylabel('Stoniness')
    matplotlib.pyplot.xlabel('Time')
    matplotlib.pyplot.ylim([1.e-2,None])
    if fig_save:
      temp = fig_name[0]+str(fig_number)+fig_name[1]
      matplotlib.pyplot.savefig(temp, dpi=fig_resolution_cambridge)
      fig_number=fig_number+1
  except:
    import traceback
    traceback.print_exc()

#
#   SYNTHETIC STONINESS SOIL PROFILES
#

  try:
    fig = matplotlib.pyplot.figure(figsize=fig_size0)
    if cambridge:
      cambridge_fig_defs()
    order = cum_obj.keys()
    order.sort()

    contents=fig.add_subplot(221)
    saves_t_const = []
    for t in saves:
      saves_t_const.append(t/float(saves[-1]))
    for dataset in order:
      if 'stoch' in dataset:
        stuff6 = contents.plot(ts_stone[dataset],saves,
                                  label=dataset,
                                  ls='--', lw=2)
      else:
        stuff6 = contents.plot(ts_stone[dataset],saves,
                                  label=dataset,
                                  ls=LINE_STYLE[k], lw=2)
    contents.set_title('Stoniness profile (const weathering), rate='+
                            str(w_rate))
    contents.legend(loc='best')
    matplotlib.pyplot.ylabel('Time')
    matplotlib.pyplot.xlabel('Stoniness')
    matplotlib.pyplot.xlim([1.e-2,None])
    matplotlib.pyplot.ylim([None,1000])

    contents=fig.add_subplot(222)
    saves_t_revexp = []
    for t in saves:
      saves_t_revexp.append(1-numpy.exp(t/float(saves[-1])-t))
    for dataset in order:
      if 'stoch' in dataset:
        stuff6 = contents.plot(ts_stone[dataset],saves_t_revexp,
                                  label=dataset,
                                  ls='--', lw=2)
      else:
        stuff6 = contents.plot(ts_stone[dataset],saves_t_revexp,
                                  label=dataset,
                                  ls=LINE_STYLE[k], lw=2)
    contents.set_title('Stoniness profile (rev exp weathering), rate='+
                            str(w_rate))
    contents.legend(loc='best')
    matplotlib.pyplot.ylabel('Time')
    matplotlib.pyplot.xlabel('Stoniness')
#    matplotlib.pyplot.xlim([1.e-2,None])
#    matplotlib.pyplot.ylim([None,1000])

    contents=fig.add_subplot(223)
    saves_t_exp = []
    for t in saves:
      saves_t_exp.append(numpy.exp(1-t/float(saves[-1])))
    for dataset in order:
      if 'stoch' in dataset:
        stuff6 = contents.plot(ts_stone[dataset],saves_t_exp,
                                  label=dataset,
                                  ls='--', lw=2)
      else:
        stuff6 = contents.plot(ts_stone[dataset],saves_t_exp,
                                  label=dataset,
                                  ls=LINE_STYLE[k], lw=2)

    contents.set_title('Stoniness profile (exp weathering), rate='+
                            str(w_rate))
    contents.legend(loc='best')
    matplotlib.pyplot.ylabel('Time')
    matplotlib.pyplot.xlabel('Stoniness')
#    matplotlib.pyplot.xlim([1.e-2,None])
#    matplotlib.pyplot.ylim([None,1000])

    if fig_save:
      temp = fig_name[0]+str(fig_number)+fig_name[1]
      matplotlib.pyplot.savefig(temp, dpi=fig_resolution_cambridge)
      fig_number=fig_number+1
  except:
    import traceback
    traceback.print_exc()

  if not calibrate_auto:
    output_stats(cal_stats)

#  print 'stone_save',stone_save
#  print 'psd_save',psd_save[0:2,0:]
#  print 'cumpsd_save',cumpsd_save[0:2,0:]

#  print 'Hit any key to finish'

#  print '### ALL figures displayed'
  matplotlib.pyplot.show()
  
#  raw_input('Hit any key to finish')
#  matplotlib.pyplot.close('all')
  


