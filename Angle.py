from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.timeseries import TimeSeries
import cdflib
import sys
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MaxNLocator, MultipleLocator, AutoMinorLocator)
import matplotlib.dates as mdates
from sunpy.timeseries import TimeSeries

trange = a.Time('2009/03/12', '2009/03/12')
dataset = a.cdaweb.Dataset('OMNI_HRO2_1MIN')
result = Fido.search(trange, dataset)
print(dataset)

print(result)

downloaded_files = Fido.fetch(result[0, 0:])
print(downloaded_files)

cdf_data = cdflib.CDF(downloaded_files[0])

#print(cdf_data.cdf_info())
#print(cdf_data.varattsget('F'))
#print(cdf_data.varget('F'))

tim = []
magne_tot = []
magne_x = []
magne_y = []
magne_z = []
speed = []
density = []
temperature = []
beta1 = []



for files in downloaded_files:
    #if filename.endswith('.cdf'):  # Check if the file is a CDF file
    print(files)
    file = cdflib.CDF(files)
    ti = file.varget('Epoch')
    ma_val = file.varget('F')
    ma_x = file.varget('BX_GSE')
    ma_y = file.varget('BY_GSE')
    ma_z = file.varget('BZ_GSE')
    spee = file.varget('flow_speed')
    densi = file.varget('proton_density')
    temper = file.varget('T')
    bet = file.varget('Beta')
    
    
    tim.append(ti)
    magne_tot.append(ma_val)
    magne_x.append(ma_x)
    magne_y.append(ma_y)
    magne_z.append(ma_z)
    speed.append(spee)
    density.append(densi)
    temperature.append(temper)
    beta1.append(bet)
    
sys.path.append(r"file path")
from functions_MVA_n import cdfdtime, var_atts, genfillval_nan, MVA, orie_GSE, combine 

time = combine(tim)
dtim = cdfdtime(time)
print(dtim)

mag_tot = combine(magne_tot)
mag_x = combine(magne_x)
mag_y = combine(magne_y)
mag_z = combine(magne_z)
spe = combine(speed)
den = combine(density)
temp =combine(temperature)
beta = combine(beta1)
mag_fillval = var_atts(cdf_data, 'F', 'FILLVAL')

spe_fillval = var_atts(cdf_data, 'flow_speed', 'FILLVAL')

den_fillval = var_atts(cdf_data, 'proton_density', 'FILLVAL') 

temp_fillval = var_atts(cdf_data, 'T', 'FILLVAL') 

beta_fillval = var_atts(cdf_data, 'Beta', 'FILLVAL') 

print(mag_fillval)

mag_tot_nan = genfillval_nan(mag_tot, mag_fillval)

mag_x_nan = genfillval_nan(mag_x, mag_fillval)

mag_y_nan = genfillval_nan(mag_y, mag_fillval)

mag_z_nan = genfillval_nan(mag_z, mag_fillval)

spe_nan = genfillval_nan(spe, spe_fillval)

den_nan = genfillval_nan(den, den_fillval)

temp_nan = genfillval_nan(temp, temp_fillval)

beta_nan = genfillval_nan(beta, beta_fillval)

print(mag_x_nan)

mc_start_time = datetime.strptime('12-03-2009 02:00:00.806000', '%d-%m-%Y %H:%M:%S.%f')
mc_end_time = datetime.strptime('12-03-2009 09:00:00.806000', '%d-%m-%Y %H:%M:%S.%f')

time_mc_requ = []

for tim in dtim:
    #print(type(tim))
    if tim>=  np.datetime64(mc_start_time) and tim <=  np.datetime64(mc_end_time):
        time_mc_requ.append(tim)

ind_mc = np.where((np.array(dtim) >= np.datetime64(mc_start_time)) & (np.array(dtim) <= np.datetime64(mc_end_time)))[0]
mag_tot_mc = mag_tot_nan[ind_mc]
#print((mag_tot_mc))

mag_x_mc = mag_x_nan[ind_mc]

mag_y_mc = mag_y_nan[ind_mc]

mag_z_mc = mag_z_nan[ind_mc]

print(mag_x_nan[34861])

B_mat = np.transpose(np.array([mag_x_mc, mag_y_mc, mag_z_mc]))

#print(B_mat)

#c, val, vec are the minimum varience matrix, eigen values in diagonal form, and are the eigen vectors respective to eigen value
# for the stereo spacecraft

c, val, vec = MVA(B_mat)

int_vec = vec[:,1]

print(int_vec)
lon, lat = orie_GSE(int_vec)

print('Longitude of the axis:',lon)
print('Latitude of the axis:',lat)

start_time = datetime.strptime('12-03-2009 01:00:00.806000', '%d-%m-%Y %H:%M:%S.%f')
end_time = datetime.strptime('13-03-2009 02:00:00.806000', '%d-%m-%Y %H:%M:%S.%f')

time_requ = []

for tim in dtim:
    if tim>= np.datetime64(start_time) and tim <=  np.datetime64(end_time):
        time_requ.append(tim)
        
#print((time_requ))

ind_requ = np.where((np.array(dtim) >= np.datetime64(start_time)) & (np.array(dtim) <= np.datetime64(end_time)))[0]

#print(len(mag_tot[ind_requ]))

mag_tot_requ = mag_tot_nan[ind_requ]
mag_x_requ = mag_x_nan[ind_requ]
mag_y_requ = mag_y_nan[ind_requ]
mag_z_requ = mag_z_nan[ind_requ]
spe_requ = spe_nan[ind_requ]
den_requ = den_nan[ind_requ]
temp_requ = temp_nan[ind_requ]
beta_requ = beta_nan[ind_requ]

print(len(mag_tot_requ))

label = 14
font =20
maj_tick = 5
min_tick = 3
#width = 0.6
row_panel = 6
col_panel = 1
alpha = 0.08
marker = ','
marksiz = 0.2074



plt.figure(figsize=(8, 16)) 
plt.rcParams['figure.dpi']=500.0


plt.subplot(row_panel,col_panel,1)
plt.scatter(time_requ, mag_tot_requ, marker = marker, s = marksiz, color = 'black')
plt.ylabel('B \n (nT)', fontsize = font)
#plt.ylim(min(mag_tot_win_nan) - 5, max(mag_tot_win_nan) +5)
#plt.title('Magnetic Field on STEREO-A', fontsize = font)


plt.subplot(row_panel,col_panel,2)
plt.scatter(time_requ, mag_x_requ, marker = marker, s = marksiz, color = 'green', label = 'X')
plt.scatter(time_requ, mag_y_requ, marker = marker, s = marksiz, color = 'darkorange', label = 'Y')
plt.scatter(time_requ, mag_z_requ, marker = marker, s = marksiz, color = 'blue', label = 'Z')
plt.ylabel('$\mathrm{B}_\mathrm{XYZ}$ \n (nT)', fontsize = font)
plt.legend(frameon = False, fontsize = 12, bbox_to_anchor=[0.0, 1], loc = 'upper left', markerscale = 10)
#plt.ylim(min(min(mag_R_win_nan),min(mag_T_win_nan),min(mag_N_win_nan)) - 10, max(max(mag_R_win_nan),max(mag_T_win_nan),max(mag_N_win_nan)) + 10)




plt.subplot(row_panel,col_panel,3)
plt.scatter(time_requ, spe_requ, marker = marker, s = marksiz, color = 'black')
plt.ylabel('$\mathrm{V}_\mathrm{p}$ \n (km s$^{-1}$)', fontsize = font)
#plt.ylim(min(vel_win_nan) - 50, max(vel_win_nan) +50)
plt.tick_params(right = True, top = True, labeltop = False,labelbottom =False, axis = 'both', which='major', length = maj_tick, direction ='in', labelsize = label, width = 1.0)
plt.tick_params(right = True,top = True, labeltop = False,labelbottom =False, axis='both', which='minor', length = min_tick, direction = 'in', labelsize = label, width = 1.0)
plt.gca().yaxis.set_major_locator(MultipleLocator(100))
plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))



#plt.plot(dtime, mag_T, linewidth = width, color = 'darkorange')
#plt.ylabel(r'$\text{B}_{\text{T}}$ (nT)', fontsize = font)



plt.subplot(row_panel,col_panel,4)
plt.scatter(time_requ, den_requ, marker = marker, s = marksiz, color = 'black')
plt.ylabel('$\mathrm{N_p}$ \n (cm$^{-3}$)', fontsize = font)
#plt.ylim(min(den_win_nan) - 5, max(den_win_nan) +5)
plt.tick_params(right = True, top = True, labeltop = False,labelbottom =False, axis = 'both', which='major', length = maj_tick, direction ='in', labelsize = label, width = 1.0)
plt.tick_params(right = True,top = True, labeltop = False,labelbottom =False, axis='both', which='minor', length = min_tick, direction = 'in', labelsize = label, width = 1.0)
plt.gca().yaxis.set_major_locator(MultipleLocator(20))
plt.gca().yaxis.set_minor_locator(AutoMinorLocator(4))



#plt.plot(dtime, mag_N, linewidth = width, color = 'blue')
#plt.ylabel(r'$\text{B}_{\text{N}}$ (nT)', fontsize = font)
#plt.ylim(min(mag_N) - 5, max(mag_N) +5)



plt.subplot(row_panel,col_panel,5)
plt.scatter(time_requ, temp_requ/10000, marker = marker, s = marksiz,color = 'black')
plt.ylabel('$\mathrm{T_p}$ \n (10$^4$ K)', fontsize = font)
#plt.ylim(min(temp_win_nan/10000) - 5, max(temp_win_nan/10000) +5)
plt.tick_params(right = True, top = True, labeltop = False,labelbottom =False, axis = 'both', which='major', length = maj_tick, direction ='in', labelsize = label, width = 1.0)
plt.tick_params(right = True,top = True, labeltop = False,labelbottom =False, axis='both', which='minor', length = min_tick, direction = 'in', labelsize = label, width = 1.0)
plt.gca().yaxis.set_major_locator(MultipleLocator(20))
plt.gca().yaxis.set_minor_locator(AutoMinorLocator(4))




plt.subplot(row_panel,col_panel,6)
plt.scatter(time_requ, beta_requ, marker = marker, s = marksiz,color = 'black')
plt.ylabel(r'$\beta$', fontsize = font)
plt.tick_params(right = True, top = True, labeltop = False,labelbottom =True, axis = 'both', which='major', length = maj_tick, direction ='in', labelsize = label, width = 1.0)
plt.tick_params(right = True,top = True, labeltop = False,labelbottom =True, axis='both', which='minor', length = min_tick, direction = 'in', labelsize = label, width = 1.0)
#plt.gca().yaxis.set_major_locator(LogLocator(base = 10))
#plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs =  'auto'))
#plt.yscale('log')


for subplot in range(1,3):
    plt.subplot(row_panel,col_panel,subplot)
    plt.tick_params(right = True, top = True, labeltop = False,labelbottom =False, axis = 'both', which='major', length = maj_tick, direction ='in', labelsize = label, width = 1.0)
    plt.tick_params(right = True,top = True, labeltop = False,labelbottom =False, axis='both', which='minor', length = min_tick, direction = 'in', labelsize = label, width = 1.0) 
    #plt.gca().yaxis.set_major_locator(MultipleLocator(20))
    #plt.gca().yaxis.set_minor_locator(AutoMinorLocator(4))


locator = mdates.HourLocator(interval =6)
for subplot in range(1, row_panel+1):
    plt.subplot(row_panel, col_panel, subplot)
    plt.gca().xaxis.set_major_locator(locator)
    date_format = mdates.DateFormatter("%d %b\n%H:%M")
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
    #plt.axvline(shock_time_st, linewidth = width, color = 'deepskyblue', linestyle = '--')
    #plt.axvline(sheath_end_time_st, linewidth = width, color = 'red', linestyle = '--')
    #plt.axvline(mc_end_time_st, linewidth = width, color = 'red', linestyle = '--')
    #plt.axvspan(shock_time_wi, sheath_end_time_wi, color='red', alpha = alpha, linewidth = 0)
    #plt.axvspan(sheath_end_time_wi, mc_end_time_wi, color='yellow', alpha = alpha, linewidth = 0)
    #plt.ylim(min(min(magn_To), min(magn_R), min(magn_T), min(magn_N)) - 5, max(max(magn_To), max(magn_R), max(magn_T), max(magn_N)) + 5)
    #plt.xlim(dtime[0], dtime[-50])



plt.xlabel('Time in {} (UT)'.format(start_time.year), fontsize=font)
plt.subplots_adjust(hspace = 0.09)

   
