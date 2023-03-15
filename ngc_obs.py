import ephem
import math
import numpy as np
import matplotlib.pyplot as plt
import sys

inputfile = sys.argv[1]

with open(inputfile,'r') as f:
	in_obs = f.readline().split(',') #reads first line for observer inputs
	in_obs[-1] = float(in_obs[-1]) #sets elevation term to float; ephem needs numeric for elevation
	in_obj = [f.readline().split(',')] #read second line/first observable line
	line = f.readline()
	if line != '': #checks if next lines are empty before adding them to observables list
		in_obj.append(line.split(','))
	else:
		f.close()

obs = ephem.Observer() #make observer
obs.date = in_obs[0] #UTC time
obs.lat, obs.lon, obs.elevation = in_obs[1:] #SPB2 coords
obs.horizon = (np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+in_obs[3])) #horizon calculated from elevation of obs

s = ephem.Sun() #make Sun
m = ephem.Moon() #make Moon

NGC = [] #start list to hold NGC objects for each observable
for x in range(len(in_obj)):
	ngc_eq = ephem.Equatorial(ephem.Galactic(in_obj[x][2], in_obj[x][3], epoch = ephem.J2000), epoch = ephem.J2000);print(in_obj[x][2], in_obj[x][3]) #converting to equatorial
	ngc_xephem_format = in_obj[x][0] + ',' + in_obj[x][1] + ',' + str(ngc_eq.ra) + ',' + str(ngc_eq.dec) + ',' + in_obj[x][4];print(ngc_xephem_format);
	NGC.append(ephem.readdb(ngc_xephem_format)) #create ncg object

days = 3 #integration period in days
step_hours = 0.08333333333 #step time in hours

for i in range(len(in_obj)):
	NGC[i].compute(obs)
	prev_alt = float(repr(NGC[i].alt)) * 180. / math.pi

azis = np.zeros([len(in_obj),int(days * 24 / step_hours)]) #azimuths for observables
alts = np.zeros([len(in_obj),int(days * 24 / step_hours)]) #altitudes for observables
altSun = np.zeros(int(days * 24 / step_hours))
phaseMoon = np.zeros(int(days * 24 / step_hours))
altMoon = np.zeros(int(days * 24 / step_hours))
times = np.empty(int(days * 24 / step_hours)) #times
dates = [] #dates

prev_alt = np.empty(len(in_obj))

for x in range(int(days * 24 / step_hours)):
	for i in range(len(in_obj)):
		NGC[i].compute(obs)
		azis[i][x] = float(repr(NGC[i].az)) * 180. / math.pi
		alts[i][x] = float(repr(NGC[i].alt)) * 180. / math.pi
	s.compute(obs)
	altSun[x] = float(repr(s.alt)) * 180 / math.pi
	m.compute(obs)
	phaseMoon[x] = float(repr(m.moon_phase))
	altMoon[x] = float(repr(m.alt)) * 180 / math.pi
	times[x] = obs.date.tuple()[2]*24 + obs.date.tuple()[3] + obs.date.tuple()[4]/60 + obs.date.tuple()[5]/3600 #add time in hours to array
	dates.append(str(obs.date))
	obs.date = obs.date + step_hours * ephem.hour #update date

altmask = (alts >= -5.0) & (alts <= 2.5)
sunmask = altSun <= ((float(repr(obs.horizon)) * 180 / math.pi) - 15.0) #set condition that sun is 
moonmask = phaseMoon <= 30.0

#print to text file
for i in range(len(in_obj)):
	with open('{}.txt'.format(in_obj[i][0]),'w') as f:
		f.write('Start Date (UTC): {}\n'.format(in_obs[0]))
		f.write('Observer Location (Lat, Long, Elevation): ({lat},{long},{elevation})\n'.format(lat=in_obs[1],long=in_obs[2],elevation=in_obs[3]))
		f.write('Format:\n')
		f.write('ALTITUDE (DEG), AZIMUTH (DEG), DATE (Y/M/D), TIME (UTC), CROSSED HORIZON?\n')
		if len((altmask[i]&sunmask&moonmask).nonzero()[0]) != 0:
			previndex = (altmask[i]&sunmask&moonmask).nonzero()[0][0]
		for j in (altmask[i]&sunmask&moonmask).nonzero()[0]:
			if j == previndex + 1:
				f.write(str(alts[i][j]) + ',' + str(azis[i][j]) + ',' + dates[j] + 'N\n')
			else:
				f.write(str(alts[i][j]) + ',' + str(azis[i][j]) + ',' + dates[j] + 'Y\n')
			previndex = j
		f.close()

for i in range(len(in_obj)):
	plt.plot(times[altmask[i]&sunmask&moonmask],azis[i][altmask[i]&sunmask&moonmask],'.',label=in_obj[i][0])
plt.xlabel('Time [hrs]')
plt.ylabel('Azimuth [$^\circ$]')
plt.legend()
plt.ylim((0,360))
plt.savefig('temp.png',dpi=1200)
