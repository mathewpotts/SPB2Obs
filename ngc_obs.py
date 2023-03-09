import ephem
import math
import numpy as np
import matplotlib.pyplot as plt

obs = ephem.Observer() #make observer

in_obs = ['2023/03/07 00:00','38.520278', '-113.2875', 0]
in_obj = [['NGC1068','f|G','172.1037082138010','-51.9337915309647','8.9']]

obs.date = in_obs[0] #UTC time
obs.lat, obs.lon, obs.elevation = in_obs[1:] #SPB2 coords
obs.horizon = 0 #horizon

NGC = []
for x in range(len(in_obj)):
	ngc_eq = ephem.Equatorial(ephem.Galactic(in_obj[x][2], in_obj[x][3], epoch = ephem.J2000), epoch = ephem.J2000) #converting to equatorial
	ngc_xephem_format = in_obj[x][0] + ',' + in_obj[x][1] + ',' + str(ngc_eq.ra) + ',' + str(ngc_eq.dec) + ',' + in_obj[x][4]#supplying fixed coord data in xephem format (https://xephem.github.io/XEphem/Site/help/xephem.html#mozTocId800642)
	NGC.append(ephem.readdb(ngc_xephem_format)) #create ncg object

days = 3 #days
step_hours = 0.08333333333 #step time in hours

print('Start Date (UTC):', obs.date)
print('Observer Location (Lat, Long):', float(repr(obs.lat)) * 180. / math.pi, float(repr(obs.lon)) * 180. / math.pi)
print("Duration (Days):", days)
print("Timestep (Hr.)", step_hours)
print()
print('FORMAT:')
print('ALTITUDE (DEG), AZIMUTH (DEG), DATE (Y/M/D), TIME (UTC), ABOVE HORIZON? (Y/N), CROSSED HORIZON?\n')

for i in range(len(in_obj)):
	NGC[i].compute(obs)
	prev_alt = float(repr(NGC[i].alt)) * 180. / math.pi

azis = np.zeros([len(in_obj),int(days * 24 / step_hours)])
times = np.empty(int(days * 24 / step_hours))

for x in range(int(days * 24 / step_hours)):
	for i in range(len(in_obj)):
		NGC[i].compute(obs) #compute horizontal coords at obs location time and date
		if ((float(repr(NGC[i].alt)) * 180. / math.pi) >= float(repr(obs.horizon))) and (prev_alt < float(repr(obs.horizon))):
			print(float(repr(NGC[i].alt)) * 180. / math.pi, float(repr(NGC[i].az)) * 180. / math.pi, obs.date, 'Y, CROSS!')
			azis[i][x] = float(repr(NGC[i].az)) * 180. / math.pi
		elif ((float(repr(NGC[i].alt)) * 180. / math.pi) <= float(repr(obs.horizon))) and (prev_alt > float(repr(obs.horizon))):
			print(float(repr(NGC[i].alt)) * 180. / math.pi, float(repr(NGC[i].az)) * 180. / math.pi, obs.date, 'N, CROSS!')
			azis[i][x] = None
		elif (float(repr(NGC[i].alt)) * 180. / math.pi) >= float(repr(obs.horizon)):
			print(float(repr(NGC[i].alt)) * 180. / math.pi, float(repr(NGC[i].az)) * 180. / math.pi, obs.date, 'Y')
			azis[i][x] = float(repr(NGC[i].az)) * 180. / math.pi
		else: #if below horizon
			print(float(repr(NGC[i].alt)) * 180. / math.pi, float(repr(NGC[i].az)) * 180. / math.pi, obs.date, 'N')
			azis[i][x] = None

		prev_alt = float(repr(NGC[i].alt)) * 180. / math.pi
	
	times[x] = obs.date.tuple()[2]*24 + obs.date.tuple()[3] + obs.date.tuple()[4]/60 + obs.date.tuple()[5]/3600
	obs.date = obs.date + step_hours * ephem.hour #update date

for i in range(len(in_obj)):
	plt.plot(times[azis[i] != None],azis[i][azis[i] != None],'.')
plt.xlabel('Time [hrs]')
plt.ylabel('Azimuth [$^\circ$]')
plt.show()