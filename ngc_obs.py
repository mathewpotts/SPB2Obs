import ephem
import math

obs = ephem.Observer() #make observer

obs.date = '2022/11/04 00:00' #UTC time
obs.lat, obs.lon = '38.520278', '-113.2875' #firsco peak coords
obs.horizon = 0 #horizon
ngc_gal_long, ngc_gal_lat = '172.1037082138010', '-51.9337915309647' #ngc galactic coords (J2000)
ngc_eq = ephem.Equatorial(ephem.Galactic(ngc_gal_long, ngc_gal_lat, epoch = ephem.J2000), epoch = ephem.J2000) #converting to equatorial

ngc_xephem_format = 'NGC1068,f|G,' + str(ngc_eq.ra) + ',' + str(ngc_eq.dec) + ',8.9,2000' #supplying fixed coord data in xephem format (https://xephem.github.io/XEphem/Site/help/xephem.html#mozTocId800642)

years = 2 #years
step_hours = 0.08333333333 #step time in hours

NGC_1068 = ephem.readdb(ngc_xephem_format) #create ncg object

print('Start Date (UTC):', obs.date)
print('Observer Location (Lat, Long):', float(repr(obs.lat)) * 180. / math.pi, float(repr(obs.lon)) * 180. / math.pi)
print("Duration (Yr.):", years)
print("Timestep (Hr.)", step_hours)
print()
print('FORMAT:')
print('ALTITUDE (DEG), AZIMUTH (DEG), DATE (Y/M/D), TIME (UTC), ABOVE HORIZON? (Y/N), CROSSED HORIZON?\n')

NGC_1068.compute(obs)
prev_alt = float(repr(NGC_1068.alt)) * 180. / math.pi

for x in range(int(365 * years * 24 / step_hours)):
	NGC_1068.compute(obs) #compute horizontal coords at obs location time and date
	
	if ((float(repr(NGC_1068.alt)) * 180. / math.pi) >= float(repr(obs.horizon))) and (prev_alt < float(repr(obs.horizon))):
		print(float(repr(NGC_1068.alt)) * 180. / math.pi, float(repr(NGC_1068.az)) * 180. / math.pi, obs.date, 'Y, CROSS!')
	elif ((float(repr(NGC_1068.alt)) * 180. / math.pi) <= float(repr(obs.horizon))) and (prev_alt > float(repr(obs.horizon))):
		print(float(repr(NGC_1068.alt)) * 180. / math.pi, float(repr(NGC_1068.az)) * 180. / math.pi, obs.date, 'N, CROSS!')
	elif (float(repr(NGC_1068.alt)) * 180. / math.pi) >= float(repr(obs.horizon)):
		print(float(repr(NGC_1068.alt)) * 180. / math.pi, float(repr(NGC_1068.az)) * 180. / math.pi, obs.date, 'Y')
	else: #if below horizon
		print(float(repr(NGC_1068.alt)) * 180. / math.pi, float(repr(NGC_1068.az)) * 180. / math.pi, obs.date, 'N')
	
	prev_alt = float(repr(NGC_1068.alt)) * 180. / math.pi
	
	obs.date = obs.date + step_hours * ephem.hour #update date
