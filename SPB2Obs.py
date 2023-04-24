#!/usr/bin/python3

# Authors: Mathew Potts, Jordan Bogdan, Andrew Wang

#Import libs
import ephem
import math
import numpy as np
#import matplotlib.pyplot as plt
import argparse
from gcn_kafka import Consumer
import tkinter as tk
from PIL import Image, ImageTk
import time
import calendar
import threading
import re
import requests
from bs4 import BeautifulSoup

# Define Mask parameters
SUN_OFFSET   = 15.0
MOON_PERCENT = 30.0

# Init alert flags for moon/sun
SUN_FLAG  = False
MOON_FLAG = False

# Read in args
def read_in_args():
    parser = argparse.ArgumentParser(description = 'SPB2Obs shows Objects of Interest (OoI) in the FoV of the CT telescope displaying the azimuth and altitude of those objects. SPB2Obs incorporates live alerts of Gamma-Ray Bursts (GRBs) from the General Coordinates Network (GCN).')
    parser.add_argument('-obj', metavar='objFile',action='store',help='Path to file containing OoI.')
    parser.add_argument('-loc', metavar='locFile',action='store',help='Path to file containing the current GPS location of SPB2 and time.')
    parser.add_argument('-test', action='store_true',default=False,help='Debugging option.')
    args = parser.parse_args()
    return args

def GUI(args):
    root = tk.Tk()
    root.geometry("920x400") # set default window size
    app = SourcesGUI(root,args)
    root.mainloop()

# Define class to observer
class SPB2Obs:
    def __init__(self,args):
        self.in_obj = self.read_file(args.obj) # read infile objects
        self.ephem_objarray = self.ephem_object_array() # Create Ephem objects

        # Set init observer
        self.obs = ephem.Observer() #make observer
        self.gpsLoc = self.read_file(args.loc) #read in gps location
        self.set_observer(self.gpsLoc) #set observer based on location

        # Init horizons, later it is updated by the horizons function
        self.default_horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.obs.elevation))))
        self.upperfov = self.default_horizon + (math.pi/180)*2.5
        self.lowerfov = self.default_horizon + (math.pi/180)*-5.0

        # Flight location Equatorial
        self.url = "https://www.csbf.nasa.gov/map/balloon10/flight728NT.htm"

        # Init masks
        self.s = ephem.Sun() #make Sun
        self.m = ephem.Moon() #make Moon

    def horizons(self, elevation):
        self.default_horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(elevation))))
        self.obs.horizon = self.default_horizon
        self.upperfov = self.default_horizon + (math.pi/180)*2.5
        self.lowerfov = self.default_horizon + (math.pi/180)*-5.0

    def update_gpsLoc(self):
        response = requests.get(self.url)

        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            rows = soup.find_all('tr')

            # extract payload time, latitude, and longitude from the second row
            data_row = rows[0]
            payload_time = list(reversed(data_row.find_all('center')[1].get_text().strip().split('Z  ')))
            payload_time = '20'+payload_time[0].split('/')[2]+'/'+payload_time[0].split('/')[0]+'/'+payload_time[0].split('/')[1]+' '+payload_time[1]
            latitude = data_row.find_all('center')[2].get_text().strip().split('  ')[1]
            longitude = data_row.find_all('center')[3].get_text().strip().split('  ')[1]
            height = data_row.find_all('center')[4].get_text().strip().split('  ')[1]
            wind = data_row.find_all('center')[5].get_text().strip()

            # Lat/Long string conversion
            latitude = self.lat_long_convert(latitude, 1)
            longitude = self.lat_long_convert(longitude, 0)

            # Update observer
            locArray = [payload_time,latitude,longitude,height]
            self.set_observer(locArray)

            # print the extracted date
            if TESTING:
                print("Payload time:", payload_time)
                print("Latitude:", latitude)
                print("Longitude:", longitude)
                print("Height:", height)
                print("Wind:", wind)

        else:
            print("Failed to retrieve webpage")

    def lat_long_convert(self, string, negFlag):
        # Extract numerical values as strings
        nums_str = re.findall(r'\d+(?:\.\d+)?', string)

        # Convert strings to floats and store in a list
        nums = [float(num_str) for num_str in nums_str]
        value = nums[0] + nums[1]/60
        return str(-1*value) if negFlag else str(value)

    @property
    def gps_loc(self):
        return [self.obs.date,self.obs.lat,self.obs.long,self.obs.elevation]

    @property
    def horizon(self):
        default_horizon = self.rad2degHMS(self.default_horizon)
        upperFoV = self.rad2degHMS(self.upperfov)
        lowerFoV = self.rad2degHMS(self.lowerfov)
        return [default_horizon, upperFoV, lowerFoV]

    def rad2degHMS(self, alpha):
        D = int(alpha * (180/math.pi))
        M = int((alpha * (180/math.pi) - D) * 60)
        S = ((alpha * (180/math.pi) - D) * 60 - M) * 60
        return "{0}:{1}:{2:.2f}".format(D,M,S)

    def check_fov(self, utctime):
        sources = []
        self.obs.date = utctime
        for i in range(len(self.ephem_objarray)):
            sources.append(self.objs_in_fov(utctime, self.ephem_objarray[i]))

        # remove objects that will never be in FoV
        filter_list = []
        for i,s in enumerate(sources):
            if s in self.ephem_objarray:
                filter_list.append(i)
        self.ephem_objarray = [obj for i,obj in enumerate(self.ephem_objarray) if i not in filter_list]
        if TESTING:
            print(self.ephem_objarray)

        sources = [s for s in sources if type(s) == str ] # filter list
        return sources

    def objs_in_fov(self, utctime, ephem_obj):
        try:
            rise = self.obs.next_rising(ephem_obj)
            sett = self.obs.next_setting(ephem_obj)
            if rise > sett: # if source is setting first
                az,alt,mask = self.masks(ephem_obj, utctime)
                if TESTING:
                    gui_str = "{0},{1},{2}".format(ephem_obj.name,az,alt)
                    print(gui_str)
                    return gui_str
                elif alt <= self.upperfov and alt >= self.lowerfov and mask:
                    gui_str = "{0},{1},{2}".format(ephem_obj.name,az,alt)
                    return gui_str
            else: # if source is rising first
                az,alt,mask = self.masks(ephem_obj, utctime)
                if TESTING:
                    gui_str = "{0},{1},{2}".format(ephem_obj.name,az,alt)
                    print(gui_str)
                    return gui_str
                elif alt <= self.upperfov and alt >= self.lowerfov and mask:
                    gui_str = "{0},{1},{2}".format(ephem_obj.name,az,alt)
                    return gui_str
        except ephem.AlwaysUpError:
            print("Warning: Object of interest {0} is always up always up, and is out of the FoV.".format(ephem_obj))
            return ephem_obj
        except ephem.NeverUpError:
            print("Warning: Object of interest {0} is never up, and is out of the FoV.".format(ephem_obj))
            return ephem_obj

    def masks(self, ephem_obj, utctime):
        global SUN_FLAG
        global MOON_FLAG
        self.obs.date = utctime # make sure obs is at current time
        ephem_obj.compute(self.obs) # compute location of object
        self.s.compute(self.obs) # compute location of sun relative to observer
        altSun = float(repr(self.s.alt)) * 180 / math.pi
        self.m.compute(self.obs) #compute location of moon relative to observer
        phaseMoon = float(repr(self.m.moon_phase))
        altMoon = float(repr(self.m.alt)) * 180 / math.pi

        sunmask = altSun <= (((self.default_horizon) * 180 / math.pi) - SUN_OFFSET) #set condition for sun
        moonmask = phaseMoon <= MOON_PERCENT
        if TESTING:
            if not sunmask:
                print(altSun,"Warning: Sun is up!")
            if not moonmask:
                print(altMoon,"Moon phase is higher.")
            SUN_FLAG = True if not sunmask else False
            MOON_FLAG = True if not moonmask else False
            print(SUN_FLAG,MOON_FLAG)
        else:
            SUN_FLAG = True if not sunmask else False
            MOON_FLAG = True if not moonmask else False

        return ephem_obj.az, ephem_obj.alt, sunmask*moonmask

    def set_observer(self, locArray):
        print(locArray)
        ls = locArray[0].split(',') if len(locArray) == 1 else locArray
        if len(locArray) == 1:
            self.obs.date = ls[0] # UTC time
        self.obs.lat = ls[1]
        self.obs.lon = ls[2]
        self.obs.elevation = float(ls[3])
        self.horizons(float(ls[3]))
        self.obs.pressure = 0 # turn off refraction
        print(self.obs)

    def read_file(self, infile):
        # Read input file
        f = open(infile,'r')
        array = f.readlines()
        f.close()
        return array

    def ephem_object_array(self):
        NGC = [] #start list to hold NGC objects for each observable
        for x in range(len(self.in_obj)):
            in_obj = self.in_obj[x].split(',')
            NGC.append(self.create_ephem_object(in_obj))
        return NGC

    def create_ephem_object(self,in_obj):
        ngc_eq = ephem.Equatorial(ephem.Galactic(in_obj[2], in_obj[3], epoch = ephem.J2000), epoch = ephem.J2000) #converting to equatorial
        ngc_xephem_format = in_obj[0] + ',' + in_obj[1] + ',' + str(ngc_eq.ra) + ',' + str(ngc_eq.dec) + ',' + in_obj[4]#supplying fixed coord data in xephem format (https://xephem.github.io/XEphem/Site/help/xephem.html#mozTocId800642)
        return ephem.readdb(ngc_xephem_format)

    def check_sun_and_moon(self, utctime):
        self.obs.date = utctime # make sure obs is at current time

        # Checking sun position
        sun = []
        self.obs.horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.obs.elevation))))- (math.pi/180)*SUN_OFFSET # define horizon for the sun calculation
        if TESTING:
            print('sun',self.obs)
        try:
            sun_rise = self.obs.next_rising(self.s) # sun rise at defined horizon
            sun_set = self.obs.next_setting(self.s) # sun set at defined horizon
        except ephem.AlwaysUpError:
            print("Warning: Sun is always up!")
            sun_rise = 'N/A\t\t'
            sun_set = 'N/A\t\t'
        except ephem.NeverUpError:
            print("Warning: Sun is never up!")
            sun_rise = 'N/A\t\t'
            sun_set = 'N/A\t\t'
        self.s.compute(self.obs) # compute the location of the sun relative to observer
        sun = [sun_rise,sun_set,self.s.az,self.s.alt]

        # Checking moon position
        moon = []
        self.obs.horizon = self.default_horizon # reset horizon back to default for moon calculation
        if TESTING:
            print('moon',self.obs)
        try:
            moon_rise = self.obs.next_rising(self.m) # moon rise at defined horizon
            moon_set = self.obs.next_setting(self.m) # moon rise at defined horizon
        except ephem.AlwaysUpError:
            print("Warning: Moon is always up!")
            moon_rise = 'N/A\t\t'
            moon_set = 'N/A\t\t'
        except ephem.NeverUpError:
            print("Warning: Moon is never up!")
            moon_rise = 'N/A\t\t'
            moon_set = 'N/A\t\t'
        self.m.compute(self.obs) # compute the location of the moon relative to observer
        moon = [moon_rise,moon_set,self.m.az,self.m.alt,self.m.moon_phase]

        return sun,moon

    def gcn_alerts(self):
        # Typical GCN alert
        #b'TITLE:           GCN/FERMI NOTICE\n
        #NOTICE_DATE:     Sat 11 Mar 23 02:01:13 UT\n
        #NOTICE_TYPE:     Fermi-GBM Test Position\n
        #RECORD_NUM:      1\n
        #TRIGGER_NUM:     99999\n
        #GRB_RA:          180.000d {+12h 00m 00s} (J2000),\n
        #180.297d {+12h 01m 11s} (current),\n
        #179.361d {+11h 57m 27s} (1950)\n
        #GRB_DEC:         -35.000d {-35d 00\' 00"} (J2000),\n
        #         -35.129d {-35d 07\' 44"} (current),\n
        #-34.722d {-34d 43\' 17"} (1950)\n
        #GRB_ERROR:       5.80 [deg radius, statistical plus systematic]\n
        #GRB_INTEN:       1000 [cnts/sec]\n
        #DATA_SIGNIF:     11.78 [sigma]\n
        #INTEG_TIME:      0.064 [sec]\n
        #GRB_DATE:        20014 TJD;    70 DOY;   23/03/11\n
        #GRB_TIME:        7265.00 SOD {02:01:05.00} UT\n
        #GRB_PHI:           0.00 [deg]\n
        #GRB_THETA:         0.00 [deg]\n
        #DATA_TIME_SCALE: 0.0000 [sec]\n
        #HARD_RATIO:      3.18\n
        #LOC_ALGORITHM:   1 (version number of)\n
        #MOST_LIKELY:      67%  Below horizon\n
        #2nd_MOST_LIKELY:  32%  Local Particles\n
        #DETECTORS:       1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,\n
        #SUN_POSTN:       351.03d {+23h 24m 07s}   -3.87d {-03d 52\' 01"}\n
        #SUN_DIST:        140.04 [deg]   Sun_angle= 11.4 [hr] (West of Sun)\n
        #MOON_POSTN:      208.95d {+13h 55m 49s}  -11.47d {-11d 28\' 29"}\n
        #MOON_DIST:        35.13 [deg]\n
        #MOON_ILLUM:      88 [%]\n
        #GAL_COORDS:      291.16, 26.69 [deg] galactic lon,lat of the burst (or transient)\n
        #ECL_COORDS:      195.56,-31.75 [deg] ecliptic lon,lat of the bust (or transient)\n
        #COMMENTS:        Fermi-GBM TEST Coordinates.  \n
        #COMMENTS:        This Notice was ground-generated -- not flight-generated.  \n'
        # I assume that all GCN alerts are of similar format
        # wait for a gcn alerts
        # Connect as a consumer.
        # Warning: don't share the client secret with others.
        consumer = Consumer(client_id='7l6biv6qbbgut8rllm527eas0e',
                            client_secret='12ev0n7c81vj0g019f5amdlhjt64vcrecnuadr9nfekr6bikc8ln')

        # Subscribe to topics and receive alerts
        consumer.subscribe(['gcn.classic.text.FERMI_GBM_FIN_POS',
                            'gcn.classic.text.FERMI_GBM_FLT_POS',
                            'gcn.classic.text.FERMI_GBM_GND_POS',
                            'gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE',
                            'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD',
                            'gcn.classic.text.ICECUBE_CASCADE',
                            'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK'])

        while True:
            with open('GCNalerts.txt','a') as f:
                for message in consumer.consume(timeout=1):
                    value = message.value().decode('ASCII')
                    print(value,file=f) # print alert to file
                    print('--------',file=f) # separate alerts

                    # create a ephem object
                    alert = value.split('\n')
                    name  = re.search(r'TITLE:\s*(.*)',alert[0]).group(1) # name of object
                    type_  = "f|G"                                        # dummy type
                    ra    = re.search(r'\d+.\d+', alert[5]).group()       # find J2000 RA
                    dec   = re.search(r'\d+.\d+', alert[8]).group()       # find J2000 DEC
                    mag   = "1.0"                                         # dummy magnitude
                    in_obj = [name,type_,ra,dec,mag]
                    obj = create_ephem_object(in_obj)
                    self.ephem_objarray.append(obj)


class SourcesGUI:
    def __init__(self, master,args):
        self.master = master

        # Get initial source list
        self.observer = SPB2Obs(args)

        # Start GCN alerts check in the background
        b = threading.Thread(name='gcn_alerts', target = self.observer.gcn_alerts) # run GCN alerts in background
        b.start()

        self.sources = [
            ["Object", "Azimuth", "Altitude"],
            ["------", "-------", "--------"]
        ]
        self.master.title("SPB2Obs")

        # Create a label for the current time
        self.time_label = tk.Label(self.master, text="", font=("Arial", 12))
        self.time_label.pack(side=tk.TOP)

        # Create a label for the gps location
        gps_loc = "Observer -    Time: {0}         Latitude: {1}         Longitude: {2}         Height: {3}\n".format(*self.observer.gps_loc)
        self.gps_loc = tk.Label(self.master, text=gps_loc, font=("Arial",12))
        self.gps_loc.pack(side=tk.TOP, anchor="w")

        # Create a label for the horizon location
        horizon = "Horizon -   Limb: {0}         Upper FoV: {1}         Lower FoV: {2}\n".format(*self.observer.horizon)
        self.horizon = tk.Label(self.master, text=horizon, font=("Arial",12))
        self.horizon.pack(side=tk.TOP, anchor="w")

        # Create a label for the Sun rise and set times
        self.sun_schedule = tk.Label(self.master, text="", font=("Arial",12))
        self.sun_schedule.pack(side=tk.TOP, anchor="w")

        # Create a label for the Moon rise and set times
        self.moon_schedule = tk.Label(self.master, text="", font=("Arial",12))
        self.moon_schedule.pack(side=tk.TOP, anchor="w")

        # Create a listbox widget and populate it with the sources
        self.listbox = tk.Listbox(self.master,font="TkFixedFont")
        for sourcerow in self.sources:
            row = "{: >20} {: >20} {: >20}".format(*sourcerow)
            self.listbox.insert(tk.END, row)

        # Create a scrollbar for the listbox
        self.scrollbar = tk.Scrollbar(self.master)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.listbox.yview)

        # Pack the listbox into the GUI
        self.listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Create a button that opens another window
        self.button = tk.Button(self.master, text="Open Window", command=self.open_window)
        self.button.pack(side=tk.BOTTOM)

        # Update the time label and sources list
        self.update_time()

    def update_time(self):
        # Get the current time and format it as a string
        #new_time = time.gmtime(calendar.timegm(time.gmtime()) + 3350)
        current_time = time.strftime("%Y/%m/%d %H:%M:%S", time.gmtime()) # in UTC

        # Update the time label
        self.time_label.config(text="Current Time: " + current_time + "\n")

        # Check for an alert every 60 seconds
        if int(time.strftime("%S")) == 0:
            b = threading.Thread(name='update_gpsLoc', target = self.observer.update_gpsLoc) # run GCN alerts in background
            b.start()

        self.check_alert(current_time)

        # Check for sources in the FoV
        self.check_sources(current_time)

        # Check the rise and set times of the Sun and Moon
        self.check_sun_and_moon(current_time)

        # Update the gps location of payload label
        self.update_gpsLoc()

        # Update the horizons
        self.update_horizons()

        # Schedule the next update in 1 second
        self.master.after(1000, self.update_time)

    def update_horizons(self):
        horizon = "Horizon -   Limb: {0}         Upper FoV: {1}         Lower FoV: {2}\n".format(*self.observer.horizon)
        self.horizon.config(text=horizon)

    def update_gpsLoc(self):
        gps_loc = "Observer -    Time: {0}         Latitude: {1}         Longitude: {2}         Height: {3:.2f}\n".format(*self.observer.gps_loc)
        self.gps_loc.config(text=gps_loc)

    def check_sources(self, current_time):
        # Check if sources are in the fov or close to it
        sources = self.observer.check_fov(current_time)

        # Update sources list
        self.sources = sources
        self.listbox.delete(2,self.listbox.size()) # Clear old list
        for source in self.sources:
            row = "{: >20} {: >20} {: >20}".format(*source.split(','))
            if TESTING:
                print(row)
            self.listbox.insert(tk.END, row)

    def check_alert(self, current_time):
        # Simulate an alert
        sun_alert = None
        moon_alert = None

        sun,moon = self.observer.check_sun_and_moon(current_time)
        current_time = str(ephem.Date(current_time))
        if current_time == str(sun[0]):
            sun_alert = "Alert: Sun is rising over the limb!"
        elif current_time == str(sun[1]):
            sun_alert = "Alert: Sun is setting over the limb!"
        elif current_time == str(moon[0]):
            moon_alert = "Alert: Moon is rising over the limb!"
        elif current_time == str(moon[1]):
            moon_alert = "Alert: Moon is setting over the limb!"
        if sun_alert:
            self.make_alert_win(sun_alert)
        if moon_alert:
            self.make_alert_win(moon_alert)

    def change_color(self, object, color):
        object.config(fg= "{0}".format(color))

    def make_alert_win(self, alert_txt): # If there is an alert, display it in a new window and add a source to the list
        # Create a new window
        top = tk.Toplevel(self.master)
        top.title("Alert!")

        # Create a label to display the alert message
        alert_label = tk.Label(top, text=alert_txt, font=("Arial", 14))
        alert_label.pack(side=tk.TOP)

        # Create an acknowledgement button
        ack_button = tk.Button(top, text="Acknowledge", command=top.destroy)
        ack_button.pack(side=tk.BOTTOM)

    def check_sun_and_moon(self, current_time):
        sun,moon = self.observer.check_sun_and_moon(current_time)
        sun_str = "Sun - \t Rise: {0} \t Set: {1} \t Azi: {2} \t Alt: {3}".format(sun[0], sun[1],sun[2],sun[3])
        self.sun_schedule.config(text=sun_str)
        if SUN_FLAG:
            self.change_color(self.sun_schedule, "red")
        else:
            self.change_color(self.sun_schedule, "black")
        moon_str = "Moon - \t Rise: {0} \t Set: {1} \t Azi: {2} \t Alt: {3} \t Phase: {4:.2f}%".format(moon[0],moon[1],moon[2],moon[3],moon[4]*100)
        self.moon_schedule.config(text=moon_str)
        if MOON_FLAG:
            self.change_color(self.moon_schedule, "red")
        else:
            self.change_color(self.moon_schedule, "black")

    def open_window(self):
        # Create a new window
        top = tk.Toplevel(self.master)
        top.title("New Window")

        # Create a label to hold the image
        image_label = tk.Label(top)
        image_label.pack()

        # Load the temporary image and resize it to fit the new window
        image = Image.open("temp.png")
        image = image.resize((700, 500))
        photo = ImageTk.PhotoImage(image)

        # Update the image label with the loaded image
        image_label.config(image=photo)
        image_label.image = photo


if __name__ == '__main__':
    args = read_in_args() # read in user input arguments
    if args.test:
        TESTING = True
    else:
        TESTING = False

    # Open the GUI
    GUI(args)
