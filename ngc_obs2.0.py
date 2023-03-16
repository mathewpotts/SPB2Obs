#!/usr/bin/python3

#Import libs
import ephem
import math
import numpy as np
import matplotlib.pyplot as plt
import argparse
from gcn_kafka import Consumer
import tkinter as tk
from PIL import Image, ImageTk
import time
import threading

# Read in args
def read_in_args():
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-obj', metavar='objFile',action='store',help='')
    parser.add_argument('-loc', metavar='locFile',action='store',help='')
    args = parser.parse_args()
    return args

def gcn_alerts(): # wait for a gcn alerts        # Connect as a consumer.
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
            for message in consumer.consume(timeout=1):
                value = message.value()
                print(value)

def GUI(sources):
    root = tk.Tk()
    app = SourcesGUI(root,list(sources))
    root.mainloop()

# Define class to observer
class SPB2Obs:
    def __init__(self,args):
        self.days = 3 #integration period in days
        self.step_hours = 0.08333333333 #step time in hours
        self.in_obj = self.read_file(args.obj) # read infile objects
        self.ephem_objarray = self.ephem_object_array() # Create Ephem objects

        # Set up observer
        self.obs = ephem.Observer() #make observer
        self.gpsLoc = self.read_file(args.loc) #read in gps location
        self.set_observer(self.gpsLoc) #set observer based on location

        # Compute objects at init location/altitude (without masks)
        for i in range(len(self.in_obj)):
            self.ephem_objarray[i].compute(self.obs)
        #print(self.obs.next_setting(self.ephem_objarray[0]),self.obs.next_rising(self.ephem_objarray[0]))

        # Init masks
        self.s = ephem.Sun() #make Sun
        self.m = ephem.Moon() #make Moon
        #self.mask_moon_and_sun() #generate masks
        self.new()

    def new(self):
        for i in range(len(self.in_obj)):
            try:
                rise = self.obs.next_rising(self.ephem_objarray[i])
                sett = self.obs.next_setting(self.ephem_objarray[i])
                if rise > sett: # if source is setting
                    self.obs.horizon = float(self.obs.horizon) + (math.pi/180)*2.5 # set to upper fov
                    upperfov = self.obs.next_setting(self.ephem_objarray[i])
                    print(self.obs)
                    print(i,upperfov)
                    self.obs.horizon = float(self.obs.horizon) + (math.pi/180)*-7.5 # set to lower fov
                    lowerfov = self.obs.next_setting(self.ephem_objarray[i])
                    print(self.obs)
                    print(i,lowerfov)
                    self.obs.horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.obs.elevation)))) #reset to horizon
                    self.fill_arrays(i,upperfov,lowerfov)
                else: # if source is rising
                    self.obs.horizon = float(self.obs.horizon) + (math.pi/180)*-5.0 # set to upper fov
                    lowerfov = self.obs.next_rising(self.ephem_objarray[i])
                    print(self.obs)
                    print(i,lowerfov)
                    self.obs.horizon = float(self.obs.horizon) + (math.pi/180)*7.5 # set to lower fov
                    upperfov = self.obs.next_rising(self.ephem_objarray[i])
                    print(self.obs)
                    print(i,upperfov)
                    self.obs.horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.obs.elevation)))) #reset to horizon
                    self.fill_arrays(i,lowerfov,upperfov)
            except ephem.AlwaysUpError:
                print("Warning: Object of interest {0} is always up always up, and is out of the FoV.".format(self.in_obj[i].split(',')[0]))
            except ephem.NeverUpError:
                print("Warning: Object of interest {0} is never up, and is out of the FoV.".format(self.in_obj[i].split(',')[0]))

    def fill_arrays(self, i, time1, time2):
        dt = (time2-time1)*86400 # crossing time in seconds
        self.obs.date = time1
        print(self.obs)
        
        # Init various arrays
        self.azis = np.zeros([len(self.in_obj),int(dt)]) #azimuths for observables
        self.alts = np.zeros([len(self.in_obj),int(dt)]) #altitudes for observables
        self.altmask = np.array(np.zeros(len(self.in_obj)), dtype=bool)
        self.altSun = np.zeros(int(dt))
        self.phaseMoon = np.zeros(int(dt))
        self.altMoon = np.zeros(int(dt))
        self.times = np.empty(int(dt)) #times
        self.dates = [] #dates
        
        for x in range(int(dt)):
            #for i in range(len(self.in_obj)):
            self.ephem_objarray[i].compute(self.obs) #compute observables for obs location time and date
            self.azis[i][x] = float(self.ephem_objarray[i].az) * 180. / math.pi #add azimuths, altitudes to arrays
            self.alts[i][x] = float(self.ephem_objarray[i].alt) * 180. / math.pi
            #print(self.obs.date,i,self.azis[i][x],self.alts[i][x])
            self.s.compute(self.obs)
            self.altSun[x] = float(repr(self.s.alt)) * 180 / math.pi
            self.m.compute(self.obs)
            self.phaseMoon[x] = float(repr(self.m.moon_phase))
            self.altMoon[x] = float(repr(self.m.alt)) * 180 / math.pi
            self.times[x] = self.obs.date.tuple()[2]*24 + self.obs.date.tuple()[3] + self.obs.date.tuple()[4]/60 + self.obs.date.tuple()[5]/3600 #add time in hours to array
            self.dates.append(str(self.obs.date))
            self.obs.date = self.obs.date + ephem.second #update date

        self.altmask[i] = True
        self.sunmask = self.altSun <= ((float(self.obs.horizon) * 180 / math.pi) - 15.0) #set condition that sun is
        self.moonmask = self.phaseMoon <= 30.0

    def set_observer(self, locArray):
        ls = locArray[0].split(',')
        self.obs.date = ls[0] #UTC time
        self.obs.lat = ls[1]
        self.obs.lon = ls[2]
        self.obs.elevation = float(ls[3])
        self.obs.horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(ls[3])))) #horizon calculated from elevation of obs
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
            ngc_eq = ephem.Equatorial(ephem.Galactic(in_obj[2], in_obj[3], epoch = ephem.J2000), epoch = ephem.J2000) #converting to equatorial
            ngc_xephem_format = in_obj[0] + ',' + in_obj[1] + ',' + str(ngc_eq.ra) + ',' + str(ngc_eq.dec) + ',' + in_obj[4]#supplying fixed coord data in xephem format (https://xephem.github.io/XEphem/Site/help/xephem.html#mozTocId800642)
            #print(ngc_xephem_format)
            NGC.append(ephem.readdb(ngc_xephem_format)) #create ncg object
        return NGC
    
    def mask_moon_and_sun(self):
        for x in range(int(self.days * 24 / self.step_hours)):
            for i in range(len(self.in_obj)):
                self.ephem_objarray[i].compute(self.obs) #compute observables for obs location time and date
                self.azis[i][x] = float(self.ephem_objarray[i].az) * 180. / math.pi #add azimuths, altitudes to arrays
                self.alts[i][x] = float(self.ephem_objarray[i].alt) * 180. / math.pi
                #print(self.obs.date,i,self.azis[i][x],self.alts[i][x])
            self.s.compute(self.obs)
            self.altSun[x] = float(repr(self.s.alt)) * 180 / math.pi
            self.m.compute(self.obs)
            self.phaseMoon[x] = float(repr(self.m.moon_phase))
            self.altMoon[x] = float(repr(self.m.alt)) * 180 / math.pi
            self.times[x] = self.obs.date.tuple()[2]*24 + self.obs.date.tuple()[3] + self.obs.date.tuple()[4]/60 + self.obs.date.tuple()[5]/3600 #add time in hours to array
            self.dates.append(str(self.obs.date))
            self.obs.date = self.obs.date + self.step_hours * ephem.hour #update date

        self.altmask = (self.alts >= (float(self.obs.horizon) * 180 / math.pi)-5.0) & (self.alts <= (float(self.obs.horizon) * 180 / math.pi)+2.5)
        #print(self.altmask)
        self.sunmask = self.altSun <= ((float(self.obs.horizon) * 180 / math.pi) - 15.0) #set condition that sun is 
        self.moonmask = self.phaseMoon <= 30.0
    
    def print_files(self):
        ls = self.gpsLoc[0].split(',')
        for i in range(len(self.in_obj)):
            in_obj = self.in_obj[i].split(',')
            f = open('{}.txt'.format(in_obj[0]),'w')
            f.write('Start Date (UTC): {}\n'.format(ls[0]))
            f.write('Observer Location (Lat, Long, Elevation): ({lat},{long},{elevation})\n'.format(lat=ls[1],long=ls[2],elevation=float(ls[3])))
            f.write('Format:\n')
            f.write('ALTITUDE (DEG), AZIMUTH (DEG), DATE (Y/M/D), TIME (UTC), CROSSED HORIZON?\n')
            if len((self.altmask[i]&self.sunmask&self.moonmask).nonzero()[0]) != 0:
                previndex = (self.altmask[i]&self.sunmask&self.moonmask).nonzero()[0][0]
            else:
                if not self.altmask[i]:
                    f.write('Warning: Not within the elevations of interest.\n')
                if not self.sunmask[i]:
                    f.write('Warning: Sun is above horizon!\n')
                if not self.moonmask[i]:
                    f.write('Warning: Moon is bright\n')
            for j in (self.altmask[i]&self.sunmask&self.moonmask).nonzero()[0]:
                if j == previndex + 1:
                    f.write(str(self.alts[i][j]) + ',' + str(self.azis[i][j]) + ',' + self.dates[j] + ' N\n')
                else:
                    f.write(str(self.alts[i][j]) + ',' + str(self.azis[i][j]) + ',' + self.dates[j] + ' Y\n')
                previndex = j
            f.close()
                            
        for i in range(len(self.in_obj)):
            in_obs = self.in_obj[i].split(',')
            plt.plot(self.times[self.altmask[i]&self.sunmask&self.moonmask],self.azis[i][self.altmask[i]&self.sunmask&self.moonmask],'.',label=in_obs[0])
            plt.xlabel('Time [hrs]')
            plt.ylabel('Azimuth [$^\circ$]')
            plt.legend()
            plt.ylim((0,360))
            plt.savefig('temp.png',dpi=1200)


class SourcesGUI:
    def __init__(self, master, sourcelist):
        self.master = master
        self.sources = sourcelist
        self.master.title("Sources")
        
        # Create a label for the current time
        self.time_label = tk.Label(self.master, text="", font=("Arial", 12))
        self.time_label.pack(side=tk.TOP)

        # Create a listbox widget and populate it with the sources
        self.listbox = tk.Listbox(self.master)
        for source in self.sources:
            self.listbox.insert(tk.END, source)

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
        current_time = time.strftime("%H:%M:%S")

        # Update the time label
        self.time_label.config(text="Current Time: " + current_time)

        # Check for an alert every 10 seconds
        if int(time.strftime("%S")) % 10 == 0:
            self.check_alert()

        # Schedule the next update in 1 second
        self.master.after(1000, self.update_time)

    def check_alert(self):
        # Simulate an alert
        alert = None
        if int(time.strftime("%M")) % 2 == 1:
            alert = "Alert: odd minute!"
            
        # If there is an alert, display it in a new window and add a source to the list
        if alert:
            # Create a new window
            top = tk.Toplevel(self.master)
            top.title("Alert!")

            # Create a label to display the alert message
            alert_label = tk.Label(top, text=alert, font=("Arial", 14))
            alert_label.pack(side=tk.TOP)
            
            # Create an acknowledgement button
            ack_button = tk.Button(top, text="Acknowledge", command=top.destroy)
            ack_button.pack(side=tk.BOTTOM)
            
            # Add a new source to the list and update the listbox
            new_source = "New Source"
            self.sources.append(new_source)
            self.listbox.insert(tk.END, new_source)

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
    observer = SPB2Obs(args)
    #observer.print_files()
    sources = ["Source 1","Source 2","Source 3"]

    #b = threading.Thread(name='gcn_alerts', target = gcn_alerts) # run GCN alerts in background
    f = threading.Thread(name='GUI', target = GUI, args = {tuple(sources)}) # run GUI in foreground
    #b.start()
    f.start()
