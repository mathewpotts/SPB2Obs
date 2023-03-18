#!/usr/bin/python3

# Authors: Mathew Potts, Jordan Bogdan, Andrew Wang

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
    parser = argparse.ArgumentParser(description = 'SPB2Obs shows Objects of Interest (OoI) in the FoV of the CT telescope displaying the azimuth and altitude of those objects. SPB2Obs incorporates live alerts of Gamma-Ray Bursts (GRBs) from the General Coordinates Network (GCN).')
    parser.add_argument('-obj', metavar='objFile',action='store',help='Path to file containing OoI.')
    parser.add_argument('-loc', metavar='locFile',action='store',help='Path to file containing the current GPS location of SPB2 and time.')
    parser.add_argument('-test', action='store_true',default=False,help='Debugging option.')
    args = parser.parse_args()
    return args

def GUI(args):
    root = tk.Tk()
    root.geometry("800x200") # set default window size
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
        self.default_horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.obs.elevation))))
        self.upperfov = self.default_horizon + (math.pi/180)*2.5
        self.lowerfov = self.default_horizon + (math.pi/180)*-5.0

        # Init masks
        self.s = ephem.Sun() #make Sun
        self.m = ephem.Moon() #make Moon

    def check_fov(self, utctime):
        sources = []
        self.obs.date = utctime
        for i in range(len(self.ephem_objarray)):
            sources.append(self.objs_in_fov(utctime, self.ephem_objarray[i]))

        # remove objects that will never be in FoV
        filter_list = []
        for i,s in enumerate(sources):
            if s is not None:
                filter_list.append(i)
        self.ephem_objarray = [obj for i,obj in enumerate(self.ephem_objarray) if i in filter_list]
        
        sources = [s for s in sources if s is not None] # filter list
        return sources

    def objs_in_fov(self, utctime, ephem_obj):
        try:
            rise = self.obs.next_rising(ephem_obj)
            sett = self.obs.next_setting(ephem_obj)
            if rise > sett: # if source is setting first
                self.obs.horizon = self.upperfov # set to upper fov
                self.obs.horizon = self.lowerfov # set to lower fov
                self.obs.horizon = self.default_horizon #reset to horizon
                az,alt,mask = self.masks(ephem_obj, utctime)
                if TESTING:
                    gui_str = "{0},{1},{2}".format(ephem_obj.name,az,alt)
                    print(gui_str)
                    return gui_str
                elif alt <= self.upperfov and alt >= self.lowerfov and mask:
                    gui_str = "{0},{1},{2}".format(ephem_obj.name,az,alt)
                    return gui_str
            else: # if source is rising firsts
                self.obs.horizon = self.upperfov # set to upper fov
                self.obs.horizon = self.lowerfov # set to lower fov
                self.obs.horizon = self.default_horizon #reset to horizon
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
        except ephem.NeverUpError:
            print("Warning: Object of interest {0} is never up, and is out of the FoV.".format(ephem_obj))

    def masks(self, ephem_obj, utctime):
        self.obs.date = utctime # make sure obs is at current time
        ephem_obj.compute(self.obs)
        self.s.compute(self.obs)
        altSun = float(repr(self.s.alt)) * 180 / math.pi
        self.m.compute(self.obs)
        phaseMoon = float(repr(self.m.moon_phase))
        altMoon = float(repr(self.m.alt)) * 180 / math.pi

        sunmask = altSun <= (((self.default_horizon) * 180 / math.pi) - 15.0) #set condition that sun is
        moonmask = phaseMoon <= 30.0
        #print(sunmask,moonmask)
        if not sunmask:
            print("Sun is up.")
        if not moonmask:
            print("Moon phase is higher.")
        return ephem_obj.az, ephem_obj.alt, sunmask*moonmask

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
            NGC.append(ephem.readdb(ngc_xephem_format)) #create ncg object
        return NGC

    def gcn_alerts(self):
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
                    value = message.value()
                    print(value,file=f)

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
        current_time = time.strftime("%Y/%m/%d %H:%M:%S",time.gmtime()) # in UTC

        # Update the time label
        self.time_label.config(text="Current Time: " + current_time)

        # Check for an alert every 10 seconds
        #if int(time.strftime("%S")) % 10 == 0:
        #    self.check_alert()

        # Check for sources in the FoV
        self.check_sources(current_time)

        # Schedule the next update in 1 second
        self.master.after(1000, self.update_time)

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
    if args.test:
        TESTING = True
    else:
        TESTING = False

    # Open the GUI
    GUI(args)
