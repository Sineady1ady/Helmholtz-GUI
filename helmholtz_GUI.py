# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 09:31:56 2017

@author: Sinéad Ryan
sineadaryan@yahoo.co.uk

"""
# -*- coding: utf-8 -*-

#       NAMING CONVENTION:
#               "l_"  is a label
#               "E_" is an entry (a box that the user can write in)
#               "b_" is a button
#               "f_" is a function
#               "ND" if a variable takes ND as it's value this means
#               we have no data for it yet, hasn't been calculated

#  Reference are at the bottom of the code


#importing the packages needed:
#visa is for serial communication
#more on this at [1]
import visa

#matlplot lib is plotting package
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

#importing mathematical tools
import numpy as np
from math import modf
from scipy.optimize import curve_fit

#Tkinter is python's GUI builder.  The syntax changes a bit with python 3 so 
#watch out for that
import Tkinter as tk
from Tkinter import *
import ttk
import tkMessageBox
from PIL import Image, ImageTk



#import time to allow us to make the code pause
import time

comp_signs=["ND","ND","ND"]

#start naming some variables that will be useful later
mag_mo_1=["ND","ND","ND"]
mag_mo_2=["ND","ND","ND"]

total_mag_mo="ND"
err_tot_mag_mo="ND"
mag_rem=["ND","ND","ND","ND","ND"]
err_mag_res="ND"
x_comp="ND"
z_comp="ND"
s_comp="ND"
err=["ND","ND","ND"]
magnet_ID="000"
field_dir="ND"
dev_angle="ND"

direction=1

# default numbers for setting menu
dfts=[500,0.30,150,5,0.2,0,0,7.068e-6,1,25]
#Number of turns of coil(number here is just a guess), radius of coil (m), 
#number of readings in sine fit, interval between readings in ms
#(actual time interval varies a lot from time set here: depends greatly on
# other keithley settings), sine fit amplitude guess, sine fit phase shift guess,
# sine fit offset guess, volume of magnet, motor velocity, motor acceleration 


# The entries of this array will change if the user adjusts the settings:
settings=[500,0.30,150,5,0.2,0,0,7.068e-6,1,25]

# t_lim is used in the sine determination section.  It sets the upper and lower
# times between which the signal will be integrated to find the component sign.
# See sign determination graphs on GUI
t_lim=[0.5,1.6]
steps_per_rot=25000  #this number is set in the ST10-S motor driver configuration.
#It is the number of steps that correspond to one rotation of the motor
number_motor_rotations=5

mu_0=4*np.pi*10**(-7) #permeability of free space

LARGE_FONT= 12 #different sized fonts will be used for different labels
SMALL_FONT=10

canvas2=0  #defining some more variables for use later
canvas=[0,0,0] 
toolbar=[0,0,0]
l_amp=[0,0,0]
l_T_av=[0,0,0]
l_A_av=[0,0,0]
l_mag_mo_1=[0,0,0]
l_mag_mo_2=[0,0,0]
b_clear=[0,0,0]
taken_data=[0,0,0]
l_comp=[[0,0,0],[0,0,0],[0,0,0]]
graph_coor=[0,0,0]
sign_data=np.zeros((3,3,2,2),dtype=object)
b_sign=[0,0,0]
b_change_fit_param=[0,0,0]
l_dev_angle=[0,0]


def f_split(lst):
    return [lst[::3], lst[1::3], lst[2::3]]
#define a function "split" that splits up lists
# the data comes out of keithley is sets of three.  voltage, time, index.
#split list to get three seperate lists of each of these quantities

def f_keithley_get_ready(argu1,argu2):
# This is a function that gets the keithley taking data
#It takes two arguments.  THe first being the number of readings (argu1) and
# the second being the time interval between them in ms (argu2)  
    keithley.write("status:measurement:enable 512; *sre 1") #the commands here come from [2] and [4]
    keithley.write("sense:voltage:range 0.3") #voltage range
    keithley.write("system:azero:state off") #turning autozero off increases the rate of sampling
    keithley.write("sense:voltage:dc:nplc 1") #the number following "nplc" dictates the integration time =number/50Hz e.g. =1/50=20ms
    keithley.write("sample:count %d" % argu1) #number of readings to take
    keithley.write("trigger:source bus")
    keithley.write("trigger:delay %f" % (argu2 / 1000.0)) #time interval between readings (converted to seconds)
    keithley.write("trace:points %d" % argu1) #number of readings to take
    keithley.write("trace:feed sense1; feed:control next")
    
    keithley.write("initiate")  #more commands from [2]
    keithley.assert_trigger()
 
def f_keithley_retrieve_data():
#after the f_keithley_get_ready function is used,this function retrieves the data    
    keithley.wait_for_srq()

    kthly_data = keithley.query("trace:data?")  #gets the data from the keithley
    keithley.query("status:measurement?")
    keithley.write("trace:clear; feed:control next") 
    
    kthly_data_list=kthly_data.split(",") #seperates items with commas
    kthly_data_stripped=[y.strip('VDCuSECSRNDG#') for y in kthly_data_list] #remove some unwanted symbols in list so only important numbers are left
    del kthly_data_stripped[-1] #delete the last entry because it never gets stripped properly (for some reason)
     
    kthly_data_floats=map(float,kthly_data_stripped)#turns the stuff in the list into floating point numbers
    length=(len(kthly_data_floats)+1)/3-1 #the entry that got deleted was just an indexing number
    kthly_data_floats.append(length) #add an extra index to make up for entry that got deleted
    
    kthly_processed=(f_split(kthly_data_floats))  #use split function defined earlier to turn list into 3 lists
    return kthly_processed #return the data
    


def f_clear(self,arg1):
#this function is used to destroy labels and graphs that are no longer wanted
#after the clear button is pressed the clear button itself is destroyed
#so long as such a button exists
            taken_data[arg1]=0    
            l_amp[arg1].destroy() #<-destroying labels i.e. text boxes
            l_T_av[arg1].destroy() #arg1 tells us which page/plane the labels are
            l_A_av[arg1].destroy() #being removed from 
            l_mag_mo_1[arg1].destroy() #arg1=0 x-s plane page. arg1=1 z-s plane
            l_mag_mo_2[arg1].destroy() #page etc.
            toolbar[arg1].forget() #<-removes graph's toolbar
            canvas[arg1]._tkcanvas.pack_forget() #<-removes graph
            b_sign[arg1].destroy()
            b_change_fit_param[arg1].destroy() #<- if there is a "change parameters" button,  get rid of it
            b_clear[arg1].destroy() #<-if there is a clear button get rid of it
            for index1 in range(3):
                for index2 in range(3):
                    try:
                        l_comp[index1][index1].destroy()
                    except:
                        pass

       
def f_clear_main_page(self):
#similar to the function above.  This function also destroys unwanted labels.
#this time the labels are from the main page.         
    l_data_1.destroy()
    l_data_2.destroy()
    l_data_x.destroy()
    l_data_z.destroy()
    l_data_s.destroy()
    l_rem_1.destroy()
    l_rem_2.destroy()
    l_rem_x.destroy()
    l_rem_z.destroy()
    l_rem_s.destroy()
    l_dev_angle[0].destroy()
    l_dev_angle[1].destroy()
    b_clear2.destroy()
    
def f_sign_determination(self, rot_no, axis):
# this function is here to figure out the magnetic field direction
# on each axis
    area=np.zeros((rot_no,2))
    global sign_data
    global t_lim
    for index1 in range(0,rot_no): #index 1 tells us what rotation number we have reached
        for index2 in range(0,2): # index 2 tells us if the reading is a 0-180 rotation (positive) or a 180-360 rotation (negative)
            f_keithley_get_ready(70,4) #70 is the number of readings, 4 is the interval between them in ms
            
            time.sleep(1)              
            motor.write("VE1")  #the velocity shouldn't be too fast here
            motor.write("DI%d" %(0.5*steps_per_rot))  #do a half turn rotation
            motor.write("FL")        #go motor!
            time.sleep(1)   
            
            k_data=f_keithley_retrieve_data() #retrieve the data from the keithley 
            sign_data[axis][index1][index2][0]=k_data[0] #save the voltages and times to the sign_data array
            sign_data[axis][index1][index2][1]=k_data[1]
            times=k_data[1]   
            voltages=k_data[0]

            #Take data within the timeframe surrounding the motor rotation (t_lim[0] to t_lim[1]).  There is some
            #leeway on either side since we don't know how long the keithley will 
            #take to start triggering            
            myIndexes = [i for i,value in enumerate(times) if (value > t_lim[0] and value < t_lim[1])]
            #selected the x and y data that falls within the time frame specified above ^
            slctd_v=[voltages[i] for i in myIndexes]
            slctd_t=[times[i] for i in myIndexes]
            #find the area using the trapezium method
            area[index1][index2]=np.trapz(slctd_v,slctd_t)

    #take the median so that if we get a "weird" outlier it will be discounted           
    median_A=np.median(area[:,0])
    median_B=np.median(area[:,1])

    
    #The list comp_signs is used to keep track of what sign each component
    #takes.  comp_signs takes entries of 1 or -1 depending on which median
    #has the greater value i.e. which rotation was more positive.
    global comp_signs
    if median_A > median_B:
        comp_signs[axis]="+"

    else:
        comp_signs[axis]="-"
    
    global l_comp
    #display the signs that have been found
    if axis==0:
        l_comp[axis][0] = tk.Label(self, text=("component x is %sve" %comp_signs[axis]),font=LARGE_FONT)
        l_comp[axis][0].pack()
    elif axis==1:
        l_comp[axis][1] = tk.Label(self, text=("component z is %sve" %comp_signs[axis]),font=LARGE_FONT)
        l_comp[axis][1].pack()
    else:
        l_comp[axis][2] = tk.Label(self, text=("component s is %sve" %comp_signs[axis]),font=LARGE_FONT)
        l_comp[axis][2].pack()
        
        
                        
def f_take_data(self,arg1,controller):

#this function is the big one.  It triggers the motor and the keithley.  
#It takes, plots, fits and displays the data.  Arg1 tells us which
#plane the data is being taken on.
#  "controller" is used to navigate between different frames
    global taken_data #variable that keeps track of whether data has been taken or not
    
    if taken_data[arg1]==0: #check whether data has already been taken, if not proceed
        taken_data[arg1]=1 #indicates that data has been taken (or in this case, will be imminently)
        l_wait = tk.Label(self, text="please wait...") #<-label to tell the use to wait while data is taken
        l_wait.pack()
        l_wait.update_idletasks() #<-makes sure label is shown before python starts other tasks
        
        global settings
        global canvas
        global toolbar
        global l_amp
        global l_T_av
        global l_A_av
        global l_mag_mo_1
        global l_mag_mo_2
        global motor_warm_up_time
        
        #get time interval and number of readings from settings
        interval_in_ms=float(settings[3]) #make sure numbers are numbers (not strings)
        number_of_readings=int(settings[2])

        motor.write("FL%d" %(5*steps_per_rot)) #5 rotations of motor (between 2 & 8 rotations was found to be the optimum number)
        time.sleep(0.2)
        
        #get the keithley ready and take the data
        f_keithley_get_ready(number_of_readings, interval_in_ms)
        data=f_keithley_retrieve_data()        
        
        #plot voltage vs. time
        times=data[1]    
        voltages=data[0]
        fig1 = matplotlib.figure.Figure(figsize=(5,5), dpi=100)
        ax = fig1.add_subplot(111)
        ax.plot(times,voltages,'bo')
        ax.set_ylabel("Voltage (V)")
        ax.set_xlabel("Time (s)")

                        
        def f_display_fit(): #this plots the sine fit found and returns values for different parameters
            def func(x, p1,p2,p3,p4):
            # defines the kind of function we want the fit to use
            #i.e. a sine fn. with 4 parameters
              return p1*np.sin(p2*x+p3)+p4
          
            global settings
            sin_a_i=float(settings[4])#make sure all the intial parameter guesses are floats            
            sin_c_i=float(settings[5]) #these can be changed on the settings menu
            sin_d_i=float(settings[6])
            
            #get the initial frequency guess for sine wave from the motor speed ("settings[8]")
            sin_b_i=2*np.pi*settings[8]
    
            
            #try to do the fit.  If not send an error message
            try:
               popt, pcov = curve_fit(func, times, voltages,p0=(sin_a_i,sin_b_i,sin_c_i,sin_d_i), bounds=((0,-np.inf,-np.pi,-np.inf),(np.inf,np.inf,np.pi,np.inf)))
               perr = np.sqrt(np.diag(pcov))
        
            except:
                result=tkMessageBox.askyesno("No fit found",
                "The data could not be fitted.  Would you like to retake the data?", icon="warning")
                if result==True:  #user has the option to retake the data
                    l_wait.pack_forget() #delete "please wait" label
                    f_take_data(self,arg1,controller) #the data will be taken again
                    
                    
                else:  #if you don't want to retake the data it will be plotted without the fit
                    ax.grid(True)
                    l_wait.pack_forget() #delete "please wait" label
                    canvas[arg1] = FigureCanvasTkAgg(fig1, self) #make a canvas for putting graph on
                    canvas[arg1].show()
                    canvas[arg1].get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        
                    toolbar[arg1] = NavigationToolbar2TkAgg(canvas[arg1], self) #add a toolbar to graph (contains save option among other icons)
                    toolbar[arg1].update()
                    canvas[arg1]._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        
            t=np.arange(0.0,times[-1],0.01)  #x values for the plotting the sine fit that is obtained
            ax.plot(t, popt[0]*np.sin(popt[1]*t+popt[2])+popt[3]) #plot the sine fit
            ax.grid(True) #make the plot show
                   
            period=2*np.pi/(popt[1])#find period from angular frequency determined in fit
        
            def delta_t(x):
                return times[x+1]-times[x]  #this function gives the time that has elapsed between consecutive readings
    
            Area=np.zeros(len(voltages)-1) #set up an array of the correct length that will soon be populated

            #this bit of the code calculates the area under the curve using the trapezium method
            #only absolute areas are considered (no negative areas).  This means 2 different calculations 
            #are needed depending on the signs of the two points being used.
            for i in range(0, len(voltages)-1):
                if (voltages[i]>=0 and voltages[i+1]>=0) or (voltages[i]<=0 and voltages[i+1]<=0):        
                    Area[i]=abs((voltages[i]+voltages[i+1])/2)*delta_t(i)                                       
                else:            
                    m=abs(voltages[i+1]-voltages[i])/delta_t(i) #gradient of line
                    L=abs(voltages[i+1])/m  #distance between voltages[i+1] and intersection
                    Area[i]=0.5*L*abs(voltages[i+1])+0.5*L*abs(voltages[i])*(delta_t(i)-L) #adding two triangles together.  The first is 1/2*distance1*height1+1/2*distance2*height2
                    
            Area_tot=np.sum(Area)   #add up all the small areas calculated using the trapezium method
        
            Area_per_period=Area_tot*period/times[-1]  #average area is total area*period/total time
            
            
            whole_periods=times[-1]/period #times[-1] is the final timestamp, i.e. total time.  Divide by period to get no. of cycles (periods)
            full_half_periods=float(int(whole_periods/0.5))  #divide by two to get no. of half periods.  We want a whole number so use int. Turn answer back into float.
         
            decimal=modf(whole_periods/0.5)[0]  #takes the decimal (i.e. fractional part) that is left over from division on above line
            
            #find out what the area represented by this additional fractional part is.
            #integral of sine is cosine.  But need conditions that take care of the case
            #where we cross the x-axis. since we want the absolute area and don't want any
            #cancellation of areas above and below the x-axis.
            #The "fraction_of_full_period" obtained is an area weighted fraction
            #of the total area of a period.
            if (np.sin(0+popt[2])>=0 and np.sin(decimal*np.pi+popt[2])>=0) or (np.sin(0+popt[2])<=0 and np.sin(decimal*np.pi+popt[2])>=0):
                fraction_of_full_period=abs((np.cos(0+popt[2])-np.cos(decimal*np.pi+popt[2]))/4)
            else:
                fraction_of_full_period=abs((np.cos(0+popt[2]))+1+abs(-1-np.cos(decimal*np.pi+popt[2])))/4
         #Calculate the area per period based on the total area and the no of period (including the area weighted fraction calculated above)
            area_per_period=Area_tot/(full_half_periods/2+fraction_of_full_period)
            corres_sine_amp=(Area_per_period/4) *(2*np.pi)/period  #reverse engineer to get the sine amplitude that would give rise the area that was calculated
        
            err_area_per_period=(((interval_in_ms / 1000.0)**2)*period/12)*popt[0]*popt[1]**2  #the error in the area found based on the error formula for the  trapezium method
                                
            R_i=float(settings[1]) #radius of coils
            N_i=int(settings[0]) #number of turns
            global mu_0
            k=((float(5)/4)**(float(2)/3)) *R_i/(mu_0*N_i)  #calulation of constant k, from [3]
            global mag_mo_1
            mag_mo_1[arg1]=abs(popt[0])*k/popt[1] #get magnetic moment from sine amplitude and k value
            global mag_mo_2
            mag_mo_2[arg1]=corres_sine_amp*k/popt[1]
            global err
            err[arg1]=perr[0]*k #the error in the magnetic moment (err) scales with error in amplitude from fit
               
            canvas[arg1] = FigureCanvasTkAgg(fig1, self)  #"canvas" for plot
            canvas[arg1].show()
            canvas[arg1].get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        
            toolbar[arg1] = NavigationToolbar2TkAgg(canvas[arg1], self)  #makes a little toolbar for graph
            toolbar[arg1].update()
            canvas[arg1]._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

            #we display the different values that have been calculated as labels
            l_amp[arg1] = tk.Label(self, text="Amplitude of sine fit is:                         %.*f +/- %.*f mV" % (3,abs(popt[0])*1000,3,perr[0]*1000),font=SMALL_FONT)
            l_T_av[arg1] = tk.Label(self, text="Average period is:                                       %0.*f +/- %0.*f s" %(3, period, 3, 2*np.pi*perr[1]/(period**2)), font=SMALL_FONT)
            l_A_av[arg1] = tk.Label(self, text="Average area per period is:                   %0.*f V*s" %(6,area_per_period), font=SMALL_FONT)
            l_mag_mo_1[arg1] = tk.Label(self, text="Magnetic moment from sine fit is:          %.*f +/- %.*f J/T" %(4,mag_mo_1[arg1],4, perr[0]*k),font=LARGE_FONT)
            l_mag_mo_2[arg1] = tk.Label(self, text="Magnetic moment from area is:             %.*f +/- %.*f J/T" %(4,mag_mo_2[arg1],4,abs(err_area_per_period)),font=LARGE_FONT)

            
            l_wait.pack_forget() #gets rid of "please wait" label
            
            global b_clear
            #create a button that will clear the data taken
            b_clear[arg1]=ttk.Button(self, text="Clear",command=lambda: f_clear(self,arg1))
                   
            #layout of labels and buttons:        
            l_amp[arg1].pack()  
            l_T_av[arg1].pack()
            l_A_av[arg1].pack()
            l_mag_mo_1[arg1].pack(pady=30,padx=10)
            l_mag_mo_2[arg1].pack(pady=30,padx=10)
            b_clear[arg1].pack()
            
            b_sign[arg1]=ttk.Button(self, text="Direction Graphs", command=lambda: controller.show_frame(Graph_Page))
            b_sign[arg1].pack()
            
            
        f_display_fit()#this function creates a sine fit and plots it
        
        
        #  The following section is for determining the signs of the field
        # There are two different sets of condition depending on whether the
        # field is horizontal or vertical        
        if field_dir=="H": #case of horizontal field
            if arg1==0:
                f_sign_determination(self, 1, 2) # "1"- 1 rotation required, "2"-measurement will give us direction of "axis 2" i.e. s
            if arg1==2: #comment this bit
                f_sign_determination(self, 3, 0) #direction along x, "3"- 3 rotations required (we should do more than one measurement for
                #the off axis components since they are more tricky, "0"-will give us direction on axis 0 i.e. x
                motor.write("FL%d" %(90*steps_per_rot/360))           
                time.sleep(0.2)
                f_sign_determination(self, 3, 1) #direction along z
                motor.write("FL%d" %(270*steps_per_rot/360))

        else:  #case of vertical field
            if arg1==0:
                f_sign_determination(self, 1, 1) #direction along z, "1"- 1 rotation required, "1"-will give us direction of z-axis
            if arg1==1: # 
                f_sign_determination(self, 3, 2)  #direction along s
                motor.write("FL%d" %(90*steps_per_rot/360))
                time.sleep(0.2)
                f_sign_determination(self, 3, 0) #direction along x
                motor.write("FL%d" %(270*steps_per_rot/360))        
        
        

                        
        def f_change_sine_fit(self,arg1):
    #this function let's us refit the data with new sine parameter guesses
            global b_change_fit_param
            b_change_fit_param[arg1].pack_forget()  #gets rid of the current "change fit parameters" button
        
            def f_quit_sine():#quitting change parameters
        #If the user presses "quit" then the entry boxes and associated labels disappear
                l_sin_a.destroy()
                E_sin_a.destroy()
                l_sin_c.destroy()
                E_sin_c.destroy()
                l_sin_d.destroy()
                E_sin_d.destroy()
                b_go_back.destroy()
                b_Reset.destroy()
                l_init_guess.destroy()
                b_quit_param.destroy()
            
            def f_apply_fit_changes(self,arg1): #when the user applies the changes the new sine parameter guesses
                    global sin_d #are taken from the entry boxes
                    sin_d=E_sin_d.get()
                    global sin_c
                    sin_c=E_sin_c.get()
                    global sin_a
                    sin_a=E_sin_a.get()
        
                    f_clear(self,arg1) #the fit is cleared to make room for new fit
                    f_display_fit() #new fit is performed
                    
             
            def f_reset_sine(): #the reset button clears the entries and resets
                    global dfts #them to their default values (as set in "dfts")
                    E_sin_a.delete(0, END)
                    E_sin_a.insert(0, dfts[4])
                    E_sin_c.delete(0, END)
                    E_sin_c.insert(0, dfts[5])
                    E_sin_d.delete(0, END)
                    E_sin_d.insert(0, dfts[6])                    
                    
            global l_sin_a  #set labels, entries, buttons etc as globals so
            global l_sin_c #they can be deleted conveniently later
            global l_sin_d
            global E_sin_a
            global E_sin_c
            global E_sin_d
            global b_Reset
            global b_go_back
            global l_init_guess
                    
            b_Reset = ttk.Button(self, text="Reset",command=lambda: f_reset_sine()) #reset button
                        
            b_quit_param=ttk.Button(self, text="Quit Sine Parameters", command=lambda: f_quit_sine()) #quit button                
            l_init_guess = tk.Label(self, text="initial guess for fit: a*sin(ωt+c)+d")
             
            #create entry boxes with labels for each parameter that can be changed.
            #Default entries are from "dfts" list
            l_sin_a = tk.Label(self, text="a")
            E_sin_a=Entry(self,bd=5)
            E_sin_a.insert(0, settings[4])
                    
            l_sin_c = tk.Label(self, text="c")        
            E_sin_c=Entry(self,bd=5)        
            E_sin_c.insert(0, settings[5])
            
            l_sin_d = tk.Label(self, text="d")
            E_sin_d=Entry(self,bd=5)
            E_sin_d.insert(0, settings[6])
            
            b_go_back=ttk.Button(self, text="Apply",command=lambda: f_apply_fit_changes(self,arg1)) #button for applying changes to sine fit
             
            #layout of sine fit parameters on page:
            b_quit_param.pack(pady=5)
            b_Reset.pack(pady=5)
            l_init_guess.pack(pady=20)
            l_sin_a.pack(side=LEFT)
            E_sin_a.pack(side=LEFT)
            l_sin_c.pack(side=LEFT)
            E_sin_c.pack(side=LEFT)
            l_sin_d.pack(side=LEFT)
            E_sin_d.pack(side=LEFT)
            b_go_back.pack()
            
        global b_change_fit_param
        b_change_fit_param[arg1]=ttk.Button(self, text="Change Sine Parameters",command=lambda: f_change_sine_fit(self,arg1)) #the button for changing the sine fit parameters
        b_change_fit_param[arg1].pack()         

        
    else:
            tkMessageBox.showwarning("Clear Data Warning","Please clear old data before taking new data")
            #if data has already been taken a warning pops up to tell the user
            #to clear the old data.
            #In this case taken_data[arg1]=1


class Helmholtz_App(tk.Tk):  #now we will create the app itself with all the frames needed
#this part of the code is based on [6]
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "Helmholtz Coils Data")
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
#the names of our different frames:
        for F in (Launch_Page, Main_Page, Magnet_Identifier, mag_meas_s_z_Page, mag_meas_x_s_Page, mag_meas_x_z_Page, Settings_Page, Motor_Page, SCL_List_Page, Graph_Page):

            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(Launch_Page)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()
    
class Launch_Page(tk.Frame): #this is the first page the user will encounter
        def __init__(self, parent, controller):
            tk.Frame.__init__(self,parent)
            
            
            def f_new_series():           
                global filename
                filename=FilenameEntry.get()
            
                try:
                    test=open("%s.txt" % filename,"r")
                    test.read()
                    tkMessageBox.showwarning(
            "Filename already in use",
           "This filename has already been used. Please choose a unique filename to start a new measurement series." 
            ) 
                except:
                    try:
                        file=open("%s.txt" % filename,"a+")
                        #  "a+" means open for readings and appending.  The file with the specified filename is created if it
                        #does not already exist
                    except:
                        tkMessageBox.showwarning(
            "Filename error",
           "This filename could not be created.  Check that the filename chosen contains allowed symbols." 
            )
                    #writes headings to textfile:
                    file.write("magnet_ID\tmain_magnetic_direction\tangle_of_deviation\tmag_mo_tot\tmag_mo_err\tmag_mo_x\tmag_mo_z\tmag_mo_s\tmag_rem_tot\tmag_rem_err\tmag_rem_x\tmag_rem_z\tmag_rem_s\tVolume\tz-s_plane\tx-s_plane\tx-z_plane\n")
                    file.close()
                    f_open_communication(Settings_Page)
            
                
            def f_old_series():
                #if the user selects "old series" these steps are carried out
                global filename
                global settings
                filename=FilenameEntry.get() #get the users filename
                
                try:
                    #attempt to read the file name
                    file=open("%s.txt" % filename,"r")
                    file.read()
                    file.close
                    try:
                        #try to extract the settings from the settings file
                        file=open("%s_settings.txt" % filename,"r") #open settings file
                        set_read=file.read() #read the settings file
                        set_split=set_read.splitlines() #split the settings file by line
                        settings_tab=set_split[1] #we want the second line as this contains the numbers
                        settings=[float(x) for x in settings_tab.split()] #take these numbers as our settings
                        f_open_communication(Magnet_Identifier) #opens hardware communication then goes to Magnet_Identifier page
                    except:
                        result=tkMessageBox.askyesno("Could not find settings",
               "No settings textfile of the correct format could be found for a file of this name.  Would you like to put the settings in by hand?" 
                , icon="warning")
                        if result==True:
                            f_open_communication(Settings_Page)
                    #user has the option to put their own settings in if the settings can't be found
                    #opens hardware communication and then takes user to settings page
#                    
                except:
                   tkMessageBox.showwarning(
            "Filename does not exist",
           "This filename does not exist.  Need an old filename to continue the data series." 
            )

            
            def f_open_communication(page):
                #this function is called once the user presses "start measurement"        
                #open communication with the motor and the keithley
                global keithley
                global motor
                global settings
                rm = visa.ResourceManager()
                
                #this section establishes communication with the keithley and 
                #then the motor.  If it is not successful it throws up error messages
                try:
                    keithley = rm.open_resource('GPIB0::19::INSTR')
                    #send some commands to initilise the keithley:
                    keithley.write("*rst; status:preset; *cls")
                    try:
                        motor = rm.open_resource("COM1")
                        
                        #once communication with motor has been established
                        #we write the chosen settings to the motor
                        motor.write("JS%f" %settings[8]) #jog speed
                        motor.write("VE%f" % settings[8]) #FL speed
                        motor.write("JA%f" % settings[9]) #jog acceleration
                        motor.write("AC%f" % settings[9]) #FL acceleration
                        motor.write("DE%f" % settings[9]) #FL deceleration
                        controller.show_frame(page) #go to whatever page was specified in f_open_communication argument
                    except:
                        tkMessageBox.showwarning(
            "Error Connecting to Instruments",
           "Could not connect to the motor.  Check that ST10-S stepper motor controller is connected to PC." 
            ) 
                except:

                   tkMessageBox.showwarning(
            "Error Connecting to Instruments",
           "Could not connect to Keithley 2700.  Check that the Keithley is switched on and connected to the PC via GPIB with a valid card.  Check that the Keithley GPIB address is set to 19.  Details on how to do this are in the manual." 
            ) 
                              

            
            FilenameEntry=Entry(self,bd=5)#box for filename              
            
            #two button, old or new series
            b_new_series = ttk.Button(self, text="New series",command=f_new_series)                                         
            b_old_series = ttk.Button(self, text="Old series",command=f_old_series) 
        
        #labels:
            l_filename = tk.Label(self, text="filename",font=LARGE_FONT)
            l_dottxt = tk.Label(self, text=".txt",font=LARGE_FONT)
            
            seperator=Frame(height=2, bd=1, relief=SUNKEN) #for layout purposes
            
            #layout:
            l_filename.pack(padx=20)
            FilenameEntry.pack()
            l_dottxt.pack()
            seperator.pack(fill=X, padx=5, pady=5)
            b_new_series.pack(pady=20)
            b_old_series.pack(pady=20)
            
            
            
class Magnet_Identifier(tk.Frame): #this is the first page the user will encounter
        def __init__(self, parent, controller):
            tk.Frame.__init__(self,parent)
            
            l_ident=tk.Label(self, text="Choose a number or other identifier for your magnet",font=SMALL_FONT)
            E_ident=Entry(self,bd=5)
            l_mag_type=tk.Label(self, text="What type of field does the magnet produce?",font=SMALL_FONT)
            
                        
            def f_ok():
            #function executed once user presses "ok" on mag ident. page
                global magnet_ID
                magnet_ID=E_ident.get() #get the magnet_ID from the textbox
                
                global filename
                #we read the text in the filename to check if the magnet ID has been used before
                file=open("%s.txt" % filename,"r")
                info=file.read()
                info_list=info.split()
                
                global field_dir
                field_dir=v.get() #get the field direction chosen by the user from the radio button ("H" or "V")
                
                #if the magnet_ID string appears anywhere else in the text file we get an error message:
                if any(magnet_ID in s for s in info_list):  
                    result=tkMessageBox.askyesno("Magnet name already used",  #asks the user if they are sure they want to quit
            "The string chosen for the filename already exists somewhere in the text file.  Do you want to continue despite this duplication?", icon="warning")
                    if result==True: #if the user doesn't mind duplication then proceed to main page              
                        controller.show_frame(Main_Page)
                else:#if name appears to be unique proceed to main page
                    controller.show_frame(Main_Page)
                
            b_ok=ttk.Button(self, text="Ok",
                                command=lambda: f_ok() )
            
            v = StringVar()  #a variable used for the field direction radio button
            
            #layout
            l_ident.pack()
            E_ident.pack()
            b_ok.pack()
            l_mag_type.pack()
            
            #radio button lets user choses between horizontal and vertical magnet types
            button_rad=Radiobutton(self, text="Horizontal (s-field)", variable=v, value="H").pack()
            button_rad=Radiobutton(self, text="Vertical (z-field)", variable=v, value="V").pack()
            v.set("H") #default is horizontal
            
            

class Main_Page(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        #labels
        l_press_btn = tk.Label(self, text="Press the buttons below to make a measurement",font=SMALL_FONT)           
        l_axis = tk.Label(self, text="Measurement Plane:",font=LARGE_FONT)
        l_process = tk.Label(self, text="Process measurements:", font=LARGE_FONT)
        
        #these buttons bring up other pages for doing measurements on specific axes
        b_x = ttk.Button(self, text="z-s",
                            command=lambda: controller.show_frame(mag_meas_s_z_Page))        
        b_z = ttk.Button(self, text="x-s",
                            command=lambda: controller.show_frame(mag_meas_x_s_Page))        
        b_s = ttk.Button(self, text="x-z",
                            command=lambda: controller.show_frame(mag_meas_x_z_Page))                   
        b_motor = ttk.Button(self, text="motor menu",
                            command=lambda: controller.show_frame(Motor_Page))
        
        #photos of the magnet being measured on each of the three planes
        #for the user's reference:        
        im1=Image.open('z-s.png')
        im1=im1.resize((150,150),Image.ANTIALIAS)
        photoIm1=ImageTk.PhotoImage(im1)
        l_img = Label(self,image=photoIm1)
        l_img.image=photoIm1
        
        
        im2=Image.open('x-s.png')
        im2=im2.resize((150,150),Image.ANTIALIAS)
        photoIm2=ImageTk.PhotoImage(im2)
        l_img2 = Label(self,image=photoIm2)
        l_img2.image=photoIm2
        
        
        im3=Image.open('x-z.png')
        im3=im3.resize((150,150),Image.ANTIALIAS)
        photoIm3=ImageTk.PhotoImage(im3)
        l_img3 = Label(self,image=photoIm3)
        l_img3.image=photoIm3
        
        #photo of the magnet with its axes marked on      
        im4=Image.open('magnet_photo.png')
        im4=im4.resize((150,150),Image.ANTIALIAS)
        photoIm4=ImageTk.PhotoImage(im4)
        l_img4 = Label(self,image=photoIm4)
        l_img4.image=photoIm4
        

        
        def f_clear_variables(self):
        #does as it sounds, resets all the variables
            global mag_mo_1
            global mag_mo_2
            global mag_rem
            global settings
            global comp_signs
            global x_comp
            global y_comp
            global z_comp
            global total_mag_mo
            global err_tot_mag_mo
            global magnet_ID
            global field_dir
            global sign_data
            global dev_angle
            mag_mo_1=["ND", "ND", "ND"] 
            mag_mo_2=["ND", "ND", "ND"]
            total_mag_mo="ND"
            err_tot_mag_mo="ND"
            mag_rem=["ND","ND","ND","ND","ND"]
            comp_signs=["ND","ND","ND"]
            x_comp="ND"
            y_comp="ND"
            z_comp="ND"
            field_dir="ND"
            dev_angle="ND"
            sign_data=np.zeros((3,3,2,2),dtype=object)
            
            #clears the canvas (i.e. graphs) on each measurement plane (if they exist):
            for i in range(0,3):
                try:
                    f_clear(self,i)
                    canvas[i]._tkcanvas.pack_forget()   
                
                except:
                    pass
                
            #clears variables on main page (if they exist):       
            try:
                f_clear_main_page(self)
            except:
                pass
                        
        def f_save(self): #this function saves and exits the measurement
            global mag_mo_1
            global mag_mo_2
            global mag_rem
            global x_comp
            global y_comp
            global z_comp
            global total_mag_mo
            global err_tot_mag_mo
            global magnet_ID
            global field_dir
            global dev_angle

            global filename #get our filename
            file=open("%s.txt" % filename,"a+")
            #  "a+" means open for readings and appending.  The file with the specified filename is created if it
            #does not already exist

            file.write("\n")
            variables_list=[magnet_ID, field_dir,dev_angle, total_mag_mo,err_tot_mag_mo, "%s%s" % (comp_signs[0],x_comp), "%s%s" % (comp_signs[1],z_comp), "%s%s" % (comp_signs[2],s_comp), mag_rem[0], mag_rem[1], "%s%s" % (comp_signs[0],mag_rem[2]),"%s%s" % (comp_signs[1],mag_rem[3]),"%s%s" % (comp_signs[2],mag_rem[4]), settings[7], mag_mo_1[0],mag_mo_1[1],mag_mo_1[2]]
            for item in variables_list:
                file.write("%s\t" % item)
            #all the items in the list of variables are written to the file through this loop ^                
            file.close()  #we are done writing so we can close the file
            
            f_clear_variables(self)  #reset the variables now that the measurement is done
            try:
                f_clear_main_page(self) #clear the main page too if it is populated
            except:
                pass
            
            for i in range(3):
                try:
                    f_clear(self,i) #clear individual measurement plane pages if possible
                except:
                    pass
                        
        def f_save_and_next_magnet(self):
            f_save(self) #saves the measurements
            controller.show_frame(Magnet_Identifier) #takes us to the magnet ID page
        
        def f_save_and_exit(self):
            f_save(self) #saves the measurement
            controller.show_frame(Launch_Page) #brings us back to start page            
                                                
        def f_quit(self): #this button is for exiting without saving
            result=tkMessageBox.askyesno("Quit Measurement",  #asks the user if they are sure they want to quit
            "Are you sure you want to exit without saving?", icon="warning")
            if result==True: #if they do then proceed
                f_clear_variables   
                f_clear_main_page(self)
                for i in range(3):
                    try:
                        f_clear(self,i) #clear individual measurement plane pages if possible
                    except:
                        pass
                controller.show_frame(Launch_Page) #takes us back to start page
       
        #three different buttons with options for what to do after a measurement is taken
        b_end=ttk.Button(self, text="save & finish", command=lambda: f_save_and_exit(self))
        b_next=ttk.Button(self, text="save & next magnet", command=lambda: f_save_and_next_magnet(self))        
        b_quit=ttk.Button(self, text="quit without saving", command=lambda: f_quit(self))
        
        #layout:               
        l_img.grid(row=2, column=3,padx=10)
        l_img2.grid(row=3, column=3)
        l_img3.grid(row=4, column=3)
        l_img4.grid(row=2, column=0)
        l_press_btn.grid(row=0,column=0)
        self.grid_rowconfigure(1, minsize=50)
        self.grid_columnconfigure(1, minsize=200)
        b_x.grid(sticky="w",row=2, column=1)
        l_axis.grid(sticky="e",row=3, column=0)        
        b_z.grid(sticky="w",row=3, column=1)
        b_s.grid(sticky="w",row=4,column=1)
        self.grid_rowconfigure(5, minsize=100)
        l_process.grid(sticky="e",row=7,column=0)
        self.grid_rowconfigure(8, minsize=50)
        self.grid_rowconfigure(14, minsize=30)
        b_quit.grid(row=15, column=1)
        b_motor.grid(row=16, column=0)
        b_next.grid(row=17, column=1)
        b_end.grid(row=16, column=1)
        

        
                         
        def f_process_data(): #this button processs all our data
           
            global err
            global mag_mo_1

            try:
                global mag_mo_1
                x=float(mag_mo_1[0])
                z=float(mag_mo_1[1])
                s=float(mag_mo_1[2])

                a=x**0.5+z**0.5+s**0.5
           #if some of the mag_mo data is still set at -1 then this indicates
           #not all the data has been taken and the calculation of "a" won't be
           #possible.  We don't need "a".  It just serves as a test.
            except:
                tkMessageBox.showwarning(
                "Process Data",
               "Need valid data on all three axes before data can be processd" 
                )
                return
            
            global total_mag_mo
            #process data from all 3 axes of rotation to get the total
            total_mag_mo=(0.5*(x**2+z**2+s**2))**0.5
            global err_tot_mag_mo
            #finding the error of the total based on the error in the components:
            err_tot_mag_mo=total_mag_mo*np.sqrt((err[0]*x)**2+(err[1]*z)**2+(err[2]*s)**2)/(x**2+z**2+s**2)
         
            global l_data_1
            global l_data_2
            global l_data_x
            global l_data_z
            global l_data_s
            global l_rem_1
            global l_rem_2
            global l_rem_x
            global l_rem_z
            global l_rem_s
            global l_dev_angle
            global settings
            
            #make a label to diplay the total magnetic moment with error
            l_data_1 = tk.Label(self, text="Magnetic moment:", font=SMALL_FONT)            
            l_data_2 = tk.Label(self, text="%.2f +/- %.2f J/T" % (total_mag_mo,err_tot_mag_mo) , font=SMALL_FONT)

            global x_comp
            global z_comp
            global s_comp
                        
            #Cartesian calculations
            #These components can't physically be negative so they are set to
            #be either zero or the calculated value.  Whichever is bigger.
            x_comp=(max((0.5*(-x**2+z**2+s**2)),0))**0.5   
            z_comp=(max((0.5*(x**2-z**2+s**2)),0))**0.5       
            s_comp=(max((0.5*(x**2+z**2-s**2)),0))**0.5
            

            global mag_rem
            #the magnetic remenant is calculated from the volume and magnetic
            #moment along with the permeability of free space.  Calculation
            #is based on [3]
            mag_rem[0]=mu_0*total_mag_mo/settings[7]
            mag_rem[1]=mu_0*err_tot_mag_mo/settings[7]
            
            #use ratios to split the total magnetic remenant into components
            mag_rem[2]=mag_rem[0]*x_comp/total_mag_mo
            mag_rem[3]=mag_rem[0]*z_comp/total_mag_mo
            mag_rem[4]=mag_rem[0]*s_comp/total_mag_mo
            
            #label to display the magnetic remenant
            l_rem_1 = tk.Label(self, text="Magnetic remanence:", font=SMALL_FONT)            
            l_rem_2 = tk.Label(self, text="%.3f +/- %.3f T" % (mag_rem[0],mag_rem[1]) , font=SMALL_FONT)
            
            #comp_signs is list of values of +1 or -1 that tell us what the sign
            #of the magnetic moment components are.  Here we convert these numbers to strings
            #of either "+" of "-" so that we can display these on the labels.
            global comp_signs
           
            #these labels diplay the individual components of the 
            #magnetic moment and remenant with +/- signs
            l_data_x = tk.Label(self, text="%s %.3f e_x" % (comp_signs[0],x_comp) , font=SMALL_FONT)            
            l_data_z = tk.Label(self, text="%s %.3f e_z" % (comp_signs[1],z_comp) , font=SMALL_FONT)            
            l_data_s= tk.Label(self, text="%s  %.3f e_s" % (comp_signs[2],s_comp) , font=SMALL_FONT)
            l_rem_x = tk.Label(self, text="%s %.3f e_x" % (comp_signs[0],mag_rem[2]) , font=SMALL_FONT)            
            l_rem_z = tk.Label(self, text=" %s %.3f e_z" % (comp_signs[1],mag_rem[3]) , font=SMALL_FONT)            
            l_rem_s= tk.Label(self, text="%s %.3f e_s" % (comp_signs[2],mag_rem[4]) , font=SMALL_FONT)
            
            #This section is for calculating the angle that the magnetic vector deviates from the main axis
            global dev_angle
            global field_dir
            #there are two slightly different calculations depending on whether magnet is verical or horizontal
            if field_dir=="H":
                dev_angle=np.arctan(s/s_comp)*360/(2*np.pi)
            else:
                dev_angle=np.arctan(z/z_comp)*360/(2*np.pi)
                        
            l_dev_angle[0]=tk.Label(self, text="Deviation from main axis:" % (dev_angle) , font=SMALL_FONT)
            l_dev_angle[1]=tk.Label(self, text="%.2f°" % (dev_angle) , font=SMALL_FONT)
            
            #create a button that can clear these new labels if needs be
            global b_clear2
            b_clear2=ttk.Button(self, text="Clear", command=lambda: f_clear_main_page(self))
                       
            #layout:
            self.grid_columnconfigure(2, minsize=200)
            self.grid_rowconfigure(10, minsize=50)            
            l_data_1.grid(row=8, column=0)
            l_data_2.grid(row=8, column=1)            
            l_data_x.grid(row=9, column=0)
            l_data_z.grid(row=9,column=1)
            l_data_s.grid(row=9,column=2)            
            l_rem_1.grid(row=10, column=0)
            l_rem_2.grid(row=10, column=1)            
            l_rem_x.grid(row=11, column=0)
            l_rem_z.grid(row=11,column=1)
            l_rem_s.grid(row=11,column=2)
            l_dev_angle[0].grid(row=12, column=0)
            l_dev_angle[1].grid(row=12, column=1)            
            b_clear2.grid(row=13, column=1)
            
          
        #the button for the above function that processs data:
        b_go=ttk.Button(self, text="go", command=lambda: f_process_data())
        #layout:
        b_go.grid(sticky="w",row=7, column=1)
           
        
class mag_meas_s_z_Page(tk.Frame):
#this is the page for measuring the magnet on the z-s plane.
#There are two buttons.  One that takes us back to the home page.
#And a second button that calls the f_take_data.
#The 2nd argument of this function ("0") tells us which axis the data
#will pertain to i.e. the magnetic moment obtained will be stored in the
#zeroth entry of the mag_mo list which correspond to rotation about x.

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        l_sz = tk.Label(self, text="Magnetic moment in z-s plane", font=LARGE_FONT)
        l_sz.pack(pady=10,padx=10)

        b_home = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(Main_Page))
        b_home.pack()
        
        b_data=ttk.Button(self, text="Take Data", command=lambda: f_take_data(self,0,controller))
        
        b_data.pack()

        
class mag_meas_x_s_Page(tk.Frame):
#Similar to the class described above.  The 2nd argument of the f_take_data being "1" means that
#an extra rotation will be done in order to determine the sign of the main
#magnetic component.  This will be seen as a condition in the f_take_data
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        l_xs = tk.Label(self, text="Magnetic moment in x-s plane", font=LARGE_FONT)
        l_xs.pack(pady=10,padx=10)

        b_home = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(Main_Page))
        b_home.pack()
            
        b_data=ttk.Button(self, text="Take Data", command=lambda: f_take_data(self,1,controller))
        b_data.pack()
        
        
class mag_meas_x_z_Page(tk.Frame):
#similar to the previous 2 classes.  The f_take_data argument of 2
#will mean extra rotations are done to find the sign of x and z since this
#is the measurement where these off axis components will be most visible.        
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        l_xz = tk.Label(self, text="Magnetic moment in x-z plane", font=LARGE_FONT)
        l_xz.pack(pady=10,padx=10)

        b_home = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(Main_Page))
        b_home.pack()
       
        b_data=ttk.Button(self, text="Take Data", command=lambda: f_take_data(self,2,controller))
        b_data.pack()
        
                
class Settings_Page(tk.Frame):
#The settings page allows to user to set a number of different parameters
#for their data taking.

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        l_settings = tk.Label(self, text="Settings", font=LARGE_FONT)
        
        
        def f_apply_new_settings():
        #this function gets all the things the user put in the settings entries
        #and sets these as global variables for using late
            fail=0 #this variable will be used to keep track of whether any of
            #the setting chosen by the user fail basic criteria

            
            global dfts
            global settings
         
            if float(E_dt.get())>=0: #time interval must be positive value
                settings[3]=float(E_dt.get()) 
            else:
                fail=1
                tkMessageBox.showwarning("Time Interval Warning",
        "The minimum time interval is 0ms.  Please select a larger time interval")
            #we "get" all the other parameters from their entries:   
                        
#the number of readings should be an integer between 0 and 10,000
#Choose 10,000 as an arbitrary upper limit.
#If these conditions are not met by the user an error box pops up.

            try:
                val=int(E_read_num.get())
                if val<=10000 and val>0:
                    settings[2] = E_read_num.get()
                else:
                    fail=1
                    tkMessageBox.showwarning("Number of Readings Warnings.",
        "The number of readings should be positive and less than 10,000.")

            except:
                fail=1
                tkMessageBox.showwarning("Number of Readings Warnings.",
        "The number of readings should be an integer.")
                
#The radius of the coil should be a number and it should be greater than 0.
            try:
                val=float(E_coil_r.get())
                if val>0:
                    settings[1]= E_coil_r.get()
                else:
                    fail=1
                    tkMessageBox.showwarning("Coil Radius Warning",
        "The radius of the coil should be a positive number.")
            except:
                fail=1
                tkMessageBox.showwarning("Coil Radius Warning",
        "The radius of the coil should be a number.")

#The number of turns of the coil should be an integer greater than zero           

            try:
                val=int(E_turns.get())
                if int(E_turns.get())>0:
                    settings[0]= E_turns.get()
                else:
                    fail=1
                    tkMessageBox.showwarning("Number of Turns Warning",
        "The number of turns must be a positive integer")   
            except:
                fail=1
                tkMessageBox.showwarning("Number of Turns Warning",
        "The number of turns must be a positive integer") 

#The values chosen for the sine parameters should be numbers            
            try:
                float(E_sin_d.get())
            #must be a number
                settings[6]=E_sin_d.get()
            except:
                fail=1
                tkMessageBox.showwarning("Sine Parameter Warning",
        "The sine fit parameter 'd' must be a number.")
                
            try:
                float(E_sin_c.get())
            #must be a number
                settings[5]=E_sin_c.get()
            except:
                fail=1
                tkMessageBox.showwarning("Sine Parameter Warning",
        "The sine fit parameter 'c' must be a number.")
                        
            global sin_a
            try:
                float(E_sin_a.get())
            #must be a number
                settings[4]=E_sin_a.get()
            except:
                fail=1
                tkMessageBox.showwarning("Sine Parameter Warning",
        "The sine fit parameter 'a' must be a number.")
              
#The user can input a number or a mathematical expression for the volume.  The 
#"eval" command evaluates any expression from the user and turns it into a number.

            Volume_form=E_mag_V.get()
            try:
                settings[7]=float(eval(Volume_form))
            except:
                fail=1
                tkMessageBox.showwarning("Volume Warning",
        "Choose a valid number or formula for V, e.g. '0.0009'or '(0.033*0.033-2*0.006*0.006)*0.00695'.")

            try:
                settings[8]= float(E_vel.get())                
                if 0.0042 <= settings[8] <= 8:
                    pass

                elif settings[8]> 8:
                    fail=1
                    tkMessageBox.showwarning(
                    "Slow Down Ya Wild Thing!",
                    "Rotational speed must be between 0.0042 and 8")
                elif settings[8]<0:
                    fail=1
                    tkMessageBox.showwarning(
                    "Stay Positive",
                    "The rotational speed must take a positive value")
                else:
                    fail=1
                    tkMessageBox.showwarning(
                    "Hurry up Mr Slowpoke!",
                    "Rotational speed must be between 0.0042 and 8")
                
            except:
                fail=1
                tkMessageBox.showwarning(
                "That's not a number",
                "Please put in a valid number for the speed")
            
#the ST10-S specifies a lower velocity limit of 0.0042.  The motor seems to 
#struggle with speeds above 8m/s.  This sets the velocity range.
# The acceleration and deceleration limits are also set by the ST10-S.
#  If appropriate values aren't chosen, various warning messages will pop up.              

                                                               
            try:
                settings[9]= float(E_ac.get())
                if 0.167 <= settings[9]<= 5451.167:
                    pass

                elif settings[9]> 5451.167:
                    fail=1
                    tkMessageBox.showwarning(
                    "Slow Down Ya Wild Thing!",
                    "Acceleration must be between 0.167 and 5451.167")
                elif settings[9]<0:
                    fail=1
                    tkMessageBox.showwarning(
                    "Stay Positive",
                    "The acceleration must take a positive value")
                else:
                    fail=1
                    tkMessageBox.showwarning(
                    "Hurry up Mr Slowpoke!",
                    "Acceleration must be between 0.167 and 5451.167")
            except:
                fail=1
                tkMessageBox.showwarning(
                "That's Not a Number",
                "Please put in a valid number for the acceleration")
                            
            
            if fail==0: #if the settings are all okay  (none of them fail)
            #then we write the new settings to the motor and the settings file
                
                motor.write("JS%f" %settings[8])
                motor.write("VE%f" % settings[8])
                
                motor.write("JA%f" % settings[9])
                motor.write("AC%f" % settings[9])
                motor.write("DE%f" % settings[9])
                
                global filename
                settings_filename=("%s_settings" %filename)
    
                file=open("%s.txt" % settings_filename,"a+")
                file.write("Number_of_turns\tRadius_of_coils\tnumber_of_readings\tinterval_in_ms\tsin_a\tsin_c\tsin_d\tVolume\tmotor_velocity\tmotor_acceleration\n")
                for item in settings:
                    file.write("%s\t" % item)
                #all the items in the list of variables are written to the file through this loop ^  
                file.write("\n")              
                file.close()
                controller.show_frame(Magnet_Identifier)

      
        def f_reset():
        #This function is there if the user wants to reset all the settings
        #entries to their default values.  Each entry is deleted and replaced
        #with the corresponding value from the defaults list (named dfts).
            E_turns.delete(0, END)
            E_turns.insert(0, dfts[0])
            E_coil_r.delete(0, END)
            E_coil_r.insert(0, dfts[1])
            E_read_num.delete(0, END)
            E_read_num.insert(0, dfts[2])
            E_dt.delete(0, END)
            E_dt.insert(0, dfts[3])
            E_sin_a.delete(0, END)
            E_sin_a.insert(0, dfts[4])
            E_sin_c.delete(0, END)
            E_sin_c.insert(0, dfts[5])
            E_sin_d.delete(0, END)
            E_sin_d.insert(0, dfts[6])
            E_mag_V.delete(0, END)
            E_mag_V.insert(0, dfts[7])
            E_vel.delete(0, END)
            E_vel.insert(0, dfts[8])
            E_ac.delete(0, END)
            E_ac.insert(0, dfts[9])

            #The default velocities and accelerations are written to the motor
            vel= float(E_vel.get())
            motor.write("JS%f" % vel)
            motor.write("VE%f" % vel)
            
            ac= float(E_ac.get())
            motor.write("JA%f" % ac)
            motor.write("AC%f" % ac)
            motor.write("DE%f" % ac)
            
            
        b_Reset = ttk.Button(self, text="Reset",command=f_reset) #reset button
        b_Enter = ttk.Button(self, text="Done",command=f_apply_new_settings) #new setting button
        
        #labels and entries for each variable:
        #(default values are inserted into entries)
        l_turns = tk.Label(self, text="No. of turns")        
        E_turns=Entry(self,bd=5)        
        E_turns.insert(0, dfts[0]) 
        
        l_coil_r = tk.Label(self, text="Radius of coils (m)")        
        E_coil_r=Entry(self,bd=5)        
        E_coil_r.insert(0, dfts[1])
        
        l_read_num = tk.Label(self, text="No. of data points")        
        E_read_num=Entry(self,bd=5)        
        E_read_num.insert(0, dfts[2])
        
        l_dt = tk.Label(self, text="Time between data points (ms)")        
        E_dt=Entry(self,bd=5)        
        E_dt.insert(0, dfts[3])
        
        l_mag_V = tk.Label(self, text="Volume of magnet (m^3) \n Type a number or formula")        
        E_mag_V=Entry(self,bd=5)        
        E_mag_V.insert(0, dfts[7])
        
        l_vel = tk.Label(self, text="Speed rev/s (range: 0.0042 to 8)")        
        E_vel=Entry(self,bd=5)
        E_vel.insert(0, dfts[8])
               
        l_ac = tk.Label(self, text="Rate of Acceleration/Deceleration rev/s^2")
        E_ac=Entry(self,bd=5)       
        E_ac.insert(0, dfts[9])
                
        l_init_guess = tk.Label(self, text="initial guess for sine fit parameters: a*sin(ωt+c)+d")        
        
        l_sin_a = tk.Label(self, text="a")
        E_sin_a=Entry(self,bd=5)        
        E_sin_a.insert(0, dfts[4])
               
        l_sin_c = tk.Label(self, text="c")
        E_sin_c=Entry(self,bd=5)        
        E_sin_c.insert(0, dfts[5])
        
        l_sin_d = tk.Label(self, text="d")
        E_sin_d=Entry(self,bd=5)        
        E_sin_d.insert(0, dfts[6])
        
        #layout:
        l_settings.grid(row=0,column=0)
        b_Reset.grid(row=1, column=2) 
        b_Enter.grid(row=6, column=4)
        self.grid_rowconfigure(1, minsize=50)
        self.grid_columnconfigure(3, minsize=50)
        l_turns.grid(sticky="w",row=2,column=1, padx=10)
        E_turns.grid(sticky="w",row=2,column=2, padx=10)
        l_coil_r.grid(sticky="w",row=3,column=1, padx=10)
        E_coil_r.grid(sticky="w",row=3,column=2, padx=10)
        l_read_num.grid(sticky="w",row=4,column=1, padx=10)
        E_read_num.grid(sticky="w",row=4,column=2, padx=10)
        l_dt.grid(sticky="w",row=5,column=1, padx=10)
        E_dt.grid(sticky="w",row=5,column=2, padx=10)
        l_mag_V.grid(sticky="w",row=6,column=1, padx=10)
        E_mag_V.grid(sticky="w",row=6,column=2, padx=10)
        l_vel.grid(sticky="w",row=7,column=1, padx=10)
        E_vel.grid(sticky="w",row=7,column=2, padx=10)
        l_ac.grid(sticky="w",row=8,column=1, padx=10)
        E_ac.grid(sticky="w",row=8,column=2, padx=10)
        l_init_guess.grid(sticky="w",row=9,column=1, padx=10)
        l_sin_a.grid(sticky="w",row=10, column=1, padx=10)       
        E_sin_a.grid(sticky="w",row=10, column=2, padx=10)
        l_sin_c.grid(sticky="w",row=11, column=1, padx=10)         
        E_sin_c.grid(sticky="w",row=11, column=2, padx=10)
        l_sin_d.grid(sticky="w",row=12, column=1, padx=10)        
        E_sin_d.grid(sticky="w",row=12, column=2, padx=10)
        
        
class Motor_Page(tk.Frame):
#this is a special page dedicated to the motor.  From here you can change the
#motor setting or make the motor move how you like.  There is also an option to send
#the motor SCL commands. A list of these can be found in [5].

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        l_mtr_cntrl = tk.Label(self, text="Motor Control", font=LARGE_FONT)
        l_mtr_cntrl.pack(pady=10,padx=10)

        b_home = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(Main_Page))
        b_home.pack(pady=5)
        
        global settings
        
        def f_motor_reset():
        #This function resets all the variables
            
            E_rotate.delete(0,END)
            E_command.delete(0,END)
            
            #We can keep track of the jogging direction of the motor
            #The direction button ("b_direc") is reset to clockwise:
            b_direc.config(text='Clockwise')            
            motor.write("CS%f" % settings[8])
            
        #reset button:    
        b_reset = ttk.Button(self, text="Reset",command=f_motor_reset)
        l_direc = tk.Label(self, text="Current direction of rotation:")
        
#When the "b_direc" button is pressed, the button toggles between different
# labels of either "clockwise" or "counterclockwise". The new direction is
# written to the motor using "CS" for jog mode. The variable "direction"
# is used to keep track of the current direction of the motor for FL mode.        
        def f_toggle_dir():
            global direction
            if b_direc.config('text')[-1] == 'Clockwise':
                b_direc.config(text='Counterclockwise')
                motor.write("CS-%f" %settings[8])
                
                direction=-1

                
            else:
                b_direc.config(text='Clockwise')                
                motor.write("CS%f" %settings[8])
                direction=1

        
        b_direc = ttk.Button(self,text="Clockwise", width=18, command=f_toggle_dir)
        l_OnOff = tk.Label(self, text="Motor Jog On/Off")
             
        def f_motor_on_off():
            #this button sets the motor jogging
            global direction
            if b_motor.config('text')[-1] == 'Jog':
                b_motor.config(text='Off') #"On" button turns into "off" button
                motor.write("CJ") #commence jogging
                
                #the motor always starts out clockwise. Must make sure directin
                #button and "direction" variable reflect this:
                b_direc.config(text='Clockwise')
                direction=1                


            else:
                #if the motor is already on we turn it off and set the button
                # to read "on" (since we now have the option to turn it on).
                motor.write("SJ")
                b_motor.config(text='Jog')
                
  
        b_motor = ttk.Button(self,text="Jog", width=10, command=f_motor_on_off)                
        l_rotate = tk.Label(self, text="Rotate motor x degrees:")        
        E_rotate=Entry(self,bd=5)
        
        #layout:
        b_reset.pack(pady=5)
        l_direc.pack()
        b_direc.pack()
        l_OnOff.pack()
        b_motor.pack(pady=5)
        l_rotate.pack()
        E_rotate.pack()
        
        def f_rotate_motor():
            #this function allows us to rotate the motor a specified number of degrees
            global direction
            global settings
            rotation=float(E_rotate.get()) #get the angle specified by the user
            motor.write("VE0.3")  #set a low speed for this rotation
            motor.write("FL%d" %(direction*rotation*steps_per_rot/360)) #rotate the specified amount
            motor.write("VE%f" % settings[8]) #return to the rotation speed chosen by the user
            
        def f_cancel_rotate():
            motor.write("SK") #this command kills the motor and can be sent to cancel a motor rotation
               
        b_rotate = ttk.Button(self, text="Rotate",command=f_rotate_motor)        
        b_cancel = ttk.Button(self, text="Cancel Rotation",command=f_cancel_rotate)


        
        
        
        
        def f_send_SCL_command():
        # This is a handy function that will send whatever text is written
        # as an SCL command to the motor.
            command=E_command.get() #gets the command from the entry
            motor.write(command) #sends it to the motor
            E_command.delete(0, END) #deletes the command from the entry box
            
        #label, entry and button for the SCl command sending function        
        l_command= tk.Label(self, text="Send other SCL command.  Type command in box:")        
        E_command=Entry(self,bd=5)   
        b_command = ttk.Button(self, text="Send",command=f_send_SCL_command)        
        
        #if the user want more info on SCL commands this button will take them to an info page
        b_info = ttk.Button(self, text="More info on SCL commands",command=lambda: controller.show_frame(SCL_List_Page))
        
        #layout
        b_rotate.pack(pady=5)
        b_cancel.pack(pady=5)
        l_command.pack()
        E_command.pack()
        b_command.pack(pady=5)
        b_info.pack(pady=5)
        
class SCL_List_Page(tk.Frame):
# this page is there to provide information on SCL commands that can be sent to the motor
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        b_back = ttk.Button(self, text="Back",command=lambda: controller.show_frame(Motor_Page))
        b_back.grid(row=0,column=0)

# There is a list of commands and a list of explanations        
        commands=["CJ:","SJ:","SK:","JS:","VE:","JA:","AC:","DE:","CS:","DI:","FL:"]
        descriptions=["Commence jogging","Stop jogging", "Kill",
                      "Jog speed (sets the speed for a CJ command)",
                      "Velocity (will set the speed for a FP, FL or FS move until another VE command is sent)",
                      "Jog acceleration (sets the acceleration and deceleration for CJ)",
                      "Sets the acceleration for FL, FS, FP, SH",
                      "Sets the deceleration for FL, FS, FP, SH",
                      "Change jog speed (on the fly)",
                      "Sets requests of move distance in steps.  Sign indicates move direction.  \n Steps per full rotation are defined using the ST configurator. \n They are currently set at the arbitrary value of 25000 ",
                      "Feed to length.  Move distance and direction of last DI command"]
        #these lists are displayed in a grid
        l_command=[0]*(len(commands))     
        for index in range(0,len(commands)):
            l_command=tk.Label(self, text=commands[index]).grid(sticky="n", row=index+1, column=0)
        for index in range(0,len(descriptions)):
            l_description=tk.Label(self, text=descriptions[index],anchor="w").grid(sticky="w",row=index+1, column=1)
            
class Graph_Page(tk.Frame):
#This page displays the graph that come from trying to find the directions
#of the different components    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        global shown
        shown=0
        
        def f_back(self):
            #when the user hits back, the canvas (i.e.graphs) and labels are destroyed
            global canvas2
            global shown
            global l_comp_signs
            global l_t_lim
            
            try:
                canvas2._tkcanvas.pack_forget()
                l_comp_signs.destroy()
                l_t_lim.destroy()
            except:
                pass
            controller.show_frame(Main_Page) #we return to the main page
            shown=0
            
        b_back = ttk.Button(self, text="Back",command=lambda: f_back(self))
        b_back.pack()
        
        
        
        def f_show_graphs(argument):
            #this function displays the direction determination graphs
            global shown
            global comp_signs
            
            if argument==1:
                tkMessageBox.showwarning("Already Showing Graphs",
                "The graphs are already on show.  To produce new graphs, go and take some more data.", icon="warning")
            else:   
                shown=1
                global sign_data
                fig = Figure(figsize=(5,8), dpi=100)
                plt=[0,0,0]
                axis_names=["x","z","s"]
                
                for i in range(0,3):
                #this loop runs through all three axes
                #it plots the 0-180 sign data in red and 180-360 sign data in blue
                #each one is plotted twice- once as a series of points (e.g. 'or')and once as a connecting line (e.g. '-r')
                    plt[i] = fig.add_subplot(3,1,1+i)
                    plt[i].plot(sign_data[i][0][0][1],sign_data[i][0][0][0],'-r', label="0-180")
                    plt[i].plot(sign_data[i][0][0][1],sign_data[i][0][0][0],'or')
                    plt[i].plot(sign_data[i][0][1][1],sign_data[i][0][1][0],'-b', label="180-360")
                    plt[i].plot(sign_data[i][0][1][1],sign_data[i][0][1][0],'ob')
                    plt[i].legend(loc="upper left")
                    plt[i].set_ylabel("Voltage (V)")
                    plt[i].set_xlabel("Time (s)")
                    plt[i].set_title("%s-axis" %axis_names[i])                        
               
                #if data from further measurements exists on any of the axes this wil also be plotted in this loop
                #(i.e. off-axis case where sign measurements are done 3 times):
                for index2 in range(1,3):
                    for index1 in range(0,3):
                        try:
                            plt[index1].plot(sign_data[index1][index2][0][1],sign_data[index1][index2][0][0],'-r')
                            plt[index1].plot(sign_data[index1][index2][0][1],sign_data[index1][index2][0][0],'or')
                            plt[index1].plot(sign_data[index1][index2][1][1],sign_data[index1][index2][1][0],'-b')
                            plt[index1].plot(sign_data[index1][index2][1][1],sign_data[index1][index2][1][0],'ob')
                            
                        except:
                            pass
                fig.tight_layout() #this line stops the graph labels from overlapping with eachother
                global canvas2
                try:
                    canvas2._tkcanvas.pack_forget()  #if the canvas already exist get rid of it
                except:
                    pass
                canvas2 = FigureCanvasTkAgg(fig, self) #make a new canvas and put the figure on it
                canvas2.show()
                canvas2.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
                
                #print some labels to supplement the graphs:
                global l_comp_signs
                l_comp_signs=tk.Label(self,text=("x=%s, z=%s, s=%s" %(comp_signs[0],comp_signs[1],comp_signs[2])), font=LARGE_FONT)
                l_comp_signs.pack()
                global l_t_lim
                l_t_lim=tk.Label(self,text=("Integration occur between %ss &%ss" %(t_lim[0],t_lim[1])), font=LARGE_FONT)
                l_t_lim.pack()
   
                
        #button that user presses to get the graphs to be displayed:
        b_show = ttk.Button(self, text="Show Graphs",command=lambda: f_show_graphs(shown))
        b_show.pack()  
        

        
app = Helmholtz_App()
app.mainloop()



##~~~~~~~~~~~~~~~~~~~~~~~REFERENCES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PyVisa documentation https://pyvisa.readthedocs.io/en/stable/
# 2. PyVisa Keithley communication example   https://pyvisa.readthedocs.io/en/stable/example.html
# 3. A NEW HELMHOLTZ COIL PERMANENT MAGNET
#    MEASUREMENT SYSTEM*
#    Joseph Z. Xu and Isaac Vasserman, ANL, Argonne, IL 60439, U.S.A
#    https://accelconf.web.cern.ch/accelconf/icalepcs2011/papers/wepkn015.pdf
# 4. Keithley 2700 communication commands from manual 
#    http://www.ee.bgu.ac.il/~acl/Equip/2700_900_01fnl.pdf
# 5. Motor control commands 
#    https://www.applied-motion.com/sites/default/files/Host-Command-Reference_920-0002N.pdf
# 6. Instructions of plotting graphs in tkinter
#    https://pythonprogramming.net/how-to-embed-matplotlib-graph-tkinter-gui/
