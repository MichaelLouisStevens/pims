import numpy as np
import math as math

from bokeh.plotting import figure
from bokeh.models.layouts import HBox
from bokeh.layouts import column
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider
from bokeh.io import curdoc


# thermal speed equivalent for a temperature given in eV 
def eV2kms(eV, Z, A):
    return np.sqrt(eV*(1/A)*(170.**2)/150.)

# speed of ion accelerated through potential a V [Volts]
def ZV2kms(V, Z, A):
    return np.sqrt(V*(Z/A)*(170.**2)/150.)

# dimensionless version of the COLD FC response on a given
# voltage window. Full response would be f * (Ze*n*Area*u)
# This is invalid for drift speed on order of or less than thermal speed
def cold_flux_function(Vlo, Vhi, ukms, TeV, Z, A):
    wkms = eV2kms(TeV, Z, A)
    vlokms = ZV2kms(Vlo, Z, A)
    vhikms = ZV2kms(Vhi, Z, A)
    ulim = (1./wkms)*(vhikms-ukms)
    llim = (1./wkms)*(vlokms-ukms)
    Qv = wkms/(ukms*np.sqrt(np.pi)) # inverse Mach number
    return 0.5*(math.erf(ulim) - math.erf(llim) -
        Qv*(np.exp(-ulim**2)-np.exp(-llim**2)))
#I THINK THERE'S AN ERROR HERE. FLUX SHOULD DECREASE WITH TEMPERATURE


# declare some constants
e0_pC = 1.6021892e-7 # electron charge in picoCoulombs
cm2km = 1e5

# declare starting values for all adjustable parameters
nwin_0 = 64     # number of window (DEFAULT)
Vlow_0 = 10.    # lowest voltage (DEFAULT)
Vhigh_0 = 8000.   # highest voltage  (DEFAULT)
ukms_0 = 100.   # corotation speed (DEFAULT)
TeV_0 = 100.    # isothermal temperature (DEFAULT)
necm3_0 = 25.   # total electron density (DEFAULT)
Xp_0 = 0.1      # relative abundance of protons (DEFAULT)    
OSratio_0 = 2.  # Oxygen-sulfur ratio (DEFAULT)
ZO2_0 = .75      # O++/O ratio (DEFAULT)
ZS2_0 = 0.6      # S++/S ratio (DEFAULT)
ZS3_0 = 0.2      # S+++/S ratio (DEFAULT)
Aeffcm2_0 = 2.5 # effective area of sensor (DEFAULT) 

# calculate the number densities for each ion species
nH_0 = Xp_0*necm3_0   
nO_0 = necm3_0 * (1.-Xp_0)/(1.+(1./OSratio_0))
nS_0 = nO_0/OSratio_0
nO2_0 = nO_0*ZO2_0
nO1_0 = nO_0 - nO2_0
nS2_0 = nS_0*ZS2_0
nS3_0 = nS_0*ZS3_0
nS1_0 = nS_0 - (nS2_0+nS3_0)

# define the voltage windows (DEFAULT)    
#Vwindows_0 = np.linspace(Vlow_0, Vhigh_0, num=nwin_0)
Vwindows_0 = np.logspace(np.log10(Vlow_0), np.log10(Vhigh_0), num=nwin_0)


# define a flux scale factor eA for conveience
# in units of pico-Coulombs/second per (km/s) per cm^-3
Aeffsw = Aeffcm2_0*cm2km # effective area in solar wind units (cm^3/km)
eA = e0_pC*Aeffsw

# calculate the currents for each species Z*eA*n*u*f() where
# f() is the flux function in the cold plasma approx
IH = np.zeros(nwin_0)
IO1 = np.zeros(nwin_0) 
IO2 = np.zeros(nwin_0)
IS1 = np.zeros(nwin_0)
IS2 = np.zeros(nwin_0)
IS3 = np.zeros(nwin_0)
for window in range(0,nwin_0-2):
    IH[window] = 1. * eA * nH_0 * ukms_0 * \
        cold_flux_function(Vwindows_0[window], Vwindows_0[window+1], 
                           ukms_0, TeV_0, 1., 1.)
    IO1[window] = 1.*eA*nO1_0*ukms_0* \
        cold_flux_function(Vwindows_0[window], Vwindows_0[window+1], 
                           ukms_0, TeV_0, 1., 16.)
    IO2[window] = 2.*eA*nO2_0*ukms_0* \
        cold_flux_function(Vwindows_0[window], Vwindows_0[window+1], 
                           ukms_0, TeV_0, 2., 16.)
    IS1[window] = 1.*eA*nS1_0*ukms_0* \
        cold_flux_function(Vwindows_0[window], Vwindows_0[window+1], 
                           ukms_0, TeV_0, 1., 32.)
    IS2[window] = 2.*eA*nS2_0*ukms_0* \
        cold_flux_function(Vwindows_0[window], Vwindows_0[window+1], 
                           ukms_0, TeV_0, 2., 32.)
    IS3[window] = 2.*eA*nS3_0*ukms_0* \
        cold_flux_function(Vwindows_0[window], Vwindows_0[window+1], 
                           ukms_0, TeV_0, 3., 32.)
 
# Set up plot object variables                         
x=Vwindows_0
yH=IH
yO1=IO1   
yO2 = IO2   
yS1 = IS1    
yS2 = IS2    
yS3 = IS3  
ytot = IH+IO1+IO2+IS1+IS2+IS3  
sourceH  = ColumnDataSource(data=dict(x=x, y=yH))
sourceO1  = ColumnDataSource(data=dict(x=x, y=yO1))
sourceO2  = ColumnDataSource(data=dict(x=x, y=yO2))
sourceS1  = ColumnDataSource(data=dict(x=x, y=yS1))
sourceS2  = ColumnDataSource(data=dict(x=x, y=yS2))
sourceS3  = ColumnDataSource(data=dict(x=x, y=yS3))
sourcetot = ColumnDataSource(data=dict(x=x, y=ytot))

# Set up plot
plot = figure(title="Faraday Cup Cold Plasma I(V)",
              #tools="crosshair,pan,reset,resize,save,wheel_zoom",
              x_range=[Vlow_0, Vhigh_0], y_range=[0.1, 100], x_axis_type="log",
              y_axis_label="current [pA]", x_axis_label = "modulator voltage [V]",
              y_axis_type = "log")

plot.line('x', 'y', source=sourcetot, line_width=1, line_alpha=1, line_color='black')
plot.line('x', 'y', source=sourceH, line_width=3, line_alpha=0.6, line_color='LimeGreen')
plot.line('x', 'y', source=sourceO1, line_width=3, line_alpha=0.6, line_color='red')
plot.line('x', 'y', source=sourceO2, line_width=3, line_alpha=0.6, line_color='DarkRed')
plot.line('x', 'y', source=sourceS1, line_width=3, line_alpha=0.6, line_color='cyan')
plot.line('x', 'y', source=sourceS2, line_width=3, line_alpha=0.6, line_color='blue')
plot.line('x', 'y', source=sourceS3, line_width=3, line_alpha=0.6, line_color='MidnightBlue')


# Set up widgets
nwin_w = Slider(title="number of [logarithmic] steps", value=nwin_0, start=0, end=128, step=4)
Vlow_w = Slider(title="minimum voltage", value=Vlow_0, start=0, end=100, step=1)
Vhigh_w = Slider(title="maximum voltage", value=Vhigh_0, start=100, end=10000, step=5)
ukms_w = Slider(title="drift speed [kms]", value=ukms_0, start=0, end=400, step=5)
TeV_w = Slider(title="ion temperature [eV]", value=TeV_0, start=1, end=1000, step=5)
necm3_w = Slider(title="e- density [cm^-3]", value=necm3_0, start=0.1, end=400, step=5)
Xp_w = Slider(title="H fraction (by number)", value=Xp_0, 
              start=0, end=1, step=0.01)   
OSratio_w = Slider(title="O:S abundance ratio", value=OSratio_0, start=0, end=5, step=0.1)
ZO2_w = Slider(title="O2+/O fraction", value=ZO2_0, start=0, end=1, step=0.01)
ZS2_w = Slider(title="S2+/S fraction", value=ZS2_0, start=0, end=1, step=0.01)
ZS3_w = Slider(title="S3+/S fraction", value=ZS3_0, start=0, end=1, step=0.01)
Aeffcm2_w = Slider(title="effective area [cm^2]", value=Aeffcm2_0, start=0, end=20, step=0.25)

# widget-driven updater
def update_data(attrname, old, new):
    e0_pC = 1.6021892e-7 # electron charge in picoCoulombs
    cm2km = 1e5
    
    # Get the current slider values
    nwin = nwin_w.value
    Vlow = Vlow_w.value
    Vhigh = Vhigh_w.value
    ukms = ukms_w.value
    TeV = TeV_w.value
    necm3 = necm3_w.value
    Xp = Xp_w.value
    OSratio = OSratio_w.value
    ZO2 = ZO2_w.value
    ZS2 = ZS2_w.value
    ZS3 = ZS3_w.value
    Aeffcm2 = Aeffcm2_w.value
    
    # calculate the number densities for each ion species
    nH = Xp*necm3   
    nO = necm3 * (1.-Xp)/(1.+(1./OSratio))
    nS = nO/OSratio
    nO2 = nO*ZO2
    nO1 = nO - nO2
    nS2 = nS*ZS2
    nS3 = nS*ZS3
    nS1 = nS - (nS2+nS3)
    
    # define the voltage windows (defaults)    
#    Vwindows = np.linspace(Vlow, Vhigh, num=nwin)
    Vwindows = np.logspace(np.log10(Vlow), np.log10(Vhigh), num=nwin)

    # define the flux scale factor eA
    # in units of pico-Coulombs/second per (km/s) per cm^-3
    Aeffsw = Aeffcm2*cm2km # effective area in solar wind units (cm^3/km)
    eA = e0_pC*Aeffsw
        
    # calculate the currents for each species Z*eA*n*u*f
    IH = np.zeros(nwin)
    IO1 = np.zeros(nwin) 
    IO2 = np.zeros(nwin)
    IS1 = np.zeros(nwin)
    IS2 = np.zeros(nwin)
    IS3 = np.zeros(nwin)
    for window in range(0,nwin-2):
        IH[window] = 1. * eA * nH * ukms * \
            cold_flux_function(Vwindows[window], Vwindows[window+1], 
                               ukms, TeV, 1., 1.)
        IO1[window] = 1.*eA*nO1*ukms* \
            cold_flux_function(Vwindows[window], Vwindows[window+1], 
                               ukms, TeV, 1., 16.)
        IO2[window] = 2.*eA*nO2*ukms* \
            cold_flux_function(Vwindows[window], Vwindows[window+1], 
                               ukms, TeV, 2., 16.)
        IS1[window] = 1.*eA*nS1*ukms* \
            cold_flux_function(Vwindows[window], Vwindows[window+1], 
                               ukms, TeV, 1., 32.)
        IS2[window] = 2.*eA*nS2*ukms* \
            cold_flux_function(Vwindows[window], Vwindows[window+1], 
                               ukms, TeV, 2., 32.)
        IS3[window] = 2.*eA*nS3*ukms* \
            cold_flux_function(Vwindows[window], Vwindows[window+1], 
                              ukms, TeV, 3., 32.)
                              
    x = Vwindows                                                                                   
    yH = IH
    yO1 = IO1    
    yO2 = IO2   
    yS1 = IS1    
    yS2 = IS2    
    yS3 = IS3  
    ytot = IH+IO1+IO2+IS1+IS2+IS3  
    sourceH.data = dict(x=x, y=yH)
    sourceO1.data = dict(x=x, y=yO1)
    sourceO2.data = dict(x=x, y=yO2)
    sourceS1.data = dict(x=x, y=yS1)
    sourceS2.data = dict(x=x, y=yS2)
    sourceS3.data = dict(x=x, y=yS3)
    sourcetot.data = dict(x=x, y=ytot)

# check all sliders for actions. When actions occur, call the updater
for this_slider in [nwin_w, Vlow_w, Vhigh_w, ukms_w, TeV_w, necm3_w, Xp_w, OSratio_w, 
          ZO2_w, ZS2_w, ZS3_w, Aeffcm2_w]:
    this_slider.on_change('value',  update_data)


# Set up layouts and add to document
inputs = column(nwin_w, Vlow_w, Vhigh_w, ukms_w, TeV_w, necm3_w, 
                            Xp_w, OSratio_w, ZO2_w, ZS2_w, ZS3_w, Aeffcm2_w)

curdoc().add_root(HBox(children=[inputs, plot], width=1200))