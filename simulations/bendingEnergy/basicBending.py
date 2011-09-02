# -*- coding: utf-8 -*-
"""
Created on Mon Jun 06 23:11:36 2011

@author: Noe
"""

from __future__ import print_function
from scipy.signal import cspline1d, cspline1d_eval
from scipy.optimize import fsolve
import matplotlib.pyplot as mpl
import numpy as np
import sys, pickle
sys.path.append('../../')
import taperedElasticaBeam as teb

##################################################################
# creates a beam and sets its properties (angle not included)
def initBeam(L,Lt,w,t):
    E = 1.75e6           # elastic modulus of beam (Pa)
#    t = 20e-6            # dimension of beam in bending direction (m)
#    w = 20e-6            # width of beam (m)
#    L = 75e-5            # length of beam (m)
#    Lt = 80e-6           # length of taper (m) at Lt beam has zero thickness
    beam=teb.taperedElasticaBeam()
    beam.setBeamDimensions(E=E, t=t, w=w, L=L, Lt=Lt)

    return beam
##################################################################
# creates the figure and axes, returns a tuple (fig, ax)
def initFig():
    figHeight = 6
    figWidth = 6
    fig = mpl.figure()
    fig.set_figheight(figHeight)
    fig.set_figwidth(figWidth)

    ax = fig.add_axes([.15,.3, .7,.6])
    ax.grid()

    return (fig, ax)
##################################################################
# takes the specified test load normal to the surface of the
# cantilever and translates it to a shear and normal load with
# respect to the beam.
# 0 is completely in the normal direction
def setLoads(beam, angle, force):
    normal = force * np.cos(angle)
    shear  = force * np.sin(angle)
    beam.setAxialLoad(axialLoad=normal)
    beam.setShearLoad(shearLoad=shear)

def evalTipError(load,beam,goal):
    beam.setEndAngle(0)
    beam.constrainEndAngle()
    beam.setShearLoad(load)
    beam.calculateSlopeFunction()
    beam.calculateDisplacements()
    return beam.yTipDisplacement() - goal
    
def dispSearch(beam,initLoad,goal,tol,right=0,left=0):
    
    myLoad = fsolve(evalTipError,initLoad,args=(beam,goal))
    # apply fsolve solution
    beam.setShearLoad(myLoad)
    beam.calculateSlopeFunction()
    beam.calculateDisplacements()
#    print ('Desired: %s, Actual: %s, returnForce: %s, actual load: %s' % (goal*scale,beam.yTipDisplacement()*scale,myLoad,beam.shearLoad))
 

def evalBeams(arrBeams,guessShearLoad,spc,debugFlag = False):
    
    lastBeam = arrBeams[-1]
    arrStrainEnergy = np.zeros(np.size(arrBeams))
    arrSurfEnergy = np.zeros(np.size(arrBeams))
    lastIndex = np.size(arrBeams)
    maxSurfEnergy = -gamma * lastBeam.w*(lastBeam.Lt - lastBeam.L)
    arrResults = np.array([])
    
    for j,beam in enumerate(arrBeams):
        if debugFlag:
            print('Beam Length: %f of %f'%(beam.L*scale,beam.Lt*scale))
        
        dispSearch(beam=beam,initLoad = guessShearLoad,goal=spc/2/scale,tol=1e-7,right=0,left=0)
        guessShearLoad = beam.shearLoad
        
        if debugFlag:
            print('Solved Beam -- TipDisp: %s Goal: %s Force: %s' % (beam.yTipDisplacement()*scale,spc/2,beam.shearLoad))
        
        arrStrainEnergy[j] = beam.calculateStrainEnergy()
        arrSurfEnergy[j] = -gamma * beam.w *(beam.Lt - beam.L)
        if arrStrainEnergy[j] >= np.abs(maxSurfEnergy): 
            # since there is more bending energy than surface energy stop computing 
            print('Super stiff beam')
            lastIndex = j
            break
    
    if lastIndex > 0:   # This ensures that we have more than one data point before trying to interpolate
        interpLens = np.linspace(arrBeamLens[0],arrBeamLens[lastIndex-1],num=100,endpoint=True) # Generate x values for which to interpolate
        csFit = cspline1d((arrStrainEnergy[0:lastIndex]+arrSurfEnergy[0:lastIndex]))    # Generate cubic spline fit to the sub dataset
        interpTotalEnergy = cspline1d_eval(csFit,interpLens,dx=(arrBeamLens[1]-arrBeamLens[0]), x0 = arrBeamLens[0])    # Generate the interpolated values from the fit and x points
        finalLen = interpLens[interpTotalEnergy.argmin()]   # find the minimum of the energy balance and grab index to choose the appropriate length
        
        if debugFlag:
            print('beamLens shape: %s arrStrain: %s'%(arrBeamLens[0:lastIndex].shape,arrStrainEnergy[0:lastIndex].shape))
            mpl.figure()
            mpl.hold(True)
            mpl.plot(arrBeamLens[0:lastIndex]*scale,arrStrainEnergy[0:lastIndex]*scale)
            mpl.plot(arrBeamLens[0:lastIndex]*scale,arrSurfEnergy[0:lastIndex]*scale)
            mpl.plot(interpLens*scale,interpTotalEnergy*scale,arrBeamLens[0:lastIndex]*scale,(arrStrainEnergy+arrSurfEnergy)[0:lastIndex]*scale,'o')
        arrResults = np.array([arrBeamLens[0:lastIndex],arrStrainEnergy[0:lastIndex]])
    else:   # since there is only one datapoint then use that as the value
        finalLen = arrBeamLens[lastIndex]
        arrResults = np.array([arrBeamLens[lastIndex],arrStrainEnergy[lastIndex]])
    
    
    
    return (finalLen,arrResults)


#fig, ax = initFig()
#fig2, ax2 = initFig()

force = 1e-9                            #test load in micronewtons
scale = 1e6
gamma = 200e-3
baseWidth = 40
t = 100
spacing = (30.0,40.0,50.0,60.0)
aspectRatios = np.linspace(0.05,0.50,num = 30, endpoint = True)
freeLenRatio = np.ones(aspectRatios.shape)

#aspectRatios = np.array([0.05])

outputFile = 'clumpingdata.pkl'
output = open(outputFile, 'wb')

data = {'keyField': 'Aspect',
        'parameters': [('gamma',gamma),('baseWidth',baseWidth),('thickness',t)]}

for i,aspect in enumerate(aspectRatios):
    print('Running aspect = %f' % aspect)
    taperLen = baseWidth/aspect/scale

    lower = 0.1
    upper = 0.99
    
    beamHt = taperLen
    arrBeamLens = np.linspace(lower*beamHt,upper*beamHt,num = 60, endpoint = True)
    arrBeamLens = arrBeamLens[::-1]
    
    beams = []
    
    for length in arrBeamLens:
        beams.append(initBeam(L=length,Lt = taperLen,w=baseWidth/scale,t=t/scale))
    
    spcFreeLen = []
    
    subDict = {'keyField':'Spacing'}
    for spc in spacing:
        arrTemp = np.array([])
        resultLen,arrTemp = evalBeams(beams,force,spc,debugFlag = False)    
        freeLenRatio[i]= resultLen/beamHt
        print('Free Ratio: %f' % (resultLen/beamHt))
        print('equil length: %f microns of %f\n' % (resultLen*scale,beamHt*scale))
        spcFreeLen.append(freeLenRatio)
        subDict[spc] = arrTemp
    data[aspect] = subDict
    
mpl.figure()
mpl.hold(True)
for freeLen in spcFreeLen:
    mpl.plot(aspectRatios,freeLen)
mpl.xlabel('$\\alpha$')
mpl.ylabel('a\L')

pickle.dump(data,output)
output.close()

#    if i>0:
#        lower = freeLenRatio[i-1] - freeLenRatio[i-1]*0.4
#        upper = freeLenRatio[i-1] + freeLenRatio[i-1]*0.4
#        if lower < 0.2:
#            lower = 0.2
#        if upper > 0.99:
#            upper = 0.99
#    else:
#        lower = 0.2
#        upper = 0.99

#def dispSearch(beam,initLoad,goal,tol,right=0,left=0):
#    # set boundary conditions    
#    beam.setEndAngle(0)
#    beam.constrainEndAngle()
#    beam.setShearLoad(initLoad)
#    # solve psi(s) function
#    beam.calculateSlopeFunction()
#    # convert psi(s) to cartesian
#    beam.calculateDisplacements()
#    
#    if (beam.yTipDisplacement() <= goal+tol and beam.yTipDisplacement() >= goal-tol):
#        return
#    else:
#
#        if beam.yTipDisplacement() > goal:
#            # found upper value
#            right = beam.shearLoad        
#        else:
#            # found lower value
#            left = beam.shearLoad
#        
#        if right < left:
#            initLoad = 2*beam.shearLoad
#        else:
#            initLoad = (right-left)/2+left
#            
#        return dispSearch(beam,initLoad,goal,tol,right,left)