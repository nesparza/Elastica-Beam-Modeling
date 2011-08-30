# -*- coding: utf-8 -*-
"""
Created on Mon Jun 06 23:11:36 2011

@author: Noe
"""

from __future__ import print_function
from scipy.signal import cspline1d, cspline1d_eval
import matplotlib.pyplot as mpl
import numpy as np
import sys
sys.path.append('../../')
import taperedElasticaBeam as teb

##################################################################
# creates a beam and sets its properties (angle not included)
def initBeam(L,Lt):
    E = 1.75e6           # elastic modulus of beam (Pa)
    t = 20e-6            # dimension of beam in bending direction (m)
    w = 20e-6            # width of beam (m)
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

def dispSearch(beam,goal,tol,right=0,left=0):
    # solve psi(s) function
    beam.calculateSlopeFunction()

    # convert psi(s) to cartesian
    beam.calculateDisplacements()

    if (beam.yTipDisplacement() <= goal+tol and beam.yTipDisplacement() >= goal-tol):
        return
    else:
        if beam.yTipDisplacement() > goal:
            # found upper value
            right = beam.shearLoad        
        else:
            # found lower value
            left = beam.shearLoad
        if right == 0:
            beam.setShearLoad(2*beam.shearLoad)
        else:
            beam.setShearLoad((right-left)/2+left)
        dispSearch(beam,goal,tol,right,left)


#fig, ax = initFig()
#fig2, ax2 = initFig()

force = 1e-6                            #test load in micronewtons
scale = 1e6
gamma = 50e-3
aspectRatios = np.linspace(0.05,0.50,num = 30, endpoint = True)
#taperRatios = np.array([0.25,0.3])
freeLenRatio = np.ones(aspectRatios.shape)

for i,aspect in enumerate(aspectRatios):
    print('Running aspect = %f' % aspect)
    taperLen = 20/aspect/scale
    if i>0:
        lower = freeLenRatio[i-1] - freeLenRatio[i-1]*0.4
        upper = freeLenRatio[i-1] + freeLenRatio[i-1]*0.4
        if lower < 0.2:
            lower = 0.2
        if upper > 0.99:
            upper = 0.99
    else:
        lower = 0.2
        upper = 0.99
    
#    if taperLen*scale > 200:
#        beamHt = 200/scale
#    else:
#        beamHt = taperLen
    beamHt = taperLen
    arrBeamLens = np.linspace(lower*beamHt,upper*beamHt,num = 10, endpoint = True)
    arrStrainEnergy = np.zeros(arrBeamLens.shape)
    arrSurfEnergy = np.zeros(arrBeamLens.shape)
    beams = []
    
    for length in arrBeamLens:
        beams.append(initBeam(L=length,Lt = taperLen))

    spacing = (20.0,)    
    spcFreeLen = []
    guessShearLoad = 1e-6
    for spc in spacing:
        for j,beam in enumerate(beams):
            beam.setEndAngle(0)
            beam.constrainEndAngle()
            beam.setShearLoad(guessShearLoad)
            dispSearch(beam=beam,goal=spc/2/scale,tol=1e-8,right=0,left=0)
            guessShearLoad = beam.shearLoad
            # plot beam position
    #        beam.plotBeam(ax,legendLabel = str(beam.L))
    #        beam.plotSlope(ax2,legendLabel = str(beam.L))
            arrStrainEnergy[j] = beam.calculateStrainEnergy()
            arrSurfEnergy[j] = -gamma * beam.w *(beam.Lt - beam.L)
        
        interpLens = np.linspace(arrBeamLens[0],arrBeamLens[-1],num=100,endpoint=True)
        csFit = cspline1d((arrStrainEnergy+arrSurfEnergy))
        interpTotalEnergy = cspline1d_eval(csFit,interpLens,dx=(arrBeamLens[1]-arrBeamLens[0]), x0 = arrBeamLens[0])
        
    #    mpl.figure()
    #    mpl.hold(True)
    #    mpl.plot(arrBeamLens*scale,arrStrainEnergy*scale)
    #    mpl.plot(arrBeamLens*scale,arrSurfEnergy*scale)
    #    mpl.plot(interpLens*scale,interpTotalEnergy*scale,arrBeamLens*scale,(arrStrainEnergy+arrSurfEnergy)*scale,'o')
        
        resultLen = interpLens[interpTotalEnergy.argmin()]
        freeLenRatio[i]= resultLen/beamHt
        print('Free Ratio: %f' % (resultLen/beamHt))
        print('equil length: %f microns of %f\n' % (resultLen*scale,beamHt*scale))
    spcFreeLen.append(freeLenRatio)
mpl.figure()
mpl.hold()
for freeLen in spcFreeLen:
    mpl.plot(aspectRatios,freeLen)
mpl.xlabel('$\\alpha$')
mpl.ylabel('a\L')