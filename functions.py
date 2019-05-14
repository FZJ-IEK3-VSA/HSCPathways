# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 12:52:20 2016

@author: m.reuss
"""
from mpl_toolkits.mplot3d.art3d import Text3D
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import xlsxwriter
import os

plt.style.use(os.path.dirname(os.path.realpath(__file__))+"\matplotlibrc.mplstyle")


from mpl_toolkits.mplot3d import proj3d
 
#%%
def orthogonalProj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    # -0.0001 added for numerical stability as suggested in:
    # http://stackoverflow.com/questions/23840756
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])




#%%    
def trisurfplotMin(demandArray,distArray,z_array,min_array,line_array,dfHSC,zmax=8):
    '''
    Function for plotting the costs as surface based on minimal costs + lines
    
    input parameter:
    X: Demand Array
    Y: Distance Array
    z_array: Array of solutions for different HSC (cutted)
    min_array: Array of minimal Costs
    line_Array: Array of intersection line
    dfHSC: DataFrame for the HSC
    '''
    demandArray=demandArray/1000
    x,y=np.meshgrid(demandArray,distArray)
    x=x.ravel()
    y=y.ravel()
    xmin=np.amin(min_array)
    xmax=np.amax(min_array)
    fig = plt.figure()
    #fig.autolayout=True
    fig.patch.set_facecolor('white')
    #ax = plt.subplot2grid((1,6), (0, 0), colspan=6, projection='3d')
    ax = fig.gca(projection='3d')
    ax.set_alpha(None)
    line=line_array[:,:].ravel()
    for i in range(len(dfHSC)): 
        z=z_array[:,:,i].ravel()
        trip=np.array([x,y,z])
        tripx=(trip.T[~np.isnan(trip.T).any(axis=1)]).T
        if len(tripx.T)>3:
            xx=tripx[0,:]
            yy=tripx[1,:]
            zz=tripx[2,:]
            xm=np.average(xx)
            ym=np.average(yy)
            zm=np.average(zz)
            ax.text(xm,ym,zm,dfHSC['General'][i],size=12,zorder=len(dfHSC)+i,ha='center',bbox=dict(facecolor=(1,1,1), alpha=0.7, edgecolor='none'))            
    
    
    trip=np.array([x,y,line])
    tripx=(trip.T[~np.isnan(trip.T).any(axis=1)]).T
    xx=tripx[0,:]
    yy=tripx[1,:]
    zz=tripx[2,:]
    
    surf=ax.plot_surface(demandArray,distArray,min_array,rstride=1,cstride=1,linewidth=0,cmap=cm.jet, antialiased=False,vmin=xmin,vmax=xmax)
    ax.plot_trisurf(x,y,line,color='black',linewidth=0,zorder=i+1)
    surf.set_alpha(None)
        
    
    ax.set_alpha(None)
    ax.set_xlabel('Demand of Hydrogen in t per day', labelpad=10, size=14)
    ax.set_ylabel('Distance in km', labelpad=10, size=14)
    ax.set_zlabel('Hydrogen Costs in €/kg')
    ax.set_zlim(0,xmax)
    proj3d.persp_transformation = orthogonalProj
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d t/day'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d km'))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%d €'))
    fig.suptitle('Wasserstoffgesamtkosten nach Elektrolyseproduktion')

    
    #Insert Colorbar
    cbaxes = fig.add_axes([0.05, 0.25, 0.02, 0.5]) 
    cb = plt.colorbar(surf, cax = cbaxes) 
    #fig.colorbar(surf,shrink=0.5, aspect=15,pad=0.2, location='left')#orientation='vertical'
#    plt.tight_layout()

    plt.show()

#%%
def plotFixDemand(data,distArray,demFix,demIndex,names):
    
    ''' plotting different technologies against a rising demand
    data=result array
    distFix=fixed demand
    demandArray=demand Array
    distIndex=Index of the fixed demand inside the results array
    names=Names of Technologies (for the legend)
    '''
    
    fig=plt.figure()
    fig.suptitle('Hydrogen Refueling Cost depending on distance')
    lineplot=plt.subplot2grid((1,6), (0, 0), colspan=4)
    for i in range (len(names)):
        lineplot.plot(distArray, data[:,demIndex,i], label=names[i])
    
    lineplot.set_xlabel('Distance between source and demand', labelpad=5)
    lineplot.set_ylabel('Cost of hydrogen at refueling station', labelpad=5)
    lineplot.xaxis.set_major_formatter(FormatStrFormatter('%d km'))
    lineplot.yaxis.set_major_formatter(FormatStrFormatter('%.2f €/kg'))
    
    lineplot.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
    
    info='Demand: ' + str(demFix/1000) + ' t/day'
    fig.text(0.2,0.8, info,
        horizontalalignment='left',
        verticalalignment='top',
        backgroundcolor=(220/255,220/255,220/255))

    plt.tight_layout()
    plt.show()
 
#%%   
def plotFixDistance(data,distFix,demandArray,distIndex,names):
    
    ''' plotting different technologies against a rising demand
    data=result array
    distFix=fixed demand
    demandArray=demand Array
    distIndex=Index of the fixed demand inside the results array
    names=Names of Technologies (for the legend)
    '''
    demandArray=demandArray.T/1000
    fig=plt.figure()
    fig.suptitle('Hydrogen Refueling Cost depending on distance')
    lineplot=plt.subplot2grid((1,6), (0, 0), colspan=4)
    for i in range (len(names)):
        lineplot.plot(demandArray, data[distIndex,:,i], label=names[i])
    
    lineplot.set_xlabel('Demand of Hydrogen', labelpad=5)
    lineplot.set_ylabel('Cost of hydrogen at refueling station', labelpad=5)
    lineplot.xaxis.set_major_formatter(FormatStrFormatter('%d t/day'))
    lineplot.yaxis.set_major_formatter(FormatStrFormatter('%.2f €/kg'))
    
    lineplot.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)    
    
    info='Distance: ' + str(distFix) + ' km'
   
    fig.text(0.5,0.8, info, horizontalalignment='left', verticalalignment='top',bbox=dict(facecolor= (220/255,220/255,220/255), alpha=0.5, edgecolor='none'))
      
    plt.ylim(6, 12)
    
#    plt.tight_layout()  
    plt.show()
    


#%%
# tornado chart example
def tornado(data, namesX, namesY ,x,dist,dem,inputVariables, kind=0,xmin=5,xmax=12):
    '''
    kind=0 --> tornado chart with fixed variables
    kind=1 --> tornado chart with fixed HSC
    '''
    xmin=np.min(data)
    xmax=np.max(data)
    
    if kind==0:
        names=namesX
        title=namesY[x]
        bases = np.array(data[x,:,1])
        base=bases[0]
        
        lows = np.array(data[x,:,0])
        
        highs = np.array(data[x,:,2])
      

        
    else:  
        names=namesY
        names[names=='electricityCostLow']='electricityCost\nLow'
        names[names=='electricityCostHigh']='electricityCost\nHigh'
        title=namesX[x]
        bases = np.array(data[:,x,1])
        base=bases[0]
        
        lows = np.array(data[:,x,0])
        
        highs = np.array(data[:,x,2])

  

    changeLows=(lows/bases-1)*100
    changeHighs=(highs/bases-1)*100
    sorts=np.absolute(highs-lows)
    sorts, lows, highs, names, bases, changeLows, changeHighs  = zip(*sorted(zip(sorts, lows, highs,  names, bases, changeLows, changeHighs), reverse=True))        

    
                
    ###############################################################################
    # The actual drawing part
    
    # The y position for each variable
    ys = range(len(highs))[::-1]  # top to bottom
    
    fig=plt.figure()
    if kind==0:
        fig.suptitle('Investigated Variable: ' + title)
    else:
        fig.suptitle('Investigated System: \n ' + title)
    axes = plt.subplot2grid((20,1), (1, 0), colspan=1, rowspan=19)
    # Plot the bars, one by one
    for y, low, high, base, changeLow, changeHigh in zip(ys, lows, highs, bases, changeLows, changeHighs):
        # The width of the 'low' and 'high' pieces
        low_width = base - low
        high_width = high - base
    
        # Each bar is a "broken" horizontal bar chart
        axes.broken_barh(
            [(low, low_width), (base, high_width)],
            (y - 0.4, 0.8),
            facecolors=[(81/255, 83/255, 86/255), (0/255, 91/255, 130/255)],  # Try different colors if you like
            edgecolors=['black', 'black'],#'none',
            linewidth=1,
        )
    
        # Display the high as text. It should be positioned in the center of
        # the 'high' bar, except if there isn't any room there, then it should be
        # next to bar instead.
#        if high<low:
#            x=-1
#        else:
        x = 1    
        if low_width > x:
            plt.text(base-low_width/2, y, str(np.around(changeLow,1))+'%', 
                     va='center', 
                     ha='center',
                     color='white')
        else:
            plt.text(low-x/2, y, str(np.around(changeLow,1))+'%', 
                     va='center', 
                     ha='center', 
                     bbox=dict(facecolor= 'white', alpha=0.5, edgecolor='none'))
            
        if high_width>x:
            plt.text(base+high_width/2, y, str(np.around(changeHigh,1))+'%',
                     va='center',
                     ha='center',
                     color='white')
        else:
            plt.text(high+x/2, y, str(np.around(changeHigh,1))+'%',
                     va='center',
                     ha='center', 
                     bbox=dict(facecolor= 'white', alpha=0.5, edgecolor='none'))
    
    # Draw a vertical line down the middle
    if kind==1: 
        axes.axvline(base, color='black')

    
            
    
    # Additional Information
    info='Distance: ' + str(dist).rstrip('0').rstrip('.') + ' km \nDemand: ' + str(dem/1000).rstrip('0').rstrip('.') + ' t/day'
#    fig.text(0.12,0.85, info,
#              backgroundcolor=(220/255,220/255,220/255)
#              )

    #info2='

    # Position the x-axis on the top, hide all the other spines (=axis lines)
    
    #axes.xaxis.set_ticks_position('top')
    axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f €'))
    axes.set_xlabel('Total Cost of Hydrogen in €/kg')
    #axes.xaxis.set_label_coords(0.5, 1.1)    
    # Make the y-axis display the variables
    plt.grid()
    plt.yticks(ys, names)
    plt.show()
    lowPatch=mpatches.Patch(color=(81/255, 83/255, 86/255), label='-50% value')
    highPatch=mpatches.Patch(color=(0/255, 91/255, 130/255), label='+50% value')
    infoPatch=mpatches.Patch(color='w', label=info)
    plt.legend(handles=[infoPatch, lowPatch, highPatch], loc=4)
    # Set the portion of the x- and y-axes to show
    plt.xlim(xmin, xmax)
    if kind == 0:
        plt.ylim(-1, len(names))

    #plt.tight_layout()
    plt.show()