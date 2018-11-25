# terminal call:
# python plot_1-D-NflexPD.py /file/path/file_name.nc plot_sediment(0/1) num_years_to_plot(counting_backwards_from_the_last)
    
from pylab import *
import netCDF4 as nc4
import netcdftime
import sys

def plot_maecs_Bpoolx2_phy():
    
    """from plot_maecs_omexdia import plot_maecs_omexdia 
    (time,z,dz,data,datanames)=plot_maecs_omexdia()"""    
    
    plottype='wc_mean' #wc_int, wc_mean,middlerow
    colmap='viridis'
    #import pdb
    if len(sys.argv) < 2: #this means no arguments were passed      
      fname='/home/onur/setups/test-BGCmodels/NflexPD/1D-40m/test/1D-40m_NflexPD'
      disp('plotting default file:'+fname)
    else:
      disp('plotting file specified:'+sys.argv[1])
      fname=sys.argv[1]
    
    if len(sys.argv)<3: #no second argument was passed
      plotsed=0
    else:
      plotsed=int(sys.argv[2])
    disp('plotsed:'+str(plotsed))
      
    if len(sys.argv)<4: #no third argument was passed
      numyears=1
    else: 
      numyears=int(sys.argv[3])
    disp('plotting last '+str(numyears)+' year of the simulation')
    
    if len(sys.argv) < 5: #this means no arguments were passed
	varnames= [ 'temp','phy_PAR','phy_NPR',
                'abio_din','abio_don','abio_detn',
                'phy_N', 'phy_Q', 'phy_mu',
                'phy_fV', 'phy_fA', 'phy_ThetaHat']
	numcol=3.0
    else: 
	varnames=sys.argv[4].split(',')
	numcol=length(varnames)
	
    if plotsed: 
        #from plot_sediment import readsed
        #sediment variables
        pickled=0
        varnames2=[]
        numsedvars=len(varnames2)
        figuresize=(15,8)
    else:
        numsedvars=0
        figuresize=(15,12) #(25,15)
    
    #pelagic variables
    nc=nc4.Dataset(fname+'.nc')
    ncv=nc.variables
    #print('available maecs variables:')        
    #disp(ncv)
    
    depth=np.squeeze(ncv['z'][:])
    if len(depth.shape)==2:
        if all(depth[0]==depth[-1]):
            depth=depth[0,:]
        else:
            depth = depth[0, :]
            #raise(Warning('plotting for varying-depths is not implemented')) #eg., it needs interpolation on a fixed grid

    #time=ncv['time'][:]/86400.

    tv = nc.variables['time']
    utime=netcdftime.utime(tv.units)
    tvec=utime.num2date(tv[:])

    #crop the data for the time period requested
    years=np.array([tvec[ti].year for ti in range(0,len(tvec))])
    years2plot=range(years[-1]+1-numyears, years[-1]+1)
    yeari=np.where((years>=years2plot[0]) * (years<=years2plot[-1]))
    tvecC=tvec[yeari[0]]

    f=figure(figsize=figuresize)
    f.subplots_adjust(top=0.95,bottom=0.05,hspace=0.5, wspace=0.5)

    numvar=len(varnames)+numsedvars

    for i,varn in enumerate(varnames):       
        print varn
        ax=subplot(ceil(numvar/numcol),numcol,i+1)

        if (varnames[i]=='skip'):
            continue
        if (not (varnames[i] in ncv)):
            ax.text(0.5,0.5,varnames[i]+'\n\n was not found',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            continue


        #pdb.set_trace()
        dat=squeeze(ncv[varnames[i]][:,:])

        #crop the data for the time period requested
        datC=dat[yeari[0],:]
        longname=ncv[varnames[i]].long_name
        shortname=longname.split('hzg_maecs ')[-1]

        if (np.max(datC)-np.min(datC)<1e-10):
            ax.text(0.5,0.5,varnames[i]+'\n\n all: %3.2f'%np.max(datC),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            continue
        else:
            units=ncv[varnames[i]].units
            if units in ['%']:
                title(shortname + ' [%s]'%units, size=10.0)
            elif units== 'Celsius':
                title('GOTM Temperature [$^oC$]', size=10.0)
            else:
                title(shortname + ' [$%s$]'%units, size=10.0)

            pcf=ax.contourf(tvecC,depth,transpose(datC),cmap=plt.get_cmap(colmap))

        #x-axis
        format_date_axis(ax,[tvecC[0], tvecC[-1]])
        ax.xaxis.grid(color='k',linestyle=':',linewidth=0.5)
        xlabel('')

        #y-axis
        #yt  = gca().get_yticks()
        #ytl = gca().get_yticklabels()
        #gca().set_yticks([yt[0],yt[-1]])
        #gca().set_yticklabels([str(yt[0]),str(yt[-1])])
        ylabel('depth [m]')

        cbar = colorbar(pcf, shrink=0.8)
        #cbar.solids.set_edgecolor("face")
        #draw()
     
        
    nc.close()
    savefig(fname+'_cont.png')
    disp('python contour plot saved in: '+fname+'_cont.png')
    #show()
        
    #if plotsed: return time,z,dz,data,datanames


def format_date_axis(ax,tspan):
    ax.set_xlim(tspan[0], tspan[1])
    if diff(tspan)[0].days<63:
        ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(byweekday=matplotlib.dates.MO) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b\n%d'))
    elif diff(tspan)[0].days<367:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=2) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=15)
    elif diff(tspan)[0].days<732:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=3) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=16)
    elif diff(tspan)[0].days<1466:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=6) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=16)
    elif diff(tspan)[0].days<3655:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
    elif diff(tspan)[0].days<9130: #25*365=9125
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=60) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
    else:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=120) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))

if __name__ == "__main__":
    # if you call this script from the command line (the shell) it will
    # run the 'main' function
    plot_maecs_Bpoolx2_phy()

