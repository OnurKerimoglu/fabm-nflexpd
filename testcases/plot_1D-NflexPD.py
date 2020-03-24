# terminal call:
# python plot_1-D-NflexPD.py /file/path/file_name.nc plot_sediment(0/1) num_years_to_plot(counting_backwards_from_the_last)
    
from pylab import *
import netCDF4 as nc4
import netcdftime
import sys
import warnings

def plot_maecs_Bpoolx2_phy():
    
    """from plot_maecs_omexdia import plot_maecs_omexdia 
    (time,z,dz,data,datanames)=plot_maecs_omexdia()"""    
    
    plottype='wc_mean' #wc_int, wc_mean,middlerow
    colmap='viridis'
    #import pdb
    if len(sys.argv) < 2: #this means no arguments were passed      
      fname='/home/onur/setups/test-BGCmodels/PMbench/1D-NS-40m/1D-40m_NflexPD.nc'
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
	varnames= [ 'temp','nuh','abio_PAR',
                'abio_din','abio_don','abio_detn',
                'phy_IOQ_PPR', 'phy_DOQ_PPR','total_PPR_calculator_result',
                'phy_IOQ_N', 'phy_IOQ_Q', 'phy_IOQ_Chl2C',
                'phy_DOQ_N', 'phy_DOQ_Q', 'phy_DOQ_Chl2C',
                'phy_IOQ_fV', 'phy_IOQ_fA', 'phy_IOQ_ThetaHat',
                'phy_DOQ_fV', 'phy_DOQ_fA', 'phy_DOQ_ThetaHat']
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
    nc=nc4.Dataset(fname)
    ncv=nc.variables
    #print('available maecs variables:')        
    #disp(ncv)
    
    z=np.squeeze(ncv['z'][:]) #depths at layer centers (fabm variables, temp, salt, etc)
    zi=np.squeeze(ncv['zi'][:]) #depths at layer interfaces (diffusivities, fluxes, etc)

    tv = nc.variables['time']
    utime=netcdftime.utime(tv.units)
    tvec=utime.num2date(list(tv[:]))

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
        if varn in ['nuh','nus']:
            depth=zi #depth at interfaces
        else:
            depth=z
            
        #if depth vector is 2-D (vary with time)
        if len(depth.shape)==2:
            # repeat the tvecC to obtain a matrix
            t=np.transpose(array([tvecC,]*depth.shape[1]))
        else:
            t=tvecC
            
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
        datC[datC<-1e10] = np.nan

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

            if len(z.shape) == 2:
                pcf = ax.contourf(t, depth, datC, cmap=plt.get_cmap(colmap))
            else:
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
    figname=fname.split('.nc')[0]+'_cont.png'
    savefig(figname)
    disp('python contour plot saved in: '+figname)
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

