# terminal call:
# python plot_1D-NflexPD.py /file/path/file_name.nc num_years_to_plot(counting_backwards_from_the_last)
    
from pylab import *
import netCDF4 as nc4
import netcdftime
import sys
import warnings

varlims={'abio_PAR_dmean':[0,35], 'temp':[5,20],
         'abio_din':[0,35], 'abio_detn':[0,6], 'abio_don':[0,8],
         'PPR':[0,20.], 'N':[0,5.0], 'Q':[0.025,0.225], 'Chl2C':[0.0,0.5],
         'fA':[0.0,1.0], 'fV':[0.0,0.5], 'ThetaHat':[0.04,0.54]}
numlevels=6
def plot_nflexpd():

    #models = ['phy_cQ','phy_IOQf', 'phy_IOQ', 'phy_DOQ', 'phy_DOQf']
    models = ['phy_cQ', 'phy_IOQ', 'phy_DOQ']
    #models = ['phy_IOQ', 'phy_DOQ']
    vars2comp = ['PPR', 'N', 'Q', 'Chl2C', 'fA', 'fV', 'ThetaHat'] #
    plottype='wc_mean' #wc_int, wc_mean,middlerow
    colmap='viridis'
    #import pdb
    if len(sys.argv) < 2: #this means no arguments were passed      
      #fname='/home/onur/setups/test-BGCmodels/nflexpd/1D-NS-40m/1D-40m_NflexPD.nc'
      fname = '/home/onur/setups/test-BGCmodels/nflexpd/1D-ideal-NA/1D-NA_mean.nc'
      disp('plotting default file:'+fname)
    else:
      disp('plotting file specified:'+sys.argv[1])
      fname=sys.argv[1]
      
    if len(sys.argv)<3: #no third argument was passed
      numyears=-1 # -1 means plot everything
    else: 
      numyears=int(sys.argv[2]) #number of years to plot (counting from the last year backwards)
    disp('plotting last '+str(numyears)+' year of the simulation')

    if len(models)==2: #show the difference between models in the 3rd column
        numcol = 3.0
        figuresize = (13, 15)  # (25,15)
        varnames = ['airt', 'temp', 'abio_PAR_dmean', #nuh
                    'abio_din', 'abio_detn', 'abio_don']
        for var in vars2comp:
            varnames.append('%s_%s' % (models[0], var))
            varnames.append('%s_%s' % (models[1], var))
            varnames.append('%s_%s-%s_%s' % (models[0], var, models[1], var))
    elif len(models)>2:
        numcol = len(models)
        if len(models)==3:
            figuresize = (1+4*len(models), 15)
            varnames = ['airt','temp','abio_PAR_dmean', #nuh
                        'abio_din', 'abio_detn', 'abio_don']
        elif len(models)==4:
            figuresize = (1+4*len(models), 15)
            varnames = ['airt','temp', 'abio_PAR_dmean', 'nuh',
                        'abio_din', 'abio_detn', 'abio_don', 'total_nitrogen_calculator_result']
        elif len(models)==5:
            figuresize = (1+4*len(models), 15)
            varnames = ['temp', 'nuh', 'abio_PAR_dmean', 'skip','skip',
                        'abio_din', 'abio_detn', 'abio_don', 'total_nitrogen_calculator_result','skip']
        for var in vars2comp:
            for i in range(len(models)):
                varnames.append('%s_%s' % (models[i], var))
    
    #pelagic variables
    nc=nc4.Dataset(fname)
    ncv=nc.variables

    tv = nc.variables['time']
    utime=netcdftime.utime(tv.units)
    tvec=utime.num2date(list(tv[:]))

    #crop the data for the time period requested
    years=np.array([tvec[ti].year for ti in range(0,len(tvec))])
    if numyears==-1: #plot all years
        numyears=years[-1]-years[0]+1
    years2plot=range(years[-1]+1-numyears, years[-1]+1)
    yeari=np.where((years>=years2plot[0]) * (years<=years2plot[-1]))
    tvecC=tvec[yeari[0]]

    z=np.squeeze(ncv['z'][ti,:]) #depths at layer centers (fabm variables, temp, salt, etc)
    zi=np.squeeze(ncv['zi'][ti,:]) #depths at layer interfaces (diffusivities, fluxes, etc)

    f=figure(figsize=figuresize)
    f.subplots_adjust(top=0.95,bottom=0.05,hspace=0.5, wspace=0.5)

    numvar=len(varnames)

    for i,varn in enumerate(varnames):       
        print varn
        varn_basic=get_basic_varname(varn,models)

        ax=subplot(ceil(numvar/numcol),numcol,i+1)

        if (varn == 'skip'):
            continue

        varfound, dat, valsat, longname, units = get_varvals(ncv, varn)

        if (not varfound):
            ax.text(0.5, 0.5, varnames[i] + '\n\n was not found',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            continue

        if valsat == 'plate':
            depth=np.array([0])
            # crop the data for the time period requested
            datC = dat[yeari[0]]
        else:
            # crop the data for the time period requested
            datC = dat[yeari[0], :]
            if valsat == 'center':
                depth = z  # depth at centers
            elif valsat == 'int':
                depth = zi  # depth at interfaces

        # if depth vector is 2-D (vary with time)
        if len(depth.shape) == 2:
            # repeat the tvecC to obtain a matrix
            t = np.transpose(array([tvecC, ] * depth.shape[1]))
        else:
            t = tvecC

        #datC[datC<-1e10] = np.nan

        if (np.max(datC)-np.min(datC)<1e-10):
            ax.text(0.5,0.5,varnames[i]+'\n\n all: %3.2f'%np.max(datC),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            continue
        else:
            if units in ['%']:
                title(longname + ' [%s]'%units, size=10.0)
            elif units== 'Celsius':
                title('GOTM Temperature [$^oC$]', size=10.0)
            else:
                title(longname + ' [$%s$]'%units, size=10.0)

            if valsat == 'plate':
                ax.plot(t,datC)
            else:
                if len(z.shape) == 2:
                    pcf = ax.contourf(t, depth, datC, cmap=plt.get_cmap(colmap),vmin=vmin,vmax=vmax)
                else:
                    if varn_basic in varlims.keys():
                        levels=linspace(varlims[varn_basic][0],varlims[varn_basic][1],numlevels)
                        pcf=ax.contourf(tvecC,depth,transpose(datC),cmap=plt.get_cmap(colmap),levels=levels)
                    else:
                        pcf = ax.contourf(tvecC, depth, transpose(datC), cmap=plt.get_cmap(colmap))

                # y-axis
                # yt  = gca().get_yticks()
                # ytl = gca().get_yticklabels()
                # gca().set_yticks([yt[0],yt[-1]])
                # gca().set_yticklabels([str(yt[0]),str(yt[-1])])
                ylabel('depth [m]')

                cbar = colorbar(pcf, shrink=0.8)
                # cbar.solids.set_edgecolor("face")
                # draw()

        #x-axis
        format_date_axis(ax,[tvecC[0], tvecC[-1]])
        ax.xaxis.grid(color='k',linestyle=':',linewidth=0.5)
        xlabel('')
        
    nc.close()
    figname=fname.split('.nc')[0]+'_cont.png'
    savefig(figname)
    disp('python contour plot saved in: '+figname)
    #show()

def get_basic_varname(varn,models):
    for model in models:
        if model in varn:
            varn_basic=varn.split(model+'_')[1]
            return varn_basic
    return varn
def get_varvals(ncv,varn0):
    if '-' in varn0:
        varn=varn0.split('-')[0]
        varn2 = varn0.split('-')[1]
        v1=squeeze(ncv[varn][:,:])
        v2 = squeeze(ncv[varn2][:,:])
        varvals=v1-v2
        longname='%s - %s'%(varn,varn2)
        units=ncv[varn].units
    else:
        if not (varn0 in ncv):
            return (False,0,0,'','')
        else:
            varn=varn0
            varvals=squeeze(ncv[varn][:])
            longname = ncv[varn].long_name
            units=ncv[varn].units

    if len(varvals.shape)==1: #if 1-dimensional variable (e.g., airt)
        valsat='plate'
    elif varn in ['nuh', 'nus']:
        valsat='int'
    else:
        valsat='center'

    return (True,varvals,valsat,longname,units)

def format_date_axis(ax,tspan):
    ax.set_xlim(tspan[0], tspan[1])
    if diff(tspan)[0].days<63:
        ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(byweekday=matplotlib.dates.MO) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b\n%d'))
    elif diff(tspan)[0].days<367:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter(''))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=2) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        #ax.xaxis.set_tick_params(which='major', pad=15)
    elif diff(tspan)[0].days<732:
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=3) )
        ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=10)
    elif diff(tspan)[0].days<1466:
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=1))
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        #ax.xaxis.set_tick_params(which='major', pad=10)
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
    plot_nflexpd()

