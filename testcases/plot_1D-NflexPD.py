# terminal call:
# python plot_1D-NflexPD.py /file/path/file_name.nc num_years_to_plot(counting_backwards_from_the_last) modname
# python plot_1D-NflexPD.py ${rootdir}/${out_pref}/${out_pref}_mean.nc 3 FS-IA-DA
# python plot_1D-NflexPD.py ${rootdir}/${out_pref}/${out_pref}_mean.nc 3 phy_FS

import matplotlib.pyplot as plt
import netCDF4 as nc4
import netcdftime
import sys
import numpy as np

varlims={'abio_PAR_dmean':[0,30], 'airt':[0,21], 'I_0':[0,250],'temp':[2,22], 'mld_surf':[-100,0],'wind':[-6,26],
         'abio_detc_sed/abio_detn_sed':[4.0,14.0],'abio_detc/abio_detn':[5.0,30.0],'abio_doc/abio_don':[5.0,30.0],
         'abio_din': [0, 30],'abio_detc':[0,80],'abio_detn':[0,6], 'abio_doc':[0,80], 'abio_don':[0,8],
         'Chl':[0,10.],'C':[0,50.0],'N':[0,7.5],'Q':[0.025,0.225],'Chl2C':[0.0,0.5],
         'PPR':[0,20.],'mu':[0,0.5],'V_N':[0,0.05],'R_N':[0,0.05],'R_Chl':[0,0.1],
         'fA':[0.0,1.0], 'fV':[0.0,0.5], 'ThetaHat':[0.04,0.54],'fN':[0,0.2],'fL':[0.,1]}
prettyunits={'abio_detc_sed/abio_detn_sed':'molC/molN','abio_detc/abio_detn':'molC/molN','abio_doc/abio_don':'molC/molN'}
prettynames={'abio_PAR_dmean':'\overline{I}','I_0':'I_{0}','mld_surf':'\mathrm{MLD}',
             'airt': 'T_{air}','temp': 'T','u10':'\mathrm{Wind \ Speed \ (-u)}','wind':'\mathrm{Wind \ Speed}',
             'abio_din':'DIN','abio_detc_sed/abio_detn_sed':'\mathrm{C:N \ of \ POM-export}','abio_detc/abio_detn':'POC:PON','abio_doc/abio_don':'DOC:DON',
             'abio_detc':'POC','abio_detn':'PON','abio_doc':'DOC','abio_don':'DON',
             'Chl':'Phy_{Chl}','C':'Phy_C','N':'Phy_N',
             'Q':'Q','V_N':'f_{DIN-Phy}','mu':'\mu',
             'R_N':'R_N','R_Chl':'R_{Chl}','fN':'L_N','fL':'L_I',
             'fA':'f_A','fV':'f_V','ThetaHat':'\hat{\Theta}'}
numlevels=6
#depth range to be shown:
prescylim=[0,0] #[0,0] means no ylimits are prescribed, full depth range will be shown'
#prescylim=[-30,0.0]

def main(fname, numyears, modname):
    
    if not '-' in modname: #i.e., single model runs
      models = [modname]
      varsets={#'abio0':['airt', 'wind', 'I_0'],
             'abio1':['abio_PAR_dmean','temp', 'mld_surf'],
             'abio2':['abio_din','abio_detc/abio_detn','abio_detc_sed/abio_detn_sed'], #'abio_doc/abio_don'],
             'abio3':['abio_detn','abio_detc','abio_don','abio_doc'],
             'phy-1':['C','N','Q'],
             'phy-2':['mu','V_N','R_N','R_Chl'],
             'phy-3':['fA','fV','fN','fL'] #, 'ThetaHat']
             }
    elif modname=='FS-IA-DA': #i.e., competition experiment
      #models = ['phy_IOQ', 'phy_DOQ']
      #models = ['phy_cQ','phy_IOQf', 'phy_IOQ', 'phy_DOQ', 'phy_DOQf']
      models = ['phy_FS', 'phy_IA', 'phy_DA']
      varsets={#'abio0':['I_0','airt', 'wind'],
             #'abio12':['temp','mld_surf','abio_din','abio_PAR_dmean',],
             #'abio3':['abio_detn','abio_detc','abio_don','abio_doc'],
             #'phy-1':['C','N','Q'],
             #'phy-2':['mu','V_N','R_N','R_Chl'],
             'phy-3':['fA','fV','fN','fL'] #, 'ThetaHat']
             }
      
    for groupname,varset in varsets.iteritems():
        #print ('%s,%s'%(groupname,varset))
        plot_nflexpd(fname,numyears,groupname,varset,models)

def plot_nflexpd(fname,numyears,groupname,varset,models):

    vars2comp = ['Chl', 'PPR', 'Q'] # 'Chl2C', 'fA', 'fV', 'ThetaHat'] #
    plottype='wc_mean' #wc_int, wc_mean,middlerow
    colmap='viridis'
    
    if 'phy' in groupname:
        if len(models)==1:
            if len(varset)>4:
                numcol = 3.0
            else:
                numcol=len(varset)
            figuresize = (1 + 4 * numcol, .5 + 1.5 * len(varset)/numcol)
            varnames = []
            for var in varset:
                varnames.append('%s_%s' % (models[0], var))
            if len(varset)/numcol==1:
               fpar = {'left':0.05, 'right':0.99, 'top': 0.85, 'bottom': 0.15, 'hspace': 0.5, 'wspace': 0.2}
            else:
                fpar = {'left':0.05, 'right':0.99, 'top': 0.93, 'bottom': 0.07, 'hspace': 0.5, 'wspace': 0.2}
        elif len(models)==2: #show the difference between models in the 3rd column
            numcol = 3.0
            figuresize = (1 + 4 * numcol, 1 + 1.5 * len(varset))
            varnames = []
            for var in varset:
                varnames.append('%s_%s' % (models[0], var))
                varnames.append('%s_%s' % (models[1], var))
                varnames.append('%s_%s-%s_%s' % (models[0], var, models[1], var))
            fpar = {'top': 0.95, 'bottom': 0.05, 'hspace': 0.5, 'wspace': 0.2}
        elif len(models)>2:
            numcol = len(models)
            figuresize = (1 + 4 * len(models), .5 + 1.5 * len(varset))
            varnames = []
            for var in varset:
                for i in range(len(models)):
                    varnames.append('%s_%s' % (models[i], var))
            if len(varset)==2:
                fpar = {'left': 0.05, 'right': 0.99, 'top': 0.93, 'bottom': 0.07, 'hspace': 0.5, 'wspace': 0.2}
            elif len(varset)>2:
                fpar = {'left':0.05, 'right':0.99, 'top': 0.95, 'bottom': 0.05, 'hspace': 0.5, 'wspace': 0.2}
            else:
                print('len(varset):%s'%(len(varset)))
    else:
        varnames = varset
        if len(varset)>4:
            numcol = 3.0
        else:
            numcol=len(varset)
        #figuresize = (1 + 4 * len(models), .5 + 1.5 * len(varset)/numcol)
        figuresize = (1 + 4 * numcol, .5 + 1.5 * len(varset)/numcol)
        if len(varset)/numcol==1:
            fpar = {'left':0.05, 'right':0.99, 'top': 0.85, 'bottom': 0.15, 'hspace': 0.5, 'wspace': 0.2}
        else:
            fpar = {'left':0.05, 'right':0.99, 'top': 0.93, 'bottom': 0.07, 'hspace': 0.5, 'wspace': 0.2}
    
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

    #assume that z does not vary and vectorize:
    z=np.squeeze(ncv['z'][0,:]) #depths at layer centers (fabm variables, temp, salt, etc)
    zi=np.squeeze(ncv['zi'][0,:]) #depths at layer interfaces (diffusivities, fluxes, etc)

    f=plt.figure(figsize=figuresize)
    f.subplots_adjust(left=fpar['left'],right=fpar['right'],top=fpar['top'],bottom=fpar['bottom'],hspace=fpar['hspace'],wspace=fpar['wspace'])

    numvar=len(varnames)

    for i,varn in enumerate(varnames):       
        print (varn)
        varn_basic,model=get_basic_varname(varn,models)

        ax=plt.subplot(np.ceil(numvar/numcol),numcol,i+1)

        if (varn == 'skip'):
            continue

        if varn=='wind':
            varfound, datU, valsat, longname, units = get_varvals(ncv, 'u10')
            varfound, datV, valsat, longname, units = get_varvals(ncv, 'v10')
            dat=np.sqrt(datU**2+datV**2)
            #restore negative velocities
            ineg=datU<0.0
            dat[ineg]=dat[ineg]*-1.0
        else:
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

        #not really plot, if all values are same
        if False: #(np.max(datC)-np.min(datC)<1e-10):
            ax.text(0.5,0.5,varnames[i]+'\n\n all: %3.2f'%np.max(datC),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            continue
        else:
            #if units in ['%']:
            #    title(longname + ' [%s]'%units, size=10.0)
            if units== 'Celsius':
                units='^oC'
            elif varn=='I_0' and units=='W/m2':
                #convert to Watts to Einstein/d
                units='E/m^2/d'
                datC=datC*4.6 * 1e-6 * 86400
            elif varn=='mld_surf':
                #convert to absolute depth
                datC=datC*-1
            if model=='':
                prettyname='$%s$'%prettynames[varn_basic]
            else:
                prettyname='(%s) $%s$'%(model.split('phy_')[1],prettynames[varn_basic])
            if varn_basic in prettyunits:
                units = prettyunits[varn_basic]
            plt.title('%s [$%s$]'%(prettyname,units), size=12.0)

            cmap = plt.get_cmap(colmap)
            if valsat == 'plate':
                ax.plot(t,datC)
                #shrink the axes width by 20% to fit that of the contour plots
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, (box.x1 - box.x0) * 0.8, box.y1 - box.y0])
                ax.set_ylim(varlims[varn][0],varlims[varn][1])
            else:
                if len(z.shape) == 2:
                    pcf = ax.contourf(t, depth, datC, cmap=cmap,vmin=vmin,vmax=vmax)
                else:
                    if varn_basic in varlims.keys():
                        levels=np.linspace(varlims[varn_basic][0],varlims[varn_basic][1],numlevels)
                        extendopt,cmap=get_extendopt(levels,datC,cmap)
                        pcf=ax.contourf(tvecC,depth,np.transpose(datC),cmap=cmap,levels=levels,extend=extendopt)
                    else:
                        pcf = ax.contourf(tvecC, depth, np.transpose(datC), cmap=cmap)
                    if not (prescylim[0]==0 and prescylim[1]==0):
                        ax.set_ylim(prescylim[0], prescylim[1])
                        ylimsuf='_ylim_%s-%s'%(prescylim[0],prescylim[1])
                    else:
                        ylimsuf=''
                plt.ylabel('depth [m]')

                cbar = plt.colorbar(pcf, shrink=0.9)
                # cbar.solids.set_edgecolor("face")
                # draw()
                if (np.max(datC)-np.min(datC)<1e-4):
                    #print mean concentration (rounded to 3rd digit)
                    meanc=np.round(np.mean(datC)*1000)/1000
                    ax.text(0.5,0.5,'constant: %3.3f'%meanc,
                            horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
        
        ax.grid(b=True, axis='y', which='major', color='0.5', linestyle='-')
        
        #x-axis
        format_date_axis(ax,[tvecC[0], tvecC[-1]])
        plt.xlabel('')

    nc.close()
    figname=fname.split('.nc')[0]+'_cont_'+groupname+ '_'+str(numyears)+'y'+ylimsuf+'.png'
    plt.savefig(figname)
    print('python contour plot saved in: '+figname)

def get_extendopt(levels, datC,cmap):

    if np.max(datC) > levels[-1]:
        if np.min(datC) < levels[0]:
            extendopt = "both"
        else:
            extendopt = "max"
    elif np.min(datC) < levels[0]:
        extendopt = "min"
    else:
        extendopt = "neither"

    if extendopt in ['max','both']:
        cmap.set_over('lightyellow')
    if extendopt in ['min','both']:
        cmap.set_under('darkslategray')
    return (extendopt,cmap)

def get_basic_varname(varn,models):
    for model in models:
        if model in varn:
            varn_basic=varn.split(model+'_')[1]
            return (varn_basic,model)
    return (varn,'')

def get_varvals(ncv,varn0):
    if '+' in varn0  or '-' in varn0 or '*' in varn0 or '/' in varn0:
        varfound,varvals,valsat,longname,units=get_varvals_op(ncv,varn0)
    else:
        if not (varn0 in ncv):
            return (False,0,0,'','')
        else:
            varn=varn0
            varvals=np.squeeze(ncv[varn][:])
            longname = ncv[varn].long_name
            units=ncv[varn].units
            if len(varvals.shape)==1: #if 1-dimensional variable (e.g., airt)
                valsat='plate'
            elif varn in ['nuh', 'nus']:
                valsat='int'
            else:
                valsat='center'
    return (True,varvals,valsat,longname,units)

def get_varvals_op(ncv,varn0):
    if '+' in varn0:
        symb='+'
    elif '-' in varn0:
        symb='-'
    elif '*' in varn0:
        symb='*'
    elif '/' in varn0:
        symb='/'
    else:
        return (False,0,0,'','')
    varn=varn0.split(symb)[0]
    varn2 = varn0.split(symb)[1]
    v1 = np.squeeze(ncv[varn][:,:])
    v2 = np.squeeze(ncv[varn2][:,:])
    longname='%s %s %s'%(varn,symb,varn2)
    if '+' in varn0:
        varvals=v1+v2
        units=ncv[varn].units
    elif '-' in varn0:
        varvals=v1-v2
        units=ncv[varn].units
    elif '*' in varn0:
        varvals=v1*v2
        units=ncv[varn].units + '*' + ncv[varn2].units
    elif '/' in varn0:
        varvals=v1/v2
        units=ncv[varn].units + '/' + ncv[varn2].units
    if len(varvals.shape)==1: #if 1-dimensional variable (e.g., airt)
        valsat='plate'
    elif varn in ['nuh', 'nus']:
        valsat='int'
    else:
        valsat='center'
    return (True,varvals,valsat,longname,units)
        
def format_date_axis(ax,tspan):
    import matplotlib.dates as mpldates
    ax.set_xlim(tspan[0], tspan[1])
    if np.diff(tspan)[0].days<63:
        ax.xaxis.set_major_locator(mpldates.WeekdayLocator(byweekday=mpldates.MO) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%b\n%d'))
    elif np.diff(tspan)[0].days<367:
        ax.xaxis.set_major_locator(mpldates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter(''))
        ax.xaxis.set_minor_locator(mpldates.MonthLocator(bymonthday=1, interval=2) )
        ax.xaxis.set_minor_formatter(mpldates.DateFormatter('%b'))
        ax.grid(b=True, axis='x', which='minor', color='0.5', linestyle='-')
        #ax.xaxis.set_tick_params(which='major', pad=15)
    elif np.diff(tspan)[0].days<732:
        ax.xaxis.set_major_locator(mpldates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mpldates.MonthLocator(bymonthday=1, interval=3) )
        ax.xaxis.set_minor_formatter(mpldates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(which='major', pad=10)
    elif np.diff(tspan)[0].days<1466:
        ax.xaxis.set_minor_locator(mpldates.MonthLocator(bymonthday=1, interval=1))
        ax.xaxis.set_major_locator(mpldates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%Y'))
        ax.grid(b=True, axis='x', which='major', color='0.5', linestyle='-')
        #ax.xaxis.set_tick_params(which='major', pad=10)
    elif np.diff(tspan)[0].days<3655:
        ax.xaxis.set_major_locator(mpldates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%Y'))
    elif np.diff(tspan)[0].days<9130: #25*365=9125
        ax.xaxis.set_major_locator(mpldates.MonthLocator(bymonthday=1, interval=60) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%Y'))
    else:
        ax.xaxis.set_major_locator(mpldates.MonthLocator(bymonthday=1, interval=120) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%Y'))

if __name__ == "__main__":
    # if you call this script from the command line (the shell) it will
    # run the 'main' function
    if len(sys.argv) < 2: #this means no arguments were passed      
      #fname='/home/onur/setups/test-BGCmodels/nflexpd/1D-NS-40m/1D-40m_NflexPD.nc'
      fname = '/home/onur/setups/test-BGCmodels/nflexpd/1D-ideal-highlat/20-05-15/Highlat-100m_wconst_FS-IA-DA/Highlat-100m_wconst_FS-IA-DA_mean.nc'
      print('plotting default file:'+fname)
    else:
      print('plotting file specified:'+sys.argv[1])
      fname=sys.argv[1]
      
    if len(sys.argv)<3: #no third argument was passed
      numyears=-1 # -1 means plot everything
    else: 
      numyears=int(sys.argv[2]) #number of years to plot (counting from the last year backwards)
    print('plotting last '+str(numyears)+' year of the simulation')
    
    if len(sys.argv)<4:
      modname='FS-IA-DA'
    else:
      modname=sys.argv[3]
    main(fname, numyears, modname)
