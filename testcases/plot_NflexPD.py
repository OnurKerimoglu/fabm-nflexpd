# terminal call:
# python plot_1D-NflexPD.py filenames idnames num_years_to_plot(counting_backwards_from_the_last) modnames
# python plot_1D-NflexPD.py file_name1.nc,file_name2.nc Nbased,Cbased 3

import matplotlib.pyplot as plt
import netCDF4 as nc4
import netcdftime
import sys
import numpy as np

varlims1D={'I_dm':[0,30], 'airt':[0,21], 'I_0':[0,250],'temp':[2,22], 'mld_surf':[-100,0],'wind':[-6,26],
         'DetC_sed/DetN_sed':[4.0,16.0],'DetC/DetN':[5.0,30.0],'DOC/DON':[5.0,30.0],
         'DIN': [1, 26],'DetC':[0,80],'DetN':[0,6], 'DOC':[0,80], 'DON':[0,6],
         'Phy-Chl':[0,10.],'Phy-C':[0,50.0],'Phy-N':[0,7.5],'Phy-Q':[0.02,0.22],'Phy-Chl2C':[0.00,0.05],
         'Phy-PPR':[0,20.],'Phy-mu':[0,0.4],'Phy-vN':[0,0.05],'Phy-f_dinphy':[0,0.5],'Phy-R_N':[0,0.04],'Phy-R_Chl':[0,0.1],
         'Phy-fA':[0.0,1.0], 'Phy-fV':[0.0,0.5], 'Phy-ThetaHat':[0.00,0.05],'Phy-fC':[0,0.2],'Phy-limfunc_L':[0.,1]}
varlims0D={'I_dm':[0,30], 'airt':[0,21], 'I_0':[0,250],'temp':[2,22], 'mld_surf':[-100,0],'wind':[-6,26],
         'DetC_sed/DetN_sed':[4.0,16.0],'DetC/DetN':[5.0,30.0],'DOC/DON':[5.0,30.0],
         'DIN': [1, 20],'DetC':[0,250],'DetN':[0,15], 'DOC':[0,250], 'DON':[0,15],
         'Phy-Chl':[0,20.],'Phy-C':[0,100.0],'Phy-N':[0,7.5],'Phy-Q':[0.02,0.22],'Phy-Chl2C':[0.00,0.1],
         'Phy-PPR':[0,20.],'Phy-mu':[0,0.4],'Phy-vN':[0,0.05],'Phy-f_dinphy':[0,0.5],'Phy-R_N':[0,0.04],'Phy-R_Chl':[0,0.1],
         'Phy-fA':[0.0,1.0], 'Phy-fV':[0.0,0.5], 'Phy-ThetaHat':[0.00,0.05],'Phy-fC':[0,0.2],'Phy-limfunc_L':[0.,1]}
namelibNbased={'I_0':'I_0','wind':'m\ s^{-1}','T':'temp','totalN':'total_nitrogen_calculator_result',
               'I-dm':'abio_PAR_dmean','I':'abio_PAR',
             'Phy-C':'phy_IA_C','Phy-N':'phy_IA_N','Phy-Q':'phy_IA_Q',
             'Phy-Chl':'phy_IA_Chl','Phy-Chl2C':'phy_IA_Chl2C',
             'mu':'phy_IA_mu','vN':'phy_IA_vN','R_N':'phy_IA_R_N','R_Chl':'phy_IA_R_Chl',
             'fA':'phy_IA_fA', 'fV':'phy_IA_fV', 'fC':'phy_IA_fC',
             'ThetaHat':'phy_IA_ThetaHat', 'limfunc_L':'phy_IA_limfunc_L',
             'DIN':'abio_din','DON':'abio_don','DetN':'abio_detn',
             'DOC':'abio_doc','DetC':'abio_detc'
            }
namelibCbased={'I_0':'I_0','wind':'m\ s^{-1}','T':'temp','totalN':'total_nitrogen_calculator_result',
               'I-dm':'phyabio_PAR_dmean','I':'phyabio_PAR',
             'Phy-C':'phyabio_C','Phy-N':'phyabio_N','Phy-Q':'phyabio_Q',
             'Phy-Chl':'phyabio_Chl','Phy-Chl2C':'phyabio_Chl2C',
             'mu':'phyabio_mu','vN':'phyabio_vN','R_N':'phyabio_R_N','R_Chl':'phyabio_R_Chl',
             'fA':'phyabio_fA', 'fV':'phyabio_fV', 'fC':'phyabio_fC',
             'ThetaHat':'phyabio_ThetaHat', 'limfunc_L':'phyabio_limfunc_L',
             'DIN':'phyabio_din','DON':'phyabio_don','DetN':'phyabio_detn',
             'DOC':'phyabio_doc','DetC':'phyabio_detc'
            }
prettyunits={'I_0':'E\ m^{-2}\ d^{-1}','wind':'m\ s^{-1}','T':'^\circ C',
             'I-dm':'E\ m^{-2}\ d^{-1}','I':'E\ m^{-2}\ d^{-1}','totalN':'mmolN\ m^{-3}',
             'Phy-C':'mmolC\ m^{-3}','Phy-N':'mmolN\ m^{-3}','Phy-Q':'molN\ molC^{-1}',
             'Phy-Chl':'mg m^{-3}','Phy-Chl2C':'gChl\ gC^{-3}',
             'DIN':'mmolN\ m^{-3}','DON':'mmolN\ m^{-3}','DetN':'mmolN\ m^{-3}',
             'DOC':'mmolC\ m^{-3}','DetC':'mmolC\ m^{-3}',
              'mu':'d^{-1}', 'vN':'molN\ molC^{-1} d^{-1}', 'R_N':'d^{-1}', 'R_Chl':'d^{-1}',
             #'ThetaHat', 'fA', 'fV', 'fC', 'limfunc_L'
             }
numlevels=6
#depth range to be shown:
prescylim=[0,0] #[0,0] means no ylimits are prescribed, full depth range will be shown'
#prescylim=[-30,0.0]

def main(fnames, sids, numyears, modnames):

  varsets={#'abio0':['airt', 'wind', 'I_0'], #I_dm
           #'abio1':['abio_PAR_dmean','temp', 'mld_surf'],
           #'abio2':['abio_din','abio_detc/abio_detn','abio_detc_sed/abio_detn_sed'],
           #'abio3':['abio_detn','abio_detc','abio_don','abio_doc'],
            'abio1':['DIN','totalN','I','I-dm',
                     'DOC','DON','DetC', 'DetN'],
            'phy-1':['Phy-C','Phy-N','Phy-Q',
                     '',     'Phy-Chl','Phy-Chl2C'],
            'phy-2': ['mu', 'vN', 'R_N', 'R_Chl'],
            'phy-3': ['fA', 'fV', 'fC',
                      'ThetaHat', 'limfunc_L',''],
            #'phy-3b': ['ThetaHat', 'fA', 'fV', 'limfunc_Nmonod', 'limfunc_L']
         }
  for groupname,varset in varsets.iteritems():
      #print ('%s,%s'%(groupname,varset))
      if len(sids)>1:
          plot_multifile(fnames, numyears, groupname, varset, sids,modnames)
      #if not (len(modnames)==2 and 'OBS' in modnames):
      #  for sidno in range(0,len(sids)):
      #      plot_singlefile(fnames[sidno], numyears, groupname, varset, sids[sidno],modnames[sidno])

def plot_multifile(fnames, numyears, groupname, varset, sids, modnames):
    cols=['darkblue','orange','green']
    linestyles=['-',':','--']
    varlims=varlims0D
    if len(varset)>3:
        if len(varset) in [3,6,9]:
            numcol=3.0
        elif len(varset) in [4,8,12]:
            numcol=4.0
    else:
        numcol = len(varset)
    numrow = np.ceil(len(varset) / numcol)
    figuresize = (1 + 4*numcol, 1. + 1.5*numrow)
    fpar = {'left': 0.05, 'right': 0.99, 'top': 0.93, 'bottom': 0.08, 'hspace': 0.5, 'wspace': 0.2}
    if numrow == 1:
        fpar['top']=0.9; fpar['bottom']=0.13

    ncL,ncvL,tL,tiL,fname_sid_str = collect_data(fnames,sids,numyears)

    f = plt.figure(figsize=figuresize)
    f.subplots_adjust(left=fpar['left'], right=fpar['right'], top=fpar['top'], bottom=fpar['bottom'],
                      hspace=fpar['hspace'], wspace=fpar['wspace'])

    for j, varn in enumerate(varset):
        print(varn)
        ax = plt.subplot(numrow, numcol, j + 1)

        if (varn == ''):
            ax.axis('off')
            continue

        units=''
        for i,sid in enumerate(sids):
            ncv=ncvL[i];t=tL[i];ti=tiL[i]
            if modnames[i]=='Nbased':
                namelib=namelibNbased
            elif modnames[i]=='Cbased':
                namelib=namelibCbased
            if varn in namelib:
                fvarn=namelib[varn]
            else:
                fvarn=varn
            varfound, dat, valsat, longname, unitsnew = get_varvals(ncv, fvarn, avg=True,zint=[-5,0])
            #update the units only if the units of the new data set exist
            if units=='' and unitsnew!='':
                units=unitsnew
            if varfound:
                datC = dat[ti[0]]
                if modnames[i]=='OBS':
                    ax.plot(t, datC, label=sid, color=cols[i], linestyle='',marker='o', markersize=4)
                else:
                    ax.plot(t, datC, label=sid, color=cols[i],linestyle=linestyles[i])
        if varn in prettyunits:
            prettyunit = prettyunits[varn]
        else:
            prettyunit = units
        if prettyunit=='':
            prettyunit='-'
        plt.title('%s [$%s$]' % (varn, prettyunit), size=12.0)

        if varn in varlims.keys():
            ax.set_ylim(varlims[varn][0], varlims[varn][1])
        if not (prescylim[0] == 0 and prescylim[1] == 0):
            ax.set_ylim(prescylim[0], prescylim[1])
            # ylimsuf='_ylim_%s-%s'%(prescylim[0],prescylim[1])

        # shrink the axes width by 20% to fit that of the contour plots, and put the legend in that space
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, (box.x1 - box.x0) * 0.8, box.y1 - box.y0])
        if j==0: #(j+1)%numcol==1 
            #ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.8),fontsize=12)
            ax.legend(loc='center left', bbox_to_anchor=(0.6, 0.7),fontsize=12)

        ax.grid(b=True, axis='y', which='major', color='0.5', linestyle='-')
        # x-axis
        format_date_axis(ax, [t[0], t[-1]])
        plt.xlabel('')

    for i in range(0,len(fnames)):
        ncL[i].close()
    figname = fname_sid_str +'_cont_' + groupname + '_' + str(numyears) + 'y.png'
    plt.savefig(figname)
    print('python line plot saved in: ' + figname)

def plot_singlefile(fname,numyears,groupname,varset,sid,modname):
    colmap='viridis'
    axgrid=True
    varlims = varlims1D
    #default:
    fpar = {'left': 0.05, 'right': 0.99, 'top': 0.85, 'bottom': 0.15, 'hspace': 0.5, 'wspace': 0.2}

    varnames = varset
    if len(varset)<3:
        numcol = len(varset)
    elif len(varset)%4 == 0:
        numcol = 4.
    elif len(varset)%3 == 0:
        numcol = 3.

    numrow=len(varset)/numcol
    #figuresize = (1 + 4 * len(models), .5 + 1.5 * len(varset)/numcol)
    figuresize = (1 + 4 * numcol, .5 + 1.5 * numrow)
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
    if 'z' in ncv.keys():
        z=np.squeeze(ncv['z'][0,:]) #depths at layer centers (fabm variables, temp, salt, etc)
        zi=np.squeeze(ncv['zi'][0,:]) #depths at layer interfaces (diffusivities, fluxes, etc)

    f=plt.figure(figsize=figuresize)
    f.subplots_adjust(left=fpar['left'],right=fpar['right'],top=fpar['top'],bottom=fpar['bottom'],hspace=fpar['hspace'],wspace=fpar['wspace'])

    numvar=len(varnames)
    #ylimsuf = ''
    for i,varn in enumerate(varnames):
        print (varn)

        ax = plt.subplot(numrow, numcol, i + 1)

        if (varn == ''):
            ax.axis('off')
            continue

        if modname == 'Nbased':
            namelib = namelibNbased
        elif modname == 'Cbased':
            namelib = namelibCbased
        if varn in namelib:
            fvarn = namelib[varn]
        else:
            fvarn = varn

        if fvarn=='wind':
            varfound, datU, valsat, longname, units = get_varvals(ncv, 'u10')
            varfound, datV, valsat, longname, units = get_varvals(ncv, 'v10')
            dat=np.sqrt(datU**2+datV**2)
            #restore negative velocities
            ineg=datU<0.0
            dat[ineg]=dat[ineg]*-1.0
        else:
            varfound, dat, valsat, longname, units = get_varvals(ncv, fvarn)

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
            t = np.transpose(np.array([tvecC, ] * depth.shape[1]))
        else:
            t = tvecC

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

        if varn in prettyunits:
            prettyunit = prettyunits[varn]
        else:
            prettyunit = units
        if prettyunit == '':
            prettyunit = '-'
        plt.title('%s [$%s$]'%(varn,prettyunit), size=12.0)

        cmap = plt.get_cmap(colmap)
        if valsat == 'plate':
            ax.plot(t,datC)
            #shrink the axes width by 20% to fit that of the contour plots
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, (box.x1 - box.x0) * 0.8, box.y1 - box.y0])
            if varn in varlims.keys():
                ax.set_ylim(varlims[varn][0],varlims[varn][1])
        else:
            if ((np.max(datC)+9.9<0.01) & (np.min(datC)+9.9<0.01)):
                datC[:,:]=np.nan
            #coloring of values exactly 0.0 becomes arbitrary. Tip them over the positive side to prevent alarming results
            datC[datC==0.0]=1e-15
            if len(z.shape) == 2:
                pcf = ax.contourf(t, depth, datC, cmap=cmap) #vmin=vmin,vmax=vmax
            else:
                if varn in varlims.keys():
                    levels=np.linspace(varlims[varn][0],varlims[varn][1],numlevels)
                    extendopt,cmap=get_extendopt(levels,datC,cmap)
                    pcf=ax.contourf(tvecC,depth,np.transpose(datC),cmap=cmap,levels=levels,extend=extendopt)
                else:
                    pcf = ax.contourf(tvecC, depth, np.transpose(datC), cmap=cmap)
                if not (prescylim[0]==0 and prescylim[1]==0):
                    ax.set_ylim(prescylim[0], prescylim[1])
                    #ylimsuf='_ylim_%s-%s'%(prescylim[0],prescylim[1])
            plt.ylabel('depth [m]')

            cbar = plt.colorbar(pcf, shrink=0.9)
            # cbar.solids.set_edgecolor("face")
            # draw()
            #not really plot, if the variable is marked as 'missing value (-9.9'
            axgrid=True
            if (np.isnan(datC).all()):
                ax.text(0.5,0.5,'N/A',transform=ax.transAxes,
                        horizontalalignment='center',verticalalignment='center')
                axgrid=False

            elif (np.max(datC)-np.min(datC)<1e-4):
                #print mean concentration (rounded to 3rd digit)
                meanc=np.round(np.mean(datC)*1000)/1000
                ax.text(0.5,0.5,'constant: %3.3f'%meanc,transform=ax.transAxes,
                        horizontalalignment='center',verticalalignment='center')
                axgrid=False
            else:
                ax.grid(b=True, axis='y', which='major', color='0.5', linestyle='-')

        #x-axis
        format_date_axis(ax,[tvecC[0], tvecC[-1]], axgrid)
        plt.xlabel('')

    nc.close()

    figname=fname.split('.nc')[0]+'_cont_'+groupname+ '_'+str(numyears)+'y.png'
    plt.savefig(figname)
    print('python contour plot saved in: '+figname)

def collect_data(fnames,sids,numyears):
    ncL=[]; ncvL=[]; tL=[]; tiL=[]
    fname_sid_str=fnames[0].split('.nc')[0]
    for fno in range(0,len(fnames)):
        fname_sid_str=fname_sid_str+ '_vs_' + sids[fno]
        nc = nc4.Dataset(fnames[fno])
        ncv = nc.variables

        #handle time
        tv = nc.variables['time']
        utime = netcdftime.utime(tv.units)
        tvec = utime.num2date(list(tv[:]))
        if fno==0: #if it's the first file, determine the time based on the # of years to plot
            # crop the data for the time period requested
            years = np.array([tvec[ti].year for ti in range(0, len(tvec))])
            if numyears == -1:  # plot all years
                numyears = years[-1] - years[0] + 1
            years2plot = range(years[-1] + 1 - numyears, years[-1] + 1)
            ti = np.where((years >= years2plot[0]) * (years <= years2plot[-1]))
        else: #if it's not the first file, determine the tie based on the t of the previous file
            tprev=tL[fno-1]
            ti=np.where((tvec>tprev[0])*(tvec<tprev[-1]))
        t = tvec[ti[0]]
        ncL.append(nc)
        ncvL.append(ncv)
        tL.append(t)
        tiL.append(ti)
    return(ncL, ncvL, tL, tiL, fname_sid_str)

def get_varvals(ncv,varn0,avg=False,zint=[-9999,9999]):

    if '+' in varn0  or '-' in varn0 or '*' in varn0 or '/' in varn0:
        varfound,varvals,longname,units=get_varvals_op(ncv,varn0)
    else:
        if not (varn0 in ncv):
            return (False,0,0,'','')
        else:
            varn=varn0
            varvals=np.squeeze(ncv[varn][:])
            try:
                longname = ncv[varn].long_name
            except:
                longname=varn
            units=ncv[varn].units
    if len(varvals.shape)==1: #if 1-dimensional variable (e.g., airt)
        valsat='plate'
    else:
        if varn0 in ['nuh', 'nus']:
            valsat = 'int'
            depth = np.squeeze(ncv['zi'][0, :])  # depths at layer interfaces (diffusivities, fluxes, etc)
        else:
            valsat = 'center'
            if 'z' in ncv:
                depth = np.squeeze(ncv['z'][0, :])  # depths at layer centers (fabm variables, temp, salt, etc)
            elif 'depth' in ncv:
                depth = np.squeeze(ncv['depth'][:])  # depths at layer centers (fabm variables, temp, salt, etc)
        if avg:
            varvals = get_varvals_avg(varvals,depth,zint)  # zint=[0,50]
            valsat='plate'

    return (True,varvals,valsat,longname,units)

def get_varvals_avg(varvals,depth,zint):
    if all(depth>0):
        depth=-1*depth
    zi=(depth>=zint[0]) * (depth<=zint[1])
    if sum(zi)==0:
        varvalsM = varvals[:, 0] * np.nan
    else:
        varvalsM=np.mean(varvals[:, zi], axis=1)

    return (varvalsM)

def get_varvals_op(ncv,varn0):
    if '+' in varn0:
        symb='+'
        varn_set = varn0.split(symb)
        varvals = np.squeeze(ncv[varn_set[0]][:])
        units = ncv[varn_set[0]].units
        longname = '%s' %(varn_set[0])
        for varn in varn_set[1:]:
            varvals = varvals + np.squeeze(ncv[varn][:])
            longname = '%s %s %s' % (longname, symb, varn)
            #if not ncv[varn].units == units:
            #    raise(Exception('Units of variables to be added do not match: %s vs %s'%(units,ncv[varn].units)))
    else:
        if '-' in varn0:
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
        if '-' in varn0:
            varvals=v1-v2
            units=ncv[varn].units
        elif '*' in varn0:
            varvals=v1*v2
            units=ncv[varn].units + '*' + ncv[varn2].units
        elif '/' in varn0:
            varvals=v1/v2
            units=ncv[varn].units + '/' + ncv[varn2].units
            if units=='mmolC/m^2/d/mmolN/m^2/d':
                units='molC/molN'
            elif units=='mmolC/m^3/mmolN/m^3':
                units='molC/molN'
    return (True,varvals,longname,units)

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
        cmap.set_under('black')
    return (extendopt,cmap)

def format_date_axis(ax,tspan,axgrid=True):
    import matplotlib.dates as mpldates
    ax.set_xlim(tspan[0], tspan[1])
    if np.diff(tspan)[0].days<63:
        ax.xaxis.set_major_locator(mpldates.WeekdayLocator(byweekday=mpldates.MO) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%m.%d'))
    elif np.diff(tspan)[0].days<367:
        ax.xaxis.set_major_locator(mpldates.MonthLocator(bymonthday=1, interval=12) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter(''))
        ax.xaxis.set_minor_locator(mpldates.MonthLocator(bymonthday=1, interval=2) )
        ax.xaxis.set_minor_formatter(mpldates.DateFormatter('%b'))
        if axgrid:
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
        if axgrid:
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
      fnames = ['/home/onur/setups/test-BGCmodels/nflexpd/ideal_highlat_NflexPD-Nbased_Cbased/0D-Highlat_wconst_dm_Nbased.nc',
                '/home/onur/setups/test-BGCmodels/nflexpd/ideal_highlat_NflexPD-Nbased_Cbased/0D-Highlat_wconst_dm_Cbased.nc']
      print('plotting default file(s):'+'.'.join(fnames))
    else:
      print('plotting file specified:'+sys.argv[1])
      fnames=sys.argv[1].split(',')

    if len(sys.argv)<3:
      #sids = ['dm', '6h']
      sids = ['Nbased','Cbased']
    else:
      sids=sys.argv[2].split(',')

    if len(sys.argv)<4: #no third argument was passed
      #numyears=-1 # -1 means plot everything
      numyears=2
    else: 
      numyears=int(sys.argv[3]) #number of years to plot (counting from the last year backwards)

    if len(sys.argv)<5: #no third argument was passed
      modnames=['Nbased','Cbased']
    else:
      modnames=sys.argv[4].split(',')

    print('plotting last ' + str(numyears) + ' year of the simulation')

    main(fnames, sids, numyears, modnames)
