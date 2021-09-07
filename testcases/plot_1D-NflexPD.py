# terminal call:
# python plot_1D-NflexPD.py /file/path/file_name.nc num_years_to_plot(counting_backwards_from_the_last) modname
# python plot_1D-NflexPD.py ${rootdir}/${out_pref}/${out_pref}_mean.nc 3 FS-IA-DA
# python plot_1D-NflexPD.py ${rootdir}/${out_pref}/${out_pref}_mean.nc 3 phy_FS

import matplotlib.pyplot as plt
import netCDF4 as nc4
import netcdftime
import matplotlib.dates as mpldates
import sys
import numpy as np

varlims_vert={'abio_PAR_dmean':[-1.5,30],'abio_din': [-1.5, 30],'abio_detc':[0,80],'abio_detn':[0,6],'abio_detc/abio_detn':[5.0,30.0],
         'Chl':[-1.1,22.],'C':[-3.0,100.0],'N':[-0.375,7.5],'Q':[0.02,0.22],'Chl2C':[-0.005,0.08],
         'PPR':[-1.1,20.],'mu':[-0.02,0.4],'vN':[-0.0025,0.05],'f_dinphy':[0,0.5],'R_N':[-0.002,0.04],'R_Chl':[-0.0125,0.25],
         'fA':[0.0,1.0], 'fV':[-0.025,0.5], 'ThetaHat':[-0.004,0.08],'fC':[-0.05,1.0],'limfunc_Nmonod':[0,1],'limfunc_L':[-0.05,1],
         'muNET':[-0.02,0.4], 'muG':[-0.02,0.4], 'muhatNET':[-0.125, 2.5], 'muhatG':[-0.125, 2.5],'Rhat_Chl':[-0.03,0.6],}
varlims={'abio_PAR_dmean':[0,30], 'airt':[0,21], 'I_0':[0,250],'temp':[2,22], 'mld_surf':[-100,0],'wind':[-6,26],
         'abio_detc_sed/abio_detn_sed':[4.0,16.0],'abio_detc/abio_detn':[5.0,30.0],'abio_doc/abio_don':[5.0,30.0],
         'abio_din': [1, 26],'abio_detc':[0,80],'abio_detn':[0,6], 'abio_doc':[0,80], 'abio_don':[0,6],
         'Chl':[0,10.],'C':[0,50.0],'N':[0,7.5],'Q':[0.02,0.22],'Chl2C':[0.00,0.05],
         'PPR':[0,20.],'mu':[0,0.4],'vN':[0,0.05],'f_dinphy':[0,0.5],'R_N':[0,0.04],'R_Chl':[0,0.1],
         'fA':[0.0,1.0], 'fV':[0.0,0.5], 'ThetaHat':[0.00,0.05],'fC':[0,0.2],'limfunc_Nmonod':[0,0.2],'limfunc_L':[0.,1]}
#prettyunits={'abio_detc_sed/abio_detn_sed':'molC/molN','abio_detc/abio_detn':'molC/molN','abio_doc/abio_don':'molC/molN'}
prettyunits={'abio_PAR_dmean':'E\ m^{-2}\ d^{-1}','I_0':'E\ m^{-2}\ d^{-1}','wind':'m\ s^{-1}',
             'abio_din':'mmolN\ m^{-3}','C':'mmolC\ m^{-3}', 'N':'mmolN\ m^{-3}',#'abio_din':'\mu M\ N','C':'\mu M\ C', 'N':'\mu M\ N',
             'Q':'molN\ molC^{-1}','abio_detc/abio_detn':'molC\ molN^{-1}', 'abio_detc_sed/abio_detn_sed':'molC\ molN^{-1}',
             'abio_detc':'mmolC\ m^{-3}','abio_doc':'mmolC\ m^{-3}','abio_detn':'mmolN\ m^{-3}','abio_don':'mmolN\ m^{-3}',
             'abio_detc_sed':'mmolC\ m^{-2}\ d^{-1}','abio_detn_sed':'mmolN\ m^{-2}\ d^{-1}',
             'Chl':'mgChl\ m^{-3}','Chl2C':'gChl\ gC^{-3}','ThetaHat':'gChl\ gC^{-3}','fC':'-', 'fV':'-',
             'mu':'d^{-1}','muG':'d^{-1}','muNET':'d^{-1}','vN':'molN\ molC^{-1}\ d^{-1}','R_N':'d^{-1}','R_Chl':'d^{-1}',
             'muhatNET':'d^{-1}','muhatG':'d^{-1}','Rhat_Chl':'d^{-1}','vNhat':'molN\ molC^{-1}\ d^{-1}',
             'PPR':'mmolC\ m^{-3}\ d^{-1}', 'NDDR':'mmolN\ m^{-3}\ d^{-1}','f_phy_detc':'mmolC\ m^{-3}\ d^{-1}', 'f_phy_detn':'mmolN\ m^{-3}\ d^{-1}'}
prettynames={'abio_PAR_dmean':'\overline{I}','I_0':'I_{0}','mld_surf':'\mathrm{MLD}',
             'airt': 'T_{air}','temp': 'T','u10':'\mathrm{Wind \ Speed \ (-u)}','wind':'\mathrm{Wind \ Speed}',
             'abio_din':'DIN','abio_detc_sed/abio_detn_sed':'\mathrm{C:N \ of \ Det_{bot}}','abio_detc/abio_detn':'Det_C:Det_N','abio_doc/abio_don':'DOC:DON',
             'abio_detc':'det_C','abio_detn':'det_N','abio_doc':'DOC','abio_don':'DON',
             'abio_detc_sed':'\mathrm{C\ export\ rate}','abio_detn_sed':'\mathrm{N\ export\ rate}',
             'Chl':'Phy_{Chl}','C':'Phy_C','N':'Phy_N','f_dinphy':'f_{DIN-Phy}',
             'Q':'Q','vN':'V','vNhat':'\hat{V}','mu':'\mu','muNET':'\mu_{net}','muG':'\mu_{g}',
             'muhatNET':'\hat{\mu}_{net}','muhatG':'\hat{\mu}_{g}','Rhat_Chl':'\hat{R}_{Chl}',
             'R_N':'R_N','R_Chl':'R_{Chl}','fC':'f_C','limfunc_Nmonod':'L_N','limfunc_L':'L_I',
             'fA':'f_A','fV':'f_V','ThetaHat':'\hat{\Theta}','Chl2C':'\Theta',
             'PPR':'\mathrm{NPP\ rate}', 'NDDR':'\mathrm{NDD\ rate}', 'f_phy_detc':'F_{PhyC-DetC}', 'f_phy_detn':'F_{PhyN-DetN}'}
numlevels=6
#depth range to be shown:
prescylim=[0,0] #[0,0] means no ylimits are prescribed, full depth range will be shown'
#prescylim=[-30,0.0]

def main(fname, numyears, modname):
    
    if not '-' in modname: #i.e., single model runs
      models = [modname]
      varsets={
             'abio0':['airt', 'wind', 'I_0'],
             'abio1':['abio_PAR_dmean','temp', 'mld_surf'],
             'abio2':['abio_din','abio_detc/abio_detn','abio_detc_sed/abio_detn_sed'],
             'abio3':['abio_detn','abio_detc','abio_don','abio_doc'],
             'abio4':['abio_detc_sed', 'abio_detn_sed','abio_detc_sed/abio_detn_sed'],
             'phy-1':['C','N','Q','Chl','Chl2C'],
             'phy-2':['mu','vN','R_N','R_Chl'],
              'phy-3':['ThetaHat', 'fA','fV','fC','limfunc_L'],
             'phy-3b':['ThetaHat','fA','fV','limfunc_Nmonod','limfunc_L'],
              'phy-avg1': ['Q_avg0-50', 'fC_avg0-50', 'C_avg0-50'],
             'phy-avg2': ['mu_avg0-50', 'vN_avg0-50', 'R_N_avg0-50', 'R_Chl_avg0-50',
                          'mu_avg50-100', 'vN_avg50-100', 'R_N_avg50-100', 'R_Chl_avg50-100'],
             'phy_int1':['PPR_int0-100', 'f_phy_detc_int0-100', 'f_phy_detn_int0-100']
             }
    elif modname=='FS-IA-DA': #i.e., competition experiment
      #models = ['phy_IOQ', 'phy_DOQ']
      #models = ['phy_cQ','phy_IOQf', 'phy_IOQ', 'phy_DOQ', 'phy_DOQf']
      models = ['phy_FS','phy_IA', 'phy_DA'] 
      varsets={'abio0':['I_0','airt', 'wind'],
             'abio12':['temp','mld_surf','abio_din','abio_PAR_dmean',],
             'abio3':['abio_detn','abio_detc','abio_don','abio_doc'],
             'phy-1':['C','N','Q','Chl','Chl2C'],
             'phy-2':['mu','vN','R_N','R_Chl'],
             'phy-3':['ThetaHat', 'fA','fV','fC','limfunc_L'],
             'phy-3b':['ThetaHat','fA','fV','limfunc_Nmonod','limfunc_L'],
             'phy-avg1': ['Q_avg0-100', 'fC_avg0-100', 'C_avg0-100'],
             'phy-avg1SB': ['Q_avg0-50', 'fC_avg0-50', 'C_avg0-50',
                            'Q_avg50-100', 'fC_avg50-100', 'C_avg50-100'],
             'phy-avg2': ['mu_avg0-100', 'vN_avg0-100', 'R_N_avg0-100', 'R_Chl_avg0-100'],
             'phy-avg2SB': ['mu_avg0-50', 'vN_avg0-50', 'R_N_avg0-50', 'R_Chl_avg0-50',
                          'mu_avg50-100', 'vN_avg50-100', 'R_N_avg50-100', 'R_Chl_avg50-100']
             }
    elif modname=='FS-IA-DA_merged': #i.e., merged single experiments
      models = ['phy_FS','phy_IA', 'phy_DA']
      varsets = {
                #'abio-avg1SB': ['din_avg0-50', 'PAR_dmean_avg0-50', 'detc/detn_avg0-50', 'detc_sed/detn_sed',
                #                 'din_avg50-100', 'PAR_dmean_avg50-100', 'detc/detn_avg50-100', 'detc_sed/detn_sed'],
                #'abio-avg1': ['din_avg0-100', 'PAR_dmean_avg0-100', 'detc/detn_avg0-100', 'detc_sed/detn_sed'],
                #'abio-int1': ['detc_sed', 'detn_sed', 'detc_sed/detn_sed'],
                #'phy-avg1SB': ['Q_avg0-50', 'fC_avg0-50', 'C_avg0-50', 'N_avg0-50',
                #               'Q_avg50-100', 'fC_avg50-100', 'C_avg50-100', 'N_avg50-100'],
                #'phy-avg1': ['Q_avg0-100', 'fC_avg0-100', 'C_avg0-100', 'N_avg0-100'],
                'phy-avg2SB': ['mu_avg0-50', 'vN_avg0-50', 'R_N_avg0-50', 'R_Chl_avg0-50',
                               'mu_avg50-100', 'vN_avg50-100', 'R_N_avg50-100', 'R_Chl_avg50-100'],
                'phy-avg2': ['mu_avg0-100', 'vN_avg0-100', 'R_N_avg0-100', 'R_Chl_avg0-100'],
                #'phy-int1': ['PPR_int0-100', 'NDDR_int0-100', 'C_int0-100'], #, 'N_int0-100'
                #'abio-vp1': ['din_vp0-100', 'PAR_dmean_vp0-100', 'detn_vp0-100', 'detc/detn_vp0-100'],
                #'phy-vp1': ['Q_vp0-100', 'Chl_vp0-100', 'C_vp0-100', 'N_vp0-100'],
                'phy-vp2': ['mu_vp0-100', 'vN_vp0-100', 'R_N_vp0-100', 'R_Chl_vp0-100'],
                #'phy-vp3': ['muhatNET_vp0-100', 'muhatG_vp0-100', 'Rhat_Chl_vp0-100', 'limfunc_L_vp0-100'],
                #'phy-vp4': ['muNET_vp0-100', 'muG_vp0-100', 'R_Chl_vp0-100', 'limfunc_Nmonod_vp0-100'],
                #'phy-vp5': ['fC_vp0-100', 'fV_vp0-100', 'ThetaHat_vp0-100', 'Chl2C_vp0-100'],
                #'phy-vp6': ['mu_vp0-100', 'C_vp0-100', 'PPR_vp0-100', 'f_phy_detc_vp0-100'],
                'phy-vp7': ['vN_vp0-100', 'C_vp0-100', 'NDDR_vp0-100', 'f_phy_detn_vp0-100']
                 }
      
    for groupname,varset in varsets.iteritems():
        #print ('%s,%s'%(groupname,varset))
        if ('avg' in groupname) or ('int' in groupname):
            plot_multivar_ts(fname, numyears, groupname, varset, models)
        elif 'vp' in groupname:
            plot_multivar_vertprof(fname, numyears, groupname, varset, models)
        else:
            plot_singlevar(fname, numyears, groupname, varset, models)

def plot_multivar_vertprof(fname, numyears, groupname, varset, models,plotcumsum=True):
    cols = ['green', 'darkblue', 'orange']
    linestyles = ['--', ':', '-']
    linews=[2,2,2]

    dates2plot = [mpldates.datetime.datetime(2008,2,1),mpldates.datetime.datetime(2008,3,1), mpldates.datetime.datetime(2008,4,1),
                  mpldates.datetime.datetime(2008,7,1),mpldates.datetime.datetime(2008,10,15),mpldates.datetime.datetime(2008,11,15)]
    numrow = len(dates2plot)
    numcol = len(varset)
    #with var names and units for each row
    #figuresize = (1.5 + 1.5 * numcol, 1. + 2.5 * numrow)
    #fpar = {'left': 0.12, 'right': 0.85, 'top': 0.95, 'bottom': 0.05, 'hspace': 0.5, 'wspace': 0.35}
    # with var names only at first row, units only at last
    figuresize = (1.5 + 1.5 * numcol, 1. + 2.1 * numrow)
    fpar = {'left': 0.12, 'right': 0.85, 'top': 0.95, 'bottom': 0.05, 'hspace': 0.2, 'wspace': 0.35}

    nc = nc4.Dataset(fname)
    ncv = nc.variables

    # handle time
    tv = nc.variables['time']
    utime = netcdftime.utime(tv.units)
    tvec = utime.num2date(list(tv[:]))

    f = plt.figure(figsize=figuresize)
    f.subplots_adjust(left=fpar['left'], right=fpar['right'], top=fpar['top'], bottom=fpar['bottom'],
                      hspace=fpar['hspace'], wspace=fpar['wspace'])

    for k, date2plot in enumerate(dates2plot):
        for j, varn_basic in enumerate(varset):
            print(varn_basic)
            ax = plt.subplot(numrow, numcol, k*len(varset)+j + 1)

            for i, model in enumerate(models):
                if 'abio' in groupname:
                    modelpref = models[i].split('phy_')[1]
                    if '/' in varn_basic:
                        varn_basic0 = 'abio_%s_%s' % (modelpref, varn_basic.split('/')[0])
                        varn_basic1 = 'abio_%s_%s' % (modelpref, varn_basic.split('/')[1])
                        varn = '%s/%s' % (varn_basic0, varn_basic1)
                    else:
                        varn = 'abio_%s_%s' % (modelpref, varn_basic)
                else:
                    varn = '%s_%s' % (models[i], varn_basic)

                if 'vp' in varn_basic:
                    varn_basic_root = varn_basic.split('_vp')[0]
                    depthintstr = varn_basic.split('_vp')[1]  # .split.('_')[0]
                else:
                    varn_basic_root = varn_basic
                    depthintstr=''
                if 'abio' in groupname:
                    if '/' in varn_basic:
                        varn_basic0 = 'abio_%s' % (varn_basic.split('/')[0])
                        varn_basic1 = 'abio_%s' % (varn_basic.split('/')[1].split('_vp')[0])
                        varn_basic_root = '%s/%s' % (varn_basic0, varn_basic1)
                    else:
                        # if models[i]=='FS' and varn_basic=='fC':
                        #    varn_basic='limfunc_Nmonod'
                        varn_basic_root = 'abio_%s' % varn_basic_root

                varfound, dat, depths, valsat, longname, units = get_varvals(ncv, varn)

                if varfound:
                    if valsat == 'plate':
                        raise (Exception('Resulting valus have no depth dimension'))
                    # crop the data for the date and depth interval
                    ti=np.where(tvec==date2plot)
                    zfull = depths[ti,:]
                    zlimstr = varn.split('_vp')[1]
                    zint=[np.float(zlimstr.split('-')[0]),np.float(zlimstr.split('-')[1])]
                    zi=(zfull>=zint[0]) * (zfull<=zint[1])
                    z=zfull[zi]
                    datC = dat[ti,:][zi]

                    if plotcumsum:
                        datCcum = np.cumsum(dat[ti, :][zi])
                        ax.plot(datCcum, z, label=model.split('phy_')[1], color=cols[i], linestyle=linestyles[i], lw=linews[i])
                    else:
                        ax.plot(datC, z, label=model.split('phy_')[1], color=cols[i], linestyle=linestyles[i], lw=linews[i])

                if models[i]=='phy_FS' and varn_basic_root=='fC':
                    varn_basic2='limfunc_Nmonod_vp'+depthintstr
                    varn2 = '%s_%s' % (models[i], varn_basic2)

                    varfound2, dat2, depths2, valsat2, longname2, units2 = get_varvals(ncv, varn2)

                    if varfound2:
                        if valsat2 == 'plate':
                            raise (Exception('Resulting valus have no depth dimension'))
                        # crop the data for the date and depth interval
                        ti = np.where(tvec == date2plot)
                        datC2 = dat2[ti, :][zi]

                        ax.plot(datC2, z, label=model.split('phy_')[1], color=cols[i], linestyle=linestyles[i],lw=linews[i],alpha=0.5) #lw=0.5) #

            if k==0:
                if 'vp' in varn_basic:
                    if plotcumsum:
                        plt.title('$\Sigma\ %s$' % (prettynames[varn_basic_root]), size=12.0)
                    else:
                        plt.title('$%s$' % (prettynames[varn_basic_root]), size=12.0)

            # ax.invert_yaxis()
            ax.set_ylim(zint[0],zint[1])
            ax.set_ylim(ax.get_ylim()[::-1]) #revert the y axis

            if j==0:
                datestr = date2plot.strftime('%d %b')
                plt.ylabel('%s\ndepth [$m^{-1}$]'%datestr)
                #pheight=(fpar['top']-fpar['bottom'])/float(numrow)
                #y_datestr=fpar['bottom']+(pheight*(1.05))*(numrow-k-1)+pheight/2.
                #plt.figtext(0.00, y_datestr, datestr, rotation=90, verticalalignment='center')
            else:
                ax.yaxis.set_ticklabels([])
            if (k==0) and (j == numcol-1):  # (j+1)%numcol==1
                ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.7), fontsize=12)
                # for the rates plot, place legend inside to avoid avorlap y axis labs
                # ax.legend(loc='center left', bbox_to_anchor=(0.4, 0.7),fontsize=12)

            ax.grid(b=True, axis='y', which='major', color='0.5', linestyle='-')
            # x-axis
            if varn_basic_root in varlims_vert.keys() and not plotcumsum:
                ax.set_xlim(varlims_vert[varn_basic_root][0], varlims_vert[varn_basic_root][1])

            if k+1==numrow:
                if varn_basic_root in prettyunits:
                    units = prettyunits[varn_basic_root]
                    plt.xlabel('[$%s$]'%units, size=12.0)

    nc.close()
    if plotcumsum:
        figname = fname.split('.nc')[0] + '_vertprof_' + groupname + '_' + str(numyears) + 'y_cumsum.png'
    else:
        figname = fname.split('.nc')[0] + '_vertprof_' + groupname + '_' + str(numyears) + 'y.png'
    plt.savefig(figname)
    print('python vertical profile plot saved in: ' + figname)

def plot_multivar_ts(fname, numyears, groupname, varset, models):
    cols=['green','darkblue','orange']
    linestyles=['--',':','-']
    linews=[2,2,2]

    if len(varset) > 3:
        if len(varset) in [3, 6, 9]:
            numcol = 3.0
        elif len(varset) in [4, 8, 12]:
            numcol = 4.0
    else:
        numcol = len(varset)
    numrow = np.ceil(len(varset) / numcol)
    figuresize = (1 + 3.5*numcol, 1. + 1.5*numrow)
    fpar = {'left': 0.04, 'right': 0.92, 'top': 0.93, 'bottom': 0.07, 'hspace': 0.5, 'wspace': 0.3}
    if numrow == 1:
        fpar['top']=0.9; fpar['bottom']=0.1

    nc = nc4.Dataset(fname)
    ncv = nc.variables

    #handle time
    tv = nc.variables['time']
    utime = netcdftime.utime(tv.units)
    tvec = utime.num2date(list(tv[:]))
    # crop the data for the time period requested
    years = np.array([tvec[ti].year for ti in range(0, len(tvec))])
    if numyears == -1:  # plot all years
        numyears = years[-1] - years[0] + 1
    years2plot = range(years[-1] + 1 - numyears, years[-1] + 1)
    yeari = np.where((years >= years2plot[0]) * (years <= years2plot[-1]))
    t= tvec[yeari[0]]

    f = plt.figure(figsize=figuresize)
    f.subplots_adjust(left=fpar['left'], right=fpar['right'], top=fpar['top'], bottom=fpar['bottom'],
                      hspace=fpar['hspace'], wspace=fpar['wspace'])

    for j, varn_basic in enumerate(varset):
        print(varn_basic)
        ax = plt.subplot(numrow, numcol, j + 1)

        for i,model in enumerate(models):
            if 'abio' in groupname:
                modelpref=models[i].split('phy_')[1]
                if '/' in varn_basic:
                    varn_basic0='abio_%s_%s' % (modelpref, varn_basic.split('/')[0])
                    varn_basic1 = 'abio_%s_%s' % (modelpref, varn_basic.split('/')[1])
                    varn='%s/%s'%(varn_basic0,varn_basic1)
                else:
                    varn = 'abio_%s_%s' % (modelpref, varn_basic)
            else:
                #if models[i]=='FS' and varn_basic=='fC':
                #    varn_basic='limfunc_Nmonod'
                varn = '%s_%s'%(models[i], varn_basic)
            if varn_basic in prettyunits:
                units = prettyunits[varn_basic]
            if 'avg' in varn_basic:
                varn_basic_root = varn_basic.split('_avg')[0]
                depthintstr = varn_basic.split('_avg')[1]
            elif 'int' in varn_basic:
                varn_basic_root = varn_basic.split('_int')[0]
                depthintstr = varn_basic.split('_int')[1]
            else:
                varn_basic_root = varn_basic
                depthintstr = ''
            if 'abio' in groupname:
                if '/' in varn_basic:
                    varn_basic0 = 'abio_%s' % (varn_basic.split('/')[0])
                    if 'avg' in varn_basic:
                        varn_basic1 = 'abio_%s' % (varn_basic.split('/')[1].split('_avg')[0])
                    elif 'int' in varn_basic:
                        varn_basic1 = 'abio_%s' % (varn_basic.split('/')[1].split('_int')[0])
                    else:
                        varn_basic1='abio_%s' % (varn_basic.split('/')[1])
                    varn_basic_root = '%s/%s' % (varn_basic0, varn_basic1)
                else:
                    varn_basic_root = 'abio_%s' % varn_basic_root

            varfound, dat, depths, valsat, longname, units = get_varvals(ncv, varn)

            if valsat == 'plate':
                # crop the data for the time period requested
                datC = dat[yeari[0]]
            else:
                raise(Exception('Resulting valus are not 1-dimensional'))
            ax.plot(t, datC, label=model.split('phy_')[1], color=cols[i],linestyle=linestyles[i], lw=linews[i])

            if models[i] == 'phy_FS' and varn_basic_root == 'fC':
                varn_basic2 = 'limfunc_Nmonod_avg'+depthintstr
                varn2 = '%s_%s' % (models[i], varn_basic2)

                varfound2, dat2, depths2, valsat2, longname2, units2 = get_varvals(ncv, varn2)

                if varfound2:
                    if valsat == 'plate':
                        # crop the data for the time period requested
                        datC2 = dat2[yeari[0]]
                    else:
                        raise (Exception('Resulting valus are not 1-dimensional'))
                    ax.plot(t, datC2, label=model.split('phy_')[1], color=cols[i], linestyle=linestyles[i], lw=linews[i], alpha=0.5)  # lw=0.5) #

        if varn_basic_root in prettyunits:
            units = prettyunits[varn_basic_root]
        if 'avg' in varn_basic:
            plt.title('%sm mean $%s$ [$%s$]' % (depthintstr, prettynames[varn_basic_root], units), size=12.0)
        elif 'int' in varn_basic:
            plt.title('%sm int. $%s$ [$%s$]' % (depthintstr, prettynames[varn_basic_root], units), size=12.0)
        else:
            plt.title('$%s$ [$%s$]' % (prettynames[varn_basic_root], units), size=12.0)

        # shrink the axes width by 20% to fit that of the contour plots, and put the legend in that space
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, (box.x1 - box.x0) * 0.8, box.y1 - box.y0])
        if varn_basic in varlims.keys():
            ax.set_ylim(varlims[varn][0], varlims[varn][1])

        if j+1==numcol: #(j+1)%numcol==1
            ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.8),fontsize=12)
            #for the rates plot, place legend inside to avoid avorlap y axis labs
            #ax.legend(loc='center left', bbox_to_anchor=(0.4, 0.7),fontsize=12)

        ax.grid(b=True, axis='y', which='major', color='0.5', linestyle='-')
        # x-axis
        format_date_axis(ax, [t[0], t[-1]])
        plt.xlabel('')

    nc.close()
    figname = fname.split('.nc')[0] + '_line_' + groupname + '_' + str(numyears) + 'y.png'
    plt.savefig(figname)
    print('python line plot saved in: ' + figname)

def plot_singlevar(fname,numyears,groupname,varset,models):
    colmap='viridis'
    axgrid=True
    #default:
    fpar = {'left': 0.05, 'right': 0.99, 'top': 0.85, 'bottom': 0.15, 'hspace': 0.5, 'wspace': 0.2}
    if 'phy' in groupname:
        if len(models)==1:
            if len(varset)>5:
                numrow = 5.0
            else:
                numrow = len(varset)
            numcol = len(varset)/numrow
            figuresize = (0.25 + 4 * numcol, .5 + 1.5 * numrow)
            varnames = []
            for var in varset:
                varnames.append('%s_%s' % (models[0], var))
            if numrow!=1:
                fpar = {'left':0.15, 'right':0.99, 'top': 0.95, 'bottom': 0.05, 'hspace': 0.5, 'wspace': 0.2}
        elif len(models)==2: #show the difference between models in the 3rd column
            numcol = 3.0
            numrow = len(varset)
            figuresize = (1 + 4 * numcol, 1 + 1.5 * numrow)
            varnames = []
            for var in varset:
                varnames.append('%s_%s' % (models[0], var))
                varnames.append('%s_%s' % (models[1], var))
                varnames.append('%s_%s-%s_%s' % (models[0], var, models[1], var))
            fpar = {'left':0.1, 'right':0.9, 'top': 0.95, 'bottom': 0.05, 'hspace': 0.5, 'wspace': 0.2}
        elif len(models)>2:
            numcol = len(models)
            numrow = len(varset)
            figuresize = (1 + 4 * numcol, .5 + 1.5 * numrow)
            varnames = []
            for var in varset:
                for i in range(len(models)):
                    varnames.append('%s_%s' % (models[i], var))
            if len(varset)==2:
                fpar = {'left': 0.05, 'right': 0.99, 'top': 0.93, 'bottom': 0.07, 'hspace': 0.5, 'wspace': 0.2}
            elif len(varset)>2:
                fpar = {'left':0.05, 'right':0.99, 'top': 0.95, 'bottom': 0.05, 'hspace': 0.5, 'wspace': 0.2}
    else:
        varnames = varset
        if len(models)==1:
            if len(varset)>5:
                numrow = 5.0
            else:
                numrow = len(varset)
            numcol = len(varset)/numrow
            figuresize = (0.25 + 4 * numcol, .5 + 1.5 * numrow)
            if numrow!=1:
                fpar = {'left':0.15, 'right':0.99, 'top': 0.95, 'bottom': 0.05, 'hspace': 0.5, 'wspace': 0.2}
        else:
            if len(varset)>4:
                numcol = 3.0
            else:
                numcol=len(varset)
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
    z=np.squeeze(ncv['z'][0,:]) #depths at layer centers (fabm variables, temp, salt, etc)
    zi=np.squeeze(ncv['zi'][0,:]) #depths at layer interfaces (diffusivities, fluxes, etc)

    f=plt.figure(figsize=figuresize)
    f.subplots_adjust(left=fpar['left'],right=fpar['right'],top=fpar['top'],bottom=fpar['bottom'],hspace=fpar['hspace'],wspace=fpar['wspace'])

    numvar=len(varnames)
    ylimsuf = ''
    for i,varn in enumerate(varnames):       
        print (varn)
        varn_basic,model=get_basic_varname(varn,models)
        
        ax=plt.subplot(numrow,numcol,i+1)

        if (varn == 'skip'):
            continue

        if varn=='wind':
            varfound, datU, depths, valsat, longname, units = get_varvals(ncv, 'u10')
            varfound, datV, depths, valsat, longname, units = get_varvals(ncv, 'v10')
            dat=np.sqrt(datU**2+datV**2)
            #restore negative velocities
            ineg=datU<0.0
            dat[ineg]=dat[ineg]*-1.0
        else:
            varfound, dat, depths, valsat, longname, units = get_varvals(ncv, varn)

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
            if varn_basic in varlims.keys():
                ax.set_ylim(varlims[varn][0],varlims[varn][1])
        else:
            if ((np.max(datC)+9.9<0.01) & (np.min(datC)+9.9<0.01)):
                datC[:,:]=np.nan
            #coloring of values exactly 0.0 becomes arbitrary. Tip them over the positive side to prevent alarming results
            datC[datC==0.0]=1e-15
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
        cmap.set_under('black')
    return (extendopt,cmap)

def get_basic_varname(varn,models):
    if len(models)>1:
        for model in models:
            if model in varn:
                varn_basic=varn.split(model+'_')[1]
                return (varn_basic, model)
        return (varn, '')
    else:
        model=models[0]
        if model in varn:
            varn_basic = varn.split(model + '_')[1]
            return (varn_basic, model)
        else:
            return (varn,model)

def get_varvals(ncv,varn0):
    if 'avg' in varn0:
        varfound, varvals,depths, longname, units = get_varvals_avg(ncv, varn0,'avg')
    elif 'int' in varn0:
        varfound, varvals,depths, longname, units = get_varvals_avg(ncv, varn0,'int')
    elif 'vp' in varn0:
        varfound, varvals,depths,longname, units = get_varvals_vp(ncv, varn0)
    elif '+' in varn0  or '-' in varn0 or '*' in varn0 or '/' in varn0:
        varfound,varvals,depths,longname,units=get_varvals_op(ncv,varn0)
    else:
        if not (varn0 in ncv):
            return (False,0,0,'','','')
        else:
            varn=varn0
            varvals=np.squeeze(ncv[varn][:])
            depths=np.squeeze(ncv['z'][:])*-1
            longname = ncv[varn].long_name
            units=ncv[varn].units
            varfound=True
    if varfound:
        if len(varvals.shape)==1: #if 1-dimensional variable (e.g., airt)
            valsat='plate'
        elif varn0 in ['nuh', 'nus']:
            valsat='int'
        else:
            valsat='center'
    else:
        valsat='NA'
    return (varfound,varvals,depths,valsat,longname,units)

def get_varvals_vp(ncv,varn0):
    varn=varn0.split('_vp')[0]
    varfound, varvals, depths, valsat,longname, units = get_varvals(ncv, varn)
    #assume that zintstr=0-X means below 0, and above Xm depth
    #zlimstr = varn0.split('_vp')[1]
    #zint=[np.float(zlimstr.split('-')[0]),np.float(zlimstr.split('-')[1])]
    #zi=(depths>=zint[0]) * (depths<=zint[1])
    #varvals=varvals[:,]
    #depths=depths[zi]

    return (varfound,varvals,depths,longname,units)

def get_varvals_avg(ncv,varn0,mode):
    if mode=='avg':
        zintstr = varn0.split('_avg')[1]
        varn = varn0.split('_avg')[0]
    elif mode=='int':
        zintstr = varn0.split('_int')[1]
        varn = varn0.split('_int')[0]
    else:
        raise(Exception('Unknown mode:%s'%mode))
    varfound, varvals, depths, valsat,longname, units = get_varvals(ncv, varn)
    #assume that zintstr=0-X means below 0, and above Xm depth
    zint=[np.float(zintstr.split('-')[0]),np.float(zintstr.split('-')[1])]
    zi=(depths>=zint[0]) * (depths<=zint[1])
    #vals=np.squeeze(ncv[varn])
    valsM=np.ma.array(varvals,mask=np.invert(zi))
    depthsM = np.ma.array(depths, mask=np.invert(zi))
    if mode in ['avg','int']:
        heights = np.squeeze(ncv['h'][:, :])  # i.e., layer thicknesses
        heightsM = np.ma.array(heights, mask=np.invert(zi))
        vals_massM = valsM * heightsM  # X/m3*m=X/m2
        vals_massM_sum = np.sum(vals_massM, axis=1)
        if mode =='avg':
            varvals = vals_massM_sum/np.sum(heightsM,axis=1) #X/m2/m = X/m3
        else:
            varvals=vals_massM_sum
            units = units.replace('m^3', 'm^2')
    depths = np.mean(depthsM, axis=1)
    return (varfound,varvals,depths,longname,units)

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
        if units=='mmolC/m^2/d/mmolN/m^2/d':
            units='molC/molN'
        elif units=='mmolC/m^3/mmolN/m^3':
            units='molC/molN'
    depths = np.squeeze(ncv['z'][:, :]) * -1
    return (True,varvals,depths,longname,units)

def format_date_axis(ax,tspan,axgrid=True):
    ax.set_xlim(tspan[0], tspan[1])
    if np.diff(tspan)[0].days<63:
        ax.xaxis.set_major_locator(mpldates.WeekdayLocator(byweekday=mpldates.MO) )
        ax.xaxis.set_major_formatter(mpldates.DateFormatter('%b\n%d'))
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
      #fname = '/home/onur/setups/test-BGCmodels/nflexpd/1D-ideal-highlat-MS/21-07-29/Highlat-100m_wconst-FS/Highlat-100m_wconst-FS_mean.nc'
      fname = '/home/onur/setups/test-BGCmodels/nflexpd/1D-ideal-highlat-MS/21-07-29/Highlat-100m_wconst_FS-IA-DA_merged/Highlat-100m_wconst_FS-IA-DA_merged_mean.nc'
      print('plotting default file:'+fname)
    else:
      print('plotting file specified:'+sys.argv[1])
      fname=sys.argv[1]
      
    if len(sys.argv)<3: #no third argument was passed
      #numyears=-1 # -1 means plot everything
      numyears=1
    else: 
      numyears=int(sys.argv[2]) #number of years to plot (counting from the last year backwards)
    print('plotting last '+str(numyears)+' year of the simulation')
    
    if len(sys.argv)<4:
      modname='FS-IA-DA_merged'
      #modname = 'phy_FS'
    else:
      modname=sys.argv[3]
    main(fname, numyears, modname)
