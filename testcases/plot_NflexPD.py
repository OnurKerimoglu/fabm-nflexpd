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
         'DIC':[0,1000],'DIN': [1, 26],'DetC':[0,80],'DetN':[0,6], 'DOC':[0,80], 'DON':[0,6],
         'Phy-Chl':[0,10.],'Phy-Q':[0.02,0.22],'Phy-Chl2C':[0.00,0.05],
         'Phy-C':[0,50.0],'Phy-C 1':[0,50.0],'Phy-C 2':[0,50.0],
         'Phy-N':[0,7.5],'Phy-N 1':[0,7.5], 'Phy-N 2':[0,7.5],
         'Phy-PPR':[0,20.],'Phy-mu':[0,0.4],'Phy-vN':[0,0.05],'Phy-f_dinphy':[0,0.5],'Phy-R_N':[0,0.04],'Phy-R_Chl':[0,0.1],
         'Phy-fA':[0.0,1.0], 'Phy-fV':[0.0,0.5], 'Phy-ThetaHat':[0.00,0.05],'Phy-fC':[0,0.2],'Phy-limfunc_L':[0.,1]}
varlims0D={'I_dm':[0,30], 'airt':[0,21], 'I_0':[0,250],'temp':[2,22], 'mld_surf':[-100,0],'wind':[-6,26],
         'DetC_sed/DetN_sed':[4.0,16.0],'DetC/DetN':[5.0,30.0],'DOC/DON':[5.0,30.0],
         'DIC': [0, 1050],'Global C':[1040.54,1040.56],'totalC':[800,1100], #,'totalC':[1040.54,1040.56],
         'DIN': [0, 25],'Global N':[21.095,21.125],'totalN':[12,22], #'totalN':[21.095,21.125], #'totalN':[21.1113,21.1115], #'totalN':[21.095,21.125], #[21.199,21.201],
         #'DIN': [0, 40],'DetC':[0,400],'DetN':[0,40], 'DOC':[0,400], 'DON':[0,40],'totalN':[36.195,36.22],
         #'Phy-C':[0,100.0],'Phy-C 1':[0,100.0],'Phy-C 2':[0,100.0],'Phy-N':[0,10],'Phy-N 1':[0,10.0], 'Phy-N 2':[0,10.0],
         'DetC':[0,300],'DetN':[0,20], 'DOC':[0,300], 'DON':[0,20],
         #'Phy-C':[0,80.0],'Phy-C 1':[0,70.0],'Phy-C 2':[0,70.0],'Phy-N':[0,8],'Phy-N 1':[0,7.0], 'Phy-N 2':[0,7.0],
         'Phy-C':[0,100.0],'Phy-C 1':[0,70.0],'Phy-C 2':[0,70.0],'Phy-N':[0,10],'Phy-N 1':[0,7.0], 'Phy-N 2':[0,7.0],
         'Phy-Chl':[0,20.],'Phy-Q':[0.02,0.22],'Phy-Chl2C':[0.00,0.1],
         'Phy-PPR':[0,20.],'Phy-mu':[0,0.4],'Phy-vN':[0,0.05],'Phy-f_dinphy':[0,0.5],'Phy-R_N':[0,0.04],'Phy-R_Chl':[0,0.1],
         'Phy-fA':[0.0,1.0], 'Phy-fV':[0.0,0.5], 'Phy-ThetaHat':[0.00,0.05],'Phy-fC':[0,0.2],'Phy-limfunc_L':[0.,1],
         'dQ_dt': [-0.01,0.01],'delQ_delt_N': [-0.01,0.01],'delQ_delt_I': [-0.01,0.01],'delQ_delt_Ld': [-0.01,0.01],'delQ_delt_T': [-0.01,0.01],
           }
namelibNbasedIA={'I_0':'I_0','wind':'m\ s^{-1}','T':'temp',
                 'totalN':'total_nitrogen_calculator_result',
               'I-dm':'abio_PAR_dmean','I':'abio_PAR',
             'Phy-C 1':'phy_IA1_C','Phy-N 1':'phy_IA1_N','PhyC-Q 1':'phy_IA1_Q',
             'Phy-C 2':'phy_IA2_C','Phy-N 2':'phy_IA2_N','PhyC-Q 2':'phy_IA2_Q',
             'Phy-C':'phy_IA1_C','Phy-N':'phy_IA1_N','PhyC-Q':'phy_IA1_Q',
             'Phy-Chl':'phy_IA1_Chl','PhyC-Chl2C':'phy_IA1_Chl2C',
             'mu':'phy_IA1_mu','vN':'phy_IA1_vN','R_N':'phy_IA1_R_N','R_Chl':'phy_IA1_R_Chl',
             'fA':'phy_IA1_fA', 'fV':'phy_IA1_fV', 'fC':'phy_IA1_fC',
             'ThetaHat':'phy_IA1_ThetaHat', 'limfunc_L':'phy_IA1_limfunc_L',
             'DIN':'abio_din','DON':'abio_don','DetN':'abio_detn',
             'DOC':'abio_doc','DetC':'abio_detc'
            }
namelibNbasedDA={'I_0':'I_0','wind':'m\ s^{-1}','T':'temp',
                 'totalN':'total_nitrogen_calculator_result',
               'I-dm':'abio_PAR_dmean','I':'abio_PAR',
             'Phy-C 1':'phy_DA1_C','Phy-N 1':'phy_DA1_N','Phy-Q 1':'phy_DA1_Q',
             'Phy-C 2':'phy_DA2_C','Phy-N 2':'phy_DA2_N','Phy-Q 2':'phy_DA2_Q',
             'Phy-C':'phy_DA1_C','Phy-N':'phy_DA1_N','Phy-Q':'phy_DA1_Q',
             'Phy-Chl':'phy_DA1_Chl','Phy-Chl2C':'phy_DA1_Chl2C',
             'mu':'phy_DA1_mu','vN':'phy_DA1_vN','R_N':'phy_DA1_R_N','R_Chl':'phy_DA1_R_Chl',
             'fA':'phy_DA1_fA', 'fV':'phy_DA1_fV', 'fC':'phy_DA1_fC',
             'ThetaHat':'phy_DA1_ThetaHat', 'limfunc_L':'phy_DA1_limfunc_L',
             'DIN':'abio_din','DON':'abio_don','DetN':'abio_detn',
             'DOC':'abio_doc','DetC':'abio_detc'
            }
namelibCbasedIA={'I_0':'I_0','wind':'m\ s^{-1}','T':'temp', 'L$_D$':'abio_Cbased_FDL',
             'Global C': 'total_carbon_calculator_result',
             'Global N': 'total_nitrogen_calculator_result',
             'totalC':'total_carbon_calculator_result-abio_Cbased_extc',
             'totalN':'total_nitrogen_calculator_result-abio_Cbased_extn',
             #'totalN':'abio_Cbased_din+abio_Cbased_don+abio_Cbased_detn+phy_Cbased_IA_C*phy_Cbased_IA_Q', #+phy_Cbased_IA_N',
             'I-dm':'abio_Cbased_PAR_dmean','I':'abio_Cbased_PAR', 'dI_dt':'phy_Cbased_IA1_dI_dt',
             'Phy-C 1':'phy_Cbased_IA1_C','Phy-N 1':'phy_Cbased_IA1_N','PhyC-Q 1':'phy_Cbased_IA1_Q',
             'Phy-C 2':'phy_Cbased_IA2_C','Phy-N 2':'phy_Cbased_IA2_N','PhyC-Q 2':'phy_Cbased_IA2_Q',
             'Phy-C':'phy_Cbased_IA1_C','Phy-N':'phy_Cbased_IA1_N','Phy-Q':'phy_Cbased_IA1_Q',
             'Phy-Chl':'phy_Cbased_IA_Chl','Phy-Chl2C':'phy_Cbased_IA_Chl2C',
             'mu':'phy_Cbased_IA_mu','vN':'phy_Cbased_IA_vN','R_N':'phy_Cbased_IA_R_N','R_Chl':'phy_Cbased_IA_R_Chl',
             'fA':'phy_Cbased_IA_fA', 'fV':'phy_Cbased_IA_fV', 'fC':'phy_Cbased_IA_fC',
             'ThetaHat':'phy_Cbased_IA_ThetaHat', 'limfunc_L':'phy_Cbased_IA_limfunc_L',
             'DIC':'abio_Cbased_dic','DIN':'abio_Cbased_din','DON':'abio_Cbased_don','DetN':'abio_Cbased_detn',
             'DOC':'abio_Cbased_doc','DetC':'abio_Cbased_detc',
             'dQ_dt':'phy_Cbased_IA1_dQdt',
             'delQ_delt_N':'phy_Cbased_IA1_delQdelt_N',
             'delQ_delt_I':'phy_Cbased_IA1_delQdelt_I',
             'delQ_delt_Ld':'phy_Cbased_IA1_delQdelt_Ld',
             'delQ_delt_T':'phy_Cbased_IA1_delQdelt_T'
            }
namelibCbasedDA={'I_0':'I_0','wind':'m\ s^{-1}','T':'temp',
             'Global C': 'total_carbon_calculator_result',
             'Global N': 'total_nitrogen_calculator_result',
            'totalC':'total_carbon_calculator_result-abio_Cbased_extc',
             'totalN':'total_nitrogen_calculator_result-abio_Cbased_extn',
             #'totalN':'abio_Cbased_din+abio_Cbased_don+abio_Cbased_detn+phy_Cbased_DA_N',
             'I-dm':'abio_Cbased_PAR_dmean','I':'abio_Cbased_PAR', 'dI_dt':'phy_Cbased_IA1_dI_dt',
             'Phy-C 1':'phy_Cbased_DA1_C','Phy-N 1':'phy_Cbased_DA1_N','Phy-Q 1':'phy_Cbased_DA1_Q',
             'Phy-C 2':'phy_Cbased_DA2_C','Phy-N 2':'phy_Cbased_DA2_N','Phy-Q 2':'phy_Cbased_DA2_Q',
             'Phy-C':'phy_Cbased_DA1_C','Phy-N':'phy_Cbased_DA1_N','Phy-Q':'phy_Cbased_DA1_Q',
             'Phy-Chl':'phy_Cbased_DA_Chl','Phy-Chl2C':'phy_Cbased_DA_Chl2C',
             'mu':'phy_Cbased_DA_mu','vN':'phy_Cbased_DA_vN','R_N':'phy_Cbased_DA_R_N','R_Chl':'phy_Cbased_DA_R_Chl',
             'fA':'phy_Cbased_DA_fA', 'fV':'phy_Cbased_DA_fV', 'fC':'phy_Cbased_DA_fC',
             'ThetaHat':'phy_Cbased_DA_ThetaHat', 'limfunc_L':'phy_Cbased_DA_limfunc_L',
             'DIC':'abio_Cbased_dic','DIN':'abio_Cbased_din','DON':'abio_Cbased_don','DetN':'abio_Cbased_detn',
             'DOC':'abio_Cbased_doc','DetC':'abio_Cbased_detc'
            }
prettyunits={'I_0':'E\ m^{-2}\ d^{-1}','wind':'m\ s^{-1}','T':'^\circ C',
             'I-dm':'E\ m^{-2}\ d^{-1}','I':'E\ m^{-2}\ d^{-1}','dI_dt':'E\ m^{-2}\ d^{-2}',
             'Phy-C':'mmolC\ m^{-3}','Phy-C 1':'mmolC\ m^{-3}', 'Phy-C 2':'mmolC\ m^{-3}',
             'Phy-N':'mmolN\ m^{-3}','Phy-N 1':'mmolN\ m^{-3}', 'Phy-N 2':'mmolN\ m^{-3}',
             'Phy-Q':'molN\ molC^{-1}','Phy-Chl':'mg m^{-3}','Phy-Chl2C':'gChl\ gC^{-3}',
             'DIN':'mmolN\ m^{-3}','DON':'mmolN\ m^{-3}','DetN':'mmolN\ m^{-3}', 'Global C':'mmolC\ m^{-3}','totalN':'mmolN\ m^{-3}',
             'DIC':'mmolC\ m^{-3}','DOC':'mmolC\ m^{-3}','DetC':'mmolC\ m^{-3}', 'Global N':'mmolN\ m^{-3}', 'totalC':'mmolC\ m^{-3}',
              'mu':'d^{-1}', 'vN':'molN\ molC^{-1} d^{-1}', 'R_N':'d^{-1}', 'R_Chl':'d^{-1}',
             'fA':'-', 'fV':'-', 'fC':'-', 'limfunc_L':'-',
             'dQ_dt':'molN\ molC^{-1}\ d^{-1}',
             'delQ_delt_N':'molN\ molC^{-1}\ d^{-1}',
             'delQ_delt_I':'molN\ molC^{-1}\ d^{-1}',
             'delQ_delt_Ld':'molN\ molC^{-1}\ d^{-1}',
             'delQ_delt_T':'molN\ molC^{-1}\ d^{-1}'
             }
prettynames={'I_0':'$I_0$','wind':'wind','T':'T',
             'I-dm':r'$\bar{I}$','I':'$I$','dI_dt':r'd$\bar{I}$/d$t$',
             'totalN':'Total N','totalC':'Total C',
             'Phy-C':r'$Phy_C$','Phy-C 1':r'$Phy_C^1$','Phy-C 2':r'$Phy_C^2$',
             'Phy-N':'$Phy_N$','Phy-N 1':r'$Phy_N^1$','Phy-N 2':r'$Phy_N^2$',
             'Phy-Q':'$Q$','Phy-Chl':'$Phy_{Chl}$','Phy-Chl2C':r'$\theta$',
             'DIN':'DIN','DON':'DON','DetN':'$Det_N$',
             'DOC':'DOC','DetC':'$Det_C$',
              'mu':'$\mu$', 'vN':'$V_N$', 'R_N':'$R_N$', 'R_Chl':'$R_{Chl}$',
             'fA':'$f_A$', 'fV':'$f_V$', 'fC':'$f_C$', 'limfunc_L':'$L_I$','ThetaHat':r'$\hat{\theta}$',
             'dQ_dt':r'$\frac{\mathrm{d}Q}{\mathrm{d}t}$',
             'delQ_delt_N':r'$\frac{\partial Q}{\partial \mathrm{DIN}} \frac{\mathrm{d} \mathrm{DIN}}{\mathrm{d} t}$',
             'delQ_delt_I':r'$\frac{\partial Q}{\partial \bar{I}} \frac{\mathrm{d} \bar{I}}{\mathrm{d} t}$',
             'delQ_delt_Ld':r'$\frac{\partial Q}{\partial \mathrm{L}_\mathrm{D}} \frac{\mathrm{d} \mathrm{L}_\mathrm{D}}{\mathrm{d} t}$',
             'delQ_delt_T':r'$\frac{\partial Q}{\partial T} \frac{\mathrm{d} \mathrm{T}}{\mathrm{d} t}$',
             }
numlevels=6
#depth range to be shown:
prescylim=[0,0] #[0,0] means no ylimits are prescribed, full depth range will be shown'
#prescylim=[-30,0.0]

def main(fnames, numyears, modnames, variants, ids):

  varsets={#'abio0':['airt', 'wind', 'I_0'], #I_dm
           #'abio1':['abio_PAR_dmean','temp', 'mld_surf'],
           #'abio2':['abio_din','abio_detc/abio_detn','abio_detc_sed/abio_detn_sed'],
           #'abio3':['abio_detn','abio_detc','abio_don','abio_doc'],
           #'abio0':['I-dm', 'dI_dt','T','L$_D$'],
           'abio1':[ 'Global C','Global N',
                     'totalC','totalN',
                     'Phy-C', 'Phy-N',
                     #'Phy-C 1', 'Phy-N 1',
                     #'Phy-C 2', 'Phy-N 2',
                     'DIC','DIN',
                     'DOC','DON',
                     'DetC', 'DetN'],
           #'phy-1':['Phy-C','Phy-N','Phy-Q',
           #          '',     'Phy-Chl','Phy-Chl2C'],
           #'phy-2': ['mu', 'vN', 'R_N', 'R_Chl'],
           #'phy-3': ['fA', 'fV', 'fC',
           #          'ThetaHat', 'limfunc_L',''],
           #'phy-4':['dQ_dt','',
           #         'delQ_delt_N','delQ_delt_I',
           #         'delQ_delt_Ld','delQ_delt_T']
         }
  for groupname,varset in varsets.iteritems():
      #print ('%s,%s'%(groupname,varset))
      #if len(variants)>1:
      plot_multifile(fnames, numyears, groupname, varset, variants,modnames,ids)
      #if not (len(modnames)==2 and 'OBS' in modnames):
      #    for variantno in range(0,len(variants)):
      #        plot_singlefile(fnames[variantno], numyears, groupname, varset, variants[variantno],modnames[variantno])

def plot_multifile(fnames, numyears, groupname, varset, variants, modnames, ids):
    cols=['darkblue','orange','green']
    linestyles=['-',':','--']
    varlims=varlims0D
    if len(varset)>3:
        if len(varset) in [4,6,8,10, 12]:
            numcol = 2.0
        elif len(varset) in [3,9]:
            numcol=3.0
        elif len(varset) in [12]:
            numcol=4.0
    else:
        numcol = len(varset)
    numrow = np.ceil(len(varset) / numcol)
    figuresize = (1 + 4*numcol, 1. + 1.5*numrow)
    fpar = {'left': 0.1, 'right': 0.99, 'top': 0.93, 'bottom': 0.08, 'hspace': 0.5, 'wspace': 0.2}
    if numrow == 1:
        fpar['top']=0.9; fpar['bottom']=0.13

    ncL,ncvL,tL,tiL,fname_str = collect_data(fnames,modnames,numyears)

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
        for i,variant in enumerate(variants):
            ncv=ncvL[i];t=tL[i];ti=tiL[i]
            if modnames[i] == 'Nbased':
                if variant == 'IA':
                    namelib = namelibNbasedIA
                elif variant == 'DA':
                    namelib = namelibNbasedDA
            elif modnames[i] == 'Cbased':
                if variant == 'IA':
                    namelib = namelibCbasedIA
                elif variant == 'DA':
                    namelib = namelibCbasedDA
            if varn in namelib:
                fvarn=namelib[varn]
            else:
                fvarn=varn
            varfound, dat, valsat, longname, unitsnew = get_varvals(ncv, fvarn, avg=True) #,zint=[-100,0])
            #update the units only if the units of the new data set exist
            if units=='' and unitsnew!='':
                units=unitsnew
            if varfound:
                datC = dat[ti[0]]
                if modnames[i]=='OBS':
                    #ax.plot(t, datC, label=variant, color=cols[i], linestyle='',marker='o', markersize=4)
                    ax.plot(t, datC, label=ids[i], color=cols[i], linestyle='', marker='o', markersize=4)
                else:
                    #ax.plot(t, datC, label=variant, color=cols[i],linestyle=linestyles[i])
                    ax.plot(t, datC, label=ids[i], color=cols[i], linestyle=linestyles[i])
                    #if varn in ['totalN','totalC']:
                    if varn in ['Global N', 'Global C']:
                        ax.text(0.6,0.6-i*0.2,r'$\delta$(%s):%.2e'%(ids[i],max(datC[1:])-min(datC[1:])), transform=ax.transAxes) #,color=cols[i]
        if varn in prettyunits:
            prettyunit = prettyunits[varn]
        else:
            prettyunit = units
        if prettyunit=='':
            prettyunit='-'
        if varn in prettynames:
            prettyname = prettynames[varn]
        else:
            prettyname = varn
        plt.title('%s [$%s$]' % (prettyname, prettyunit), size=12.0)

        if varn in varlims.keys():
            ax.set_ylim(varlims[varn][0], varlims[varn][1])
        if not (prescylim[0] == 0 and prescylim[1] == 0):
            ax.set_ylim(prescylim[0], prescylim[1])
            # ylimsuf='_ylim_%s-%s'%(prescylim[0],prescylim[1])

        ax.ticklabel_format(axis='y', useOffset=False) #style='plain',) #

        # shrink the axes width by 20% to fit that of the contour plots, and put the legend in that space
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, (box.x1 - box.x0) * 0.8, box.y1 - box.y0])
        if j==0: #(j+1)%numcol==1
            #ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.8),fontsize=12)
            ax.legend(loc='center left', bbox_to_anchor=(0.05, 0.7),fontsize=12)

        ax.grid(b=True, axis='y', which='major', color='0.5', linestyle='-')
        # x-axis
        format_date_axis(ax, [t[0], t[-1]])
        plt.xlabel('')

    for i in range(0,len(fnames)):
        ncL[i].close()
    figname = fname_str +'_cont_' + groupname + '_' + str(numyears) + 'y.png'
    plt.savefig(figname)
    print('python line plot saved in: ' + figname)

def plot_singlefile(fname,numyears,groupname,varset,variant,modname):
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
    elif len(varset)%5 == 0:
        numcol = 2.

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
            if variant=='IA':
                namelib = namelibNbasedIA
            elif variant=='DA':
                namelib = namelibNbasedDA
        elif modname == 'Cbased':
            if variant == 'IA':
                namelib = namelibCbasedIA
            elif variant == 'DA':
                namelib = namelibCbasedDA
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
        if varn in prettynames:
            prettyname = prettynames[varn]
        else:
            prettyname = varn
        plt.title('%s [$%s$]' % (prettyname, prettyunit), size=12.0)

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

def collect_data(fnames,variants,numyears):
    ncL=[]; ncvL=[]; tL=[]; tiL=[]
    for fno in range(0,len(fnames)):
        if fno==0:
            fname_str=fnames[fno].split('/')[-1].split('.nc')[0]
        else:
            fname_str=fname_str+ '_vs_' + fnames[fno].split('/')[-1].split('.nc')[0]
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
    return(ncL, ncvL, tL, tiL, fname_str)

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
            if '*' in varn or '/' in varn:
                varfound, varvalsOP, longnameOP, unitsOP = get_varvals_op(ncv, varn)
                varvals = varvals + varvalsOP
                longname = '%s %s %s' % (longname, symb, varn)
            else:
                varvals = varvals + np.squeeze(ncv[varn][:])
                longname = '%s %s %s' % (longname, symb, varn)
            # if not ncv[varn].units == units:
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
      #fnames = ['/home/onur/setups/test-BGCmodels/nflexpd/ideal_highlat_NflexPD-Nbased_Cbased/0D-Highlat_wconst_dm_NbasedDA.nc',
               # '/home/onur/setups/test-BGCmodels/nflexpd/ideal_highlat_NflexPD-Nbased_Cbased/0D-Highlat_wconst_dm_CbasedDA.nc']
      fnames = ['/home/onur/setups/test-BGCmodels/nflexpd/ideal_highlat_NflexPD-Nbased_Cbased/0D-Highlat_wconst_lext_Ldvar_Tvar_2P_CbasedIA_modular_24h.nc',
                '/home/onur/setups/test-BGCmodels/nflexpd/ideal_highlat_NflexPD-Nbased_Cbased/0D-Highlat_wconst_lext_Ldvar_Tvar_2P_CbasedIA_modular_24h.nc']
                #'/home/onur/setups/test-BGCmodels/nflexpd/ideal_highlat_NflexPD-Nbased_Cbased/0D-Highlat_wconst_lint_CbasedIA_dm.nc']
      print('plotting default file(s):'+'; '.join(fnames))
    else:
      print('plotting file specified:'+sys.argv[1])
      fnames=sys.argv[1].split(',')

    if len(sys.argv)<3:
      #variants = ['dm', '6h']
      #variants = ['IA','DA']
      variants = ['IA', 'IA']
    else:
      variants=sys.argv[2].split(',')

    if len(sys.argv)<4: #no third argument was passed
      #modnames=['Nbased','Cbased']
      modnames=['Cbased','Cbased']#
    else:
      modnames=sys.argv[3].split(',')

    if len(sys.argv)<5: #no third argument was passed
      #ids=['sim']
      #ids=['Nbased','Cbased']
      #ids = ['IA', 'DA']
      ids=['PAR:N','PAR:A']
    else:
      ids=sys.argv[4].split(',')

    if len(sys.argv)<6: #no third argument was passed
      #numyears=-1 # -1 means plot everything
      numyears=1
    else:
      numyears=int(sys.argv[5]) #number of years to plot (counting from the last year backwards)

    print('plotting last ' + str(numyears) + ' year of the simulations: ' + ','.join(ids))

    main(fnames, numyears, modnames, variants, ids)
