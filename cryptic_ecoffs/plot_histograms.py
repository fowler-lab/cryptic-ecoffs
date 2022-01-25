#! /usr/bin/env python

import pandas
import copy
import numpy
import math
import pathlib
from scipy.stats import norm

import cryptic_ecoffs
import matplotlib.pyplot as plt

def plot_histograms(    plate_design,\
                        plate_design_lookup,\
                        df,\
                        output_stem1,\
                        output_stem2,\
                        exclude_site=None,\
                        print_histogram=False,\
                        save_histogram=False,\
                        print_direct=False,\
                        save_ecoff=False,
                        plot_fitted_normal=None,\
                        plot_components=1,\
                        plot_type=None,\
                        ecoff_percentile=0.99,\
                        bar_colour='lightgrey',\
                        show_graph=False        ):

    PLATE_LAYOUT=pandas.read_pickle("tables/PLATE_LAYOUT.pkl.gz")
    PLATE_LAYOUT.reset_index(inplace=True)
    PLATE_LAYOUT.replace(to_replace="<=",value="≤",inplace=True,regex=True)

    direct_marker_lookup={0.95:'#a1d99b',0.975:'#41ab5d',0.99:'#005a32',0.999:'#f768a1'}

    color_lookup={'direct99':'#41ab5d','visualraw':'#377eb8','ecoffinderraw':'#253494','ecoffinderSS*S':'#984ea3','intregraw':'#1f78b4','tentativeraw':'#e41a1c','ecoffinderSSSSnotR':'#984ea3','intregSS*S':'#e41a1c','intregSSSSnotR':'#e41a1c','intregSSSSnotRclean':'#e41a1c','ecoffinderSSSS':'#984ea3','intregSSSS':'#e41a1c'}

    marker_lookup={'visualraw':'<','ecoffinderraw':'<','ecoffinderSS*S':'<','intregSSSSnotR':'<','intregSS*S':'<','tentativeraw':'<','direct99':'<','intregSSSSnotRclean':'<'}

    label_drug_y_lookup={'LZD':1.8,'EMB':1.8,'MXF':7.8,'LEV':7.8,'BDQ':6.8,'ETH':6.8,'INH':6.8,"KAN":5.8,"AMI":6.8,"CFZ":6.8,"RIF":3.8,"RFB":5.8,"DLM":4.8}

    # define the exact layout of the final graph
    # drug_list=[['LZD','EMB','MXF','LEV','BDQ','ETH','INH'],['KAN','AMI','CFZ','RIF','RFB','DLM','PAS']]
    drug_list=[["INH",'RIF','EMB','MXF','LEV','KAN','AMI'],['ETH','RFB','CFZ','LZD','DLM','BDQ','PAS']]

    if len(df)>1:
        stacked_bars=True
    else:
        stacked_bars=False
        df=df[0]


    layout={}
    for drug_name in plate_design_lookup.keys():
        layout[drug_name]=copy.deepcopy(PLATE_LAYOUT.loc[(PLATE_LAYOUT.PLATEDESIGN.isin(plate_design)) & (PLATE_LAYOUT.DRUG==drug_name)] )
        layout[drug_name].sort_values('DILUTION',inplace=True)

    # find out the largest number of dilutions for this plate design
    y_max=PLATE_LAYOUT.loc[PLATE_LAYOUT.PLATEDESIGN.isin(plate_design)].DILUTION.max()+0.5

    data={}
    if stacked_bars:
        data1={}

    for drug_name in plate_design_lookup.keys():
        if stacked_bars:
            data[drug_name]=copy.deepcopy(df[0].loc[(df[0].PLATEDESIGN.isin(plate_design_lookup[drug_name]))  & (df[0].DRUG==drug_name)])
            data1[drug_name]=copy.deepcopy(df[1].loc[(df[1].PLATEDESIGN.isin(plate_design_lookup[drug_name]))  & (df[1].DRUG==drug_name)])
        else:
            data[drug_name]=copy.deepcopy(df.loc[(df.PLATEDESIGN.isin(plate_design_lookup[drug_name]))  & (df.DRUG==drug_name)])

    # find out, dynamically, the largest number of samples in any one bin
    x_max=0
    bins=[i for i in range(1,13)]
    direct_ecoff={}
    histogram_lines=''

    for row in range(2):
        for col in range(7):
            drug_name=drug_list[row][col]
            direct_ecoff[drug_name]={}

            dilutions=copy.deepcopy(data[drug_name].DILUTION)

            if not dilutions.empty:

                (values,dil)=numpy.histogram(dilutions,bins)

                hist,edges=numpy.histogram(dilutions, density=True, bins=bins )
                hist=numpy.cumsum(hist)/numpy.cumsum(hist).max()


                if print_histogram or save_histogram:
                    c=list(layout[drug_name].CONC)
                    for i in dil:
                        if i <= len(c) and i<=len(values):
                            if print_histogram:
                                print("%s,%s,%s,%i,%s,%i" % (drug_name, output_stem2, output_stem1, i, c[i-1],values[i-1]))
                            if save_histogram:
                                histogram_lines+="%s,%s,%s,%i,%s,%i\n" % (drug_name, output_stem2, output_stem1, i, c[i-1],values[i-1])

                for k in [0.95,0.975,0.99,0.999]:
                    direct_ecoff[drug_name][k]=None
                    for i,j in zip(edges,hist):
                        if j>=k and direct_ecoff[drug_name][k] is None:
                            direct_ecoff[drug_name][k]=i
                            break


                    if print_direct:
                        print("%s %20s %6.1f %10s %s %7i" % (drug_name,'direct',i*100,output_stem2,output_stem1, direct_ecoff[drug_name][i]))

                if max(values) > x_max:
                    x_max=max(values)
    x_max=int(1.02*x_max)


    if output_stem2=='h37rv':
        label_values=True
        normed=False
        label_drug_y_lookup={'LZD':6.8,'EMB':6.8,'MXF':6.8,'LEV':6.8,'BDQ':6.8,'ETH':6.8,'INH':6.8,"KAN":5.8,"AMI":6.8,"CFZ":6.8,"RIF":6.8,"RFB":5.8,"DLM":6.8}
        x_max=int(1.7*x_max)
    else:
        label_values=False
        normed=True
        x_max=1

    # create a new figure object with an embedded axes objects with 14 elements
    fig, axes = plt.subplots(2,7, sharey=False,figsize=(14,(12/9.5)*y_max))

    # subtly adjust the spacing to use as much white space as possible
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    if exclude_site is not None:
        graph_filename="graphs/exclude"+exclude_site
    else:
        graph_filename="graphs/"

    pathlib.Path(graph_filename).mkdir(parents=True, exist_ok=True)

    graph_filename+="/hist-all_drugs-mic-"+output_stem1+"-"+output_stem2

    if plot_fitted_normal is not None:
        graph_filename+="-"+plot_fitted_normal
        graph_filename+="-"+str(plot_components)

    if exclude_site is not None:
        graph_filename+="-exclude"+exclude_site

    if plot_fitted_normal=="intreg":
        ecoff_dataframe=pandas.read_csv("tables/STATA_INTERVAL_REGRESSION_RESULTS.csv")
        graph_filename+='-ecoff-'+str(ecoff_percentile)
    elif plot_fitted_normal=="ecoffinder":
        foo=pandas.read_csv("tables/ECOFFINDER-RESULTS.csv")
        ecoff_dataframe=foo.loc[foo.dataset==output_stem1]
        graph_filename+='-ecoff-'+str(ecoff_percentile)
    elif plot_fitted_normal in ["summary"]:
        ecoff_dataframe=pandas.read_csv("tables/ECOFFS_SUMMARY.csv")

    if save_ecoff:
        OUTPUT_ECOFF=open(graph_filename+".txt",'w')

    if save_histogram:
        OUTPUT = open(graph_filename+'.csv','w')
        OUTPUT.write(histogram_lines)
        OUTPUT.close()

    # now iterate over the graph structure
    for row in range(2):
        for col in range(7):

            drug_name=drug_list[row][col]

            # add an asterisk to the drug if it has the same range on both UKMYC5 and UKMYC6
            drug_label_name=drug_name
            if drug_name in ['KAN','MXF','LEV','ETH','RFB']:
                drug_label_name+="*"

            # check that there are data to plot (e.g. there won't be for PAS on UKMYC6!)
            if not data[drug_name].empty and drug_name!="PAS":
                if not stacked_bars:
                    axes[row][col]=cryptic_ecoffs.plot_mic_histogram(axes[row][col],
                                                                  data[drug_name].DILUTION,
                                                                  drug_name=drug_label_name,
                                                                  conc_table=layout[drug_name],
                                                                  colour=bar_colour,
                                                                  label_values=label_values, #True
                                                                  label_drug=True,
                                                                  label_drug_y=label_drug_y_lookup[drug_name],
                                                                  label_N=True,
                                                                  label_yaxis=True,
                                                                  label_xaxis=False,
                                                                  max_label_value=2000,
                                                                  x_max=x_max,
                                                                  y_max=y_max,
                                                                  normed=normed,
                                                                  axis_fontsize=12,
                                                                  label_fontsize=14) # remove
                else:
                    axes[row][col]=cryptic_ecoffs.plot_mic_histogram(axes[row][col],
                                                                  [data[drug_name].DILUTION,data1[drug_name].DILUTION],
                                                                  drug_name=drug_label_name,
                                                                  conc_table=layout[drug_name],
                                                                  colour=['lightgrey','white'],
                                                                  label_values=label_values, #True
                                                                  label_drug=True,
                                                                  label_drug_y=label_drug_y_lookup[drug_name],
                                                                  label_N=False,
                                                                  label_yaxis=True,
                                                                  label_xaxis=False,
                                                                  max_label_value=2000,
                                                                  x_max=1, #chgange back to x_max
                                                                  y_max=y_max,
                                                                  normed=True,
                                                                  axis_fontsize=12,
                                                                  label_fontsize=14) # remove

                if plot_fitted_normal in ['visual','tentative']:

                    lowest_conc=float(PLATE_LAYOUT.loc[(PLATE_LAYOUT.DRUG==drug_name) & (PLATE_LAYOUT.PLATEDESIGN.isin(plate_design)) & (PLATE_LAYOUT.DILUTION==1)].CONC.str[1:])
                    lowest_log2conc=math.log(lowest_conc,2)
                    shift_to_dilution=1-lowest_log2conc
                    if exclude_site is None:
                        df=ecoff_dataframe.loc[(ecoff_dataframe.drug==drug_name) & (ecoff_dataframe.type==plot_type)]
                    else:
                        df=ecoff_dataframe.loc[(ecoff_dataframe.drug==drug_name) & (ecoff_dataframe.type==plot_type+"-exclude"+exclude_site)]
                    if not numpy.isnan(float(df.comp1_mic_99)):
                        ecoff_conc=float(df.comp1_mic_99)
                        log2ecoff=math.log(ecoff_conc,2)
                        ecoff_dilution=log2ecoff+shift_to_dilution
                        axes[row][col].plot(0.1,ecoff_dilution,marker='<',markersize=14.0,markerfacecolor=color_lookup[plot_fitted_normal+plot_type],markeredgewidth=0)
                        axes[row][col].text(0.2,ecoff_dilution,str("%.2f" % ecoff_conc),verticalalignment='center',color=color_lookup[plot_fitted_normal+plot_type],size=14,weight='bold')

                elif plot_fitted_normal in ['summary']:

                    lowest_conc=float(PLATE_LAYOUT.loc[(PLATE_LAYOUT.DRUG==drug_name) & (PLATE_LAYOUT.PLATEDESIGN.isin(plate_design)) & (PLATE_LAYOUT.DILUTION==1)].CONC.str[1:])
                    lowest_log2conc=math.log(lowest_conc,2)
                    shift_to_dilution=1-lowest_log2conc
                    position=0.2
                    df=ecoff_dataframe.loc[(ecoff_dataframe.drug==drug_name) & (ecoff_dataframe.type==plot_type) & (ecoff_dataframe.plate_design==plate_design[0])]
                    for type in ['ecoffinderraw','direct99','intregSSSSnotRclean']:
                        if type in df.columns and ~numpy.isnan(float(df[type])):
                            ecoff_conc=float(df[type])
                            log2ecoff=math.log(ecoff_conc,2)
                            ecoff_dilution=log2ecoff+shift_to_dilution+0.5
                            axes[row][col].plot(position,ecoff_dilution,marker=marker_lookup[type],markersize=18.0,markerfacecolor=color_lookup[type],markeredgewidth=0)
                            position+=0.12

                elif plot_fitted_normal == 'direct':

                    prev_ecoff=0
                    x=0.1
                    for i in [0.95,0.975,0.99]:
                        ecoff=direct_ecoff[drug_name][i]
                        if ecoff==prev_ecoff:
                            x+=0.1
                        if ecoff is not None:
                            axes[row][col].plot(x,ecoff+0.5,marker='<',markersize=14.0,markerfacecolor=direct_marker_lookup[i],markeredgewidth=0)
                        prev_ecoff=ecoff


                elif plot_fitted_normal is not None:

                    lowest_conc=float(PLATE_LAYOUT.loc[(PLATE_LAYOUT.DRUG==drug_name) & (PLATE_LAYOUT.PLATEDESIGN.isin(plate_design)) & (PLATE_LAYOUT.DILUTION==1)].CONC.str[1:])
                    lowest_log2conc=math.log(lowest_conc,2)
                    shift_to_dilution=1-lowest_log2conc
                    bins=numpy.arange(1,13,0.1)
                    bins2=numpy.arange(1,13,1)
                    assert plot_type in ['raw','SSUS','SSSSnotR','SSSS','SSSSnotRclean'], "incorrect plot type specified! "+plot_type
                    if exclude_site is None:
                        df=ecoff_dataframe.loc[(ecoff_dataframe.drug==drug_name) & (ecoff_dataframe.type==plot_type)]
                    else:
                        df=ecoff_dataframe.loc[(ecoff_dataframe.drug==drug_name) & (ecoff_dataframe.type==plot_type+"-exclude"+exclude_site)]

                    if plot_components==1:
                        if 'method' in ecoff_dataframe.columns:
                            df=df.loc[df.method.isin(["intreg_all","intreg_1pop"])]

                        if not df.empty:
                            mu=float(df.comp1_b1)
                            sigma=float(df.comp1_v1)

                            if sigma>0:
                                estimated_normal=numpy.random.normal(mu,sigma,1000000)+shift_to_dilution
                                hist,edges=numpy.histogram(estimated_normal,density=True,bins=bins)
                                axes[row][col].plot(hist,edges[:-1],color=color_lookup[plot_fitted_normal+plot_type],linewidth=2.0)

                                # work out the middle percentage from a percentile e.g. a percentile of 99% means the middle 98% of the distribution
                                alpha=(2*ecoff_percentile)-1

                                # this will return the endpoints of the middle percentage, so take the higher (second) value
                                ecoff=norm.interval(alpha,loc=mu+shift_to_dilution,scale=sigma)[1]

                                # convert to MIC
                                ecoff_conc=2**(ecoff-shift_to_dilution)

                                if save_ecoff:
                                    OUTPUT_ECOFF.write("%s %20s %10s %7.2f %7.3f\n" % (drug_name,plot_fitted_normal,plot_type,ecoff,ecoff_conc))

                                axes[row][col].plot(0.1,ecoff+0.5,marker='<',markersize=14.0,markerfacecolor=color_lookup[plot_fitted_normal+plot_type],markeredgewidth=0)
                                axes[row][col].text(0.2,ecoff+0.5,str("%.2f" % ecoff_conc),verticalalignment='center',color=color_lookup[plot_fitted_normal+plot_type],size=14,weight='bold')

                    elif plot_components==2:

                        mu1=None
                        mu2=None

                        if 'method' in ecoff_dataframe.columns:
                            df=df.loc[df.method.isin(["intreg_all"])]
                        if not df.empty:

                            # fig2=plt.figure(2,8)
                            # axes2=fig2.gca()
                            perc1=float(df.comp2_perc1)/100.
                            perc2=float(df.comp2_perc2)/100.

                            if ~numpy.isnan(perc1) and ~numpy.isnan(perc2):

                                # # swap the coefficients if, like it did for LZD, get them the wrong way around
                                # if float(df.comp2_b2)<float(df.comp2_b1):
                                #
                                #     mu2=float(df.comp2_b1)
                                #     mu1=float(df.comp2_b2)
                                #
                                #     sigma2=float(df.comp2_v1)
                                #     perc2=float(df.comp2_perc1)/100.
                                #     sigma1=float(df.comp2_v2)
                                #     perc1=float(df.comp2_perc2)/100.
                                #
                                # else:

                                mu1=float(df.comp2_b1)
                                mu2=float(df.comp2_b2)

                                sigma1=float(df.comp2_v1)
                                sigma2=float(df.comp2_v2)

                                if perc2<0.02 or sigma2>3:
                                    # print(drug_name,perc2, sigma2)
                                    mu2=None
                                    sigma2=None
                                    perc2=None

                            elif ~numpy.isnan(perc1):

                                    mu1=float(df.comp2_b1)
                                    mu2=None
                                    sigma1=float(df.comp2_v1)


                            if mu1 is not None:
                                estimated_normal1=numpy.random.normal(mu1,sigma1,1000000)+shift_to_dilution
                                hist,edges=numpy.histogram(estimated_normal1,density=True,bins=bins)
                                axes[row][col].plot(perc1*hist,edges[:-1],color='#bd0026')
                                hist1,edges1=numpy.histogram(estimated_normal1,density=True,bins=bins2)
                                # print(perc1*hist1,edges1)

                                # work out the middle percentage from a percentile e.g. a percentile of 99% means the middle 98% of the distribution
                                alpha=(2*ecoff_percentile)-1

                                # this will return the endpoints of the middle percentage, so take the higher (second) value
                                ecoff=norm.interval(alpha,loc=mu1+shift_to_dilution,scale=sigma1)[1]

                                # convert to MIC
                                ecoff_conc=2**(ecoff-shift_to_dilution)

                                # print("%s %10s %7.2f %7.3f" % (drug_name,plot_type,ecoff,ecoff_conc))
                                if ecoff_conc!=0.0 and ecoff!=1.1:
                                    axes[row][col].plot(0.1,ecoff+0.5,marker='<',markersize=14.0,markerfacecolor="#bd0026",markeredgewidth=0)
                                    axes[row][col].text(0.2,ecoff+0.5,str("%.2f" % ecoff_conc),verticalalignment='center',color="#bd0026",size=14,weight='bold')
                                ecoff1=ecoff


                            if mu2 is not None:
                                estimated_normal2=numpy.random.normal(mu2,sigma2,1000000)+shift_to_dilution
                                hist,edges=numpy.histogram(estimated_normal2,density=True,bins=bins)
                                axes[row][col].plot(perc2*hist,edges[:-1],color="#fd8d3c",linewidth=2.0)
                                hist2,edges2=numpy.histogram(estimated_normal2,density=True,bins=bins2)
                                # print("2",perc2*hist2,edges2)
                                R=perc2*hist2
                                S=perc1*hist1
                                d=1
                                # for p_s,p_r in zip(S,R):
                                #     if p_s==0:
                                #         print("%s %3s %2i %7.1f" % (plate_design[0],drug_name,d,100))
                                #     elif p_r==0:
                                #         print("%s %3s %2i %7.1f" % (plate_design[0],drug_name,d,0))
                                #     else:
                                #         print("%s %3s %2i %7.1f" % (plate_design[0],drug_name,d,100*p_r/(p_s+p_r)))
                                #     d+=1

                                # work out the middle percentage from a percentile e.g. a percentile of 99% means the middle 98% of the distribution
                                alpha=(2*ecoff_percentile)-1

                                # this will return the endpoints of the middle percentage, so take the lower (first) value
                                ecoff=norm.interval(alpha,loc=mu2+shift_to_dilution,scale=sigma2)[0]

                                # convert to MIC
                                ecoff_conc=2**(ecoff-shift_to_dilution)

                                # print("%s %10s %7.2f %7.3e" % (drug_name,plot_type,ecoff,ecoff_conc))
                                if ecoff>-1.0 and ecoff<10.0 and ecoff_conc>lowest_conc:
                                    axes[row][col].plot(0.1,ecoff+0.5,marker='<',markersize=14.0,markerfacecolor="#fd8d3c",markeredgewidth=0)
                                    if abs(ecoff1-ecoff)<0.5:
                                        axes[row][col].text(0.2,ecoff+0.5,str("%.2f" % ecoff_conc),verticalalignment='center',color="#fd8d3c",size=14,weight='bold')
                                    else:
                                        axes[row][col].text(0.2,ecoff+0.5,str("%.2f" % ecoff_conc),verticalalignment='center',color="#fd8d3c",size=14,weight='bold')

                            # extent=axes[row][col].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                            # fig.savefig(graph_filename+"-"+drug_name+".pdf",bbox_inches=extent.expanded(1.02, 1.02))

                    else:
                        raise ValueError("components can only be 1 or 2, instead was "+str(plot_components))

            else:
                axes[row][col].axis("off")

    fig.savefig(graph_filename+".pdf",bbox_inches='tight',transparent=True)

    if not show_graph:
        plt.close()

    if save_ecoff:
        OUTPUT_ECOFF.close()
