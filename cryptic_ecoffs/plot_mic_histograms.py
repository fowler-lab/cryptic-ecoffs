import numpy
import matplotlib.pyplot as plt


def plot_mic_histogram( axis,\
                        dilutions_series,\
                        drug_name=None,\
                        mutation_name=None,\
                        conc_table=None,\
                        colour="grey",\
                        label_yaxis=True,\
                        label_xaxis=False,\
                        label_drug=True,\
                        label_drug_y=3.2,
                        label_values=True,\
                        max_label_value=100,\
                        label_N=True,\
                        label_mutation=False,\
                        linewidth=0.5,\
                        x_max=None,\
                        y_max=9.0,\
                        normed=False,\
                        label_fontsize=9,\
                        axis_fontsize=7,\
                        ecoff=False):

    # set the limits of the graph
    # these are passed as arguments so they can be dynamically determined using the original data
    if x_max is not None and x_max>0:
        axis.set_xlim([0,x_max])
    axis.set_ylim([0.4,y_max])

    if ecoff:
        ecoff=conc_table.loc[conc_table.BINARY_PHENOTYPE=="S"].DILUTION.max()
        axis.plot([-0.5,x_max],[ecoff+0.5,ecoff+0.5],color="black",linewidth=1)
        atu=conc_table.loc[conc_table.BINARY_PHENOTYPE=="I"].DILUTION.max()
        if ~numpy.isnan(atu):
            axis.plot([-0.5,x_max],[atu+0.5,atu+0.5],color="grey",linewidth=1,linestyle="-")

    # set the ticks to zero length
    axis.tick_params(axis=u'both', which=u'both',length=0)

    if not label_xaxis:
        axis.axes.get_xaxis().set_visible(False)
        axis.spines['bottom'].set_visible(False)

    # hide most of the frame
    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)

    # if specified, label the drug
    if label_mutation:
        cols=mutation_name.split("@")
        label_text=cols[0]+"\n"+cols[1]
    elif label_drug:
        label_text=drug_name

    if label_mutation or label_drug:
        axis.text((0.45*x_max),label_drug_y,label_text,horizontalalignment='center',size=label_fontsize,color='black')

        # if also specified, label the number of samples
        if label_N:
            # text="n=%i" % dilutions_series.count()
            text="n={:,}".format(dilutions_series.count())
            axis.text((0.45*x_max),label_drug_y-0.4,text,horizontalalignment='center',weight='bold',size=label_fontsize,color='black')

    # if specified, hide the y-tics and labels
    if not label_yaxis:
        axis.axes.get_yaxis().set_visible(False)
        axis.spines['left'].set_visible(False)

    # find out the highest dilution
    max_dilution=conc_table.DILUTION.max()

    # thereby create a list of dilutions
    dilutions=[i for i in range(1,int(max_dilution)+2)]

    # and a list of concentrations, including ≤ and >
    labels=list(conc_table.CONC)

    axis.set_yticks(dilutions[:-1])
    axis.set_yticklabels(labels,fontdict={'size':axis_fontsize})
    axis.spines['left'].set_bounds(0.5, max_dilution+0.5)

    if len(dilutions_series)>1:
        axis.hist(  dilutions_series,
                    bins=dilutions,
                    density=normed,
                    histtype='barstacked',
                    align='left',
                    orientation='horizontal',
                    alpha=0.8,
                    color=colour,
                    edgecolor='black',
                    linewidth=linewidth)
    else:
        axis.hist(  dilutions_series,
                    bins=dilutions,
                    density=normed,
                    histtype='stepfilled',
                    align='left',
                    orientation='horizontal',
                    alpha=0.8,
                    color=colour,
                    edgecolor='black',
                    linewidth=linewidth)

    # if specified, label the number of samples in small bins, as defined by max_label_value
    if label_values:
        # if normed:
        #     x_max = axis.get_xlim()[1]
        #
        # print(x_max,normed)
        # have to deal with being given a list of series, or just a single series
        hist=[]
        if isinstance(dilutions_series, list):
            for i in dilutions_series:
                hist+=list(i.values)
        else:
            hist=list(dilutions_series.values)

        if normed:
            (values,dil)=numpy.histogram(hist,dilutions,density=True)
        else:
            (values,dil)=numpy.histogram(hist,dilutions)
        for v,d in zip(values,dil):
            if v<max_label_value and v>0:
                x_value=int(x_max*(v/x_max))+int(0.03*x_max)
                axis.text(x_value,d,v,horizontalalignment='left',verticalalignment='center',color='black',size=label_fontsize)

    return(axis)
