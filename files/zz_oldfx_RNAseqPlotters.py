
"""juanPlot

    Collection of plotting and formatting functions for Angueyra et al., 2021
    Created: September 2021
    Updated: October 2021
        To build wheel run ```python setup.py bdist_wheel --universal;```
        Then ```cp ./dist/juanPlot-0a2-py2.py3-none-any.whl ~/Documents/Repositories/jupyterLiteDemo/content```
"""
# import required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from scipy.stats import zscore

def defaultFonts(ax=None):
    # standarizes font sizes across plots
    fontTicks = font_manager.FontProperties(size=24)
    fontLabels = font_manager.FontProperties(size=28)
    fontTitle = font_manager.FontProperties(size=28)
    ax.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 24)
    ax.yaxis.offsetText.set_fontsize(24)
    return fontTicks, fontLabels, fontTitle

"""bar plots"""

def plotBars(barData, geneSymbol, ax=None, pC=None):
    """Creates a bar plot for a single gene
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,7) # Rods
    n = np.append(n, 7 + np.arange(1,6)) # UV
    n = np.append(n, 13 + np.arange(1,7)) # S
    n = np.append(n, 20 + np.arange(1,8)) # M
    n = np.append(n, 28 + np.arange(1,7)) # L
    h_start = 7
    h_end = 37
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    barColors = [
        pC['r'],pC['r'],pC['r'],pC['r'],pC['r'],pC['r'],
        pC['u'],pC['u'],pC['u'],pC['u'],pC['u'],
        pC['s'],pC['s'],pC['s'],pC['s'],pC['s'],pC['s'],
        pC['m'],pC['m'],pC['m'],pC['m'],pC['m'],pC['m'],pC['m'],
        pC['l'],pC['l'],pC['l'],pC['l'],pC['l'],pC['l']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot(geneSymbol, ax=ax)
    return pH

def formatBarPlot(geneSymbol, ax=None):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks([3.5,10,16.5,24,31.5])
    ax.set_xticklabels(['Rods','UV','S','M','L']);
    ax.set_ylabel('FPKM', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

def plotBars_Ogawa2021(barData, geneSymbol, ax=None, pC=None, pctPlot=False):
    """Creates a bar plot for a single gene for data from Ogawa et al. (2021) (https://doi.org/10.1038/s41598-021-96837-z)
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,9)
    delta = 0
    if pctPlot:
        delta = 9
    h_start = 2 + delta
    h_end = 10 + delta
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    barColors = [
        pC['r'],pC['u'],pC['s'],
        pC['m'],pC['l'],pC['m4'],
        pC['onBC'],pC['offBC']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Ogawa2021(geneSymbol, ax=ax, pctPlot=pctPlot)
    return pH

def formatBarPlot_Ogawa2021(geneSymbol, ax=None, pctPlot=False):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks(np.arange(1,9))
    ax.set_xticklabels(['Rods','UV','S','M','L', 'M4','BC$_{on}$','BC$_{off}$']);
    ax.xaxis.set_tick_params(rotation=45)
    ax.set_ylabel('avg. counts', fontproperties=fontLabels)
    if pctPlot:
        ax.set_ylabel('% expressing', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

def plotBars_Hoang2020(barData, geneSymbol, ax=None, pC=None, pctPlot=False):
    """Creates a bar plot for a single gene for data from Hoang et al. (2020) (https://doi.org/10.1126/science.abb8598)
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,8)
    delta = 0
    if pctPlot:
        delta = 8
    h_start = 2 + delta
    h_end = 9 + delta
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    barColors = [
        pC['r'],pC['u'],pC['s'],
        pC['m'],pC['m'],pC['m4'],
        pC['l']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Hoang2020(geneSymbol, ax=ax, pctPlot=pctPlot)
    return pH

def formatBarPlot_Hoang2020(geneSymbol, ax=None, pctPlot=False):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks(np.arange(1,8))
    ax.set_xticklabels(['Rods','UV','S','M1','M3', 'M4','L']);
    ax.xaxis.set_tick_params(rotation=45)
    ax.set_ylabel('avg. counts', fontproperties=fontLabels)
    if pctPlot:
        ax.set_ylabel('% expressing', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

def plotBars_Sun2018(barData, geneSymbol, ax=None, pC=None):
    """Creates a bar plot for a single gene for data from Sun, Galicia and Stenkamp (2018) (https://doi.org/10.1186/s12864-018-4499-y)
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,5) # Rods (GFP+)
    n = np.append(n, 5 + np.arange(1,5)) # Not Rods (GFP-)
    h_start = 7
    h_end = 16
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    barColors = [
        pC['r'],pC['r'],pC['r'],pC['r'],
        pC['m4'],pC['m4'],pC['m4'],pC['m4']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Sun2018(geneSymbol, ax=ax)
    return pH

def formatBarPlot_Sun2018(geneSymbol, ax=None):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks([2.5,7.5])
    ax.set_xticklabels(['Rods','notRods']);
    ax.set_ylabel('cpm', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

def plotBars_Hoang2020_Ret(barData, geneSymbol, ax=None, pC=None, pctPlot=False):
    """Creates a bar plot for a single gene for data from Hoang et al. (2020) (https://doi.org/10.1126/science.abb8598)
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,18)
    delta = 0
    if pctPlot:
        delta = 19
    h_start = 2 + delta
    h_end = 19 + delta
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {
            'RPC' : '#DADADA', # Retinal progenitor cell
            'PRPC' : '#dfdac8', # Photoreceptor progenitor cell
            'Cones_larval' : '#dcc360', #
            'Cones_adult' : '#ffd429', #
            'Rods' : '#7d7d7d', #
            'HC' : '#FC7715', # Horizontal cells
            'BC_larval' : '#ccf2ff', # Bipolar cell (developing)
            'BC_adult' : '#663d00', # Bipolar cell (mature)
            'AC_larval' : '#3DF591', # Amacrine cell (developing)
            'ACgaba' : '#3DF5C3', #
            'ACgly' : '#56F53D', #
            'RGC_larval' : '#F53D59', # Retinal Ganglion cell (developing)
            'RGC_adult' : '#BB0622', # Retinal Ganglion cell (mature)
            'MGi' : '#EA9D81', # Muller glia (immature)
            'MG1' : '#A2644E', # Muller glia (mature)
            'MG2' : '#7E4835', # Muller glia (mature)
            'MG3' : '#613728', # Muller glia (mature)
        }
    barColors = [
        pC['RPC'],pC['PRPC'],
        pC['Cones_larval'],pC['Cones_adult'],pC['Rods'],
        pC['HC'],
        pC['BC_larval'],pC['BC_adult'],
        pC['AC_larval'],pC['ACgaba'],pC['ACgly'],
        pC['RGC_larval'],pC['RGC_adult'],
        pC['MGi'],pC['MG1'],pC['MG2'],pC['MG3']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Hoang2020_Ret(geneSymbol, ax=ax, pctPlot=pctPlot)
    return pH

def formatBarPlot_Hoang2020_Ret(geneSymbol, ax=None, pctPlot=False):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks(np.arange(1,18))
    ax.set_xticklabels(['RPC','PRPC','C$_{larval}$','C$_{adult}$','R$_{ods}$','HC','BC$_{larval}$','BC$_{adult}$','AC$_{larval}$','AC$_{GABA}$','AC$_{Gly}$','RGC$_{larval}$','RGC$_{adult}$','MGi','MG1','MG2','MG3']);
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", va="center",rotation_mode="anchor")
    ax.set_ylabel('avg. counts', fontproperties=fontLabels)
    if pctPlot:
        ax.set_ylabel('% expressing', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

def plotBars_Hoang2020_PRDev(barData, geneSymbol, ax=None, pC=None, pctPlot=False):
    """Creates a bar plot for a single gene for data from Hoang et al. (2020) (https://doi.org/10.1126/science.abb8598)
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,12)
    delta = 12
    if pctPlot:
        delta = 12
    h_start = 2 + delta
    h_end = 13 + delta
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {
            'PRP' : "#dfdac8",
            'eslPR' : '#dacd9a',
            'mslPR' : '#dcc360',
            'lslPR' : '#cca819',
            'adPR' : '#ffd429',
            'lslR' : '#a3a3a3',
            'r' : '#7d7d7d',
            'u' : '#B540B7',
            's' : '#4669F2',
            'm' : '#04CD22',
            'l' : '#CC2C2A',
        }
    barColors = [
        pC['PRP'],
        pC['eslPR'],pC['mslPR'],pC['lslPR'],
        pC['adPR'],pC['lslR'],
        pC['r'],pC['u'],pC['s'],pC['m'],pC['l']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Hoang2020_PRDev(geneSymbol, ax=ax, pctPlot=pctPlot)
    return pH

def formatBarPlot_Hoang2020_PRDev(geneSymbol, ax=None, pctPlot=False):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks(np.arange(1,12))
    ax.set_xticklabels(['PRPC','PR$_{larval-early}$','PR$_{larval-mid}$','PR$_{larval-late}$','PR$_{adult}$','Rod$_{larval-late}$','Rod$_{adult}$','UV$_{adult}$','S$_{adult}$','M$_{adult}$','L$_{adult}$']);
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", va="center",rotation_mode="anchor")
    ax.set_ylabel('avg. counts', fontproperties=fontLabels)
    if pctPlot:
        ax.set_ylabel('% expressing', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

def plotBars_Nerli2022(barData, geneSymbol, ax=None, pC=None):
    """Creates a bar plot for a single gene
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,6) # Retinal progenitors
    n = np.append(n, 6 + np.arange(1,6)) # Photoreceptors
    n = np.append(n, 11 + np.arange(1,6)) # Amacrine/Horizontal cells
    n = np.append(n, 16 + np.arange(1,6)) # Retinal ganglion cells
    h_start = 1
    h_end = 21
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {'RPC' : '#DADADA', 'PR' : '#dcc360', 'HC_AC' : '#3DF591', 'RGC' : '#F53D59'}
    barColors = [
        pC['RPC'],pC['RPC'],pC['RPC'],pC['RPC'],pC['RPC'],
        pC['PR'],pC['PR'],pC['PR'],pC['PR'],pC['PR'],
        pC['HC_AC'],pC['HC_AC'],pC['HC_AC'],pC['HC_AC'],pC['HC_AC'],
        pC['RGC'],pC['RGC'],pC['RGC'],pC['RGC'],pC['RGC'],
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Nerli2022(geneSymbol, ax=ax)
    return pH

def plotBars_Nerli2022(barData, geneSymbol, ax=None, pC=None):
    """Creates a bar plot for a single gene
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(0,5) # Retinal progenitors
    n = np.append(n, 4.5 + np.arange(1,6)) # Photoreceptors
    n = np.append(n, 10 + np.arange(1,6)) # Amacrine/Horizontal cells
    n = np.append(n, 15.5 + np.arange(1,6)) # Retinal ganglion cells
    h_start = 1
    h_end = 21
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {'RPC' : '#DADADA', 'PR' : '#dcc360', 'HC_AC' : '#3DF591', 'RGC' : '#F53D59'}
    barColors = [
        pC['RPC'],pC['RPC'],pC['RPC'],pC['RPC'],pC['RPC'],
        pC['PR'],pC['PR'],pC['PR'],pC['PR'],pC['PR'],
        pC['HC_AC'],pC['HC_AC'],pC['HC_AC'],pC['HC_AC'],pC['HC_AC'],
        pC['RGC'],pC['RGC'],pC['RGC'],pC['RGC'],pC['RGC'],
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Nerli2022(geneSymbol, ax=ax)
    return pH

def formatBarPlot_Nerli2022(geneSymbol, ax=None):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks([2,7.5,13,18.5])
    ax.set_xticklabels(['RPC','Photo','HC/AC','RGC']);
    ax.set_xticks([0,1,2,3,4,5.5,6.5,7.5,8.5,9.5,11,12,13,14,15,16.5,17.5,18.5,19.5,20.5], minor=True)
    ax.set_ylabel('counts (norm.)', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

"""heatmaps"""

def heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """Creates a heatmap for a list of genes
    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels for the columns (overriding for custom color)
        groupsN    : number of each subtype
        groupsColors: colors for each subtype
        groupsLabels: label for each subtype
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """
    fontTicks = font_manager.FontProperties(size=36)
    fontLabels = font_manager.FontProperties(size=22)
    fontTitle = font_manager.FontProperties(size=28)

    if data.shape[0]==0:
        data = np.ones([2,data.shape[1]])
        row_labels = np.array(['not found', 'not found'])
    if not ax:
        ax = plt.gca()
    # Plot the heatmap
    # perceptually responsible colormaps are: inferno, viridis, plasma, magma, cividis
    im = ax.imshow(data, cmap = "bone", **kwargs)
#     im = ax.imshow(data, cmap = "inferno", **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', shrink=.75, pad=0.05, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, fontproperties=fontLabels, rotation=0, ha="right", va="center",rotation_mode="anchor")
    cbar.ax.tick_params(labelsize=22)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontproperties=fontLabels)
    ax.set_yticklabels(row_labels, fontproperties=fontLabels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,labeltop=True, labelbottom=False)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_yticklabels(), rotation=30, ha="right", va="center",rotation_mode="anchor")
    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_linewidth(.5)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(which="minor", bottom=False, left=False)
    # Custom grid according to photoreceptor subtype
    for h in np.arange(-.5,data.shape[0]+.5):
        ax.axhline(y = h, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')
    for v in np.arange(-.5,data.shape[1]+.5):
        ax.axvline(x = v, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')

    ax.axvline(x = -0.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = np.sum(groupsN)-.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    for i in np.arange(groupsN.shape[0]):
        ax.axvline(x = np.sum(groupsN[:i+1])-.5, color = 'white', linewidth = 3, alpha = 1, solid_capstyle='butt')
        ax.plot([np.sum(groupsN[:i])-.5,np.sum(groupsN[:i+1])-.5], [-.5,-.5], '-', lw=8, color = groupsColors[i], solid_capstyle='butt')
        ax.plot([np.sum(groupsN[:i])-.5,np.sum(groupsN[:i+1])-.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=8, color = groupsColors[i], solid_capstyle='butt')
        ax.text(((np.sum(groupsN[:i])+np.sum(groupsN[:i+1]))/2)-.5, -1.0, groupsLabels[i], color = groupsColors[i], horizontalalignment='center', fontproperties=fontTicks)
    return im, cbar

def heatmap(heatmapData, ax=None, pC=None, norm=False):
    """Main call for heatmap for data from Angueyra et al. (2021)
    Arguments:
        heatmapData : pandas dataframe containing expression data to be plotted
        pC : dict with custom photoreceptor colors
        norm : boolean that determines if plotting is raw data or row-normalized
    Returns:
        hmH : heatmap handle
        cbH : colorbar handle
    """
    genenames = heatmapData['symbol'].values
    data = heatmapData.iloc[0:,7:37].values #in FPKM
    cbarlabel = "FPKM"
    if norm:
        data = heatmapData.iloc[0:,7:37].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. FPKM"
    groupsN = np.array([6,5,6,7,6])
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    groupsColors = np.array([pC['r'],pC['u'],pC['s'],pC['m'],pC['l']])
    groupsLabels = np.array(['Rods','UV','S','M','L'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH


def heatmap_Ogawa2021(heatmapData, ax=None, pC=None, pctPlot=False, norm=False):
    """Main call for heatmap for reanalyzed data from Ogawa et al. (2021) (https://doi.org/10.1038/s41598-021-96837-z)
    Arguments:
        heatmapData : pandas dataframe containing expression data to be plotted
        pC : dict with custom photoreceptor colors
        norm : boolean that determines if plotting is raw data or row-normalized
    Returns:
        hmH : heatmap handle
        cbH : colorbar handle
    """
    genenames = heatmapData['symbol'].values
    cbarlabel = "avg."
    delta = 0
    if pctPlot:
        delta = 9
        cbarlabel = "%"
    data = heatmapData.iloc[0:,2+delta:10+delta].values #avg. counts or percent expression
    if norm:
        data = heatmapData.iloc[0:,2+delta:10+delta].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. " + cbarlabel
    groupsN = np.array([1,1,1,1,1,1,1,1])
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    groupsColors = np.array([pC['r'],pC['u'],pC['s'],pC['m'],pC['l'],pC['m4'],pC['onBC'],pC['offBC']])
    groupsLabels = np.array(['Rods','UV','S','M','L', 'M4','B$_{on}$','B$_{off}$'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH

def heatmap_Hoang2020(heatmapData, ax=None, pC=None, pctPlot=False, norm=False):
    """Main call for heatmap for reanalyzed data from Hoang et al. (2020)
    Arguments:
        heatmapData : pandas dataframe containing expression data to be plotted
        pC : dict with custom photoreceptor colors
        norm : boolean that determines if plotting is raw data or row-normalized
    Returns:
        hmH : heatmap handle
        cbH : colorbar handle
    """
    genenames = heatmapData['symbol'].values
    cbarlabel = "avg."
    delta = 0
    if pctPlot:
        delta = 8
        cbarlabel = "%"
    data = heatmapData.iloc[0:,2+delta:9+delta].values #avg. counts or percent expression
    if norm:
        data = heatmapData.iloc[0:,2+delta:9+delta].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. " + cbarlabel
    groupsN = np.array([1,1,1,1,1,1,1])
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    groupsColors = np.array([pC['r'],pC['u'],pC['s'],pC['m'],pC['m'],pC['m4'],pC['l']])
    groupsLabels = np.array(['Rods','UV','S','M','M3', 'M4','L'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH

def heatmap_Sun2018(heatmapData, ax=None, pC=None, norm=False):
    """Main call for heatmap for data from Sun, Galicia and Stenkamp (2018)
    Arguments:
        heatmapData : pandas dataframe containing expression data to be plotted
        pC : dict with custom photoreceptor colors
        norm : boolean that determines if plotting is raw data or row-normalized
    Returns:
        hmH : heatmap handle
        cbH : colorbar handle
    """
    genenames = heatmapData['symbol'].values
    data = heatmapData.iloc[0:,7:16].values #in cpm
    cbarlabel = "cpm"
    if norm:
        data = heatmapData.iloc[0:,7:16].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. cpm"
    groupsN = np.array([4,4])
    if not pC:
        pC = {'r' : '#747474','u' : '#B540B7','s' : '#4669F2','m' : '#04CD22','l' : '#CC2C2A',
        'm4': '#cdcd04','onBC': '#ccf2ff','offBC': '#663d00'}
    groupsColors = np.array([pC['r'],pC['m4']])
    groupsLabels = np.array(['Rods','notRods'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH

def heatmap_Nerli2022(heatmapData, ax=None, pC=None, norm=False):
    """Main call for heatmap for data from Nerli et al. (2022)
    Arguments:
        heatmapData : pandas dataframe containing expression data to be plotted
        pC : dict with custom photoreceptor colors
        norm : boolean that determines if plotting is raw data or row-normalized
    Returns:
        hmH : heatmap handle
        cbH : colorbar handle
    """
    genenames = heatmapData['symbol'].values
    data = heatmapData.iloc[0:,1:21].values #in FPKM
    cbarlabel = "Counts (norm.)"
    if norm:
        data = heatmapData.iloc[0:,1:21].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. FPKM"
    groupsN = np.array([5,5,5,5])
    if not pC:
        pC = {'RPC' : '#DADADA', 'PR' : '#dcc360', 'HC_AC' : '#3DF591', 'RGC' : '#F53D59'}
    groupsColors = np.array([pC['RPC'],pC['PR'],pC['HC_AC'],pC['RGC']])
    groupsLabels = np.array(['RPC','Photo','HC/AC','RGC'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH
