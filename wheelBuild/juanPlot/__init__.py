"""juanPlot

    Collection of plotting and formatting functions for Angueyra et al., 2022
    Created: September 2021 (Angueyra)
    Updated: October 2022 (Angueyra)
        - Added plotting routines for _Dev notebook
    Updated: November 2023 (Angueyra)
        - Added plotting routines for Yoshimatsu et al, acute zone data and notebook
    Updated: June 2024 (Angueyra)
        - Added plotting routines for Xu et al, 2020 retinal development data and added to Dev notebook
        
        To build wheel:
        copy .py function as "__init__.py"
        create directory "juanPlot"
        create "test.py" inside "juanPlot"
            from juanPlot.__init__ import *

            def func_test():
                print("Successfully Imported test.py file")
        create "setup.py" inside "juanPlot"
            from setuptools import setup

            setup(
                name='juanPlot',
                version='0a3',
                packages=['juanPlot'],
                options={"bdist_wheel": {"universal": True}},
            )
        run
            ```python setup.py bdist_wheel;```
        Then
            ```cp ./dist/juanPlot-0a3-py2.py3-none-any.whl ~/Documents/Repositories/drRNAseq/content```
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
            'BC_adult' : '#1B98C3', # Bipolar cell (mature)
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
    
def plotBars_Hoang2020_PhotoDev(barData, geneSymbol, ax=None, pC=None, pctPlot=False):
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
        pC = {
            'RPC' : "#DADADA",
            'PRP' : "#dfdac8",
            'Cle' : '#dacd9a',
            'Clm' : '#dcc360',
            'Cll' : '#cca819',
            'C' : '#ffd429',
            'R' : '#7d7d7d',
        }
    barColors = [
        pC['RPC'],pC['PRP'],
        pC['Cle'],pC['Clm'],pC['Cll'],
        pC['C'],pC['R']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Hoang2020_PhotoDev(geneSymbol, ax=ax, pctPlot=pctPlot)
    return pH

def formatBarPlot_Hoang2020_PhotoDev(geneSymbol, ax=None, pctPlot=False):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks(np.arange(1,8))
    ax.set_xticklabels(['RPC','PRPC','Cone$_{larval-early}$','Cone$_{larval-mid}$','Cone$_{larval-late}$','Cone','Rod']);
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", va="center",rotation_mode="anchor")
    ax.set_ylabel('avg. counts', fontproperties=fontLabels)
    if pctPlot:
        ax.set_ylabel('% expressing', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

def plotBars_Hoang2020_PhotoDev_withcontaminatedlarvalrods(barData, geneSymbol, ax=None, pC=None, pctPlot=False):
    """Creates a bar plot for a single gene for data from Hoang et al. (2020) (https://doi.org/10.1126/science.abb8598)
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
        pC = {
            'RPC' : "#DADADA",
            'PRP' : "#dfdac8",
            'Cle' : '#dacd9a',
            'Clm' : '#dcc360',
            'Cll' : '#cca819',
            'C' : '#ffd429',
            'Rll' : '#a3a3a3',
            'R' : '#7d7d7d',
        }
    barColors = [
        pC['RPC'],pC['PRP'],
        pC['Cle'],pC['Clm'],pC['Cll'],
        pC['C'],
        pC['Rll'],pC['R']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Hoang2020_PhotoDev(geneSymbol, ax=ax, pctPlot=pctPlot)
    return pH

def formatBarPlot_Hoang2020_PhotoDev_withcontaminatedlarvalrods(geneSymbol, ax=None, pctPlot=False):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks(np.arange(1,9))
    ax.set_xticklabels(['RPC','PRPC','Cone$_{larval-early}$','Cone$_{larval-mid}$','Cone$_{larval-late}$','Cone$','Rod$_{larval-late}$','Rod']);
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", va="center",rotation_mode="anchor")
    ax.set_ylabel('avg. counts', fontproperties=fontLabels)
    if pctPlot:
        ax.set_ylabel('% expressing', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)


def plotBars_Xu2020_RetDev(barData, geneSymbol, ax=None, pC=None, pctPlot=False):
    """Creates a bar plot for a single gene for data from Xu et al. (2020) (https://doi.org/10.1242/dev.185660)
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,81)
    delta = 0
    if pctPlot:
        delta = 82
    h_start = 2 + delta
    h_end = 82 + delta
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {
            'RPC24h_prim' : '#E9E9E9',
            'RPC24h' : '#D2D2D2',
            'RPC24h_D' : '#D2B6AF',
            'RPC24h_V' : '#AFAFD2',
            'RPC24h_T' : '#AFD2B0',
            'RPC24h_N' : '#D2D1AF',
            'RPC24h_neuro' : '#F5E5E8',
            'RPC36h_prim' : '#E9E9E9',
            'RPC36h' : '#D2D2D2',
            'RPC36h_D' : '#D2B6AF',
            'RPC36h_V' : '#AFAFD2',
            'RPC36h_T' : '#AFD2B0',
            'RPC36h_N' : '#D2D1AF',
            'RPC36h_hc' : '#C5B39D',
            'RPC36h_neuro' : '#F5E5E8',
            'RPC48h_prim' : '#E9E9E9',
            'RPC48h' : '#D2D2D2',
            'RPC48h_neuro' : '#F5E5E8',
            'RPC48h_neuroMix' : '#F5E5E8',
            'RPC48h_rgc' : '#EA9D81',
            'RPC48h_MG' : '#F5CCD2',
            'PRPC48h' : '#dfdac8',
            'Photo48h' : '#dacd9a',
            'HC48h' : '#FCCAA5',
            'BC48h' : '#7DAAAF',
            'BC48h_on' : '#8398B1',
            'BC48h_off' : '#B16D8A',
            'AC48h' : '#95F5C1',
            'AC48h_gaba' : '#91F5DA',
            'AC48h_sac' : '#78AD93',
            'RGC48h' : '#F53D59',
            'RPC72h_prim' : '#E9E9E9',
            'RPC72h_neuro' : '#F5E5E8',
            'RPC72h_bc' : '#CCF2FF',
            'RPC72h_rgc' : '#F5CCD2',
            'PRPC72h' : '#DFDAC8',
            'Rod72h_early' : '#A3A3A3',
            'Rod72h' : '#7D7D7D',
            'Cone72h_early' : '#DACD9A',
            'Cone72h_earlyL' : '#D69F9E',
            'Cone72h' : '#DCC360',
            'Cone72h_UV' : '#B778B9',
            'Cone72h_S' : '#7183CC',
            'Cone72h_M' : '#57CB69',
            'Cone72h_L' : '#C67271',
            'Cone72h_lateL' : '#CE524F',
            'HC72h_early' : '#FCCBA8',
            'HC72h' : '#FCA668',
            'HC72h_H1' : '#C59965',
            'HC72h_H2H3' : '#F7909C',
            'BC72h_on' : '#8398B1',
            'AC72h' : '#3DF591',
            'AC72h_gly' : '#91F5DA',
            'AC72h_gaba' : '#95F5C1',
            'AC72h_sac' : '#78AD93',
            'RGC72h' : '#F53D59',
            'MG72h' : '#EA9D81',
            'RPC336h_prim' : '#E9E9E9',
            'RPC336h_neuro' : '#F5E5E8',
            'RPC336h_bc' : '#CCF2FF',
            'RPC336h_ac' : '#B5F5D2',
            'RPC336h_rgc' : '#F5CCD2',
            'PRPC336h' : '#DFDAC8',
            'Rod336h' : '#7D7D7D',
            'Cone336h' : '#FFD429',
            'Cone336h_U' : '#B778B9',
            'Cone336h_S' : '#7183CC',
            'Cone336h_M' : '#57CB69',
            'Cone336h_L' : '#C67271',
            'Cone336h_lateL' : '#CE524F',
            'HC336h' : '#FCA668',
            'BC336h_on' : '#8398B1',
            'BC336h_off' : '#B16D8A',
            'AC336h_gaba' : '#3DF5C3',
            'AC336h_gly' : '#56F53D',
            'AC336h_ngng' : '#89AC25',
            'AC336h_sac' : '#78AD93',
            'AC336h_dac' : '#A4D1D8',
            'RGC336h' : '#F53D59',
            'MG336h' : '#EA9D81',
        }
    barColors = [
        pC['RPC24h_prim'],pC['RPC24h'],pC['RPC24h_D'],pC['RPC24h_V'],pC['RPC24h_T'],pC['RPC24h_N'],pC['RPC24h_neuro'],
        pC['RPC36h_prim'],pC['RPC36h'],pC['RPC36h_D'],pC['RPC36h_V'],pC['RPC36h_T'],pC['RPC36h_N'],pC['RPC36h_hc'],pC['RPC36h_neuro'],
        pC['RPC48h_prim'],pC['RPC48h'],pC['RPC48h_neuro'],pC['RPC48h_neuroMix'],pC['RPC48h_rgc'],pC['RPC48h_MG'],pC['PRPC48h'],pC['Photo48h'],
        pC['HC48h'],pC['BC48h'],pC['BC48h_on'],pC['BC48h_off'],pC['AC48h'],pC['AC48h_gaba'],pC['AC48h_sac'],pC['RGC48h'],
        pC['RPC72h_prim'],pC['RPC72h_neuro'],pC['RPC72h_bc'],pC['RPC72h_rgc'],pC['PRPC72h'],pC['Rod72h_early'],pC['Rod72h'],pC['Cone72h_early'],
        pC['Cone72h_earlyL'],pC['Cone72h'],pC['Cone72h_UV'],pC['Cone72h_S'],pC['Cone72h_M'],pC['Cone72h_L'],pC['Cone72h_lateL'],pC['HC72h_early'],
        pC['HC72h'],pC['HC72h_H1'],pC['HC72h_H2H3'],pC['BC72h_on'],pC['AC72h'],pC['AC72h_gly'],pC['AC72h_gaba'],pC['AC72h_sac'],pC['RGC72h'],pC['MG72h'],
        pC['RPC336h_prim'],pC['RPC336h_neuro'],pC['RPC336h_bc'],pC['RPC336h_ac'],pC['RPC336h_rgc'],pC['PRPC336h'],pC['Rod336h'],pC['Cone336h'],
        pC['Cone336h_U'],pC['Cone336h_S'],pC['Cone336h_M'],pC['Cone336h_L'],pC['Cone336h_lateL'],pC['HC336h'],pC['BC336h_on'],pC['BC336h_off'],
        pC['AC336h_gaba'],pC['AC336h_gly'],pC['AC336h_ngng'],pC['AC336h_sac'],pC['AC336h_dac'],pC['RGC336h'],pC['MG336h']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Xu2020_RetDev(geneSymbol, ax=ax, pctPlot=pctPlot)
    return pH

def formatBarPlot_Xu2020_RetDev(geneSymbol, ax=None, pctPlot=False):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks(np.arange(1,81))
    ax.set_xticklabels(["RPC24h_prim","RPC24h","RPC24h_D","RPC24h_V","RPC24h_T","RPC24h_N","RPC24h_neuro","RPC36h_prim","RPC36h","RPC36h_D","RPC36h_V","RPC36h_T","RPC36h_N","RPC36h_hc","RPC36h_neuro","RPC48h_prim","RPC48h","RPC48h_neuro","RPC48h_neuroMix","RPC48h_rgc","RPC48h_MG","PRPC48h","Photo48h","HC48h","BC48h","BC48h_on","BC48h_off","AC48h","AC48h_gaba","AC48h_sac","RGC48h","RPC72h_prim","RPC72h_neuro","RPC72h_bc","RPC72h_rgc","PRPC72h","Rod72h_early","Rod72h","Cone72h_early","Cone72h_earlyL","Cone72h","Cone72h_UV","Cone72h_S","Cone72h_M","Cone72h_L","Cone72h_lateL","HC72h_early","HC72h","HC72h_H1","HC72h_H2H3","BC72h_on","AC72h","AC72h_gly","AC72h_gaba","AC72h_sac","RGC72h","MG72h","RPC336h_prim","RPC336h_neuro","RPC336h_bc","RPC336h_ac","RPC336h_rgc","PRPC336h","Rod336h","Cone336h","Cone336h_U","Cone336h_S","Cone336h_M","Cone336h_L","Cone336h_lateL","HC336h","BC336h_on","BC336h_off","AC336h_gaba","AC336h_gly","AC336h_ngng","AC336h_sac","AC336h_dac","RGC336h","MG336h"]);
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

def plotBars_Yoshimatsu2020(barData, geneSymbol, ax=None, pC=None):
    """Creates a bar plot for a single gene for data from Yoshimatsu et al. (2020) (https://doi.org/10.1016/j.neuron.2020.04.021)
    Arguments:
        barData         : a 1D numpy array
        geneSymbol      : gene Symbol for plot title
        ax              : pyplot axis handle
        pC              : photoreceptor colors for plotting
    """
    n = np.arange(1,5) # non-strike zone UV cones
    n = np.append(n, 5 + np.arange(1,5)) # strike zone UV cones
    h_start = 1
    h_end = 9
    h = barData.iloc[0,h_start:h_end].to_numpy()
    # color array for bar plot
    if not pC:
        pC = {'u' : '#B540B7','u2' : '#B789A7'}
    barColors = [
        pC['u2'],pC['u2'],pC['u2'],pC['u2'],
        pC['u'],pC['u'],pC['u'],pC['u']
    ]
    if not ax:
        ax = plt.gca()
    pH = ax.bar(n, h, width=0.8, bottom=None, align='center', data=None, color=barColors)
    formatBarPlot_Yoshimatsu2020(geneSymbol, ax=ax)
    return pH

def formatBarPlot_Yoshimatsu2020(geneSymbol, ax=None):
    if not ax:
        ax = plt.gca()
    [fontTicks, fontLabels, fontTitle] = defaultFonts(ax = ax);
    ax.set_xticks([2.5,7.5])
    ax.set_xticklabels(['UV$_{non-sz}$','UV$_{sz}$']);
    ax.set_ylabel('cpm', fontproperties=fontLabels)
    ax.set_title(geneSymbol, fontproperties=fontTitle)

"""heatmaps"""

def heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=None,
            cbar_kw={}, cbarlabel="", groupRotation=0, **kwargs):
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
        if groupRotation==0:
            ax.text(((np.sum(groupsN[:i])+np.sum(groupsN[:i+1]))/2)-.5, -1.0, groupsLabels[i], color = groupsColors[i], ha='center', fontproperties=fontTicks, rotation=groupRotation)
        else:
            ax.text(((np.sum(groupsN[:i])+np.sum(groupsN[:i+1]))/2)-.5, -1.0, groupsLabels[i], color = groupsColors[i], ha='left', va='bottom', rotation_mode="anchor", fontproperties=fontTicks, rotation=groupRotation)
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

def heatmap_Yoshimatsu2020(heatmapData, ax=None, pC=None, norm=False):
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
    data = heatmapData.iloc[0:,1:9].values #in cpm
    cbarlabel = "cpm"
    if norm:
        data = heatmapData.iloc[0:,1:9].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. cpm"
    groupsN = np.array([4,4])
    if not pC:
        pC = {'u' : '#B540B7','u2' : '#B789A7'}
    groupsColors = np.array([pC['u2'],pC['u']])
    groupsLabels = np.array(['UV$_{non-sz}$','UV$_{sz}$'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH

def heatmap_Hoang2020_Ret(heatmapData, ax=None, pC=None, pctPlot=False, norm=False):
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
        delta = 19
        cbarlabel = "%"
    data = heatmapData.iloc[0:,2+delta:19+delta].values #avg. counts or percent expression
    if norm:
        data = heatmapData.iloc[0:,2+delta:19+delta].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. " + cbarlabel
    groupsN = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    if not pC:
        pC = {'RPC' : '#DADADA', 'PRPC' : '#dfdac8', 'Cones_larval' : '#dcc360', 'Cones_adult' : '#ffd429', 'Rods' : '#7d7d7d',
              'HC' : '#FC7715', 'BC_larval' : '#ccf2ff', 'BC_adult' : '#663d00', 'AC_larval' : '#3DF591', 'ACgaba' : '#3DF5C3', 'ACgly' : '#56F53D',
              'RGC_larval' : '#F53D59', 'RGC_adult' : '#BB0622', 'MGi' : '#EA9D81', 'MG1' : '#A2644E', 'MG2' : '#7E4835', 'MG3' : '#613728'}
    groupsColors = np.array([pC['RPC'],pC['PRPC'],pC['Cones_larval'],pC['Cones_adult'],pC['Rods'],
                             pC['HC'],pC['BC_larval'],pC['BC_adult'],pC['AC_larval'],pC['ACgaba'],pC['ACgly'],
                             pC['RGC_larval'],pC['RGC_adult'],pC['MGi'],pC['MG1'],pC['MG2'],pC['MG3'],])
    groupsLabels = np.array( ['RPC','PRPC','C$_{larval}$','C$_{adult}$','R$_{ods}$',
                              'HC','BC$_{larval}$','BC$_{adult}$','AC$_{larval}$','AC$_{GABA}$','AC$_{Gly}$',
                              'RGC$_{larval}$','RGC$_{adult}$','MGi','MG1','MG2','MG3'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, groupRotation = 45, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH

def heatmap_Hoang2020_PhotoDev(heatmapData, ax=None, pC=None, pctPlot=False, norm=False):
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
        pC = {
            'RPC' : "#DADADA",
            'PRP' : "#dfdac8",
            'Cle' : '#dacd9a',
            'Clm' : '#dcc360',
            'Cll' : '#cca819',
            'C' : '#ffd429',
            'R' : '#7d7d7d',
         }
    groupsColors = np.array([
        pC['RPC'],pC['PRP'],
        pC['Cle'],pC['Clm'],pC['Cll'],
        pC['C'],pC['R']
    ])
    groupsLabels = np.array(['RPC','PRP','Cone$_{larval-early}$','Cone$_{larval-mid}$','Cone$_{larval-late}$','Cone','Rod'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, groupRotation = 45, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH

def heatmap_Hoang2020_PhotoDev_withcontaminatedlarvalrods(heatmapData, ax=None, pC=None, pctPlot=False, norm=False):
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
        delta = 9
        cbarlabel = "%"
    data = heatmapData.iloc[0:,2+delta:10+delta].values #avg. counts or percent expression
    if norm:
        data = heatmapData.iloc[0:,2+delta:10+delta].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. " + cbarlabel
    groupsN = np.array([1,1,1,1,1,1,1,1])

    if not pC:
        pC = {
            'RPC' : "#DADADA",
            'PRP' : "#dfdac8",
            'Cle' : '#dacd9a',
            'Clm' : '#dcc360',
            'Cll' : '#cca819',
            'C' : '#ffd429',
            'Rll' : '#a3a3a3',
            'R' : '#7d7d7d',
         }
    groupsColors = np.array([
        pC['RPC'],pC['PRP'],
        pC['Cle'],pC['Clm'],pC['Cll'],
        pC['C'],
        pC['Rll'],pC['R']
    ])
    groupsLabels = np.array(['RPC','PRPC','Cone$_{larval-early}$','Cone$_{larval-mid}$','Cone$_{larval-late}$','Cone','Rod$_{larval-late}$','Rod'])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, groupRotation = 45, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH

def heatmap_Xu2020_RetDev(heatmapData, ax=None, pC=None, pctPlot=False, norm=False):
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
        delta = 81
        cbarlabel = "%"
    data = heatmapData.iloc[0:,2+delta:82+delta].values #avg. counts or percent expression
    if norm:
        data = heatmapData.iloc[0:,2+delta:82+delta].apply(lambda x: x/x.max(), axis=1).values #normalized by max
        cbarlabel = "norm. " + cbarlabel
    groupsN = np.array([
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        ])

    if not pC:
        pC = {
            'RPC24h_prim' : '#E9E9E9',
            'RPC24h' : '#D2D2D2',
            'RPC24h_D' : '#D2B6AF',
            'RPC24h_V' : '#AFAFD2',
            'RPC24h_T' : '#AFD2B0',
            'RPC24h_N' : '#D2D1AF',
            'RPC24h_neuro' : '#F5E5E8',
            'RPC36h_prim' : '#E9E9E9',
            'RPC36h' : '#D2D2D2',
            'RPC36h_D' : '#D2B6AF',
            'RPC36h_V' : '#AFAFD2',
            'RPC36h_T' : '#AFD2B0',
            'RPC36h_N' : '#D2D1AF',
            'RPC36h_hc' : '#C5B39D',
            'RPC36h_neuro' : '#F5E5E8',
            'RPC48h_prim' : '#E9E9E9',
            'RPC48h' : '#D2D2D2',
            'RPC48h_neuro' : '#F5E5E8',
            'RPC48h_neuroMix' : '#F5E5E8',
            'RPC48h_rgc' : '#EA9D81',
            'RPC48h_MG' : '#F5CCD2',
            'PRPC48h' : '#dfdac8',
            'Photo48h' : '#dacd9a',
            'HC48h' : '#FCCAA5',
            'BC48h' : '#7DAAAF',
            'BC48h_on' : '#8398B1',
            'BC48h_off' : '#B16D8A',
            'AC48h' : '#95F5C1',
            'AC48h_gaba' : '#91F5DA',
            'AC48h_sac' : '#78AD93',
            'RGC48h' : '#F53D59',
            'RPC72h_prim' : '#E9E9E9',
            'RPC72h_neuro' : '#F5E5E8',
            'RPC72h_bc' : '#CCF2FF',
            'RPC72h_rgc' : '#F5CCD2',
            'PRPC72h' : '#DFDAC8',
            'Rod72h_early' : '#A3A3A3',
            'Rod72h' : '#7D7D7D',
            'Cone72h_early' : '#DACD9A',
            'Cone72h_earlyL' : '#D69F9E',
            'Cone72h' : '#DCC360',
            'Cone72h_UV' : '#B778B9',
            'Cone72h_S' : '#7183CC',
            'Cone72h_M' : '#57CB69',
            'Cone72h_L' : '#C67271',
            'Cone72h_lateL' : '#CE524F',
            'HC72h_early' : '#FCCBA8',
            'HC72h' : '#FCA668',
            'HC72h_H1' : '#C59965',
            'HC72h_H2H3' : '#F7909C',
            'BC72h_on' : '#8398B1',
            'AC72h' : '#3DF591',
            'AC72h_gly' : '#91F5DA',
            'AC72h_gaba' : '#95F5C1',
            'AC72h_sac' : '#78AD93',
            'RGC72h' : '#F53D59',
            'MG72h' : '#EA9D81',
            'RPC336h_prim' : '#E9E9E9',
            'RPC336h_neuro' : '#F5E5E8',
            'RPC336h_bc' : '#CCF2FF',
            'RPC336h_ac' : '#B5F5D2',
            'RPC336h_rgc' : '#F5CCD2',
            'PRPC336h' : '#DFDAC8',
            'Rod336h' : '#7D7D7D',
            'Cone336h' : '#FFD429',
            'Cone336h_U' : '#B778B9',
            'Cone336h_S' : '#7183CC',
            'Cone336h_M' : '#57CB69',
            'Cone336h_L' : '#C67271',
            'Cone336h_lateL' : '#CE524F',
            'HC336h' : '#FCA668',
            'BC336h_on' : '#8398B1',
            'BC336h_off' : '#B16D8A',
            'AC336h_gaba' : '#3DF5C3',
            'AC336h_gly' : '#56F53D',
            'AC336h_ngng' : '#89AC25',
            'AC336h_sac' : '#78AD93',
            'AC336h_dac' : '#A4D1D8',
            'RGC336h' : '#F53D59',
            'MG336h' : '#EA9D81',
        }
    groupsColors = np.array([
        pC['RPC24h_prim'],pC['RPC24h'],pC['RPC24h_D'],pC['RPC24h_V'],pC['RPC24h_T'],pC['RPC24h_N'],pC['RPC24h_neuro'],
        pC['RPC36h_prim'],pC['RPC36h'],pC['RPC36h_D'],pC['RPC36h_V'],pC['RPC36h_T'],pC['RPC36h_N'],pC['RPC36h_hc'],pC['RPC36h_neuro'],
        pC['RPC48h_prim'],pC['RPC48h'],pC['RPC48h_neuro'],pC['RPC48h_neuroMix'],pC['RPC48h_rgc'],pC['RPC48h_MG'],pC['PRPC48h'],pC['Photo48h'],
        pC['HC48h'],pC['BC48h'],pC['BC48h_on'],pC['BC48h_off'],pC['AC48h'],pC['AC48h_gaba'],pC['AC48h_sac'],pC['RGC48h'],
        pC['RPC72h_prim'],pC['RPC72h_neuro'],pC['RPC72h_bc'],pC['RPC72h_rgc'],pC['PRPC72h'],pC['Rod72h_early'],pC['Rod72h'],pC['Cone72h_early'],
        pC['Cone72h_earlyL'],pC['Cone72h'],pC['Cone72h_UV'],pC['Cone72h_S'],pC['Cone72h_M'],pC['Cone72h_L'],pC['Cone72h_lateL'],pC['HC72h_early'],
        pC['HC72h'],pC['HC72h_H1'],pC['HC72h_H2H3'],pC['BC72h_on'],pC['AC72h'],pC['AC72h_gly'],pC['AC72h_gaba'],pC['AC72h_sac'],pC['RGC72h'],pC['MG72h'],
        pC['RPC336h_prim'],pC['RPC336h_neuro'],pC['RPC336h_bc'],pC['RPC336h_ac'],pC['RPC336h_rgc'],pC['PRPC336h'],pC['Rod336h'],pC['Cone336h'],
        pC['Cone336h_U'],pC['Cone336h_S'],pC['Cone336h_M'],pC['Cone336h_L'],pC['Cone336h_lateL'],pC['HC336h'],pC['BC336h_on'],pC['BC336h_off'],
        pC['AC336h_gaba'],pC['AC336h_gly'],pC['AC336h_ngng'],pC['AC336h_sac'],pC['AC336h_dac'],pC['RGC336h'],pC['MG336h']
    ])
    groupsLabels = np.array(["RPC24h_prim","RPC24h","RPC24h_D","RPC24h_V","RPC24h_T","RPC24h_N","RPC24h_neuro","RPC36h_prim","RPC36h","RPC36h_D","RPC36h_V","RPC36h_T","RPC36h_N","RPC36h_hc","RPC36h_neuro","RPC48h_prim","RPC48h","RPC48h_neuro","RPC48h_neuroMix","RPC48h_rgc","RPC48h_MG","PRPC48h","Photo48h","HC48h","BC48h","BC48h_on","BC48h_off","AC48h","AC48h_gaba","AC48h_sac","RGC48h","RPC72h_prim","RPC72h_neuro","RPC72h_bc","RPC72h_rgc","PRPC72h","Rod72h_early","Rod72h","Cone72h_early","Cone72h_earlyL","Cone72h","Cone72h_UV","Cone72h_S","Cone72h_M","Cone72h_L","Cone72h_lateL","HC72h_early","HC72h","HC72h_H1","HC72h_H2H3","BC72h_on","AC72h","AC72h_gly","AC72h_gaba","AC72h_sac","RGC72h","MG72h","RPC336h_prim","RPC336h_neuro","RPC336h_bc","RPC336h_ac","RPC336h_rgc","PRPC336h","Rod336h","Cone336h","Cone336h_U","Cone336h_S","Cone336h_M","Cone336h_L","Cone336h_lateL","HC336h","BC336h_on","BC336h_off","AC336h_gaba","AC336h_gly","AC336h_ngng","AC336h_sac","AC336h_dac","RGC336h","MG336h"])
    if not ax:
        ax = plt.gca()
    hmH, cbH = heatmap_general(data, genenames, [], groupsN, groupsColors, groupsLabels, groupRotation = 45, ax=ax, cbarlabel=cbarlabel)
    return hmH, cbH