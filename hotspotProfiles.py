#!/usr/bin/env python3
# script: hotspotProfiles.py
# author: Daniel Desiro'
script_usage="""
usage
    hotspotProfiles.py -vst <vRNAsite_table> -pfx <out_prefix> [options]

version
    hotspotProfiles.py 0.0.1 (alpha)

dependencies
    Python v3.9.7, Scipy v1.8.0, NumPy v1.22.2, Matplotlib v3.5.1

description
    Evaluates vRNAsite like interaction tables and creates hotspot profiles for 
    each individual sequence, showing regions with high interaction counts. 
    Example call: python hotspotProfiles.py -pfx example -vst example.tsv -ovr

################################################################

--prefix,-pfx
    output directory and prefix for result files

--vRNAsite,-vst
    input vRNAsite table

--candidatePeak,-cdp
    define minimum peak MFE value (default: -10.0)

--candidateMFE,-cdm
    define minimum MFE value (default: -10.0)

--candidateSPLASH,-cds
    define the minimum SPLASH count value (default: 0)

--candidateLength,-cdl
    minimum length for an interaction (default: 5)

--intra,-tra
    do intra molecular interactions from vRNAsite predictions (default: False)

--plotSize,-pls
    defines the plot size modifier (default: 1.0)

-xTickScale,-xts
    set x tick scale (default: 50)

-yTickScale,-yts
    set y tick scale (default: 0.10)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--namingExtension,-nex
    use this parameter as an extension for naming; possible are peak, mfe, 
    length, splash (default: "")

--loadMatrix,-mat
    load vRNAsite matrix data (default: )

--loadShapeMap,-lsh
    load shape map data (default: )

--quantilePeakMatrix,-qtp
    define threshold for peak matrix (default: 0.01)

--quantileSplashMatrix,-qts
    define threshold for splash matrix (default: 0.01)

--quantileDistribution,-pqd
    plot quantile distribution (default: False)

--quantileDistributionList,-pql
    set quantile distribution list (default: 0.1,0.05,0.01,0.005,0.001,0.0005,0.0001)

--rollingSHAPE,-rls
    use sliding window to calculate median SHAPE data (default: 1)

--rollingSHAPEIteration,-rli
    number of rolling iterations (default: 1)

--nameSHAPE,-nsh
    use alternative SHAPE name for plots (default: shape)

--subsetSequences,-sub
    only use a subset of sequences for analysis; these should match the aSeq and 
    bSeq columns and should be divided by commas (default: )

--createScatter,-csc
    creates a scatter plot for every segment (default: False)

--noPeakSmoothing,-nps
    deactivate peak smoothing (default: False)

--MFEasPeak,-mpk
    use mfe values instead of peak values (default: False)

--setTitle,-stt
    set title for plots (default: "")

--fixedLogMax,-lmx
    use fixed maximum log for plotting (default: 0)

reference
    D. Desirò, A. Borodavka and M. Marz.
    "HotspotProfiles: Prediction of hotspot regions capable of mediating the selective but flexible genome packaging of influenza A viruss."
    In Preparation, 2022.
    https://github.com/desiro/hotspotProfiles
"""

import argparse as ap
import sys
import os
import re
import time
import pickle
from operator import attrgetter, itemgetter
from itertools import product
from math import ceil, floor, log, isnan, log10
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from numpy import arange, mean, zeros, quantile, concatenate, median, std, cov, corrcoef, array, linspace, polyfit
from scipy.stats.stats import pearsonr, spearmanr, kendalltau
from scipy import stats
from scipy.signal import savgol_filter
from random import seed, randint

try:
    import multiprocessingPickle4
    import multiprocessing as mp
    ctxx = mp.get_context()
    ctxx.reducer = multiprocessingPickle4.MultiprocessingPickle4()
except:
    import multiprocessing as mp



################################################################################
## main
################################################################################

def main(opt):
    ############################################################################
    ## load fasta file
    opt["var_rvp"] = False
    time_s = getTime()
    print(f"Status: Read tables ...")
    tlist, segments, genomes, header, opt = readTable(**opt)
    time_s = getTime(time_s, f"Read table")
    ############################################################################
    ## create output folder
    opt = makeDir(**opt)
    ############################################################################
    ## load matrices
    if opt["var_mat"]:
        print(f"Status: Load matrices ...")
        mat_dict = loadMatrix(**opt)
        time_s = getTime(time_s, f"Load matrices")
    else:
        mat_dict = dict()
    ############################################################################
    ## load shape map
    if opt["var_lsh"]:
        print(opt["var_lsh"])
        print(f"Status: Load shape ...")
        shape_dict = loadShape(genomes, **opt)
        time_s = getTime(time_s, f"Load shape")
    else:
        shape_dict = dict()
    ############################################################################
    ## inter and intra
    #if opt["var_als"] == "both": analysis = ["inter","intra"]
    #else:                        analysis = [opt["var_als"]]
    genome_dict = dict()
    #for als in analysis:
    #als = opt["var_als"]
    if not opt["var_tra"]:
        als = "inter"
        als_tlist = removeIntra(tlist)
    else:
        als = "intra"
        als_tlist = removeInter(tlist)
        if not als_tlist:
            print("Error: There are no intra interactions in the dataset.")
            sys.exit()
    ############################################################################
    ## read viral titer table
    print(f"Status: Get {als} stats ...")
    getStats(als_tlist, als, **opt)
    time_s = getTime(time_s, f"Get stats")
    ############################################################################
    ## get link hotspots
    print(f"Status: Get {als} link hotspots ...")
    genome_dict = analyzeLinkHotSpots(als_tlist, segments, genomes, als, genome_dict, shape_dict, **opt)
    time_s = getTime(time_s, f"Get {als} link hotspots")
    if opt["var_mat"]:
        ############################################################################
        ## get matrix hotspots
        print(f"Status: Get {als} matrix hotspots ...")
        genome_dict = analyzeMatrixHotSpots(genomes, mat_dict, shape_dict, als, genome_dict, **opt)
        time_s = getTime(time_s, f"Get {als} matrix hotspots")
    ############################################################################
    ## extract hot-spot interactions
    ## for every two hot-spots get all interactions intersecting
    ## rate by biggest hot-spot
    if not opt["var_pqd"]:
        ########################################################################
        ## get hotspots
        print(f"Status: Do line plots ...")
        doPlots(genome_dict, genomes, segments, **opt)
        time_s = getTime(time_s, f"Do line plots")
    ############################################################################
    return opt["var_pfx"]




################################################################################
## functions
################################################################################

def readTable(**opt):
    ## read table file
    #--orientation,-ort
    #only allow specific orientations of the vRNPs; A = arbitrary, U = upright, R = reverse (default: A) (choices: A,U,R,UR,AU,AR,AUR)
    tlist, segments = list(), dict()
    with open(opt["var_vst"], "r") as tabin:
        header = next(tabin)
        header = header.strip().split()
        for i,line in enumerate(tabin):
            lnk = dict()
            line = line.strip().split()
            if len(header) != len(line):
                print(f"Error: Not the same number of header and list elements!")
                sys.exit()
            for name,item in zip(header,line):
                lnk[name] = item
            ## check for orientation
            if "aType" not in lnk.keys():
                asplit = lnk["aSeq"].split("-")
                if len(asplit) == 2:
                    lnk["aSeq"] = asplit[0]
                    lnk["aType"] = asplit[1]
                else:
                    lnk["aSeq"] = asplit[0]
                    lnk["aType"] = "-".join(asplit[1:2])
            if "bType" not in lnk.keys():
                bsplit = lnk["bSeq"].split("-")
                if len(bsplit) == 2:
                    lnk["bSeq"] = bsplit[0]
                    lnk["bType"] = bsplit[1]
                else:
                    lnk["bSeq"] = bsplit[0]
                    lnk["bType"] = "-".join(bsplit[1:2])
            if len(lnk["aType"].split("-")) == 1: aSequence = f"{lnk['aType']}-A"
            if len(lnk["bType"].split("-")) == 1: bSequence = f"{lnk['bType']}-A"
            aOrt = aSequence.split("-")[1]
            bOrt = bSequence.split("-")[1]
            #if opt["var_ort"]:
            #    if aOrt not in list(opt["var_ort"]) or bOrt not in list(opt["var_ort"]): continue
            if aOrt != "A" or bOrt != "A": continue
            ## create class and remove items
            lk = links(**lnk)
            if lk.aj < lk.ai:
                opt["var_rvp"] = True
                lk.ai, lk.aj, lk.bi, lk.bj = lk.alen-lk.ai, lk.alen-lk.aj+1, lk.blen-lk.bi, lk.blen-lk.bj+1 # removed -1
            else:
                lk.ai, lk.aj, lk.bi, lk.bj = lk.ai-1, lk.aj, lk.bi-1, lk.bj
            if opt["var_sub"]:
                if lk.aSeq not in opt["var_sub"].split(",") or lk.bSeq not in opt["var_sub"].split(","): continue
            sta = segments.get(lk.aSeq,(0,set()))[1]
            stb = segments.get(lk.bSeq,(0,set()))[1]
            sta.add(lk.aType)
            stb.add(lk.bType)
            segments[lk.aSeq] = (lk.alen,sta)
            segments[lk.bSeq] = (lk.blen,stb)
            if type(lk.peak) is int:
                lk.peak = float(lk.peak)
            if type(lk.mfe) is int:
                lk.mfe = float(lk.mfe)
            if lk.aj - lk.ai < opt["var_cdl"]: continue
            if lk.bj - lk.bi < opt["var_cdl"]: continue
            if type(lk.peak) is float:
                if lk.peak > opt["var_cdp"]: continue
            if type(lk.mfe) is float:
                if lk.mfe > opt["var_cdm"]: continue
            if "sp_mean" in lk.__dict__.keys():
                if type(lk.sp_mean) is int:
                    lk.sp_mean = float(lk.sp_mean)
                if type(lk.sp_mean) is float:
                    if lk.sp_mean < opt["var_cds"]: continue
            else:
                lk.sp_mean = float("nan")
            if opt["var_mpk"]:
                lk.peak, lk.mfe = lk.mfe, lk.peak
            tlist.append(lk)
            ##########################
            #TODO: remove optional interactions, only take one with highest peak
            ##########################
    #genomes = product(*[{(s,x) for x in t} for s,(l,t) in segments])
    #segments = {s:l for s,(l,t) in segments}
    tlist = sorted(tlist, key=attrgetter("peak"))
    if not opt["var_tra"] and len(segments) == 0:
        print(f"Error: Not enough sequences for inter molecular interactions!")
        sys.exit()
    segments, genomes = getGenomes(segments, **opt)
    return tlist, segments, genomes, header, opt

def getGenomes(segments, **opt):
    ## get all complete genomes
    genomes = list()
    #segments["SC35M_NA"] = (1461, {"WT","mutx"})
    segments = sorted(segments.items(), key=itemgetter(1,0), reverse=True)
    itr = [{(s,x) for x in t} for s,(l,t) in segments]
    #itr = [t for s,(l,t) in segments]
    segments = {s:l for s,(l,t) in segments}
    genomes = list(product(*itr))
    return segments, genomes

class links(object):
    def __init__(self, **data):
        self.__dict__.update((k,trans(v)) for k,v in data.items())
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    def plot(self, sep):
        ldat = sep.join([f"{var}" for key,var in vars(self).items()])
        return ldat

def makeDir(**opt):
    ## create directory
    dir_name, dir_base = opt["var_pfx"], opt["var_pfx"]
    if not opt["var_ovr"]:
        i = 1
        while os.path.isdir(dir_name):
            dir_name = f"{dir_base}_{i}"
            i += 1
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    opt["var_pfx"] = dir_name
    return opt

def removeIntra(tlist):
    ## remove intra interactions
    return [tx for tx in tlist if tx.aSeq != tx.bSeq]

def removeInter(tlist):
    ## remove inter interactions
    return [tx for tx in tlist if tx.aSeq == tx.bSeq and tx.aj < tx.bi]

def getTime(time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        print(f"Status: {name} finished at {time_c} in {time_e}")
    return time.time()

def getStats(tlist, als, **opt):
    ## plot peaks
    p_min = min([t.peak for t in tlist if not isnan(t.peak)])
    p_max = max([t.peak for t in tlist if not isnan(t.peak)])
    plotBarplot(tlist, p_min, p_max, f"peak_{als}", "average free energy in kcal/mol", "", **opt)
    ## plot average len
    all_len = [mean([lk.aj-lk.ai,lk.bj-lk.bi]) for lk in tlist]
    l_min = min(all_len)
    l_max = max(all_len)
    plotBarplot(tlist, l_min, l_max, f"len_{als}", "interaction length", "", **opt)

def plotList(tlist, name, **opt):
    # plot list of interactions
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    out_write = os.path.join(opt["var_pfx"], f"{out_name}{name}{getExt(**opt)}.tsv")
    with open(out_write, "w") as outfile:
        for i,lk in enumerate(tlist):
            if i == 0:
                header = lk.__dict__.keys()
                outfile.write("\t".join(header)+"\n")
            outfile.write(lk.plot("\t")+"\n")

def plotBarplot(links, lk_min, lk_max, dtype, xname, fname, **opt):
    ## pyplot matrix plot
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outfile = os.path.join(opt["var_pfx"], f"{out_name}_{dtype}")
    pd = getDistribution(links, lk_min, lk_max, dtype[:-6])
    s = opt["var_pls"]
    pz = plt.figure(figsize=(len(pd.keys())*0.7*s,4*s))
    x = arange(len(pd.keys()))
    h = [p for k,p in pd.items()]
    plt.bar(x, height=h, color="#ef6548")
    l = [f"{k}" for k in pd.keys()]
    plt.xticks(x, l)
    label = [f"{i}" for i in h]
    for i in range(len(pd.keys())):
        plt.text(x=x[i] , y=h[i], s=label[i], ha="center", va="bottom")
    plt.xlabel(xname)
    plt.ylabel("number of interacions")
    plt.subplots_adjust(bottom=0.2, top=1.2)
    pz.savefig(f"{outfile}{fname}{getExt(**opt)}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{outfile}{fname}{getExt(**opt)}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)

def plotDotplot(data_list, name, **opt):
    ## pyplot matrix plot
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outfile = os.path.join(opt["var_pfx"], f"{out_name}_{name}")
    pd = [quantile(data_list, float(i)) for i in opt["var_pql"].split(",")]
    s = opt["var_pls"]
    pz = plt.figure(figsize=(12*s,4*s))
    plt.plot(pd, color="#d53e4f", linewidth=5)
    plt.xlabel("quantiles")
    plt.ylabel(name)
    plt.xticks(arange(len(pd)), opt["var_pql"].split(","))
    plt.subplots_adjust(bottom=0.2, top=1.2)
    pz.savefig(f"{outfile}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{outfile}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)

def linePlot(data_list, name, **opt):
    ## pyplot matrix plot
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outfile = os.path.join(opt["var_pfx"], f"{out_name}_{name}")
    s = opt["var_pls"]
    pz = plt.figure(figsize=(ceil(len(data_list)/100)*s,4*s))
    plt.plot(data_list, color="#d53e4f", linewidth=5)
    plt.xlabel("nucleotides")
    plt.ylabel(name)
    #plt.xticks(arange(len(data_list)), opt["var_pql"].split(","))
    plt.subplots_adjust(bottom=0.2, top=1.2)
    #pz.savefig(f"{outfile}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{outfile}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)

def addArea(hots_dict, hots_minima, hots_maxima, **opt):
    ## add area for peak minima
    area_dict = dict()
    for seg,typ in hots_dict.keys():
        hots_min = sorted(hots_minima[(seg,typ)], key=itemgetter(1))
        hots_max = sorted(hots_maxima[(seg,typ)], key=itemgetter(1))
        data_dict = hots_dict.get((seg,typ), [])
        ## profile data
        p_list = [i for i in data_dict.get("peak_sum", [])]
        s_list = [i for i in data_dict.get("splash_sum", [])]
        ## area data
        if hots_min[0][1] < hots_max[0][1]:
            hots_max.insert(0,(p_list[0],0,0))
        if hots_min[-1][1] > hots_max[-1][1]:
            p = len(p_list)-1
            hots_max.append((p_list[-1],p,0))
        #test = [(pos1,hots_min[i+n][1],pos2,sum([p_list[pos] for pos in range(pos1,pos2+1)])) for i,((peak1,pos1,count1),(peak2,pos2,count2)) in enumerate(zip(hots_max[:-1],hots_max[1:]))]
        #print(test[:100])
        p_area = {hots_min[i][1]:sum([p_list[pos] for pos in range(pos1,pos2+1)]) for i,((peak1,pos1,count1),(peak2,pos2,count2)) in enumerate(zip(hots_max[:-1],hots_max[1:]))}
        s_area = {hots_min[i][1]:sum([s_list[pos] for pos in range(pos1,pos2+1)]) for i,((peak1,pos1,count1),(peak2,pos2,count2)) in enumerate(zip(hots_max[:-1],hots_max[1:]))}
        area_dict[(seg,typ)] = (p_area,s_area)
    return area_dict

def addOuterExtrema(hots_dict, hots_minima, hots_maxima):
    ## adds outer extrama
    for seg,typ in hots_dict.keys():
        hots_min = sorted(hots_minima[(seg,typ)], key=itemgetter(1))
        hots_max = sorted(hots_maxima[(seg,typ)], key=itemgetter(1))
        data_dict = hots_dict.get((seg,typ), [])
        ## profile data
        p_list = [i for i in data_dict.get("peak_sum", [])]
        ## area data
        if hots_min[0][1] < hots_max[0][1]:
            hots_max.insert(0,(p_list[0],0,0))
        if hots_min[-1][1] > hots_max[-1][1]:
            p = len(p_list)-1
            hots_max.append((p_list[-1],p,0))
        #print([b for a,b,c in hots_min])
        #print([b for a,b,c in hots_max])
        hots_maxima[(seg,typ)] = hots_max
    return hots_maxima

def plotScatter(hots_minima, hots_maxima, hots_dict, segments, name, **opt):
    var_pfx = opt["var_pfx"]
    opt["var_pfx"] = os.path.join(var_pfx, "correlation")
    makeDir(**opt)
    ## pyplot matrix plot
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    ## create colors
    colors = [(141,211,199),(255,255,179),(190,186,218),(251,128,114),(128,177,211),(253,180,98),(179,222,105),(252,205,229),(217,217,217),(188,128,189),(204,235,197),(255,237,111)]
    if len(segments) > 12:
        colors = addColors(len(segments), colors)
    else:
        colors = colors[:len(segments)]
    colors = {list(segments.keys())[i]:(r/255.0,g/255.0,b/255.0,0.5) for i,(r,g,b) in enumerate(colors)}
    outname = os.path.join(opt["var_pfx"], f"{out_name}_{name}")
    peak_list, splash_list, color_list = list(), list(), list()
    peak_minima, splash_minima, color_minima = list(), list(), list()
    peak_maxima, splash_maxima, color_maxima = list(), list(), list()
    peak_area, splash_area, color_area = list(), list(), list()
    plot_list, minima_list, maxima_list, area_list = list(), list(), list(), list()
    prev_pos = 0
    for seg,typ in hots_dict.keys():
        hots_min = sorted(hots_minima[(seg,typ)], key=itemgetter(1))
        hots_max = sorted(hots_maxima[(seg,typ)], key=itemgetter(1))
        data_dict = hots_dict.get((seg,typ), [])
        ## profile data
        p_list = [i for i in data_dict.get("peak_sum", [])]
        s_list = [i for i in data_dict.get("splash_sum", [])]
        c_list = [colors[seg] for i in range(len(data_dict["peak_sum"]))]
        peak_list += p_list
        splash_list += s_list
        color_list += c_list
        ## minima data
        p_minima = [p_list[pos] for peak,pos,count in hots_min]
        s_minima = [s_list[pos] for peak,pos,count in hots_min]
        c_minima = [colors[seg] for i in range(len(hots_min))]
        peak_minima   += p_minima
        splash_minima += s_minima
        color_minima  += c_minima
        ## minima data
        p_maxima = [p_list[pos] for peak,pos,count in hots_max]
        s_maxima = [s_list[pos] for peak,pos,count in hots_max]
        c_maxima = [colors[seg] for i in range(len(hots_max))]
        peak_maxima   += p_maxima
        splash_maxima += s_maxima
        color_maxima  += c_maxima
        ## area data
        p_area = [sum([p_list[pos] for pos in range(pos1,pos2+1)]) for (peak1,pos1,count1),(peak2,pos2,count2) in zip(hots_max[:-1],hots_max[1:])]
        s_area = [sum([s_list[pos] for pos in range(pos1,pos2+1)]) for (peak1,pos1,count1),(peak2,pos2,count2) in zip(hots_max[:-1],hots_max[1:])]
        c_area = [colors[seg] for i in range(len(hots_max)-1)]
        peak_area   += p_area
        splash_area += s_area
        color_area  += c_area
        ## append plot lists
        plot_list.append((p_list,s_list,c_list,seg,typ,"profile"))
        minima_list.append((p_minima,s_minima,c_minima,seg,typ,"minima"))
        maxima_list.append((p_maxima,s_maxima,c_maxima,seg,typ,"maxima"))
        area_list.append((p_area,s_area,c_area,seg,typ,"area"))
    plot_list.append((peak_list,splash_list,color_list,"all","-","profile"))
    minima_list.append((peak_minima,splash_minima,color_minima,"all","-","minima"))
    maxima_list.append((peak_maxima,splash_maxima,color_maxima,"all","-","maxima"))
    area_list.append((peak_area,splash_area,color_area,"all","-","area"))
    mpeak, msplash = min(peak_list), max(splash_list)
    apeak, asplash = min(p_area), max(s_area)
    s = opt["var_pls"]
    for l in ["raw", "log10"]:
        with open(f"{outname}_{l}.tsv", "w") as outcorr:
            #outcorr.write(f"segment\ttype\tdata\tpearson\tpearson-p-value\tspearman\tspearman-p-value\tkendalltau\tkendalltau-p-value\n")
            outcorr.write(f"segment\ttype\tdata\tspearman\tp-value\n")
            for p1,p2,p3,p4 in zip(plot_list,minima_list,maxima_list,area_list):
                doScatter(outname, "area", l, s, apeak, asplash, outcorr, [p4], colors)
                #doScatter(outname, "profile", l, s, mpeak, msplash, outcorr, [p1,p3,p2], colors)
    opt["var_pfx"] = var_pfx

def doScatter(outname, stype, l, s, mpeak, msplash, outcorr, pdata, colors):
    ## create plots
    pz = plt.figure(figsize=(4*s,4*s))
    for peak,splash,color,seg,typ,dname in pdata:
        if l == "log10":
            peak = [log10(-i) if i != 0 else 0.0 for i in peak]
            splash = [log10(i) if i != 0 else 0.0 for i in splash]
        elif l == "normalized":
            peak = [i/mpeak for i in peak]
            splash = [i/msplash for i in splash]
        pear = pearsonr(peak,splash)
        spear = spearmanr(peak,splash)
        kendl = kendalltau(peak,splash)
        #outcorr.write(f"{seg}\t{typ}\t{dname}\t{pear[0]}\t{pear[1]}\t{spear[0]}\t{spear[1]}\t{kendl[0]}\t{kendl[1]}\n")
        outcorr.write(f"{seg}\t{typ}\t{dname}\t{spear[0]}\t{spear[1]}\n")
        #print(corrcoef([i/mpeak for i in peak_extrema],[i/msplash for i in splash_extrema]))
        #plt.scatter(peak_list, splash_list, s=3, c=color_list, alpha=1.0)
        if dname == "minima":
            s1 = plt.scatter(peak, splash, s=3, c=color, alpha=1.0, edgecolors="#000000", linewidth=0.5)
        elif dname == "maxima":
            s2 = plt.scatter(peak, splash, s=3, c=color, alpha=1.0, edgecolors="#737373", linewidth=0.5)
        elif dname == "area":
            s3 = plt.scatter(peak, splash, s=3, c=color, alpha=1.0, edgecolors="#000000", linewidth=0.3, marker="X")
        else:
            s4 = plt.scatter(peak, splash, s=3, c=color, alpha=1.0)#, edgecolors=(0,0,0,1.0), linewidth=0.5)
        axes = plt.gca()
        x, y = polyfit(peak, splash, 1)
        cor_line = linspace(axes.get_xlim()[0],axes.get_xlim()[1],100)
        plt.plot(cor_line, x*cor_line + y, linestyle="-", color="black", linewidth=0.8)
        if  opt["var_stt"]: stt = opt["var_stt"]+"_"
        else: stt = ""
        outProfile = f"{outname}_{stt}{l}_{stype}_{seg}-{typ}"
    if l == "log10":
        plt.xlim([-0.2,0.2+ceil(max(log10(-mpeak),log10(msplash),log10(opt["var_lmx"])))])
        plt.ylim([-0.2,0.2+ceil(max(log10(-mpeak),log10(msplash),log10(opt["var_lmx"])))])
    if dname == "area":
        plt.xlabel("peak area sum")#, fontsize=12)
        plt.ylabel("RPM area sum")#, fontsize=12)
    else:
        plt.xlabel("peak sum")#, fontsize=12)
        plt.ylabel("RPM sum")#, fontsize=12)
    if stype == "profile":
        plt.legend((s1, s2, s4), ("minima", "maxima", "profile"), scatterpoints=1, loc='lower left', ncol=1, fontsize=8)#, title=list(colors.keys())[0].split("_")[0])
    else:
        pop = [mpatches.Patch(color=i, label=k) for k,i in colors.items()]
        L = plt.legend(handles=pop, scatterpoints=1, loc='lower right' if l == "log10" else 'lower left', ncol=1, fontsize=8)#, title=list(colors.keys())[0].split("_")[0])
        for ltex in L.get_texts():
            lx = ltex.get_text().split("_")[1] if len(ltex.get_text().split("_")) > 1 else ltex.get_text()
            ltex.set_text(lx)
    #if l == "log10":
    #    plt.xscale("log")
    #    plt.yscale("log")
    if l == "log10":
        tk_list = list(arange(0,ceil(max(log10(-mpeak),log10(msplash),log10(opt["var_lmx"]))),step=1))
        tk_names = [f"$10^{i}$" for i in tk_list]
        plt.xticks(tk_list,labels=tk_names)
        plt.yticks(tk_list,labels=tk_names)
    if  opt["var_stt"]: plt.title(opt["var_stt"].replace("_"," "))#, fontsize=16)
    plt.subplots_adjust(bottom=0.2, left=0.2)
    #outProfile = f"{outname}_{l}_{dname}_{seg}-{typ}"
    pz.savefig(f"{outProfile}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{outProfile}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)

def addColors(num_colors, colors):
    ## random color generator
    seed(1337)
    for i in range(13,num_colors+1):
        colors += (randint(0,255),randint(0,255),randint(0,255))
    return colors

def getDistribution(links, lk_min, lk_max, dtype):
    ## calculate distribution
    if dtype == "peak": pd = {m:0 for m in arange(int(round(lk_max,0)),int(round(lk_min,0))-1,-1)}
    else:               pd = {m:0 for m in arange(int(round(lk_min/2,0)*2),int(round(lk_max/2,0)*2)+1,2)}
    for vk in links:
        if dtype == "peak": 
            if not isnan(vk.peak): pd[int(round(vk.peak,0))] += 1
        else:               pd[int(round(mean([vk.aj-vk.ai,vk.bj-vk.bi])/2,0)*2)] += 1
    return pd

def getExt(**opt):
    ## get naming extension
    ext = ""
    extensions = opt["var_nex"].split(",")
    if extensions:
        for e in extensions:
            if e in ["peak","mfe","length","splash"]:
                ext += f"_{e}{opt[f'var_cd{e[0]}']}"
    return ext

def loadMatrix(**opt):
    ## load vRNAsite matrices
    heat_dict, splash_dict = loadData(**opt)
    ## create hotspot data
    mat_dict = dict()
    for ((aSeq,aType),(bSeq,bType)),peak_mat in heat_dict.items():
        ## remove orientation data
        if len(aType.split("-")) > 1 or len(bType.split("-")) > 1: continue
        splash_mat = splash_dict.get(((aSeq,aType),(bSeq,bType)), zeros(shape=peak_mat.shape))
        mat_dict[((aSeq,aType),(bSeq,bType))] = (peak_mat,splash_mat)
    return mat_dict

def loadData(**opt):
    ## load data with pickle
    with open(opt["var_mat"], "r+b") as mtin:
        heat_dict, splash_dict = pickle.load(mtin)
    return heat_dict, splash_dict

def loadShape(genomes, **opt):
    ## load shape data
    if not os.path.isdir(opt["var_lsh"]):
        print(f"Error: Not a valid directory!")
        sys.exit()
    shape_files = os.listdir(opt["var_lsh"])
    total = len(shape_files)
    shape_dict = dict()
    genomes_list = [g for gen in genomes for g in gen]
    for it,sfile in enumerate(shape_files):
        # initiate shape list
        print(f"Status: Reading shape data {it+1} of {total} ...             ", end="\r")
        aSeq,aType  = os.path.splitext(sfile)[0].split("-")
        if (aSeq,aType) not in genomes_list:
            print(f"Warning: No entry for \">{aSeq} {aType}\"!")
            continue
        bound_thresh = 0.85
        free_thresh = 0.4
        max_thresh = 10.0
        shape_low, shape_medium, shape_high, shape_all = list(), list(), list(), list()
        with open(os.path.join(opt["var_lsh"],sfile), "r") as sin:
            for line in sin:
                i,s = line.strip().split()
                s = float(s)
                if s <= -999.00:
                    shape_low.append(float("nan"))
                    shape_medium.append(float("nan"))
                    shape_high.append(float("nan"))
                    shape_all.append(float("nan"))#append(0.0)
                elif s < -0.5:
                    shape_low.append(0.0)
                    shape_medium.append(float("nan"))
                    shape_high.append(float("nan"))
                    shape_all.append(s)#append(0.0)
                elif s < free_thresh:
                    shape_low.append(s+0.5)
                    shape_medium.append(float("nan"))
                    shape_high.append(float("nan"))
                    shape_all.append(s)#append(s+0.5)
                elif s < bound_thresh:
                    shape_low.append(float("nan"))
                    shape_medium.append(s+0.5)
                    shape_high.append(float("nan"))
                    shape_all.append(s)#append(s+0.5)
                elif s < max_thresh:
                    shape_low.append(float("nan"))
                    shape_medium.append(float("nan"))
                    shape_high.append(s+0.5)
                    shape_all.append(s)#append(s+0.5)
                else:
                    shape_low.append(float("nan"))
                    shape_medium.append(float("nan"))
                    shape_high.append(max_thresh)
                    shape_all.append(s)#append(max_thresh)
        old_flat = shape_all
        s, t = floor(opt["var_rls"]/2),ceil(opt["var_rls"]/2)
        for j in range(opt["var_rli"]):
            new_flat = list()
            for i in range(len(old_flat)):
                nonan = [h for h in old_flat[max(0,i-s):min(i+t,len(old_flat))] if not isnan(h)]
                if nonan: x = median(nonan)
                else: x = float("nan")
                new_flat.append(x)
            old_flat = new_flat
        gnonan = [h for h in shape_all if not isnan(h)]
        if gnonan: global_median = [median(gnonan)]*len(shape_all)
        else:      global_median = [float("nan") ]*len(shape_all)
        shape_dict[(aSeq,aType)] = (shape_low, shape_medium, shape_high, old_flat, global_median)
    return shape_dict

def analyzeLinkHotSpots(tlist, segments, genomes, als, genome_dict, shape_dict, **opt):
    ## get hotspots
    lk_min = min([t.peak for t in tlist if not isnan(t.peak)])
    lk_max = max([t.peak for t in tlist if not isnan(t.peak)])
    # genomes = [(("S1,"WT"),("S2,"WT"),("S3,"WT")),(("S1,"WT"),("S2,"mut"),("S3,"WT")),...]
    for gen in genomes:
        names = "".join([t for s,t in list(gen) if t != "WT"])
        print(f"Status: Do {als} {names} links ...")
        if names: names = "_" + names
        links = [lk for lk in tlist if (lk.aSeq,lk.aType) in gen and (lk.bSeq,lk.bType) in gen]
        hots_dict, hots_minima, hots_maxima = identifyHotSpots(links, gen, segments, shape_dict, **opt)
        hots_maxima = addOuterExtrema(hots_dict, hots_minima, hots_maxima)
        plotBarplot(links, lk_min, lk_max, f"peak_{als}", "average free energy in kcal/mol", f"{names}", **opt)
        writeHotSpots(hots_minima, hots_maxima, hots_dict, f"extrema_links_{als}{names}", "links", **opt)
        if opt["var_csc"]:
            plotScatter(hots_minima, hots_maxima, hots_dict, segments, f"{als}{names}", **opt)
        genome_dict[(gen, "links", als)] = hots_dict
    return genome_dict

def analyzeMatrixHotSpots(genomes, mat_dict, shape_dict, als, genome_dict, **opt):
    # genomes = [(("S1,"WT"),("S2,"WT"),("S3,"WT")),(("S1,"WT"),("S2,"mut"),("S3,"WT")),...]
    for gen in genomes:
        print(f"Status: Do {als} matrices ...")
        names = "".join([t for s,t in list(gen) if t != "WT"])
        if names: names = "_" + names
        mat_plots = {((aSeq,aType),(bSeq,bType)):(peak_mat,splash_mat) for ((aSeq,aType),(bSeq,bType)),(peak_mat,splash_mat) in mat_dict.items() if (aSeq,aType) in gen and (bSeq,bType) in gen}
        p_quant, s_quant = getQuantiles(mat_plots, **opt)
        if not opt["var_pqd"]:
            hots_dict, hots_minima, hots_maxima = matHotSpots(mat_plots, shape_dict, als, p_quant, s_quant, **opt)
            hots_maxima = addOuterExtrema(hots_dict, hots_minima, hots_maxima)
            writeHotSpots(hots_minima, hots_maxima, hots_dict, f"extrema_matrix_{als}{names}", "matrix", **opt)
            genome_dict[(gen, "matrix", als)] = hots_dict
    return genome_dict

def createDifference(i, j, hot_d1, vtt_d1, vtt_b1, hot_d2, vtt_d2, vtt_b2, total1, total2, opt):
    ## multi analysis
    print(f"Status: Analysing {i+1}/{total1} x {j+1}/{total2}                    ", end="\r")
    segments = opt["var_seg"].split(",")
    diff_dict = dict()
    for i,seg in enumerate(hot_d1.keys()):
        nt_l1, nt_l2 = hot_d1[seg][0], hot_d2[seg][0]
        ct_l1, ct_l2 = hot_d1[seg][1], hot_d2[seg][1]
        nt_diff, ct_diff = list(), list()
        for i in range(len(nt_l1)):
            nt_diff.append(nt_l1[i]-nt_l2[i])
            ct_diff.append(ct_l1[i]-ct_l2[i])
        diff_dict[seg] = (nt_diff, ct_diff)
    ## extrema[seg] = extrema_list = (extrema,position,count)
    analysis_minima = getLocalExtrema(diff_dict, "min", **opt)
    analysis_maxima = getLocalExtrema(diff_dict, "max", **opt)
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outfile = os.path.join(opt["var_pfx"], f"{out_name}_hot-diff-{vtt_b1}-{vtt_b2}-differences.tsv")
    with open(outfile, "w") as outf:
        outf.write(f"segment\tstrain\tbinary\tminima\tmaxima\n")
        for i,seg in enumerate(segments):
            analysis_minima[seg]
            analysis_maxima[seg]
            nt_min = ";".join([f"{n},{p+1},{c}" for n,p,c in analysis_minima[seg]])
            nt_max = ";".join([f"{n},{p+1},{c}" for n,p,c in analysis_maxima[seg]])
            outf.write(f"{seg}\t{vtt_d1[seg]}-{vtt_d2[seg]}\t{vtt_b1[i]}-{vtt_b2[i]}\t{nt_min}\t{nt_max}\n")
    return (diff_dict, vtt_d1, vtt_b1, vtt_d2, vtt_b2, analysis_minima, analysis_maxima)

def writeHotSpots(hots_minima, hots_maxima, hots_dict, fname, htype, **opt):
    ## save hot spots to file
    #local_list = hots_local[(seg,typ)]
    #round(ntc,2), int(round((start+1+end)/2,0)), ctl[int(round((start+1+end)/2,0))] = local_list
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outfile = os.path.join(opt["var_pfx"], f"{out_name}_{fname}")
    ## define output list columns
    header = f"sequence\ttype\tposition"
    if htype == "matrix":
        key_list = ["peak_sum",   "peak_count",   "peak_mean",   "peak_m1std",   "peak_m2std",   "peak_m3std",   "peak_min", 
                    "splash_sum", "splash_count", "splash_mean", "splash_m1std", "splash_m2std", "splash_m3std", "splash_max"]
        header += f"\tpeak_sum\tpeak_count\tpeak_mean\tpeak_m1std\tpeak_m2std\tpeak_m3std\tpeak_min"
        header += f"\tsp_sum\tsp_count\tsp_mean\tsp_m1std\tsp_m2std\tsp_m3std\tsp_max"
    else:
        key_list = ["link_11", "link_13", "link_15", "link_17", "peak_sum", "peak_count", "splash_sum"]
        header += f"\tpeak>-12\tpeak>-14\tpeak>-16\tpeak>-18\tpeak<-18\tpeak_sum\tpeak_count\tsp_mean"
    ## add shape-map data
    if opt["var_lsh"]:
        key_list += ["shape_low", "shape_medium", "shape_high", "shape_median", "shape_gmedian"]
        header += "\tshape_low\tshape_medium\tshape_high\tshape_median\tshape_gmedian"
    #header += "\tpeak_area\tsp_area"
    with open(f"{outfile}_minima.tsv", "w") as outf:
        outf.write(f"{header}\n")
        #extrema\t
        for (seg,typ),extrema in hots_minima.items():
            data_dict = hots_dict[(seg,typ)]
            #p_area,s_area = area_dict[(seg,typ)]
            for peak,pos,count in extrema:
                outf.write(f"{seg}\t{typ}\t{pos}")
                #\tminima")
                last_link = 0.0
                for key in key_list:
                    data_list = data_dict.get(key, [])
                    if data_list:
                        if key in ["link_11", "link_13", "link_15", "link_17"] or key == "peak_sum":
                            if last_link == data_list[pos]:
                                outf.write(f"\t0.0")
                            else:
                                outf.write(f"\t{data_list[pos]-last_link:.2f}")
                                last_link = data_list[pos]
                            if key == "peak_sum":
                                outf.write(f"\t{data_list[pos]:.2f}")
                        else:
                            if key == "peak_count":
                                outf.write(f"\t{-data_list[pos]:.2f}")
                            else:
                                outf.write(f"\t{data_list[pos]:.2f}")
                    else:
                        outf.write(f"\tnan")
                #outf.write(f"\t{p_area[pos]}\t{s_area[pos]}")
                outf.write("\n")
    with open(f"{outfile}_maxima.tsv", "w") as outf:
        outf.write(f"{header}\n")
        #extrema\t
        for (seg,typ),extrema in hots_maxima.items():
            data_dict = hots_dict[(seg,typ)]
            for peak,pos,count in extrema:
                outf.write(f"{seg}\t{typ}\t{pos}")
                #\tminima")
                last_link = 0.0
                for key in key_list:
                    data_list = data_dict.get(key, [])
                    if data_list:
                        if key in ["link_11", "link_13", "link_15", "link_17"] or key == "peak_sum":
                            if last_link == data_list[pos]:
                                outf.write(f"\t0.0")
                            else:
                                outf.write(f"\t{data_list[pos]-last_link:.2f}")
                                last_link = data_list[pos]
                            if key == "peak_sum":
                                outf.write(f"\t{data_list[pos]:.2f}")
                        else:
                            if key == "peak_count":
                                outf.write(f"\t{-data_list[pos]:.2f}")
                            else:
                                outf.write(f"\t{data_list[pos]:.2f}")
                    else:
                        outf.write(f"\tnan")
                outf.write("\n")

def getQuantiles(mat_plots, **opt):
    ## calculate quantules
    peak_list = zeros(shape=(0,0)).flatten()
    splash_list = zeros(shape=(0,0)).flatten()
    for peak_mat,splash_mat in mat_plots.values():
        p_flat = peak_mat.flatten()
        s_flat = splash_mat.flatten()
        peak_list = concatenate((peak_list,p_flat))
        splash_list = concatenate((splash_list,s_flat))
    ## quantile distribution
    if opt["var_pqd"]:
        print(f"Status: Plot peak quantile distribution ...")
        plotDotplot(peak_list, f"quantile_peaks", **opt)
        print(f"Status: Plot RPM quantile distribution ...")
        plotDotplot(splash_list, f"quantile_splash", **opt)
    ## filter peaks
    p_quant = quantile(peak_list, opt["var_qtp"])
    pq_list = peak_list[peak_list <= p_quant]
    # filter splash
    s_quant = quantile(splash_list, 1-opt["var_qts"])
    sq_list = splash_list[splash_list >= s_quant]
    plotQuantiles(p_quant, s_quant, peak_list, splash_list, pq_list, sq_list, "quantiles", **opt)
    return p_quant, s_quant

def matHotSpots(mat_plots, shape_dict, als, p_quant, s_quant, **opt):
    ## identify matrix hotspots
    hots_dict = dict()
    total = len(mat_plots.keys())
    for k,(((aSeq,aType),(bSeq,bType)),(peak_mat,splash_mat)) in enumerate(mat_plots.items()):
        if als == "intra":
            if aSeq != bSeq:
                continue
        elif als == "inter":
            if aSeq == bSeq:
                continue
        print(f"Status: Calculate hot-spots for matrix {k+1} of {total}.", end="\r")
        peak_mat_quant = peak_mat.copy()
        splash_mat_quant = splash_mat.copy()
        peak_mat_quant[peak_mat_quant > p_quant] = 0.0
        splash_mat_quant[splash_mat_quant < s_quant] = 0.0
        hots_dict[(aSeq,aType)], hots_dict[(bSeq,bType)] = addMat(hots_dict, aSeq, aType, bSeq, bType, peak_mat, splash_mat, peak_mat_quant, splash_mat_quant, **opt)
    ## shl(shape low -0.5 to 0.4), shm (shape medium 0.4 to 0.85), shh (shape high 0.85 4.0)
    for (aSeq,aType) in hots_dict.keys():
        data_dict = dict()
        peak_sum, peak_count, splash_sum, splash_count, peak_data, splash_data = hots_dict[(aSeq,aType)]
        m = len(peak_sum)
        peak_mean, peak_m1std, peak_m2std, peak_m3std, peak_min = [0.0]*m, [0.0]*m, [0.0]*m, [0.0]*m, [0.0]*m
        splash_mean, splash_m1std, splash_m2std, splash_m3std, splash_max = [0.0]*m, [0.0]*m, [0.0]*m, [0.0]*m, [0.0]*m
        for i in range(m):
            pdata, sdata = peak_data[i], splash_data[i]
            if not pdata:
                pdata = [0.0]
            if not sdata:
                sdata = [0.0]
            pm, ps, pmn, sm, ss, smx = mean(pdata), std(pdata), min(pdata), mean(sdata), std(sdata), max(sdata)
            peak_mean[i], peak_m1std[i], peak_m2std[i], peak_m3std[i], peak_min[i] = pm, pm-ps, pm-2*ps, pm-3*ps, pmn
            splash_mean[i], splash_m1std[i], splash_m2std[i], splash_m3std[i], splash_max[i] = sm, sm+ss, sm+2*ss, sm+3*ss, smx
        data_dict["peak_mean"], data_dict["peak_m1std"], data_dict["peak_m2std"], data_dict["peak_m3std"], data_dict["peak_min"] = peak_mean, peak_m1std, peak_m2std, peak_m3std, peak_min
        data_dict["splash_mean"], data_dict["splash_m1std"], data_dict["splash_m2std"], data_dict["splash_m3std"], data_dict["splash_max"] = splash_mean, splash_m1std, splash_m2std, splash_m3std, splash_max
        data_dict["peak_sum"], data_dict["peak_count"], data_dict["splash_sum"], data_dict["splash_count"] = peak_sum, peak_count, splash_sum, splash_count
        ## add shape data
        data_dict["shape_low"], data_dict["shape_medium"], data_dict["shape_high"], data_dict["shape_median"], data_dict["shape_gmedian"] = shape_dict.get((aSeq,aType), ([0.0]*m,[0.0]*m,[0.0]*m,[0.0]*m,[0.0]*m))
        hots_dict[(aSeq,aType)] = data_dict
    hots_minima = getLocalExtrema(hots_dict, "min", **opt)
    hots_maxima = getLocalExtrema(hots_dict, "max", **opt)
    return hots_dict, hots_minima, hots_maxima

def plotQuantiles(p_quant, s_quant, peak_list, splash_list, pq_list, sq_list, name, **opt):
    ## plot quantile data
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outfile = os.path.join(opt["var_pfx"], f"{out_name}_{name}_stats.txt")
    with open(outfile, "w") as out:
        out.write(f"Status: Peak {opt['var_qtp']} quantile (q) {p_quant}\n")
        out.write(f"Status: Peak min {peak_list.min():.2f} q-min {pq_list.min():.2f}\n")
        out.write(f"Status: Peak max {peak_list.max():.2f} q-max {pq_list.max():.2f}\n")
        out.write(f"Status: Peak mean {peak_list.mean():.2f} q-mean {pq_list.mean():.2f}\n")
        out.write(f"Status: Peak median {median(peak_list):.2f} q-median {median(pq_list):.2f}\n")
        out.write(f"Status: RPM {opt['var_qts']} quantile (q) {s_quant}\n")
        out.write(f"Status: RPM min {splash_list.min():.2f} q-min {sq_list.min():.2f}\n")
        out.write(f"Status: RPM max {splash_list.max():.2f} q-max {sq_list.max():.2f}\n")
        out.write(f"Status: RPM mean {splash_list.mean():.2f} q-mean {sq_list.mean():.2f}\n")
        out.write(f"Status: RPM median {median(splash_list):.2f} q-median {median(sq_list):.2f}\n")

def addMat(hots_dict, aSeq, aType, bSeq, bType, peak_mat, splash_mat, peak_mat_quant, splash_mat_quant, **opt):
    ## add links to hotspot
    m, n = peak_mat.shape
    peak_sum_a, peak_count_a, splash_sum_a, splash_count_a, peak_data_a, splash_data_a = hots_dict.get((aSeq,aType), ([0.0]*m,[0]*m,[0.0]*m,[0]*m,[[] for i in range(m)], [[] for i in range(m)]))
    peak_sum_b, peak_count_b, splash_sum_b, splash_count_b, peak_data_b, splash_data_b = hots_dict.get((bSeq,bType), ([0.0]*n,[0]*n,[0.0]*n,[0]*n,[[] for i in range(n)], [[] for i in range(n)]))
    # peaks
    apSum = list(peak_mat_quant.sum(axis=1))
    bpSum = list(peak_mat_quant.sum(axis=0))
    # splash
    asSum = list(splash_mat_quant.sum(axis=1)) if splash_mat_quant.any() else [0.0]*m
    bsSum = list(splash_mat_quant.sum(axis=0)) if splash_mat_quant.any() else [0.0]*n
    # peak count
    peak_count_mat = zeros(shape=peak_mat_quant.shape)
    peak_count_mat[peak_mat_quant < 0.0] = 1.0
    akSum = list(peak_count_mat.sum(axis=1)) if peak_count_mat.any() else [0.0]*m
    bkSum = list(peak_count_mat.sum(axis=0)) if peak_count_mat.any() else [0.0]*n
    # splash count
    splash_count_mat = zeros(shape=splash_mat_quant.shape)
    splash_count_mat[splash_mat_quant > 0.0] = 1.0
    acSum = list(splash_count_mat.sum(axis=1)) if splash_count_mat.any() else [0.0]*m
    bcSum = list(splash_count_mat.sum(axis=0)) if splash_count_mat.any() else [0.0]*n
    # add
    for i in range(m):
        peak_sum_a[i]     += apSum[i]
        peak_count_a[i]   -= akSum[i]
        splash_sum_a[i]   += asSum[i]
        splash_count_a[i] += acSum[i]
        #print([f"{r:.0f}" for r in list(peak_mat_quant[i,:])])
        #exit()
        peak_data_a[i]    += [r for r in list(peak_mat_quant[i,:])]# if r < 0.0]
        #if aSeq == "SC35M_PB2":
        #    print(",".join([r for r in peak_mat[i,:] if r < 0.0]))
        #    exit()
        splash_data_a[i]  += [r for r in list(splash_mat_quant[i,:])]# if r > 0.0]
    for j in range(n):
        peak_sum_b[j]     += bpSum[j]
        peak_count_b[j]   -= bkSum[j]
        splash_sum_b[j]   += bsSum[j]
        splash_count_b[j] += bcSum[j]
        #print([f"{r:.0f}" for r in list(peak_mat_quant[:,j])])
        #exit()
        #if bSeq == "SC35M_PB2":
        #    print(j, str(len([f"{r:.0f}" for r in list(peak_mat[:,j])])), end=", ")
        peak_data_b[j]    += [r for r in list(peak_mat_quant[:,j])]# if r < 0.0]
        splash_data_b[j]  += [r for r in list(splash_mat_quant[:,j])]# if r > 0.0]
        #if bSeq == "SC35M_PB2":
        #    print(",".join([f"{r}" for r in peak_mat[:,j] if r < 0.0]))
        #    exit()
    return (peak_sum_a, peak_count_a, splash_sum_a, splash_count_a, peak_data_a, splash_data_a), (peak_sum_b, peak_count_b, splash_sum_b, splash_count_b, peak_data_b, splash_data_b)

def identifyHotSpots(links, gen, segments, shape_dict, **opt):
    ## identify hot-spots
    hots_dict = dict()
    ## initialize data
    for g in gen:
        hots_dict[g] = initializeHotSpots(segments[g[0]], g, shape_dict)
    for i,lk in enumerate(links):
        hots_dict[(lk.aSeq,lk.aType)] = addLinks(hots_dict, lk.aSeq, lk.aType, lk.peak, lk.sp_mean, range(lk.ai,lk.aj))
        hots_dict[(lk.bSeq,lk.bType)] = addLinks(hots_dict, lk.bSeq, lk.bType, lk.peak, lk.sp_mean, range(lk.bi,lk.bj))
    hots_minima = getLocalExtrema(hots_dict, "min", **opt)
    hots_maxima = getLocalExtrema(hots_dict, "max", **opt)
    return hots_dict, hots_minima, hots_maxima

def applyFilter(ntl):
    ## apply data filter for smoothing
    gtl = list(savgol_filter(ntl, 51, 6)) # window size 51, polynomial order 6
    return gtl

def initializeHotSpots(l, g, shape_dict):
    ## initialize hot spots
    s0,s1,s2,s3,s4 = shape_dict.get(g, ([0.0]*l,[0.0]*l,[0.0]*l,[0.0]*l,[0.0]*l))
    data_dict = {"peak_count":  [0]*l, "peak_sum"    :[0.0]*l, "link_17"   :[0.0]*l, "link_15"     :[0.0]*l,
                 "link_13"   :[0.0]*l, "link_11"     :[0.0]*l, "splash_sum":[0.0]*l, "splash_count":  [0]*l,
                 "shape_low" :     s0, "shape_medium":     s1, "shape_high":     s2, "shape_median"  :     s3, "shape_gmedian"  :     s4}
    return data_dict

def addLinks(hots_dict, seg, typ, peak, sp_mean, pos):
    ## add links to hotspot
    data_dict = hots_dict[(seg,typ)]
    peak_sum, peak_count, splash_sum, splash_count = data_dict["peak_sum"], data_dict["peak_count"], data_dict["splash_sum"], data_dict["splash_count"]
    link_11 , link_13   , link_15   , link_17      = data_dict["link_11"] , data_dict["link_13"]   , data_dict["link_15"]   , data_dict["link_17"]
    for j in pos:
        peak_sum[j]     += peak    if not isnan(peak)                     else 0.0
        peak_count[j]   -= 1       if not isnan(peak)    and peak < 0.0   else 0
        splash_sum[j]   += sp_mean if not isnan(sp_mean)                  else 0
        splash_count[j] += 1       if not isnan(sp_mean) and sp_mean > 0  else 0
        link_11[j]      += peak    if not isnan(peak)    and peak > -12.0 else 0.0
        link_13[j]      += peak    if not isnan(peak)    and peak > -14.0 else 0.0
        link_15[j]      += peak    if not isnan(peak)    and peak > -16.0 else 0.0
        link_17[j]      += peak    if not isnan(peak)    and peak > -18.0 else 0.0
    #data_dict["link_11"], data_dict["link_13"], data_dict["link_15"], data_dict["link_17"] = link_11, link_13, link_15, link_17
    data_dict["peak_sum"], data_dict["peak_count"], data_dict["splash_sum"], data_dict["splash_count"] = peak_sum, peak_count, splash_sum, splash_count
    data_dict["link_11"] , data_dict["link_13"]   , data_dict["link_15"]   , data_dict["link_17"]      = link_11 , link_13   , link_15   , link_17
    return data_dict

def getLocalExtrema(hots_dict, local, **opt):
    ## calculate local minima
    hots_local = dict()
    for (seg,typ),data_dict in hots_dict.items():
        #(ntl,ktl,ctl,st,shl,shm,shh,shf)
        #peak_sum, peak_count, splash_count, splash_sum, shape_low, shape_medium, sha
        ntl, ctl = data_dict["peak_sum"], data_dict["peak_count"]
        if opt["var_nps"]: gtl = ntl
        else:              gtl = applyFilter(ntl)
        local_list = list()
        direction, current, start, end = "up", 0.0, 0, 0
        for i,nt in enumerate(gtl):
            if local == "max": nt = -nt
            if nt < current:
                direction, current, start, end = "down", nt, i, i
            elif nt == current:
                end = i
            elif nt > current and direction == "down":
                if local == "max": ntc = -ntl[int(round((start+1+end)/2,2))]
                else: ntc = ntl[int(round((start+1+end)/2,2))]
                local_list.append((ntc,int(floor((start+1+end)/2)),ctl[int(round((start+1+end)/2,0))]))
                direction, current, start, end = "up", nt, i, i
            else:
                current, start, end = nt, i, i
        local_list = sorted(local_list, key=itemgetter(0))
        hots_local[(seg,typ)] = local_list
    return hots_local

def doPlots(genome_dict, genomes, segments, **opt):
    ## do all lien plots (data_type = links/matrix; als = inter/intra)
    if opt["var_tra"]:
        als = "intra"
    else:
        als = "inter"
    for gen in genomes:
        names = "".join([t for s,t in list(gen) if t != "WT"])
        if names: names = "_" + names
        data_types = ["links"]
        if opt["var_mat"]: data_types.append("matrix")
        for dtype in data_types:
            #if als == "both":
            #    plotLineplot(genome_dict[(gen, dtype, "inter")], segments, f"hot-spot_both_{dtype}{names}", 0.0, 1.0, genome_dict[(gen, dtype, "intra")], **opt)
            #else:
            plotLineplot(genome_dict[(gen, dtype, als)], segments, gen, f"hot-spot_{als}_{dtype}{names}", 0.0, 1.0, {}, **opt)

def plotLineplot(hots_dict, segments, gen, fname, y_min, y_max, hots_dict2, **opt):
    ## pyplot matrix plot vtt_dict, vtt_bin
    y_min, y_max = -1.0, 1.0
    color_dict = {"peak_sum"     :"#238b45", "peak_count"  :"#00441b", "peak_mean"   :"#deebf7",
                  "peak_m1std"   :"#c6dbef", "peak_m2std"  :"#9ecae1", "peak_m3std"  :"#6baed6",
                  "peak_min"     :"#4292c6", "splash_sum"  :"#cb181d", "splash_count":"#67000d",
                  "splash_mean"  :"#efedf5", "splash_m1std":"#dadaeb", "splash_m2std":"#bcbddc",
                  "splash_m3std" :"#9e9ac8", "splash_max"  :"#807dba", "shape_low"   :"#000000",
                  "shape_medium" :"#ffa500", "shape_high"  :"#ff0000", "shape_median":"#919191",
                  "shape_gmedian":"#666666", "link_11"     :"#2c7bb6", "link_13"     :"#abd9e9",
                  "link_15"      :"#fdae61", "link_17"     :"#d7191c"}#41ab5d
    line_dict  = {"peak_sum"     :1.0, "peak_count"  :0.1, "peak_mean"   :0.5,
                  "peak_m1std"   :0.5, "peak_m2std"  :0.5, "peak_m3std"  :0.5,
                  "peak_min"     :0.5, "splash_sum"  :1.0, "splash_count":0.1,
                  "splash_mean"  :0.5, "splash_m1std":0.5, "splash_m2std":0.5,
                  "splash_m3std" :0.5, "splash_max"  :0.5, "shape_low"   :1.0,
                  "shape_medium" :1.0, "shape_high"  :1.0, "shape_median":0.5,
                  "shape_gmedian":0.5, "link_11"     :1.0, "link_13"     :1.0,
                  "link_15"      :1.0, "link_17"     :1.0}
    style_dict = {"peak_sum"     :'-', "peak_count"  :'-', "peak_mean"   :'-',
                  "peak_m1std"   :'-', "peak_m2std"  :'-', "peak_m3std"  :'-',
                  "peak_min"     :'-', "splash_sum"  :'-', "splash_count":'-',
                  "splash_mean"  :'-', "splash_m1std":'-', "splash_m2std":'-',
                  "splash_m3std" :'-', "splash_max"  :'-', "shape_low"   :'-',
                  "shape_medium" :'-', "shape_high"  :'-', "shape_median":'-',
                  "shape_gmedian":'--', "link_11"     :'-', "link_13"     :'-',
                  "link_15"      :'-', "link_17"     :'-'}
    labl_dict  = {"peak_sum"    :"peak $\Sigma$",           "peak_count"  :"peak #",
                  "peak_mean"   :"peak $\mu$",              "peak_m1std"  :"peak $\mu+\sigma$",
                  "peak_m2std"  :"peak $\mu+2\cdot\sigma$", "peak_m3std"  :"peak $\mu+3\cdot\sigma$",
                  "peak_min"    :"peak min",                "splash_sum"  :"RPM $\Sigma$",
                  "splash_count":"RPM #",                   "splash_mean" :"RPM $\mu$",
                  "splash_m1std":"RPM $\mu+\sigma$",        "splash_m2std":"RPM $\mu+2\cdot\sigma$",
                  "splash_m3std":"RPM $\mu+3\cdot\sigma$",  "splash_max"  :"RPM max",
                  "shape_low"   :f"{opt['var_nsh']} < 0.4",             "shape_medium":f"{opt['var_nsh']} < 0.85",
                  "shape_high"  :f"{opt['var_nsh']} > 0.85",            "shape_median":f"{opt['var_nsh']} median",
                  "shape_gmedian":f"{opt['var_nsh']} global median",    "link_11"     :"peak > -12.0 kcal/mol",
                  "link_13"     :"peak > -14.0 kcal/mol",   "link_15"     :"peak > -16.0 kcal/mol",
                  "link_17"     :"peak > -18.0 kcal/mol"}
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outfile = os.path.join(opt["var_pfx"], f"{out_name}_{fname}")
    s = opt["var_pls"]
    min_y, max_y = floor(y_min/opt["var_yts"])*opt["var_yts"], ceil(y_max/opt["var_yts"])*opt["var_yts"]
    keyx = {s:t for s,t in gen}
    #for i,((seg,typ),(nt,ct,st)) in enumerate(hots_dict.items()):
    div_list = {k:0.0 for k in hots_dict[list(hots_dict.keys())[0]].keys()}
    #nts, kts, cts, sts, shp, fhp = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    shp_list = ["shape_low", "shape_medium", "shape_high", "shape_median", "shape_gmedian"]
    neg_list = ["peak_sum", "peak_count", "peak_mean", "peak_m1std", "peak_m2std", "peak_m3std", "peak_min", "link_11", "link_13", "link_15", "link_17"]
    pos_list = ["splash_sum", "splash_count", "splash_mean", "splash_m1std", "splash_m2std", "splash_m3std", "splash_max"]
    for i,seg in enumerate(segments.keys()):
        typ = keyx[seg]
        data_dict = hots_dict.get((seg,typ),dict())
        for k,d in data_dict.items():
            if not opt["var_nps"] and k in ["link_11", "link_13", "link_15", "link_17", "peak_sum"]:
                d = applyFilter(d)
            if k in shp_list:
                if k == "shape_median":
                    htest = [h for h in d if not isnan(h)]
                    if htest:
                        if max(htest) > div_list[k]:
                            div_list[k] = max(htest)
                            div_list["shape_gmedian"] = max(htest)
                elif k == "shape_gmedian": continue
                elif abs(sum([s for s in d if not isnan(s)])) > 0.0: div_list[k] = 10.0
                else: div_list[k] = 0.0
            elif k in neg_list:
                if min(d) < div_list[k]: div_list[k] = min(d)
                if k == "peak_min":
                    div_list["peak_mean"], div_list["peak_m1std"]  = div_list["peak_min"], div_list["peak_min"]
                    div_list["peak_m2std"], div_list["peak_m3std"] = div_list["peak_min"], div_list["peak_min"]
                if k == "peak_sum":
                    div_list["link_11"], div_list["link_13"] = div_list["peak_sum"], div_list["peak_sum"]
                    div_list["link_15"], div_list["link_17"] = div_list["peak_sum"], div_list["peak_sum"]
            elif k in pos_list:
                if max(d) > div_list[k]: div_list[k] = max(d)
                if k == "splash_max":
                    div_list["splash_mean"], div_list["splash_m1std"] = div_list["splash_max"], div_list["splash_max"]
                    div_list["splash_m2std"], div_list["splash_m3std"] = div_list["splash_max"], div_list["splash_max"]
    ## start plotting
    nseg = len(segments.keys())
    max_len = ceil(max([len(data_dict["peak_sum"]) for (seg,typ),data_dict in hots_dict.items()])/opt["var_xts"])*opt["var_xts"]
    pz = plt.figure(figsize=(0.015*s*max_len,6*s*nseg))
    axs = pz.subplots(nseg, sharex=True, sharey=True)
    splash_test = div_list.get("splash_sum", 0.0)
    shape_test  = div_list.get("shape_low", 0.0)
    if splash_test == 0.0 and shape_test == 0.0 :
        max_y = 0.0
        y_max = 0.0
    plt.setp(axs, xticks=arange(0,max_len+1,opt["var_xts"]), yticks=arange(min_y,max_y+1,opt["var_yts"]))
    for i,seg in enumerate(segments.keys()):
        typ = keyx[seg]
        data_dict = hots_dict.get((seg,typ),dict())
        x_min, x_max = 1, segments[seg]+1
        x = arange(1,x_max,1)
        link_ntn_last = [0.0]*(x_max-1)
        links = False
        if nseg > 1:
            axs[i].plot([0], [y_min], color="#ffffff", linewidth=0, alpha=0.0)
            axs[i].plot([0], [y_max], color="#ffffff", linewidth=0, alpha=0.0)
        else:
            axs.plot([0], [y_min], color="#ffffff", linewidth=0, alpha=0.0)
            axs.plot([0], [y_max], color="#ffffff", linewidth=0, alpha=0.0)
        for k,d in data_dict.items():
            if abs(div_list[k]) == 0.0:
                continue
            ## test plotong with savgol
            if not opt["var_nps"] and k in ["link_11", "link_13", "link_15", "link_17", "peak_sum"]:
                d = applyFilter(d)
            if opt["var_rvp"]:
                d = d[::-1]
            ntn = [float(i)/abs(div_list[k]) if div_list[k] else 0.0 for i in d]
            if nseg > 1:
                axs[i].plot(x, ntn, color=color_dict[k], linewidth=line_dict[k], label=labl_dict[k], linestyle=style_dict[k])
            else:
                axs.plot(x, ntn, color=color_dict[k], linewidth=line_dict[k], label=labl_dict[k], linestyle=style_dict[k])
            #if k[:4] == "link" or (k == "peak_sum" and links):
            #    links = True
            #    axs[i].fill_between(x, ntn, link_ntn_last, color=color_dict[k], alpha=0.5)
            #    link_ntn_last = ntn
            #else:
            #    pass
            #    axs[i].plot(x, ntn, color=color_dict[k], linewidth=line_dict[k], label=labl_dict[k])
        links = False
        # manage dualplot
        #nt, kt, ct, st, shl, shm, shh, shf = 
        #if opt["var_rvp"]:
        #    nt, kt, ct, st, shl, shm, shh, shf = nt[::-1], kt[::-1], ct[::-1], st[::-1], shl[::-1], shm[::-1], shh[::-1], shf[::-1]
        #x_min, x_max = 1, len(nt)+1
        #x = arange(1,x_max,1)
        #axs[i].plot([0], [y_min], color="#ffffff", linewidth=0, alpha=0.0)
        #axs[i].plot([0], [y_max], color="#ffffff", linewidth=0, alpha=0.0)
        #ntn = [float(i)/nts if nts else 0.0 for i in nt]
        #ktn = [float(i)/kts if kts else 0.0 for i in kt]
        #axs[i].plot(x, ntn, color="#3288bd", linewidth=1) # vsite blue
        #axs[i].plot(x, ktn, color="#21258f", linewidth=0.5) # non 0 counts dark blue
        #if not opt["var_lsh"]:
        #    stn = [float(i)/sts if sts else 0.0 for i in st]
        #    ctn = [float(i)/cts if cts else 0.0 for i in ct]
        #    axs[i].plot(x, stn, color="#d53e4f", linewidth=1) # splash red
        #    axs[i].plot(x, ctn, color="#b50d0d", linewidth=0.5) # non 0 counts dark red
        ## add second intra hots_dict peaks and splash if available
        #if hots_dict2:
        #    nt, kt, ct, st, dmy, dmy, dmy, dmy = hots_dict2[(seg,typ)]
        #    if opt["var_rvp"]:
        #        nt, kt, ct, st = nt[::-1], kt[::-1], ct[::-1], st[::-1]
        #    ntn = [float(i)/nts if nts else 0.0 for i in nt]
        #    ktn = [float(i)/kts if kts else 0.0 for i in kt]
        #    axs[i].plot(x, ntn, color="#a6d7f5", linewidth=1) # vsite blue
        #    axs[i].plot(x, ktn, color="#464882", linewidth=0.5) # non 0 counts dark blue
        #    if not opt["var_lsh"]:
        #        stn = [float(i)/sts if sts else 0.0 for i in st]
        #        ctn = [float(i)/cts if cts else 0.0 for i in ct]
        #        axs[i].plot(x, stn, color="#d48790", linewidth=1) # splash red
        #        axs[i].plot(x, ctn, color="#993f3f", linewidth=0.5) # non 0 counts dark red
        #if opt["var_lsh"]:
        #    # set shp to max plot shape of 10.0 (100%) to be consistent with all figures
        #    shp = 10.0
        #    shln = [float(i)/shp for i in shl]
        #    shmn = [float(i)/shp for i in shm]
        #    shhn = [float(i)/shp for i in shh]
        #    fhhn = [float(i)/shp for i in shf]
        #    axs[i].plot(x, fhhn, color="#919191", linewidth=0.5) # shape green rolling
        #    axs[i].plot(x, shln, color="#000000", linewidth=1) # shape black low
        #    axs[i].plot(x, shmn, color="#ffa500", linewidth=1) # shape orange medium
        #    axs[i].plot(x, shhn, color="#ff0000", linewidth=1) # shape red high
        if nseg > 1:
            axs[i].set_xlabel("nucleotide position", fontsize=20)
            axs[i].set_ylabel("normalized interactions", fontsize=20)
        else:
            axs.set_xlabel("nucleotide position", fontsize=20)
            axs.set_ylabel("normalized interactions", fontsize=20)
        if i == len(segments.keys())-1:
            if nseg > 1:
                axs[i].legend()
            else:
                axs.legend()
        if nseg > 1:
            name = " ".join(seg.split("_"))
            axs[i].set_title(f"{name} {typ}", fontsize=24)
            axs[i].tick_params(axis="x",reset=True)
            axs[i].grid()
        else:
            name = " ".join(seg.split("_"))
            axs.set_title(f"{name} {typ}", fontsize=24)
            axs.tick_params(axis="x",reset=True)
            axs.grid()
    plt.subplots_adjust(hspace=0.3)
    pz.savefig(f"{outfile}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{outfile}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)
    ## plot data
    header_list = list()
    for i,((seg,typ),data_dict) in enumerate(hots_dict.items()):
        keys = list(data_dict.keys())
        if not header_list or len(header_list) < len(keys): header_list = keys
    with open(f"{outfile}.tsv", "w") as outf:
        header = "\t".join(header_list)
        outf.write(f"segment\ttype\tposition\t{header}\n")
        for i,((seg,typ),data_dict) in enumerate(hots_dict.items()):
            for j in range(len(data_dict["peak_sum"])):
                outf.write(f"{seg}\t{typ}\t{j+1}")
                for k in header_list:
                    d = data_dict.get(k,[])
                    if d: outf.write(f"\t{d[j]}")
                    else: outf.write(f"\t0")
                outf.write(f"\n")
    ## create scatter plots
    #pz = plt.figure(figsize=(5*s,5*s))



################################################################################
## parser
################################################################################

if __name__ == "__main__":

    ############################################################################
    ## get time and save call
    sscript = sys.argv[0]
    start_time = time.time()
    current_time = time.strftime('%x %X')
    scall = " ".join(sys.argv[1:])
    with open(f"{sscript}.log", "a") as calllog:
        calllog.write(f"Start : {current_time}\n")
        calllog.write(f"Script: {sscript}\n")
        calllog.write(f"Call  : {scall}\n")
    print(f"Call: {scall}")
    print(f"Status: Started at {current_time}")
    ############################################################################
    ## transform string into int, float, bool if possible
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    ############################################################################
    ## save documentation
    rx_text = re.compile(r"\n^(.+?)\n((?:.+\n)+)",re.MULTILINE)
    rx_oneline = re.compile(r"\n+")
    rx_options = re.compile(r"\((.+?)\:(.+?)\)")
    help_dict, type_dict, text_dict, mand_list = {}, {}, {}, []
    for match in rx_text.finditer(script_usage):
        argument = match.groups()[0].strip()
        text = " ".join(rx_oneline.sub("",match.groups()[1].strip()).split())
        argopts = {"action":"store", "help":None, "default":None, "choices":None}
        for option in rx_options.finditer(text):
            key = option.group(1).strip()
            var = option.group(2).strip()
            if var == "False": argopts["action"] = "store_true"
            if var == "True": argopts["action"] = "store_false"
            if key == "choices": var = [vs.strip() for vs in var.split(",")]
            if key == "default": var = trans(var)
            argopts[key] = var
        if argopts["default"]: add_default = f" (default: {str(argopts['default'])})"
        else: add_default = ""
        argopts["help"] = rx_options.sub("",text).strip()+add_default
        argnames = argument.split(",")
        if len(argnames) > 1:
            if argopts["default"] == None:
                mand_list.append(f"var_{argnames[1][1:]}")
            type_dict[f"var_{argnames[1][1:]}"] = argopts["default"]
            argopts["argshort"] = argnames[1]
            help_dict[argnames[0]] = argopts
        else:
            text_dict[argnames[0]] = argopts["help"]
    ############################################################################
    ## get arguments
    if text_dict["dependencies"]:
        desc = f"{text_dict['description']} (dependencies: {text_dict['dependencies']})"
    else:
        desc = text_dict['description']
    p = ap.ArgumentParser(prog=sscript, prefix_chars="-", usage=text_dict["usage"],
                          description=desc, epilog=text_dict["reference"])
    p.add_argument("-v", "--version", action="version", version=text_dict["version"])
    for argname,argopts in help_dict.items():
        argshort = argopts["argshort"]
        if argopts["choices"]:
            p.add_argument(argshort, argname,            dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"],   choices=argopts["choices"])
        else:
            p.add_argument(argopts["argshort"], argname, dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"])
    p._optionals.title = "arguments"
    opt = vars(p.parse_args())
    ############################################################################
    ## validate arguments
    if None in [opt[mand] for mand in mand_list]:
        print("Error: Mandatory arguments missing!")
        print(f"Usage: {text_dict['usage']} use -h or --help for more information.")
        sys.exit()
    for key,var in opt.items():
        if key not in mand_list:
            arg_req, arg_in = type_dict[key], trans(var)
            if type(arg_req) == type(arg_in):
                opt[key] = arg_in
            else:
                print(f"Error: Argument {key} is not of type {type(arg_req)}!")
                sys.exit()
    ############################################################################
    ## call main function
    try:
        #saved = main(**opt)
        saved = main(opt)
    except KeyboardInterrupt:
        print("Error: Interrupted by user!")
        sys.exit()
    except SystemExit:
        print("Error: System exit!")
        sys.exit()
    except Exception:
        print("Error: Script exception!")
        traceback.print_exc(file=sys.stderr)
        sys.exit()
    ############################################################################
    ## finish
    started_time = current_time
    elapsed_time = time.time()-start_time
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    current_time = time.strftime('%x %X')
    if saved:
        with open(f"{sscript}.log", "a") as calllog,\
             open(os.path.join(saved,f"call.log"), "a") as dirlog:
            calllog.write(f"Save  : {os.path.abspath(saved)}\n")
            calllog.write(f"Finish: {current_time} in {elapsed_time}\n")
            ## dirlog
            dirlog.write(f"Start : {started_time}\n")
            dirlog.write(f"Script: {sscript}\n")
            dirlog.write(f"Call  : {scall}\n")
            dirlog.write(f"Save  : {os.path.abspath(saved)}\n")
            dirlog.write(f"Finish: {current_time} in {elapsed_time}\n")
    print(f"Status: Saved at {saved}")
    print(f"Status: Finished at {current_time} in {elapsed_time}")
    sys.exit(0)
