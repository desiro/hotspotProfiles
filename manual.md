# <samp>hotspotProfiles</samp> - manual

***

## Usage
```
hotspotProfiles.py -fsa <in_fasta> -pfx <out_prefix> [options]
```

## Version
```
hotspotProfiles.py 0.0.1 (alpha)
```

## Dependencies
```Python v3.9.7```, ```Scipy v1.8.0```, ```NumPy v1.22.2```, ```Matplotlib v3.5.1```

## Description
<samp>hotspotProfiles</samp> evaluates vRNAsite like interaction tables and creates hotspot profiles for each individual sequence, showing regions with high interaction counts. Example call: python hotspotProfiles.py -pfx example -vst example.tsv -ovr

## Options
```
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
```

## References
```
D. Desirò, A. Borodavka and M. Marz.
"HotspotProfiles: Prediction of hotspot regions capable of mediating the selective but flexible genome packaging of influenza A viruss."
In Preparation, 2022.
https://github.com/desiro/hotspotProfiles
```