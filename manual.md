# flexibleRegions - manual

## Usage
```
flexibleRegions.py -fsa <in_fasta> -pfx <out_prefix> [options]
```

## Version
```
flexibleRegions.py 0.0.1 (alpha)
```

## Dependencies
```python v3.7.1```, ```ViennaRNA v2.4.13```

## Description
```flexibleRegions``` reduces the interaction strength between a given area of the first sequence in regard to all other sequences. The tool is written in ```Python 3.7.1``` and the calculations are performed with the ```RNAcofold``` python site-package of the ```ViennaRNA Package 2.4.13```.

## Options
```
--prefix,-pfx
    output directory and prefix for result files

--fasta,-fsa
    fasta with all sequences; the first sequence has to be the short sequence 
    snipped which has to be deoptimized with all following sequences

--positionStart,-pss
    define the starting position of the sequence part which should be 
    deoptimized (default: 0)

--positionEnd,-pse
    define the ending position of the sequence part which should be deoptimized;
    a 0 will be interpreted as the end of the sequence (default: 0)

--seed,-sed
    set the seed for random (default: 0)

--reverse,-rev
    creates reverse of each strain if set (default: False)

--complement,-cmp
    creates complements of each strain if set (default: False)

--threads,-thr
    number of threads to use for RNAcofold (default: 1)

--slice,-slc
    window size for standard window search or minimum length for a seed to be
    accepted if the seed option is specified (default: 20)

--candidateIncrease,-cdi
    increases the energy of the top candidate if there is no possible mutation
    candidate (default: 1.0)

--candidateStop,-cds
    stopping energy for deoptimization (default: -10.0)

--candidateIteration,-cdt
    tries to increase the cdi value at high iteration numbers (default: 1000)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--dangles,-dng
    use dangling ends for foldings (default: 2) (choices: 0,1,2,3)

--noLP,-nlp
    disable lonely pairs for RNAcofold (default: False)
```
