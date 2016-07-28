# SAMPL5-DC

https://www.drugdesigndata.org/about/sampl5

The script empiricalSES.py predicts log D as in submission #49. It estimates:

1. free energy of hydration
2. free energy of solvation in cyclohexane
3. cyclohexane/water partition coefficient

### Requirements / dependencies
* Works on Linux
* Python 
* openbabel         (v2.3.2 available through apt-get)
* python-openbabel  (available through apt-get)
* python-numpy      (available through apt-get)
* MSMS              (http://mgltools.scripps.edu/downloads)

### Setup
1. Download the script (empiricalSES.py)
2. Install the dependencies
3. Edit line 4 in empiricalSES.py to specify the path to MSMS binary

### Usage examples
```bash
python /path/to/empiricalSES.py --help
python /path/to/empiricalSES.py molecule.sdf
python /path/to/empiricalSES.py mymols/*.mol2
```

### Output example
```bash
               Molecule name, dGsol_wat, dGsol_chex, chex/wat logD
SAMPL5_074_intra_H-bond.pdb ,   -25.357,     -8.482,       -12.584
SAMPL5_074_original_conf.pdb,   -27.567,     -8.613,       -14.135
```
