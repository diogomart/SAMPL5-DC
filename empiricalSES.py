#!/usr/bin/env python

# Please define the full path to the MSMS executable!
MSMS = '/home/users/diogom/Software/msms/msms.x86_64Linux2.2.6.1'
TMP = '/tmp/empiricalSES'

import sys
import subprocess
from glob import glob
import os
import openbabel, pybel
import argparse

MYVDW = {
  1:1.20,  6:1.70,  7:1.55,  8:1.52, 9:1.47, 15:1.80,
 16:1.80, 17:1.75, 35:1.85, 53:1.98
}

atype_order = [
    'H','HD','C','A','N','NHA','NH','NA','O1',
    'F','P','S','Cl','Br', 'I']

# element, hbond-acc, polar-hydrogen, aromatic
atype_dict = {
    (1,  False,False,False) : 'H',
    (1,  False, True,False) : 'HD',
    (6,  False,False,False) : 'C',
    (6,  False,False, True) : 'A',
    (7,  False,False,False) : 'N',
    (7,   True,False, True) : 'NHA',
    (7,   True,False,False) : 'NH',
    (7,  False,False, True) : 'NA',
    (8,   True,False,False) : 'O1',
    (8,   True,False, True) : 'O1',
    (9,   True,False,False) : 'F',
    (15, False,False,False) : 'P',
    (16, False,False, True) : 'S',
    (16, False,False,False) : 'S',
    (17, False, False, False) : 'Cl',
    (35, False, False, False) : 'Br',
    (53, False, False, False) : 'I'
}

flags = [
    'ANALYTICAL SURFACE AREA',
    'TRIANGULATION',
    'NUMERICAL VOLUMES AND AREA',
    'MSMS terminated'
]

class SurfComponent():
    def __init__(self, ses_ana, sas_ana, ses_num, vol_num, ses_ana_byatom,
        sas_ana_byatom):
        self.ses_ana = ses_ana
        self.sas_ana = sas_ana
        self.ses_num = ses_num
        self.vol_num = vol_num
        self.ses_ana_byatom = ses_ana_byatom
        self.sas_ana_byatom = sas_ana_byatom

class RunMSMS():
    """ This is a wrapper to run MSMS entirely from python
        mslib should be prefered otherwise """

    def __init__(self, mol_fname, readcharges):

        self.name = os.path.splitext(os.path.basename(mol_fname))[0]

        # define filenames
        os.mkdir(TMP)
        msmsexe = MSMS
        af                  = '%s/mymsms_%s.area' % (TMP, self.name)
        of                  = '%s/mymsms_%s.of'   % (TMP, self.name)
        stdname             = '%s/mymsms_%s.std'  % (TMP, self.name)
        self.xyzrn_fname    = '%s/mymsms_%s.xyzrn'% (TMP, self.name)

        # create xyzrn 
        self._writexyzrn(mol_fname)
        self.atypes, self.chrg = self._get_atomtypes(
            mol_fname, atype_dict, readcharges)

        # run MSMS
        cmd = '%s -if %s -af %s -of %s' % (msmsexe, self.xyzrn_fname, af, of)
        stdout = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)
        self.stdout = stdout.communicate()[0]
        self.af = open(af, 'r').read()

        # parse data
        self._parse()#MSMSout3(self.stdout, self.af)

        # remove tmp files
        if True:
            [os.remove(of_) for of_ in glob('%s*' % of)]
            os.remove(af)
            os.remove(self.xyzrn_fname)
            os.rmdir(TMP)

        # perform analysis
        self.pos_comps = [i for i,c in enumerate(self.msms) if c.vol_num > 0]
        self.neg_comps = [i for i,c in enumerate(self.msms) if c.vol_num < 0]


    def _writexyzrn(self, fname):
        # open temporary file
        tmp = open(self.xyzrn_fname, 'w')
        # read molecule
        extension = os.path.splitext(fname)[1].strip('.')
        obmol = pybel.readfile(extension, fname).next().OBMol
        for atom in openbabel.OBMolAtomIter(obmol):
            # residue info for label
            r = atom.GetResidue()
            name = r.GetAtomID(atom).strip()
            resname = r.GetName().strip()
            resnum = r.GetNum()
            atomlabel = '%s_%s_%d' % (name, resname, resnum)
            # radius
            radius = MYVDW[atom.GetAtomicNum()]
            # coordinates
            xyz =  '%11.6f %11.6f %11.6f ' % (
                atom.GetX(), atom.GetY(), atom.GetZ())
            # put it all together
            stuff = '%8.6f %d %14s' % (radius, 1, atomlabel)
            tmp.write(xyz + stuff + '\n')
        tmp.close()

    def _get_atomtypes(self, fname, atype_dict, charges_from_file=False):
        extension = os.path.splitext(fname)[1].strip('.')
        obmol = pybel.readfile(extension, fname).next().OBMol
        atypes = []
        partial_charges = []
        
        if not charges_from_file:
            obmol.UnsetPartialChargesPerceived()

        for atom in openbabel.OBMolAtomIter(obmol):
            key = (
                atom.GetAtomicNum(),
                atom.IsHbondAcceptor(),
                atom.IsPolarHydrogen(),
                atom.IsAromatic()
            )
            atypes.append(atype_dict[key])
            partial_charges.append(atom.GetPartialCharge())
        return atypes, partial_charges

    def _parse(self):
        (anases,anasas), (numvol, numses) = self._parse_stdout(self.stdout)
        atomlabels, a_ses, a_sas = self._parse_af(self.af)
        n_comp = len(a_ses)
        self.msms = []
        for i in range(n_comp):
            self.msms.append(SurfComponent(
                ses_ana = anases[i],
                sas_ana = anasas[i],
                ses_num = numses[i],
                vol_num = numvol[i],
                ses_ana_byatom = a_ses[i],
                sas_ana_byatom = a_sas[i]
            ))

    def _parse_stdout(self, stdout_text):
        txt = stdout_text.split('\n')

        # define text section
        pointers = {}
        for i, line in enumerate(txt):
            for j, flag in enumerate(flags): # let's see if flags is found
                if line.startswith(flag):
                    pointers[j] = i

        # parse analytical text
        a = pointers[flags.index('ANALYTICAL SURFACE AREA')]
        b = pointers[flags.index('TRIANGULATION')]
        ana = self._parse_analytical(txt[a+2:b])

        #parse numerical text
        a = pointers[flags.index('NUMERICAL VOLUMES AND AREA')]
        b = pointers[flags.index('MSMS terminated')]
        nume = self._parse_numeric(txt[a+2:b-1])

        return (ana, nume)

    def _parse_analytical(self, txt):
        ses, sas = [], []
        for line in txt:
            fields = line.split()
            ses.append(float(fields[5]))
            sas.append(float(fields[6]))
        return (ses, sas)

    def _parse_numeric(self, txt):
        ses, vol = [], []
        for line in txt:
            fields = line.split()
            vol.append(float(fields[2]))
            ses.append(float(fields[3]))
        return (vol, ses)

    def _parse_af(self, af_text):
        lines = af_text.split('\n')
        fields = lines[0].split()
        n_comp = (len(fields)-1)/2
        ses, sas = [[] for _ in range(n_comp)], [[] for _ in range(n_comp)]
        atomlabel = []

        for line in lines[1:-1]:
            fields = line.split()
            atomlabel.append(fields[3])
            for i in range(n_comp): 
                ses[i].append(float(fields[1+3*i]))
                sas[i].append(float(fields[2+3*i]))
        return (atomlabel, ses, sas)

    def get_ses_by_atype(self):
        byatype = {}
        for key in atype_order:
            byatype[key] = []
        chrg = 0.0
        for surfcomp in self.msms:
            [byatype[key].append(0.0) for key in byatype]
            for i,area in enumerate(surfcomp.ses_ana_byatom):
                byatype[self.atypes[i]][-1] += area
                chrg += area*abs(self.chrg[i])
        return byatype, chrg

    def get_areas(self):
        areas = []
        byatype_dict, chrg = self.get_ses_by_atype()
        for atype in atype_order:
            areas.append(sum(byatype_dict[atype]))
        areas.append(chrg)
        return areas

# atomic solvation parameters for water
w_dict = {
    "H"     :  0.011220,
    "HD"    : -0.193517,
    "A"     : -0.012638,
    "NHA"    : -0.185590,
    "NH"    : -0.128712,
    "NA"    : -0.626621,
    "O1"    : -0.042211,
    "F"     :  0.072437,
    "Cl"    :  0.013931,
    "charge": -0.246878,
    "C":0, "N":0, "O2":0, "P":0, "S":0, "Br":0, "I":0
}
w_vec = [w_dict[atype] for atype in atype_order + ['charge']]

# atomic solvation parameter for cyclohexane
o_vec = [-0.036 for atype in atype_order ]
o_vec += [0.0] # charge


def get_args():
    parser = argparse.ArgumentParser(description="""
        Estimates:
            (1) free energy of hydration
            (2) free energy of solvation in cyclohexane
            (3) cyclohexane / water partition coefficient
        Typical usage:
            python /path/to/empiricalSES.py *.mol2
        """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('molecules', nargs='+',
        help='any number of molecules')
    parser.add_argument('--retrospective3',
        help=' Use parameters from retrospective experiment 3 in the paper,\n\
the default is to use parameters from submission 49 to SAMPL5\n',
            default=False,
            action='store_true')
    parser.add_argument('--readcharges',
            help=' Read pre-calculated charges from file, the default\n\
is to compute Gasteiger-Marsili charges with openbabel;\n\
this would require recalibrating the parameters\n',
            default=False,
            action='store_true')
    args = parser.parse_args()
    for fname in args.molecules:
        if not os.path.exists(fname):
            sys.stderr.write('%s does not exist. Aborting.\n' % fname)
            sys.exit(2)
    return args

def main():

    # check temporary folder
    if os.path.exists(TMP):
        sys.stderr.write('Unable to create temporary folder %s\n' % TMP)
        sys.stderr.write('Most likely it has been created by an unsuccessful\
run of empiricalSES.py.\n\
You can delete it if there is nothing important in it,\n\
or set a different temporary folder at line 5 of this script\n')
        sys.exit(2)

    # temperature and other stuff
    T = 293
    R = 8.3144598 / 1000.0
    R /= 4.184
    RT = R*T
    K = -2.303 * RT
    
    args = get_args() 
    
    # print header
    largest_fname = max([len(fname) for fname in args.molecules])
    largest_fname = max(largest_fname, 14)
    header = '%' + '%d' % largest_fname + 's, dGsol_wat, dGsol_chex, chex/wat logD'
    print header % 'Molecule name'
    line = '%-' + '%d' % largest_fname + 's,%10.3f,%11.3f,%14.3f'
    
    for mol_fname in args.molecules:
        x = RunMSMS(mol_fname, args.readcharges)
        areas = x.get_areas()
    
        wat_dGsol = np.dot(w_vec, areas)
        chex_dGsol = np.dot(o_vec, areas)
    
        logD = (chex_dGsol - wat_dGsol)/K
    
        print line % (mol_fname, wat_dGsol, chex_dGsol, logD)

if __name__ == '__main__':
    main()

