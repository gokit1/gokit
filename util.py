import os
from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser

"""
    GO-Kit: A toolkit for generating flavours of Go models. 
    Copyright (C) <2018>  <Sridhar Neelamraju>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""



class Utils(object):
    def file_exists(self, file_name):
        if os.path.exists(file_name):
            # print ('Found file:',file_name)
            return True
        else:
            print('file not found: ', file_name)
            return False

    def download_pdb_from_id(self, pdbid):
        pdb1 = PDBList()
        pdb1.retrieve_pdb_file(pdbid, pdir='PDB')

    def amino_acid_radius_dict(self):
        # Radius of side-chains for two-bead model
        # http://www.pnas.org/cgi/doi/10.1073/pnas.1019500108. Table S2, SI.
        d = {'CYS': 2.74, 'ASP': 2.79, 'SER': 2.59, 'GLN': 3.01, 'LYS': 3.18,
             'ILE': 3.09, 'PRO': 2.78, 'THR': 2.81, 'PHE': 3.18, 'ASN': 2.84,
             'GLY': 2.25, 'HIS': 3.04, 'LEU': 3.09, 'ARG': 3.28, 'TRP': 3.39,
             'ALA': 2.52, 'VAL': 2.93, 'GLU': 2.96, 'TYR': 3.23, 'MET': 3.09}
        return d

    def get_dict_from_csv(self, csvfile):
        # read parameters as table and save as list
        import csv
        with open(csvfile, 'rb') as f:
            # print ('Reading',csvfile)
            reader = csv.reader(f, delimiter=' ')
            with open('coors_new.csv', mode='w') as outfile:
                mydict = {rows[0]: rows[1] for rows in reader}
        f.close()
        return mydict
    def get_three_letters(self,string):
        #return three letter string for residues from one letter string.
        #e.g for tleap input.
        d=self.amino_acid_1_3()
        #string=list(string)
        three=list((d[i] for i in string))
        return three
    def write_tleap_input(self,sequence,ff):
        #generate a tleap input file.
        print ('source /home/sridhar/wales/AMBERTOOLS/dat/leap/cmd/leaprc.ff03.r1')
        seq=self.get_three_letters(sequence)
        seq=" ".join(str(x) for x in seq)
        print ('mol = sequence {',seq,'}')
        print ('check mol')
        print ('saveamberparm mol coords.prmtop coords.inpcrd ')
        print ('savepdb mol coords.pdb')
        if self.file_exists('coords.prmtop'):
            return True
        else:
            print ('Fatal: Check libraries,path,etc.')
    @property
    def amino_acid_3_1(self):
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        return d
    @property
    def amino_acid_1_3(self):
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        return dict([[v, k] for k, v in d.items()])

    def add_H(self, pdbfile):
        # add hydrogens by running reduce from terminal or see python module!.
        return False

    def remove_duplicates(self,a):
        seen = set()
        for i in list(a):
            if i not in seen:
                seen.add(i)
        return seen

    def print_preamble(self):
        import textwrap
        print('-----------------------------------------')
        gpl=""" GO-Kit: A toolkit for generating flavours of Go models. 
            Copyright (C) <2018>  <Sridhar Neelamraju>
            This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
               
            This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
            along with this program.  If not, see <http://www.gnu.org/licenses/>.
        
            Author: Sridhar Neelamraju
        
                """
        print (textwrap.fill(gpl))
        print ('-----------------------------------------\n\n\n\n\n\n')
        print ('---------------GO-KIT--------------------')

    def print_conclusion(self):
        import textwrap
        print('-----------------------------------------------------------------------------------------------')
        gpl = """Don't forget to generate the table files with tables.py. So Long, and Thanks for All the Fish."""
        print(textwrap.fill(gpl))

    def split_chain(self,pdbfile):
        #split domain-swapped dimer into individual chains.
        chains = ['A', 'B']
        p = PDBParser(PERMISSIVE=1)
        structure = p.get_structure('test', pdbfile)
        for chain in chains:
            pdb_chain_file = 'chain_{}.pdb'.format(chain)
            io_w_no_h = PDBIO()
            io_w_no_h.set_structure(structure)
            io_w_no_h.save('{}'.format(pdb_chain_file), ChainSelect(chain))

    def check_file_extension(self,filename,extension):
        return True if filename.endswith(extension) else False
    def make_dir(self,dirname):
        print ('>>in make_dir:',dirname)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        return
    def make_dir_struc(self,pathdir,mddir):
        import shutil
        mdfiles=['gromacs.gro','gromacs.top'];pathfiles=['SBM.INP','HYDROPHOBIC','IDENTICAL','odata']
        for i in mdfiles:
            try:
                shutil.copy2(i,mddir)
            except:
                pass
        for j in pathfiles:
            try:
                shutil.copy2(j,pathdir)
            except:
                pass

    def dsb_check(self,dsbflag):
        if dsbflag:
            os.remove(PATH/SBM.INP)
            print ("No dsb+hphobic implementation in PATHSAMPLE.")

        # shutil.move(os.path.join('./','SBM.INP'), os.path.join(pathdir, 'SBM.INP'))
        # shutil.move(os.path.join('./','gromacs.gro'), os.path.join(mddir, 'gromacs.gro'))
        # shutil.move(os.path.join('./','gromacs.top'), os.path.join(mddir, 'gromacs.top'))
        # self.make_dir(pathfiles)
        # self.make_dir(mdfiles)
        # os.rename("./SBM.INP", pathfiles+'/SBM.INP');shutil.move("./SBM.INP", pathfiles+'/SBM.INP')
        # os.rename("./gromacs.gro", mdfiles+'/gromacs.gro');shutil.move("./gromacs.gro", mdfiles+'/gromacs.gro')



        return


    def get_xyz_from_minfile(self,natoms,extractmin):
        f1 = open('pathdata', "w+")
        f1.write('%s %d\n' % ('NATOMS', natoms))
        f1.write('%s %d\n' % ('EXTRACTMIN',extractmin))
        f1.close()
        return True


    def fatal_errors(self,errorid):
        d = {1: "Fatal: Check file extension",
             2: "Fatal:Columns in file not equal",
             3: "Fatal: Provide list of pairs through --ext_pair option.",
             4: "Fatal: Need --traj ***.xtc or **.xyz argument",
             5: "Fatal: pairfile must be specified with --ext_pair argument.",
             6: "Fatal: Give native pdb file with --native option.eg.--xyz2pdb 1.pdb --native native_ca.pdb",
             7: "Fatal: Atomtype not implemented. Try --attype 1 or 2",
             8: "Fatal: Dihedral function type must be 1",
             9: "Fatal: Found Glycine with 0 hydrogen atoms. Add hydrogens, or set skip_glycine==True",
             10: "Fatal: Found HETATM in pdbfile.",
             11: "Fatal: --attype must be 1 or 2."}
        print (d[errorid])
        exit()
        return True



class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain
    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:
            return 0
def main():
    X=Utils()
    #X.make_dir_struc('PATH','MD')
    #X.make_dir('peaceout')
    #X.print_preamble()
    #X.make_copy()
    #seq='GSGQVRTIWVGGTPEELKKLKEEAKKANIRVTFWGD'
    #X.write_tleap_input(seq,'ff03')
    G98='TTYKLILNLKQAKEEAIKELVDAGTAEKYFKLIANAKTVEGVWTLKDEIKTFTVTE' #3ALPHA
    G98='TTYKLILNLKQAKEEAIKELVDAGTAEKYFKLIANAKTVEGVWTYKDEIKTFTVTE'
    #X.make_copy('1k50.pdb')
    #print (len(G98))
    #print list(X.get_three_letters(G98))



main()