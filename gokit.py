#!/usr/bin/env python

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
PyeSBM:
usage: python esbm.py --options

Author: Sridhar Neelamraju)

"""


from __future__ import print_function
import numpy as np
from conmaps import conmaps
from Bio.PDB.PDBParser import PDBParser
import collections
import mdtraj as md
import os
from util import Utils
from table import tables

#Some conventions:
#First atom is always a C-alpha. This helps find non-bonded atoms easily.

class esbm(object):
    #Write GROMACS topology and SBM.INP files for two-bead models.
    ##!Residue numbering starts from 0 within the python code. Always add 1 to be compatible with GROMACS, OPTIM,etc.!

    #Make separate classes for each potential viz. SOPSC and Cheung-Thirumalai.
    #! CONTACTTYPES for OPTIM
    #1 is 6 - 12
    #2 is 10 - 12
    #5 as gaussian, no ex - vol
    #6 as gaussian, w ex-vol
    #7 is dual gaussian
    #Currently implemented potentials are: 1) Cheung-Thirumalai and 2) SOP-SC with BT params. See parser.
    # For help type python eSBM.py --help


        def __init__(self):
            return

        def globals(self,Ka1,Kb1,Kd1,CA_rad1,skip_glycine1,sopc1,dswap1,btparams1,CAcom1,hphobic1,hpstrength1,hpdist1,dsb1,mjmap1,btmap1,CBfar1,casep1,cbsep1,cabsep1,scaling1,gauss1):
            #Give default values to global variables.
            global Ka,Kb,Kd,CA_rad,hpdist,hpstrength,casep,cbsep,cabsep,scaling
            Ka=Ka1;Kb=Kb1;Kd=Kd1;CA_rad=CA_rad1;hpdist=float(hpdist1);hpstrength=float(hpstrength1);casep=casep1;cbsep=cbsep1;cabsep=cabsep1
            global skip_glycine,sopc,dswap,btparams,CAcom,hphobic,dsb,btmap,mjmap,CBfar,gauss
            skip_glycine=skip_glycine1;dswap=dswap1;btparams=btparams1;sopc=sopc1;CAcom=CAcom1;hphobic=hphobic1;dsb=dsb1;CBfar=CBfar1;scaling=scaling1
            mjmap=mjmap1;btmap=btmap1;gauss=gauss1
            #print ('Global logicals are as follows:\n--------------------------')
            #print (' skip_glycine=',skip_glycine,'\n','dswap=',dswap,'\n','bt=',btparams,'\n','sopc=',sopc,'\n','CAcom=',CAcom,'\n','hphobic=',hphobic)
            global radtodeg,kcaltokj
            #write to SBM.INP
            global w_sbm,w_gro
            w_sbm=True;w_gro=True
            radtodeg = 57.295779513;kcaltokj=4.184
            #print (dsb)
            return
        def two_lists_to_dict(self,keys,values):
            return dict(zip(keys, values))
        def get_dict_from_csv(self,csvfile):
            #read parameters as table and save as list
            import csv
            with open(csvfile, 'rb') as f:
                print ('Reading',csvfile)
                reader=csv.reader(f, delimiter=' ')
                with open('coors_new.csv', mode='w') as outfile:
                    mydict = {rows[0]: rows[1] for rows in reader}
            return mydict
        def amino_acid_dict(self):
            d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
            inv_d = {v: k for k, v in d.iteritems()}
            return inv_d
        def amino_acid_dict2(self):
            d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
            return d
        def amino_acid_radius_dict(self):
            #Radius of side-chains for two-bead model
            #Change here for custom radii.
            #http://www.pnas.org/cgi/doi/10.1073/pnas.1019500108. Table S2, SI.
             d = {'CYS': 2.74, 'ASP': 2.79, 'SER': 2.59, 'GLN': 3.01, 'LYS': 3.18,
            'ILE': 3.09, 'PRO': 2.78, 'THR': 2.81, 'PHE': 3.18, 'ASN': 2.84,
            'GLY': 2.25, 'HIS': 3.04, 'LEU': 3.09, 'ARG': 3.28, 'TRP': 3.39,
            'ALA': 2.52, 'VAL': 2.93, 'GLU': 2.96, 'TYR': 3.23, 'MET': 3.09}
             return d
        def read_external_radius(self,radiusfile):
            print ('>>reading radii from file:',radiusfile)
            r= self.get_dict_from_csv(radiusfile)
            #print (r)
            return r
            return r
        def potential_dict(self,potential):
            #three potentials only implemented for now.
            d = {'lj612': 1, 'lj1012': 2,'dsb': 8}
            return
        def get_CA(self,pdbfile):
            #Returns CA atomcoords and residue number as numpy array
            p = PDBParser(PERMISSIVE=0)
            structure = p.get_structure('test', pdbfile)
            CA = ['CA']
            coord=[]
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            if atom.name in CA:
                                coord.append(atom.get_coord())
            coord=np.array(coord)
            return coord
        def get_side_chain_CB(self,pdbfile):
            Y=conmaps()
            #print (len(Y.get_CB_coordinates(pdbfile)))
            CB_coord=np.asarray(Y.get_CB_coordinates(pdbfile))
            index=Y.get_CB_index(pdbfile)
            topology=md.load(pdbfile).topology
            CB1=[]
            count=0
            for i in index:
                CB1+=[[CB_coord[count]],Y.get_residue_number(topology,i)]
                count=count+1
            CB1 = np.array(CB1).reshape(len(CB1) / 2, 2)
            print (CB1)
            return CB1


        def get_back_bone_COM(self,pdbfile):
            print (">> in get_back_bone_COM",pdbfile)
            p = PDBParser(PERMISSIVE=0)
            structure = p.get_structure('test', pdbfile)
            backbone = ['N', 'CA', 'C', 'O']
            COM1 = []
            for model in structure:
                for chain in model:
                    count=0
                    for residue in chain:
                        coord = []
                        mass = []
                        for atom in residue:
                            if atom.name in backbone:
                                #get coordinates and atomic mass.
                                coord.append(atom.get_coord())
                                mass.append(atom._assign_atom_mass())
            #                    print atom._assign_atom_mass()
                                assert len(mass)!=0, "Zero length mass array! CB atom missing"
                                #print atom.get_parent()
                            # get COM for every
                        coord = np.array(coord)
                        mass = np.array(mass)
                        #res = residue.get_resname().strip()
                        assert len(coord) == len(mass)
                        mass_sum = sum(mass)
                        # print mass
                        COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]), coord)], [count]]
             #           print COM
                        COM1 += COM
                        count = count + 1
                COM1 = np.array(COM1).reshape(len(COM1) / 2, 2)
            return COM1

        def get_side_chain_COM(self,pdbfile,skip_glycine):
            print (">> in get_side_chain_COM:",pdbfile,sopc)
            #returns COM of each side-chain as numppy array
            #SOPC puts CB on glycine hydrogen!!!!
            #GLYCINES are skipped.
            #assert (sopc!=skip_glycine)
            #print(sopc,skip_glycine)
            if not skip_glycine:
                glyname= 'GLY111'
                print (">>Adding CB to Glycines! Make sure Hydrogen atoms exist in pdbfile", pdbfile)
            else:
                glyname='GLY'
            # res numbering starts from 1
            p = PDBParser(PERMISSIVE=0)
            structure = p.get_structure('test', pdbfile)
            backbone=['N','CA','C','O']
            COM1=[]
            count = 0
            for model in structure:
                for chain in model:
                    print (chain)
                    for residue in chain:
                        coord = []
                        mass = []
                        for atom in residue:
                            #get side chain only
                            if atom.name not in backbone:
                                #get coordinates and atomic mass.
                                coord.append(atom.get_coord())
                                #print len(coord)
                                mass.append(atom._assign_atom_mass())
                                assert len(mass)!=0, "Zero length mass array!"
                                #print atom.get_parent()
                        #get COM for every
                        coord=np.array(coord)
                        mass=np.array(mass)
                        res=residue.get_resname().strip()
                        #Make sure GLYCINE is excluded.
                        if res==glyname and skip_glycine:
                            print ("Excluding glycine, res = ",count+1)
                            #print atom_mass
                        else:
                            assert len(coord)==len(mass)
                            mass_sum=sum(mass)
                            #print mass
                            COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]),coord)],[count]]
                            #print COM
                            COM1+=COM
                        count = count + 1
                COM1=np.array(COM1).reshape(len(COM1)/2,2)
                #coordinates,resnum
                # print (COM1)
                # print (len(COM1))
            return COM1

        def write_CB_to_native(self,pdbfile,sopc):
            #Get CB COM         #print ("skip_glycine,sopc",skip_glycine,sopc)

            #print ("Writing native files in native_ca.pdb and native_cb.pdb.\n")
            if not skip_glycine:
                 Y=conmaps()
                 Y.check_glycineh(pdbfile, skip_glycine)
            else:
                sopc=False
             # if not CAcom:
            #     print ("C-alpha at CA atom position. Not COM of backbone.", CAcom)
            Y=conmaps()
            if CBfar:
                CB=self.get_farthest_CB_coords(pdbfile)
            else:
                CB=self.get_side_chain_COM(pdbfile,skip_glycine)
            CA=self.get_back_bone_COM(pdbfile)
            if CAcom and not skip_glycine:
                #print len(CA),len(CB)
                assert len(CB)==len(CA)
            #print CB[:,1]
            #if no center of mass.
            #CB=self.get_side_chain_CB(pdbfile)
            CB_coords = CB[:, 0].tolist()
            CA_coords = CA[:, 0].tolist()
            res_num = CB[:, 1].tolist()
            #print (CB[:,1])
            #print (res_num)
            res_num = [x + 1 for x in res_num]
            #print (res_num)
            CB_dict = self.two_lists_to_dict(res_num, CB_coords)
            #print (CB_dict[88])
            #CA_dict = self.two_lists_to_dict(res_num,CA_coords)
            #print CA_dict
            #print res_num
            f1 = open('native_cb.pdb', "w+"); f2 = open('native_ca.pdb', "w+")
            #Put each CB under the corresponding CA residue
            #Read CA from pdbfile
            print (">>in write_CB_to_native:",pdbfile,sopc)
            with open(pdbfile) as f:
                line=f.readlines()
            line=[x.strip() for x in line]
            count_CA=0
            count_CB=0
            count=0
            for i in line:
                if i[12:16].strip()=="CA":
                    count=count+1
                    k= i.split()
                    k[0] = k[0].ljust(6)  # atom#6s
                    # j[1] = j[1].rjust(5)#aomnum#5d
                    k[1] = str(count).rjust(5)
                    k[2] = 'CA'.center(4)  # atomname$#4s
                    k[3] = k[3].ljust(3)  # resname#1s
                 #   print k[3]
                    k[4] = k[4].rjust(1)  # Astring
                    k[5] = k[5].rjust(4)  # resnum
                    if CAcom:
                        k[6] = str('%8.3f' % (float(CA_coords[count_CA][0]))).rjust(8)  # x
                        k[7] = str('%8.3f' % (float(CA_coords[count_CA][1]))).rjust(8)  # y
                        k[8] = str('%8.3f' % (float(CA_coords[count_CA][2]))).rjust(8)  # z
                    else:
                        k[6] = str('%8.3f' % (float(k[6]))).rjust(8)  # x
                        k[7] = str('%8.3f' % (float(k[7]))).rjust(8)  # y
                        k[8] = str('%8.3f' % (float(k[8]))).rjust(8)  # z
                    k[9] = str('%6.2f' % (float(k[9]))).rjust(6)  # occ
                    k[10] = str('%6.2f' % (float(k[10]))).ljust(6)  # temp
                    k[11] = k[11].rjust(12)  # elname
#                    print("%s%s %s %s %s%s    %s%s%s%s%s%s" % (
#                    k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10], k[11]))
                    f1.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n" % (
                        k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10], k[11]))
                    f2.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n" % (
                        k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10], k[11]))
                    #groline = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"
                    #pdbline = "ATOM  %5i %-3s %3s%2s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                    #final_pdb.append(i)
                    j=i.split()
                    count_CA = count_CA + 1
                    #if not GLYCINE
                    #print (CB_dict.keys())
                    if int(j[5]) in CB_dict.keys():
                        #print (int(j[5]))
                        count=count+1
                        j[0] = j[0].ljust(6)#atom#6s
                        #j[1] = j[1].rjust(5)#aomnum#5d
                        j[1] = str(count).rjust(5)
                        j[2] = 'CB'.center(4)#atomname$#4s
                        j[3] = j[3].ljust(3)#resname#1s
                        j[4] = j[4].rjust(1) #Astring
                        j[5] = j[5].rjust(4) #resnum
                     #  print (str('%8.3f' % (float(CB_coords[count_CB][0]))).rjust(8))
                        j[6] = str('%8.3f' % (float(CB_coords[count_CB][0]))).rjust(8) #x
                        j[7] = str('%8.3f' % (float(CB_coords[count_CB][1]))).rjust(8)#y
                        j[8] = str('%8.3f' % (float(CB_coords[count_CB][2]))).rjust(8) #z\
                        j[9] =str('%6.2f'%(float(j[9]))).rjust(6)#occ
                        j[10]=str('%6.2f'%(float(j[10]))).ljust(6)#temp
                        j[11]=j[11].rjust(12)#elname
                        count_CB=count_CB+1
                        f1.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n"%(j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10],j[11]))
                        #print("%s%s %s %s %s%s    %s%s%s%s%s%s"%(j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10],j[11]))
                        #print j[0]+"    "+j[1]+"  "+j[2]+"  "+j[3]+" "+j[4]+"  "+j[5]+"       "+j[6]+"  "+j[7]
                    else:
                        continue
            print ('See file: native_cb.pdb with ',count_CB+count_CA,'atoms')
            print ('CB=',count_CB,'CA=', count_CA)
            f.close()
            f1.close()
            f2.close()
            return 1

        def get_farthest_atom_in_residue(self,pdbfile,index):
            #returns coordinates of atom farthest from CA in residue.
            Y=conmaps();traj=md.load_pdb(pdbfile)
            topology=traj.topology; nres=int(topology.n_residues)
            heavy_atoms=list(topology.select_atom_indices(selection='all'))
            sc_atoms=Y.get_sc_heavy_atoms_in_residue(topology,index)
            atomlist=[]
            for a in sc_atoms:
             if a in heavy_atoms:
                    atomlist.append(a)
            ca_atom=Y.get_CA_atoms_in_residue(topology,index)
            pairs=[]
            for i in atomlist:
                pairs+=[([i,ca_atom[0]])]
            pairs=np.asarray(pairs,dtype=int)
            pairs=np.reshape(pairs, (len(pairs), 2))
            dist=md.compute_distances(traj,pairs,periodic=False, opt=True)[0]
            farthest=pairs[np.argmax(dist)][0]

            return traj.xyz[0][farthest]

        def get_farthest_CB_coords(self,pdbfile):
            CB=[]
            Y=conmaps()

            n_res=int(Y.get_total_number_of_residues(pdbfile))
            for i in range(n_res):
               CB.append([self.get_farthest_atom_in_residue(pdbfile,i),i])

            #COM1 = np.array(COM1).reshape(len(COM1) / 2, 2)
            CB=np.array(CB).reshape(len(CB),2)
            print (CB)
            return (CB)


        def get_CB_radius(self,atomtypes,nativefile,sopc,scaling,CBradii):
            print ('>> in get_CB_radius',atomtypes,nativefile,sopc,scaling,CBradii)
            #excluded volume for each cb bead is distance from its C-alpha.
            #sopc gets side-chain radii from BT params!
            scaling=float(scaling)
            Y = conmaps();U=Utils()
            if skip_glycine:
                sopc=False
            if CBradii:
                d=self.read_external_radius('radii.dat')
                assert (len(d)==20),'Must have 20 residues in file!'
            else:
                d = self.amino_acid_radius_dict()
            ncb = Y.get_sidechain_atoms(nativefile)
            #set BT parameters
            #from greddy JPCB Supplementary.(https://pubs.acs.org/doi/suppl/10.1021/acs.jpcb.6b13100/suppl_file/jp6b13100_si_001.pdf)
            assert (atomtypes==2),'Atomtypes not set to two.'
            radius=[]
            for i in ncb:
                radius.append(d[Y.get_residue_name(nativefile,i)[0].strip()])
            radius=np.asarray(radius).reshape(1,len(ncb))
            #radius=(0.8*(radius+3.8))/2
            #print (radius)
            return radius

            # pdbfile='native_cb.pdb'
            # Y=conmaps()
            # CA_indices=Y.get_CA_index(pdbfile)
            # natoms=len(CA_indices)
            # pairs=[]
            # #See Cheung paper JPCB
            # #scaling = 0.7
            # traj=md.load(pdbfile)
            # topology=traj.topology
            # for i in xrange(0,natoms):
            #     #atoms in same residue
            #     s=Y.get_atoms_in_residue(topology,i)
            #     if len(s)>1:
            #         pairs.append(s)
            #     # else:
            #     #     print i
            # print('Found',len(pairs),'C-betas. Each will be assigned a separate atom-type.')
            # #print Y.get_distances(pdbfile, pairs).shape
            # radius=Y.get_distances(traj,pairs)*10
            #return radius


        def get_atom_types(self,pdbfile,atomtypes,sopc):
            print ('>> in get_atom_types:')
            Y=conmaps()
            seq = Y.get_sequence(pdbfile)
            #d = self.amino_acid_dict()
            glyname='G'
            gly = [pos + 1 for pos, char in enumerate(seq) if char == 'G']
            ncb = len(seq) - len(gly)
            if sopc:
                glyname='G111'
                ncb=len(seq)
            natoms = 1 + ncb #one for CA
            count_cb=0
            atomname=[]
            atomname.append('CA')
            for i in seq:
                if i!=glyname:
                    x = 'CB' + str(count_cb)
                    count_cb=count_cb+1
                    atomname.append(x)
            #make dict
            l1=np.arange(0,natoms)
            #print (len(l1),len(atomname))
            assert len(l1)==len(atomname)
            d=self.two_lists_to_dict(atomname,l1)
            return collections.OrderedDict(sorted(d.items()))

        def get_atom_names(self,pdbfile,atomtype,sopc):
            Y=conmaps()
            seq = Y.get_sequence(pdbfile)
            d = self.amino_acid_dict()
            glyname='G'
            if atomtype==1:
                atomname=[]
                for i in seq:
                    atomname.append('CA')
                return atomname
            elif atomtype==2:
                seq = Y.get_sequence(pdbfile)
                d = self.amino_acid_dict()
                if sopc:
                    glyname='G111'
                gly = [pos + 1 for pos, char in enumerate(seq) if char == glyname]
                nca = len(seq)

                ncb = len(seq) - len(gly)
                natoms = len(seq) + ncb
                count_cb=0
                atomname=[]
                for i in seq:
                    atomname.append('CA')
                    if i!=glyname:
                        x = 'CB' + str(count_cb)
                        count_cb=count_cb+1
                        atomname.append(x)
                #print atomname,natoms

                assert len(atomname)==natoms
            #make dict
                l1=np.arange(0,natoms)
                d=self.two_lists_to_dict(l1,atomname)
                #print d
                return d

        def write_header_SBM(self):
            f = open('SBM.INP', "w+")
            f.write('%s\n' % ('; Topology file generated from Go-kit. https://github.org/gokit1/gokit/'))
            f.write('%s\n' % ('Debye-Huckel Parameters: PREFACTOR, Dielectric Constant, Monovalent Ion Concentration, DH switching distance, DH cutoff distance'))
            f.write('%s\n' % (' 332.000   80.000    0.100   12.000   15.000'))

            f.write('%s\n' % ('nonbonded switching distance, truncate distance'))
            f.write('%s\n' % ('7.999999 8.999999'))
            f.close()

        def write_atomtypes_section(self, pdbfile, atomtype,sopc,CA_rad,CBcom,CBradii):
            CA_rad=float(CA_rad)
            print ('>> in write_atomtypes_section',pdbfile,atomtype,sopc,CA_rad)
            assert atomtype<=2
            Y=conmaps()
            f = open('SBM.INP', "a")
            if skip_glycine:
                sopc=False
            seq = Y.get_sequence(pdbfile);nca=len(seq);d=self.amino_acid_dict();glyname='G'
            if atomtype==2:
                #gly = [pos + 1 for pos, char in enumerate(seq) if char == 'G']
                ncb = len(Y.get_CB_index('native_cb.pdb'))
                if skip_glycine:
                    glyname='G'
                natoms = nca + ncb
                #print (nca,ncb)
                assert(natoms==Y.get_total_number_of_atoms('native_cb.pdb'))
            elif atomtype==1:
                natoms=len(seq)


            natomtypes=len(self.get_atom_types(pdbfile,atomtype,sopc))
            if atomtype==1:
                natomtypes=1
            if w_sbm:
             f.write('%d %s\n' % (natomtypes,' atomtypes'))
            atomname=[]
            atommass=np.ones(natoms,dtype=np.float32)
            atomcharge=np.zeros(natoms,dtype=np.float32)
            atomptype=['A']*natoms
            CA_rad = np.ones(nca, dtype=np.float32) * CA_rad
            #default CB radius
            #print len(CB_rad[0]),len(CA_rad),len(gly)
                    #assert len(CB_rad[0]) + len(CA_rad) == natoms
            count_ca=0
            c6=0.000
            if w_sbm:
             f.write('%d\t%8.5f\t%s\n' % (1, CA_rad[count_ca], '1.00'))

            if atomtype==1:
                atomname.append(['CA', atommass[0], atomcharge[0], atomptype[0], c6, CA_rad[0]])
                #print (atomname)
                return atomname

            elif atomtype==2:
                #add CA first.
                atomname.append(['CA', atommass[0], atomcharge[0], atomptype[0], c6, CA_rad[0]])
                CB_rad = self.get_CB_radius(atomtype, 'native_cb.pdb', sopc, 1, CBradii)
                print (len(CB_rad[0]))
                assert(len(CB_rad[0])==ncb)
                count_cb=0
                #loop over CA!
                for i in seq:
                    if i != glyname:
                        x = 'CB' + str(count_cb)
                        #print (count_cb)
                        #print (CB_rad[0][count_cb],count_cb,i,len(seq ))
                        #print (x,ncb,i,count_cb,glyname)
                        #print (atomcharge[count], atomptype[count], c6, CB_rad[0][count_cb], d[i])
                        #print (x, '1', '0', atomptype[count_cb], c6, CB_rad[0][count_cb],)
                        atomname.append([x,atommass[count_cb],atomcharge[count_cb],atomptype[count_cb],c6,CB_rad[0][count_cb],d[i]])
                        if w_sbm:
                            f.write('%d\t%8.5f\t%s\n' % (count_cb+2, CB_rad[0][count_cb], '1.00'))
                        count_cb = count_cb + 1
                #print (count_cb)
                return atomname



        def write_atoms_section(self,pdbfile,atomtype,skip_glycine):
            print ('>> in write_atoms_section\t',atomtype,pdbfile,atomtype,sopc)
            assert atomtype<=2
            Y=conmaps()
            seq = Y.get_sequence(pdbfile)
            d=self.amino_acid_dict()
            f = open('SBM.INP', "a")
            nca=len(seq)
            gly= [pos+1 for pos,char in enumerate(seq) if char== 'G']
            print('Glycine in residues:', gly)
            ncb=len(seq) - len(gly)
            glyname='GLY'
            if skip_glycine:
                glyname='GLY'
            str1= 'atoms (atomnum, atomtype, resnum, resname, atomname, charge, mass)'
            g=self.get_atom_names(pdbfile,atomtype,skip_glycine)
            h=self.get_atom_types(pdbfile,atomtype,skip_glycine)
            atomtype_ca=1
            atomtype_cb=2
            charge=0.00000
            mass=1.00000
            if atomtype==2:
                natoms = nca + ncb
                if w_sbm:
                    f.write('\n%d %s\n' % (natoms, str1))
            if atomtype==1 and not dswap:
                natoms=nca
                if w_sbm:
                    f.write('\n%d %s\n' % (nca, str1))
            if atomtype==1 and dswap:
                natoms=nca
                if w_sbm:
                    f.write('\n%d %s\n' % (2*nca, str1))

            #write Calphas first
            rescount=0
            count=0
            all=[]
            for i in seq:
                all.append([count+1, atomtype_ca,rescount+1,d[i],'CA',charge,mass])
                if w_sbm:
                 f.write(' %d\t%d\t%d\t%s\t%s\t%f\t%f\n' % (count+1, h[g[count].strip()]+1,rescount+1,d[i],g[count],charge,mass))
                count=count+1
                if atomtype==2:
                    #print count
                    if d[i].strip()!=glyname:
                        all.append([count + 1, atomtype_cb, rescount + 1, d[i], 'CB', charge, mass])
                        if w_sbm:
                         f.write(' %d\t%d\t%d\t%s\t%s\t%f\t%f\n' % (count + 1, h[g[count].strip()]+1, rescount + 1, d[i], g[count], charge, mass))
                        count=count+1
                rescount=rescount+1
            print (count,natoms,glyname)
            assert count==natoms

            return all

        def split_chains(self,pdbfile):
            U=Utils()
            U.split_chain(pdbfile)

        def get_pair_coefficients(self,pottype):
            return


        def write_contacts_section(self,pdbfile,nativefile,btfile,cutoff,atomtype,sopc,btparams,dswap):
            assert atomtype <= 2
            U=Utils()
            print ('>>writing contacts section,atomtype=',atomtype,'SOPC=',sopc)
            cutoff=float(cutoff)
            Y = conmaps()
            d = self.amino_acid_dict2()
            f = open('SBM.INP', "a")
            f1= open('contacts.txt', "w+")
            #f2 = open('backbone.txt', "w+")
            #f3 = open('sidechain.txt', "w+")
            #f4 = open('bb_sc.txt', "w+")
            traj=md.load(nativefile)
            contacttype=2
            strength_CA=float(1.0) #eps bb-bb
            strength_CB=float(1.0) #eps bb-sc
            # if U.file_exists(btfile) and btparams:
            #     btmap = self.get_dict_from_csv(btfile)
            #get CA-CA contacts.(CA atomtype=1,CB atomtype=2)
            #dict_seq = self.two_lists_to_dict(np.arange(0, len(seq)).tolist(), list(seq))
            #nc=native,sc=side-chain
            if atomtype==1:
                self.write_CB_to_native(pdbfile,skip_glycine)
                nc = Y.all_atom_contacts(pdbfile, cutoff,1.2)
                dist_ca = nc[:, 3];
                #print (dist_ca)
                contacts_ca = nc[:, (0, 1)]
                l1 = Y.get_CA_index(nativefile)  # CA from native
                num_atoms = len(l1)
                l2 = np.arange(0, len(l1))
                assert len(l1) == len(l2)
                dict1 = self.two_lists_to_dict(l2, l1)
             #   print (dict1)
                num_contacts=len(contacts_ca)
                if hphobic:
                    hp_pairs=Y.write_hydrophobic_contacts('contacts.txt',pdbfile)
                    num_contacts=num_contacts+len(hp_pairs)
                if w_sbm:
                 f.write('\n%d\t%s\n' % (num_contacts, 'contacts'))
                contacts = []
                count=0
                #print ("CA_contacts\n")
                if dsb:
                    contacttype=8;strength_CA=1.00000e+00
                for i in contacts_ca:
                    if dsb:
                        contacts.append([dict1[int(i[0])] + 1, dict1[int(i[1])] + 1, contacttype, dist_ca[count] * 10, strength_CA])
                    contacts.append([dict1[int(i[0])]+1,dict1[int(i[1])]+1,contacttype,dist_ca[count]*10,strength_CA])
             #       f2.write(' %s\t%d\t%s\t%d\n' % ('1', dict1[int(i[0])] + 1, '1', dict1[int(i[1])] + 1))
                    if w_sbm:
                        f.write('%s\t%s\t%d\t%f\t%12.9e\n'%(dict1[int(i[0])]+1,dict1[int(i[1])]+1,contacttype, dist_ca[count]*10,strength_CA))
                    count=count+1
                if hphobic:
                    print("hydrophobic contacts\n")
                    print ("Adding hydrophobic contacts: Distance =",hpdist, "Angstroms")
                    #hp_pairs=Y.write_hydrophobic_contacts('contacts.txt',pdbfile)
                    p=hp_pairs-1
                    #dist_hp = md.compute_distances(traj_hp, p, periodic=False, opt=True)[0][0]
                    counthp=0
                    for i in p:
                        #print (i[0],i[1],counthp)
                        contacts.append([dict1[int(i[0])]+1,dict1[int(i[1])]+1,contacttype,hpdist,hpstrength])
                        #f2.write(' %s\t%d\t%s\t%d\n' % ('1', dict1[int(i[0])] + 1, '1', dict1[int(i[1])] + 1))
                        if w_sbm:
                            f.write('%s\t%s\t%d\t%f\t%e\n' % (
                                dict1[int(i[0])] + 1, dict1[int(i[1])] + 1, contacttype, hpdist, hpstrength))
                        counthp = counthp + 1
                f.close()
                return contacts
            elif atomtype==2:
                contacts = []
                #scaling = 1.2;separation=3
                if skip_glycine:
                    sopc=False
                if sopc:
                   nc = Y.get_bb_contacts_SOPC(pdbfile,nativefile,4,1.0,4)
                   sc = Y.get_ss_contacts_SOPC(pdbfile,nativefile,4,1.0,2)
                   nc_sc=Y.get_bs_contacts_SOPC(pdbfile,nativefile,4,1.0,2)
                   nn_angles = Y.get_SOPC_non_native_angles(pdbfile, nativefile, 3.8, 1.0)
                else:
                    #separation is always >=
                    #cacasep=4;cacbsep=3;cbcbsep=3
                    nc = Y.get_backbone_contacts(pdbfile, cutoff,scaling,casep)
                    sc = Y.get_side_chain_contacts(pdbfile, cutoff,scaling,3)
                    nc_sc = Y.get_bb_sc_contacts(pdbfile, cutoff,scaling,3)
                topology=md.load(nativefile).topology
                #f.write('\n%d\t%s\n' % (ncontacts, 'contacts'))
                print ('Found', len(nc), 'backbone contacts')
                print ('Found', len(sc), 'side-chain contacts')
                print ('Found', len(nc_sc), 'bb-sc contacts')
                CA_indices = Y.get_CA_index(nativefile)
                num_ca = len(CA_indices)
                #residue_list=np.arange(0,num_ca)
                #d_ca=self.two_lists_to_dict(residue_list,CA_indices)

                #write native-contacts first.
                for i in xrange(0, len(nc)):
                    res1 = int(nc[i][0]); res2 = int(nc[i][1])
                    #print nc[i][0],nc[i][1]
                    res1 = Y.get_atoms_in_residue(topology, res1);res2 = Y.get_atoms_in_residue(topology, res2)
                    assert len(res1) >= 1;assert len(res2) >= 1
                    res1 = res1[0];res2 = res2[0]
                    pair = ([res1, res2])
                    pair = np.reshape(pair, (1, 2))
                    dist = md.compute_distances(traj, pair, periodic=False, opt=True)[0][0]*10
                    #f2.write(' %s\t%d\t%s\t%d\n' % ('1', int(res1)+1, '1', int(res2)+1))
                    contacts.append([res1+1, res2+1, int(contacttype), dist, strength_CA])

                # write side-chain contacts next
                sc_contacts=[]
                for i in xrange(0, len(sc)):
                    ##map CB in residue to CB in two bead file.
                    res1 = int(sc[i][0]);res2=int(sc[i][1])
                    #print (res1,res2,nativefile)
                    res1=Y.get_atoms_in_residue(topology,res1);res2=Y.get_atoms_in_residue(topology,res2)
                    #print (res1,res2)
                    assert len(res1)==2;assert len(res2)==2
                    #second atom is always CB. Put assertion.
                    res1=res1[1];res2=res2[1]
                    #print res1,res2,int(sc[i][0]),int(sc[i][1])
                    dist = sc[i][3]
                    #print (dist)
                    #add two-bead potentials here!!
                    #btparams=True
                    if btparams:
                        #print (" --interactions flag set to true. External interaction matrix will be read from file p.interaction.dat")
                        boltzmann_kcal_mol=float(0.0019872041)
                        Kb=1*boltzmann_kcal_mol
                        #print d[Y.get_residue_name(pdbfile,int(res1))[0]i]
                        res1bt = d[Y.get_residue_name(nativefile, int(res1))[0]]
                        res2bt = d[Y.get_residue_name(nativefile, int(res2))[0]]
                        strength_CB=float(Y.get_btmap_val(res1bt,res2bt,'interaction.dat'))
                        strength_CB=0.5*(0.7- float(strength_CB))*300*Kb #kcal/mol
                        contacts.append([res1 + 1, res2 + 1, int(contacttype), dist,strength_CB])
                        #print ([res1 + 1, res2 + 1, int(contacttype), dist, strength_CB, 'SCSCbt'])
                    else:
                        contacts.append([res1+1, res2+1, int(contacttype), dist, strength_CB])
                        #print ([res1 + 1, res2 + 1, int(contacttype), dist, strength_CB,'SCSC'])
                        #f3.write(' %s\t%d\t%s\t%d\n' % ('1', int(res1) + 1, '1', int(res2) + 1))

                #write bb-sc contacts next
                for i in xrange(0,len(nc_sc)):
                    res1 = int(nc_sc[i][0]);res2 = int(nc_sc[i][1]);
                    res1 = Y.get_atoms_in_residue(topology, res1); res2 = Y.get_atoms_in_residue(topology, res2)
                    assert len(res1)>=1 and len(res2)>=1
                    #res1 = CA, res2=CB
                    res1=res1[0];res2=res2[1]
                    #print (res1, res2)
                    #exit()
                    pair = ([res1, res2])
                    pair = np.reshape(pair, (1, 2))
                    dist = md.compute_distances(traj, pair, periodic=False, opt=True)[0][0] * 10
                    contacts.append([res1 + 1, res2 + 1, int(contacttype), dist, strength_CA])
                    #print ([res1 + 1, res2 + 1, int(contacttype), dist, strength_CA, 'bbsc'])
                    #f4.write(' %s\t%d\t%s\t%d\n' % ('1', int(res1) + 1, '1', int(res2) + 1))
                #write angle-angle repulsion term for sop-sc model only!
                if sopc:
                    epsl=1 #kcal/mol
                    for i in xrange(0,len(nn_angles)):
                        contacttype=1
                        #Use the 6 term from LJ potential for repulsion!(supply negative to make repulsive)
                        # fudgeQQ,qi,qj,V,W for gromacs
                        res1=int(nn_angles[i][0]);res2=int(nn_angles[i][1])
                        #print res1,res2
                        pair = ([res1, res2]);  pair = np.reshape(pair, (1, 2))
                        #print pair
                        dist = md.compute_distances(traj, pair, periodic=False, opt=True)[0][0] * 10
                        contacts.append([res1 + 1, res2 + 1, int(contacttype),dist,epsl])


            else:
                U.fatal_errors(7)

            # check if contacts are being repeated.
            contacts=np.asarray(contacts)
            contacts= contacts[np.argsort(contacts[:,0])]
            test=[]
            for i in contacts:
                j=tuple(i)
                print (j)
                test.append(j)
            #remove repeat contacts.
            a=set(test)
            print ("Found",len(a)-len(contacts),"duplicates in list. Removing...")
            a=list(a)
            f.write('\n%d\t%s\n' % (len(contacts), 'contacts'))
            for i in a:
                #print i[0],i[1],i[2],i[3],i[4]
                #"%.10lf\n
                f1.write(' %s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('1',i[0],'1',i[1], i[2], i[3], i[4]))
                if w_sbm:
                    f.write(' %d\t%d\t%d\t%.12f\t%10.3f\n' % (i[0],i[1],i[2],i[3],i[4]))
            f.close()
            f1.close();#f2.close()#;f3.close();f4.close()
            return contacts

        def write_exclusions_section(self):
            print (">>Writing exclusions section.")
            f = open('SBM.INP', "a")
            #contacts=self.write_contacts_section(pdbfile,nativefile,btfile,cutoff,atomtype,sopc,btparams,dswap)
            f.write('\n%s\n' % ('0 exclusions'))
            f.write('\n%s\n' % ('0 position restraints'))
            f.close()
            return


        def write_bonds_section(self,nativefile,Kb,atomtype):
            Y=conmaps()
            print (">>writing bonds_section: Native= ",nativefile)
            f = open('SBM.INP', "a")
            #nativefile='native_cb.pdb'
            CA_indices=Y.get_CA_index(nativefile)
            num_ca=len(CA_indices)
            CB_indices=Y.get_CB_index(nativefile)
            num_cb=len(CB_indices)
            indices=Y.get_all_atom_indices(nativefile)
            num_atoms=len(indices)
            assert (num_atoms==num_ca+num_cb)
            pairs_CA=[]
            pairs_CB=[]
            top=md.load(nativefile).topology
            for i in xrange(len(CA_indices)-1):
                assert (Y.get_atom_name(top, CA_indices[i])[0].strip() == 'CA')
                pairs_CA.append([CA_indices[i],CA_indices[i+1]])
            # print (Y.get_residue_number(nativefile, CA_indices[i]))
            traj=md.load_pdb(nativefile)
            CA_dist = Y.get_distances(traj, pairs_CA) * 10
            traj = md.load(nativefile)
            topology = traj.topology
            for i in xrange(num_ca):
                #if two atoms per residue then make pair, else GLYCINE
                j=Y.get_atoms_in_residue(topology,i)
                if len(j)==2:
                    pairs_CB.append(j)
            num_CAB=0
            CAB_dist=np.array([])
            if pairs_CB:
                CAB_dist=Y.get_distances(traj,pairs_CB) * 10
                num_CAB=CAB_dist.size
            num_bonds=num_CAB + CA_dist.size
            if w_sbm:
             f.write('\n%d\t%s\n' % (num_bonds, 'bonds'))
            num_bonds=num_CAB + CA_dist.size
            bondtype=1
            #read externally.
            K_bond = np.array(Kb, dtype=float)
            assert len(CA_dist.T)==len(pairs_CA)
            pairs=np.asarray(pairs_CA+pairs_CB) #cast lists together
            #print CAB_dist.size
            if CAB_dist.size > 0 :
                dists=np.asarray(CA_dist.T.tolist()+CAB_dist.T.tolist())
            else :
                dists=np.asarray(CA_dist.T.tolist())
            btype = np.ones(num_bonds).T*bondtype
            #print len(pairs),len(dists),len(btype)
            assert len(pairs) == len(dists) == len(btype)
            all=np.hstack((pairs,dists))
            all=all[np.argsort(all[:, 0])]
            count=0
            all1=[]
            for i in all:
                #print Y.get_atom_name(nativefile,i[0]),Y.get_atom_name(nativefile,i[1]),dists[count]
                all1.append([i[0]+1, i[1]+1, bondtype, i[2], K_bond])
                if w_sbm:
                    f.write(' %d\t%d\t%d\t%8.5f\t%e\n' % (i[0]+1, i[1]+1, bondtype, i[2], K_bond))
                    #print (dsb_pairs)
            f.close()
            return all1

        def write_angles_section(self,nativefile,Ka):
            #dswap=False
            print ('>>writing angles_section',nativefile,Ka)
            Y=conmaps()
            coords = np.asarray(Y.get_pdb_coords(nativefile) )
            num_atoms=len(coords)
            triplets=[]
            for i in xrange(1,len(coords)-1):
                triplets.append([i-1, i, i+1])
            #triplets=list(combinations(CA_indices,3))
            num_angles=len(coords)-2
            f = open('SBM.INP', "a")
            if w_sbm:
             f.write('\n%d\t%s\n' % (num_angles, 'angles'))
            CA_angles=Y.get_angles(nativefile,triplets)
            K_angle=float(Ka)
            #write CA-CA angles.
            count=0
            all=[]
            count=0
            for i in xrange(len(triplets)):
                #print triplets[i],CA_angles[0][i]
                all.append([triplets[i][0]+1,triplets[i][1]+1,triplets[i][2]+1, CA_angles[0][i],K_angle])
                if w_sbm:
                 f.write(' %s\t%s\t%d\t%f\t%e\n' % (triplets[i][0]+1,triplets[i][1]+1,triplets[i][2]+1, CA_angles[0][i],K_angle))
                count=count+1
            if dswap:
                 count=0
                 for i in xrange(len(triplets)):
                     all.append([triplets[i][0]+1+num_atoms,triplets[i][1]+1+num_atoms,triplets[i][2]+1+num_atoms, CA_angles[0][i], K_angle])
                     if w_sbm:
                      f.write(' %s\t%s\t%d\t%f\t%e\n' % (triplets[i][0]+1+num_atoms,triplets[i][1]+1+num_atoms,triplets[i][2]+1+num_atoms, CA_angles[0][i],K_angle))
                     count=count+1
            f.close()
            return  all

        def write_dihedrals_section(self,nativefile,Kd):
            #Only over C-alpha atoms for Cheung-Thirumalai model.
            #dswap=False
            print ('>>writing dihedrals section',nativefile,Kd)
            Y = conmaps()
            f = open('SBM.INP', "a")
            #Select only C-alpha atoms.
            indices = Y.get_CA_index(nativefile)
            num_atoms=len(indices)
            #print indices
            quadruplets = []
            K_dihedral=float(Kd)
            for i in xrange(2, len(indices) - 1):
                quadruplets.append([indices[i-2],indices[i-1], indices[i], indices[i+1]])
            if dswap:
                num_dihedrals=2*len(quadruplets)
            else:
                num_dihedrals = len(quadruplets)
            if w_sbm:
             f.write('\n%d\t%s\n' % (num_dihedrals, 'dihedrals'))
            dihedrals=Y.get_dihedrals(nativefile,quadruplets)
            count=0
            all=[]

            count=0
            for i in xrange(len(quadruplets)):
                #print triplets[i],CA_angles[0][i]
                all.append([quadruplets[i][0],quadruplets[i][1],quadruplets[i][2],quadruplets[i][3], '1 ', dihedrals[0][i],K_dihedral])
                if w_sbm:
                 f.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\n'% (quadruplets[i][0]+1,quadruplets[i][1]+1,quadruplets[i][2]+1,quadruplets[i][3]+1, '1 ', dihedrals[0][i],K_dihedral))
                count=count+1
            if dswap:
                for i in xrange(len(quadruplets)):
                    all.append([quadruplets[i][0]+num_atoms, quadruplets[i][1]+num_atoms, quadruplets[i][2]+num_atoms, quadruplets[i][3]+num_atoms, '1 ',
                                dihedrals[0][i], K_dihedral])
                    f.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\n' % (
                    quadruplets[i][0]+1+num_atoms, quadruplets[i][1]+1+num_atoms, quadruplets[i][2]+1+num_atoms, quadruplets[i][3]+1+num_atoms, '1 ', dihedrals[0][i], K_dihedral))
            all=np.asarray(all)
            f.close()
            return all

########GROMACS topology section. Make separate class later.##########
        def write_gro_moleculetype(self,topfilename):
            print ('>>writing GROMACS moleculetypesection', topfilename)
            f = open(topfilename, "a")
            if w_gro:
             f.write('%s\n' % ('[ moleculetype ]'))
             f.write('%s\n' % ('; name            nrexcl'))
             f.write('%s\n' % ('  Macromolecule   3'))
            f.close()
            return 1
        def write_gro_header(self,topfilename,atomtypes):
            print ('>>writing GROMACS header section', topfilename)
            f = open(topfilename, "w+")
            if w_gro:
             f.write('%s\n' % (';'))
             f.write('%s\n' % ('; Topology file generated from Go-kit. '))
             f.write('%s\n' % ('; https://github.org/gokit1/gokit/wiki/Home'))
             f.write('\n%s\n' % ('[ defaults  ]'))
             f.write('%s\n' % ('; nbfunc comb-rule gen-pairs'))
             f.write('%s\n' % ('  1      1         no   \n'))
            return





        def write_gro_nonbondparams(self,topfilename,sopc):
            assert (sopc)
            if not w_gro:
                return
            f = open(topfilename, "w+")
            if w_gro:
             f.write('%s\n' % ('[ nonbond_params ]'))
            #add non-bonded r6 term in sop-sc model for all non-bonded non-native interactions.
             f.write('%s\n' % ('; i j  func sigma(c6)    eps(c12)'))
            f.close()
            return 1
        def write_gro_atomtypes(self,topfilename,atomtypes,pdbfile,sopc,CA_rad,CBcom,CBradii):
            if not w_gro:
                return
            print (">> writing GROMACS atomtypes section",topfilename,atomtypes,pdbfile,sopc,CA_rad)
            #CA_rad=float(3.8)
            f = open(topfilename, "a")
            assert atomtypes<=2
            a=self.write_atomtypes_section(pdbfile,atomtypes,sopc,CA_rad,CBcom,CBradii)
            assert a[0][0]=='CA'
            b=a[0]
            CA_rad=(CA_rad/10)**12
            f.write('%s\n'%('[ atomtypes ]'))
            f.write('%s\n' % ('; name mass  charge ptype c6    c12'))
            f.write('%s %8.3f %8.3f %s\t%e\t%e \n' % (b[0], b[1], b[2], b[3], b[4], CA_rad))
            for i in (a):
                if i[0]!='CA':
                    r=(i[5]/10)**12 #C12 term
                    i[0]=i[0].rjust(5)
                    f.write('%s %8.3f %8.3f %s\t%e\t%e \n'%(i[0],i[1],i[2],i[3],i[4],r))
            f.close()

        def write_gro_atoms(self,topfilename,atomtypes,pdbfile,nativefile,sopc):
            if not w_gro:
                return
            print (">> writing GROMACS atom section",topfilename,atomtypes,pdbfile,nativefile,sopc)
            f = open(topfilename, "a")
            assert atomtypes<=2
            atoms = self.write_atoms_section(pdbfile,atomtypes,sopc)
            num_atoms=len(atoms)
            #print len(atoms)
            f.write('\n%s\n' % ('[ atoms ]'))
            f.write('%s\n' % (';nr  type  resnr residue atom  cgnr'))
            d = self.get_atom_names(pdbfile,atomtypes,sopc)
            e=self.get_atom_types(pdbfile,atomtypes,sopc)
            for i in xrange(0,len(atoms)):
                type=int(e[d[i].strip()])
                d[i]=d[i].rjust(4)
                f.write('\t%d\t%s\t%d\t%s\t%s\t%d \n' %(atoms[i][0],d[i],atoms[i][2],atoms[i][3],d[i].strip()[:2],atoms[i][0]))
            if dswap:
                for i in xrange(0, len(atoms)):
                     type = int(e[d[i].strip()])
                     d[i] = d[i].rjust(4)
                     f.write('\t%d\t%s\t%d\t%s\t%s\t%d \n' % (
                     atoms[i][0]+num_atoms, d[i], atoms[i][2]+num_atoms,atoms[i][3], d[i].strip()[:2], atoms[i][0]))
            f.close()
            return

        def write_gro_bonds(self,topfilename,atomtypes,nativefile,Kb,ptype,dsb):
            if not w_gro:
                return
            #ptype=1 #READ POTENTIAL FROM DICTIONARY HERE IF NEEDED.
            print (">> writing GROMACS bonds section", topfilename, atomtypes,nativefile,Kb)
            Kb = float(Kb * 100)
            #GROMACS 4.5.2 : FENE=7 AND HARMONIC=1
            allowed_pots=[1,7,9];assert ptype in allowed_pots;assert atomtypes <= 2
            Y=conmaps()
            num_atoms=Y.get_total_number_of_atoms(nativefile)
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ bonds ]'))
            f.write('%s\t%s\t%s\t%s\t%s\n' % (';ai', 'aj', 'func', 'r0(nm)', 'Kb'))
            bonds =self.write_bonds_section(nativefile,Kb,atomtypes)
            assert ((ptype in allowed_pots)),'Make sure potential is included in allowed_pots.'
            if dsb:
                ptype=9
            if ptype==1:
                for i in xrange(0,len(bonds)):
                    f.write('\t%d\t%d\t%s %12.9e %12.9e\n' %(bonds[i][0],bonds[i][1], str(ptype),bonds[i][3]/10,Kb))
                f.close()
                #print (bonds)
                return bonds
            elif ptype==7:
                #FENE potential.
                kcalAtokjA=418.4 #kcal/mol/A2 to kcal/mol/nm2
                Kb=float(20*kcalAtokjA) #Gromacs units
                R=float(0.2) #nm
                for i in xrange(0,len(bonds)):
                    #not sure of format for top file. (where does distance go?)
                    f.write('\t%d\t%d\t%s %12.9e %12.9e %12.9e\n' %(bonds[i][0]+1,bonds[i][1]+1, str(ptype),bonds[i][2]/10,R,Kb))
                f.close()
            elif ptype==9:
                for i in xrange(0,len(bonds)):
                    f.write('\t%d\t%d\t%s %12.9e %12.9e\n' %(bonds[i][0],bonds[i][1], str(ptype),bonds[i][3]/10,Kb))
                #contacts written in bonds section.
                #get_contacts from contacts.txt
                dsb_contacts=Y.get_pairs_ext('contacts.txt')
                print (dsb_contacts[0][0])
                count_dsb=0
                for i in xrange(0,len(dsb_contacts)):
                    f.write('\t%d\t%d\t%s %s %s\n' % (dsb_contacts[i][0],dsb_contacts[i][1], str(ptype),str(count_dsb), '1'))
                    count_dsb=count_dsb+1
                #generate table_files here
                Z=tables()
                if hphobic:
                    Z.gen_many_table_file('contacts.txt', 'native_ca.pdb', 0.06, 350, 0.02,True)
                else:
                    Z.gen_many_table_file('contacts.txt','native_ca.pdb',0.06,350,0.02,False)
                return bonds




        def write_gro_angles(self,topfilename,atomtypes,nativefile,Ka):
            if not w_gro:
                return
            print (">> writing GROMACS angle section", topfilename, atomtypes, nativefile, Ka)
            assert atomtypes <= 2
            Y=conmaps()
            Ka=float(Ka)
            num_atoms=Y.get_total_number_of_atoms(nativefile)
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ angles ]'))
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (';ai', 'aj', 'ak','func', 'th0(deg)', 'Ka'))
            angles=self.write_angles_section(nativefile,Ka)
            for i in xrange(0, len(angles)):
                f.write('\t%d\t%d\t%d\t%s %12.9e %12.9e\n' % (angles[i][0],angles[i][1],angles[i][2],'1',angles[i][3]*radtodeg, Ka))
            f.close()
            return

        def write_gro_dihedrals(self,topfilename,atomtypes,nativefile,Kd):
            if not w_gro:
                return
            print (">> writing GROMACS dihedrals section", topfilename, atomtypes, nativefile, Kd)
            Y=conmaps()
            num_atoms=Y.get_total_number_of_atoms(nativefile)
            pi=np.pi
            Kd=float(Kd)
            assert atomtypes <= 2
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ dihedrals ]'))
            d=self.write_dihedrals_section(nativefile,Kd)
            f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (';ai','aj','ak','al','func','phi0(deg)','Kd','mult'))
            for i in xrange(0, len(d)):
                #GROMAC dihedrals in degrees. For dihedraltype 1 add 180.
                l1= float(d[i][5])*radtodeg+180
                d1=int(d[i][0])
                d2 = int(d[i][1])
                d3 = int(d[i][2])
                d4 = int(d[i][3])
                d5 = int(d[i][4])
                if d5!=1:
                    U.fatal_errors(8)
                #print l1,d1,d2,d3,d4
                f.write('\t%d\t%d\t%d\t%d %d %12.9e %12.9e %s\n' % (d1+1, d2+1, d3+1, d4+1,d5,l1,Kd,'1'))
                f.write('\t%d\t%d\t%d\t%d %d %12.9e %12.9e %s\n' % (d1+1, d2+1, d3+1, d4+1,d5,l1*3,Kd/2,'3'))
            f.close()
            return

        def write_gro_pairs(self,topfilename,atomtypes,nativefile,pdbfile,contacttype,cutoff,sopc,btparams):
            if not w_gro:
                return
            print (">> writing GROMACS pairs sections",topfilename,atomtypes,nativefile,pdbfile,contacttype,cutoff,sopc,btparams)
            contacts_allowed=[1,2]
            assert atomtypes <= 2
            if skip_glycine:
                sopc=False
            #OPTIM works in kcal, GROMACS in kJ
            Y=conmaps()
            #num_atoms=Y.get_total_number_of_atoms(nativefile)
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ pairs ]'))
            if dsb:
                return
            f.write('\t%s\n' % ('; ai aj type A B'))
            btfile='interaction.dat'
            #cutoff=4.5
            contacts=self.write_contacts_section(pdbfile,nativefile,btfile,cutoff,atomtypes,sopc,btparams,dswap)
            assert (contacttype in contacts_allowed),'Only the 10-12 potential contacttype=2 is implemented'
            assert atomtypes<=2
            count_angles=0
            for i in xrange(0,len(contacts)):
                contacttype=int(contacts[i][2])
                #print contacttype
                if contacttype==2:
                    #SBM contactype=2 => GROMACS contacttype=1
                    #A=5*epsilon*(sigmaij)^12
                    # gromacs contact type is 1 always for table potentials.
                    # optim contact types differ.
                    epsilonij=float(contacts[i][4])
                    #print epsilonij
                    #12-10 potential.
                    # Add all sorts of potentials here. Change.
                    #default is 12-10 with table files.
                    #Add 12-6 LJ here!
                    B=((contacts[i][3]*0.1)**12)*5*epsilonij
                    A=((contacts[i][3]*0.1)**10)*6*epsilonij
                    #i-j
                    f.write('\t%d\t%d\t%s\t%12.9e\t%12.9e\n' % (contacts[i][0],contacts[i][1],'1',A,B))
                # elif contacttype==1 and sopc:
                #     count_angles=count_angles+1
                #     A=float(0.0)
                #     B=-(contacts[i][3]*0.1)**6 #sigmabb*6 and (sigma_bs)**6 (minus for repilsive)
                #     #print contacts[i][3]
                #     q1=0;q2=0;fudgeqq=1
                #     #Fudgeqq,qi,qj,V,W
                #     f.write('\t%d\t%d\t%s\t%d\t%e\t%e\t%e\t%e\n' % (contacts[i][0], contacts[i][1],'2',fudgeqq,q1,q2,B,A))
                #     print ("Writing", count_angles, "additional entries in pairs section. SOPC =",sopc)
            #write exlusions
            f.write('\n\t%s\n' % ('[ exclusions ]'))
            f.write('\t%s\n' % ('; ai aj'))
            for i in xrange(0,len(contacts)):
                f.write('\t%d\t%d\n'%(contacts[i][0],contacts[i][1]))
            f.close()
            return

        def write_gro_tail(self,topfilename):
            if not w_gro:
                return
            print (">> writing GROMACS tail section",topfilename)
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ system ]'))
            f.write('%s\n' % (';name'))
            f.write('  %s\n' % ('Macromolecule'))
            f.write('\n%s\n' % ('[ molecules ]'))
            f.write('%s\n' % (';name    #molec'))
            f.write('%s\n' % ('Macromolecule     1'))
            f.close()
            return
        def write_odata_header(self):
            f1=open('odata',"w+")
            f1.write("%s\n" % ('STEPS 10000'))
            f1.write("%s\n" % ('BFGSMIN 1.0D-3'))
            f1.write("%s\n" % ('POINTS'))
            f1.close()

        def write_gro_gro(self,grofilename,pdbfile,atomtypes,skip_glycine):
            global w_sbm;w_sbm=False
            self.write_CB_to_native(pdbfile,skip_glycine)
            sopc=True
            if skip_glycine:
                sopc=False
            print ('>> write_gro_gro',grofilename,pdbfile,atomtypes,sopc)
            Y=conmaps()
            f = open(grofilename, "w+")
            self.write_odata_header()
            f1= open('odata', "a")
            f.write('%s\n'%('Generated with Go-kit'))
            a = self.write_atoms_section(pdbfile,atomtypes,sopc)
            if atomtypes==1:
                c=Y.get_coordinates('native_ca.pdb')
                d = self.get_atom_names(pdbfile, atomtypes,sopc)
                #print d
            elif atomtypes==2:
                c=Y.get_coordinates('native_cb.pdb')
            else:
                U.fatal_errors[11]
                #print ('atomtypes must be 1 or 2')
            f.write('%d\n'%(len(a)))
            #print ((len(c),len(a)))
            #assert len(c)==len(a)
            d=self.get_atom_names(pdbfile,atomtypes,sopc)
            count=0
            for i in xrange(0,len(a)):
                #print a[i][2],a[i][3],a[i][4],d[count],c[i][0],c[i][1],c[i][2],a[i][0]
                f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (a[i][2],a[i][3],d[count][:2],a[i][0],c[i][0],c[i][1],c[i][2]))
                f1.write("%-5s%8.3f%8.3f%8.3f\n" % ('SB',c[i][0],c[i][1],c[i][2]))
                count=count+1
            f.write('%8.4f %8.4f %8.4f'%(50.000,50.000,50.0000))
            f.close()
            f1.close()
            return
        def set_false_global_wsbm(self):
            global w_sbm
            w_sbm=False
            return


        def write_gromacs_top(self,topfilename,atomtypes,pdbfile,nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii):
            #Order for SBM file
            #global w_sbm, w_gro
            assert atomtypes <= 2
            print (w_sbm)
            #SBM.INP
            self.write_header_SBM()
            self.write_atomtypes_section(pdbfile, atomtypes, sopc, CA_rad, CBcom, CBradii)
            self.write_atoms_section(pdbfile,atomtypes,skip_glycine)
            btfile='interaction.dat'
            self.write_contacts_section(pdbfile,nativefile,btfile,cutoff,atomtypes,sopc,btparams,dswap)
            self.write_bonds_section(nativefile, Kb, atomtypes)
            self.write_angles_section(nativefile,Ka)
            self.write_dihedrals_section(nativefile,Kd)
            self.write_exclusions_section()
            contacttype=2;bond_type=1
            #assert contacttype==2
            #write gromacs format files.
            #SBM order: atoms, contacts, bonds, angles, dihedrals.
            #GROMACS order: atomtypes, moleculetype,atoms,bonds, angles, dihedrals,pairs,exclusions,system,
            #global w_sbm
            #w_sbm = False
            #print (w_sbm)
            self.set_false_global_wsbm()
            #print (w_sbm)
            self.write_gro_header(topfilename,atomtypes)
            self.write_gro_atomtypes(topfilename,atomtypes,pdbfile,sopc,CA_rad,CBcom,CBradii)
            self.write_gro_moleculetype(topfilename)
            self.write_gro_atoms(topfilename, atomtypes, pdbfile, nativefile, sopc)
            self.write_gro_bonds(topfilename, atomtypes, nativefile, Kb, bond_type,dsb)
            self.write_gro_angles(topfilename,atomtypes,nativefile,Ka)
            self.write_gro_dihedrals(topfilename,atomtypes,nativefile,Kd)
            self.write_gro_pairs(topfilename, atomtypes, nativefile, pdbfile, contacttype, cutoff, sopc, btparams)
#            global w_sbm,w_gro
#            w_sbm=True;w_gro=False
 #           self.write_atomtypes_section(pdbfile, atomtypes, sopc, CA_rad, CBcom, CBradii)


            #for dswap write original protein files again.
            #self.write_CB_to_native(pdbfile,sopc)
#            bond_ptype = 1
            #self.write_gro_atoms(topfilename,atomtypes,pdbfile,nativefile)
 #           if not sopc:
 #               self.write_gro_angles(topfilename,atomtypes,nativefile,Ka)
 #               self.write_gro_dihedrals(topfilename,atomtypes,nativefile,Kd)

            self.write_gro_tail(topfilename)
            self.write_gro_gro('gromacs.gro',pdbfile,atomtypes,skip_glycine)
            print ('See file:',topfilename,'for GROMACS topology.')
            return

        def get_gro_from_xyz(self, pdbfile,nativefile, xyzfile, grofilename):
                # convert .xyz to .gro
                Y=conmaps()
                traj = Y.read_xyz_as_traj(xyzfile, nativefile)
                xyz=traj.xyz
                atomnames=self.get_atom_names(pdbfile,1)
                residues=Y.get_sequence(pdbfile)
                #print (residues)
                l1=np.arange(0,len(residues))
                d=self.two_lists_to_dict(l1,residues)
                d1=self.amino_acid_dict()
                f = open(grofilename, "w+")
                f.write('%s\n' % ('Grofile generated from Go-kit'))
                f.write('%d\n' % (len(atomnames)))
                for i in xrange(0,len(residues)):
                    #print i,d1[d[i]],atomnames[i],xyz[0][i][1]
                    f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (
                    i+1,d1[d[i]],atomnames[i],i+1,xyz[0][i][0],xyz[0][i][1],xyz[0][i][2]))

                f.write('%8.4f %8.4f %8.4f' % (5.0000, 5.0000, 5.0000))
                f.close()
                return 1

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate GROMACS and OPTIM potential files for enhanced SBM models.")
    parser.add_argument("--CA_rad","-CA_rad", help="Radius for C-alpha atom. Default=4.0")
    parser.add_argument("--CAcom","-CAcom",action='store_true',help="Place C-alpha at COM of backbone")
    parser.add_argument("--Cb_rad","-CB_rad", help="Statistically derived values by default for attype 2.")
    parser.add_argument("--Kb","-Kb", help="Kbond")
    parser.add_argument("--Ka","-Ka", help="Kangle")
    parser.add_argument("--Kd","-Kd", help="Kdihedral")
    parser.add_argument("--cutoff","-cutoff", help="Cut-off for contact-map generation")
    parser.add_argument("--scaling","-scaling", help="Scaling for mapping to all-atom contact-map.")
    parser.add_argument("--attype", "-attype",help="Number of atom types. E.g. 1 for CA, 2 for CA and CB")
    parser.add_argument("--matrix","-matrix",action='store_true', default=False, help='User defined interactions in file interaction.dat.')
    parser.add_argument("--CBcom","-CBcom", action='store_true', default=False,help='Put CB at center of mass of side-chain (no hydrogens)')
    parser.add_argument("--CBfar", "-CBfar", action='store_true', help="Place C-beta on farthest non-hydrogen atom.")
    parser.add_argument("--dsb", "-dsb",action='store_true', help="Use desolvation barrier potential for contacts.")
    parser.add_argument("--native_ca","-native_ca", help='Native file with only C-alphas. Just grep pdb. ')
    parser.add_argument("--grotop","-grotop",help='Gromacs topology file output name.')
    parser.add_argument("--aa_pdb","-aa_pdb", help='all-atom pdbfile e.g. 1qys.pdb')
    parser.add_argument("--pdbgro","-pdbgro", help='Name for .gro file.')
    parser.add_argument("--w_native","-w_native", help='Write native files, CA and CA-CB from all atom PDB file.')
    parser.add_argument("--pl_map","-pl_map", action='store_true', default=False, help='Plot contact map for two bead model')
    parser.add_argument("--skip_glycine","-skip_glycine", action='store_true', default=False, help='Skip putting Cbeta on glycine')
    parser.add_argument("--dswap","-dswap", action='store_true', default=False, help='For domain swapping runs. Symmetrised SBM is generated.')
    parser.add_argument('--hphobic',"-hphobic",action='store_true',help='Generate hydrophobic contacts.')
    parser.add_argument('--hpdist', "-hpdist", help='Equilibrium distance for hydrophobic contacts.')
    parser.add_argument('--btmap',"-btmap", action='store_true', help='Use Betancourt-Thirumalai interaction matrix.')
    parser.add_argument('--mjmap',"-mjmap", action='store_true', help='Use Miyazawa-Jernighan interaction matrix.')
    parser.add_argument('--hpstrength',"-hpstrength",help='Strength with which hydrophobic contacts interact.')
    parser.add_argument('--ext_conmap',"-ext_conmap",help='External contact map in format chain res chain res')
    parser.add_argument('--CB_radii',"-CB_radii",action='store_true', help='External contact map in format chain res chain res')
    parser.add_argument('--CA_sep',"-CA_sep", help='Backbone separation. Default is >=4')
    parser.add_argument('--CB_sep',"-CB_sep", help='Side-chain separation. Default is >=2')
    parser.add_argument('--CAB_sep',"-CAB_sep", help='Backbone-sidechain separation. Default is >=2')
    parser.add_argument('--gauss',"-gauss", help="Gaussian interactions for contacts.")


    args = parser.parse_args()
    X = esbm()
    Y = conmaps()
    U = Utils()
    #U.print_preamble()
    #X.MJparams('miyazawa-jernighan.dat')
    #Y.write_hydrophobic_contacts('contacts.txt','1ris.pdb')
    # import json
    # d=X.amino_acid_radius_dict()
    # with open('file.txt', 'w') as file:
    #     file.write(json.dumps(d))


    #Test here
    #X.write_gro_gro('test.gro','1ake.pdb',1)
    #i_back_bone_COM('1ubqH.pdb')
    #X.get_farthest_atom_in_residue('1qys_H.pdb',91)
    #X.get_farthest_CB_coords('1qys_H.pdb')
    #X.read_external_radius('radii.dat')
    #exit()


    #Set default parameters
    sopc=True
    btparams=False
    mjmap=False
    btmap=False
    dswap=False
    skip_glycine=False
    Kb = 200
    Ka = 40
    Kd = 1
    cutoff = 4.5
    CA_rad = 4.0
    CAcom = False
    hphobic=False
    dsb=False
    CBfar=False
    CBcom=True
    casep=4;cbsep=2;cabsep=2
    scaling=1.2
    gauss=False
    if args.Kb:
        Kb=float(args.Kb)
    if args.Kd:
        Kd=float(args.Kd)
    if args.Ka:
        Ka=float(args.Ka)
    if args.cutoff:
        cutoff=float(args.cutoff)
    if args.matrix:
        U.file_exists('interaction.dat')
        btparams=True
    else:
        btparams=False
    if args.CB_radii:
        U.file_exists('radii.dat')
        CBradii=True
    else:
        CBradii=False
    if args.skip_glycine:
        skip_glycine=True
        #assert(args.attype==2)
    if args.CAcom:
        CAcom=True
    if args.dswap:
        dswap=True
        print ('This one assumes both chains are identical.')
    if args.CA_rad:
        CA_rad=float(args.CA_rad)
        print ("Setting CA_rad to ",CA_rad)
    if args.hphobic:
        hphobic=True
    if args.dsb:
        dsb=True
    if args.hpdist:
        hpdist=args.hpdist
    else:
        hpdist=5.5
    if args.hpstrength:
        hpstrength = args.hpstrength
    else:
        hpstrength=1
    #set global variables
    if args.CBfar:
        CBcom=False
        CBfar=True
    if args.btmap:
        assert (int(args.attype)==2),'Attype set to 1!'
        import shutil
        shutil.copy2('btmap.dat','interaction.dat')
        btparams=True;btmap=True
    if args.mjmap:
        assert (int(args.attype)==2),'Attype set to 1!'
        import shutil
        shutil.copy2('mjmap.dat', 'interaction.dat')
        btparams=True;mjmap=True
    if args.CA_sep:
        casep=int(args.CA_sep)
    if args.CB_sep:
        cbsep=int(args.CB_sep)
    if args.CAB_sep:
        cabsep=int(args.CAB_sep)
    if args.scaling:
        scaling=float(scaling)
    if args.gauss:
        gauss=True
    X.globals(Ka,Kb,Kd,CA_rad,skip_glycine,sopc,dswap,btparams,CAcom,hphobic,hpstrength,hpdist,dsb,mjmap,btmap,CBfar,casep,cbsep,cabsep,scaling,gauss)

    if not args.ext_conmap and args.aa_pdb:
        pdbfile=args.aa_pdb
        Y.all_atom_contacts(pdbfile,cutoff,scaling)

    if args.w_native:
        pdbfile=args.w_native
        Y.check_hetatm(pdbfile)
        if not skip_glycine:
            Y.check_glycineh(pdbfile,False)
        X.write_CB_to_native(pdbfile,False)
    if args.aa_pdb:
        pdbfile=args.aa_pdb
        X.write_CB_to_native(pdbfile,sopc)
    if args.attype and args.aa_pdb:
        topfile = 'gromacs.top'
        atomtypes = int(args.attype)
        if int(args.attype)==2:
            assert not hphobic, '--hphobic flag applicable only to single-bead models.'
            nativefile='native_cb.pdb'
            X.write_gromacs_top(topfile,int(atomtypes),pdbfile,nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii)
            # else:
            #     print ('need --CBcom argument.')
        if int(args.attype)==1:
            U=Utils()
            if U.file_exists('native_ca.pdb'):
                nativefile='native_ca.pdb'
                X.write_gromacs_top(topfile, int(atomtypes), pdbfile, nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii)
            else:
                print ('Need native_ca.pdb file. Just grep "CA" pdbfile.')
        if args.pdbgro:
            #nativefile='native_cb.pdb'
            grofile='gromacs.gro'
            X.write_gro_gro(grofile,pdbfile,atomtypes,sopc)
        if args.pl_map:
            l=np.loadtxt('contacts.txt')
            xi=l[:,1]
            yi=l[:,3]
            Y.plot_map(xi,yi,'conmap','Res2','Res1')
        U.make_dir('PATH');U.make_dir('MD')
        U.make_dir_struc('PATH','MD')
if __name__ == '__main__':
    main()