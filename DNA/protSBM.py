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
usage: python protsbm.py --options

Author: Sridhar Neelamraju)

"""


from __future__ import print_function
import numpy as np
import math as m
from gr import conmaps
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO,Select
import collections
import os
from util import Utils
from table import tables
#Some conventions:
#First atom is always a C-alpha. This helps find non-bonded atoms easily.

class Measure():
    def distances(self,file,pairs):
        #loading co-ordinates
        fin = open(file)
        xyz = list()
        for i in fin:
            if i.startswith("ATOM"):
                xyz.append(np.float_([i[30:38],i[38:46],i[46:54]]))
        xyz = np.array(xyz)
        fin.close()

        #calcularing bond lengths
        r = list()
        for i in pairs:
            r.append(sum((xyz[i[1]]-xyz[i[0]])**2)**0.5)
        
        #storing bond lengths along with bond indexes
        dist = list()
        dist.append(r)
        for i in pairs:
            dist.append(i)
        
        return dist
    def angles(self,file,triplets):
        #loading co-ordinates
        fin = open(file)
        xyz = list()
        for i in fin:
            if i.startswith("ATOM"):
                xyz.append(np.float_([i[30:38],i[38:46],i[46:54]]))
        xyz = np.array(xyz)
        fin.close()
        
        #calculating angles
        theta = list()
        for t in triplets:
            #angle xzy0,xyz1,xyz2
            AB = xyz[t[0]]-xyz[t[1]]; BC = xyz[t[2]]-xyz[t[1]]
            AB_BC = (sum(AB**2)**0.5)*(sum(BC**2)**0.5)
            theta.append(m.acos(sum(AB*BC)/AB_BC))  #   in rads
        
        #returning angles along with teir triplets
        angles = list()
        angles.append(theta)
        for i in triplets:
            angles.append(i)
        return angles
    def dihedrals(self,file,quads):
        #loading co-ordinates
        fin = open(file)
        xyz = list()
        for i in fin:
            if i.startswith("ATOM"):
                xyz.append(np.float_([i[30:38],i[38:46],i[46:54]]))
        xyz = np.array(xyz)
        fin.close()

        #calculating torsions
        phi = list()
        for q in quads:
            #vectors
            #let the four atoms be A,B,C,D
            BA = xyz[q[0]]-xyz[q[1]]; BC = xyz[q[2]]-xyz[q[1]]
            CD = xyz[q[3]]-xyz[q[2]]; CB = BC


            #Calculating cross product()
            #normals to plane 1 and 2
            #n1 = np.array([BA[1]*BC[2]-BA[2]*BC[1],BA[2]*BC[0]-BA[0]*BC[2],BA[0]*BC[1]-BA[1]*BC[0]])
            #n2 = np.array([CD[1]*CB[2]-CD[2]*CB[1],CD[2]*CB[0]-CD[0]*CB[2],CD[0]*CB[1]-CD[1]*CB[0]])
            n1 = np.cross(BA,BC)
            n2 = np.cross(CD,BC)

            #checking the sign (+/-) of the angle
                #Cos-1(x) usually has 2 solutions (+phi) and (-phi)
                #By default acos (from math lib) or arccos (from numpy) 
                #will return only +phi as they are calculated in the limits
                # if 0<phi<pi. This will give a mirrored dihedrals
                #for places where phi was actually negative
            #a vector prependicular to both n1 and m2 will be
            #parallel or anti-parallel to BC

            direction = np.cross(n1,n2)
            cosine_d = sum(direction*BC)/((sum(direction**2)**0.5)*(sum(BC**2)**0.5))

            #if cosine_d is -1, the n1Xn2 is antiparallel to BC and the dihedral is negative
            #if cosine_d is  1, the n1Xn2 is parallel to BC and the dihedral is positive         
            sign = int(np.round(cosine_d))

            #making dihedtals are negative, they lie b/w pi to 2*pi positive
            if sign == 1:
                angle_normalizer = 0
            elif sign == -1:
                angle_normalizer = 2*np.pi

            #calculating dihedral
            n1_n2 = (sum(n1**2)**0.5)*(sum(n2**2)**0.5)
            phi.append(angle_normalizer + sign*np.arccos(sum(n1*n2)/n1_n2))    #   in rads
            

        #returning dihedrals along with teir quadruplets
        dihedrals = list()
        dihedrals.append(phi)
        for i in quads:
            dihedrals.append(i)
        return dihedrals

class pdb_kit():
    def get_backbone_contacts(self,pdbfile,cutoff,scaling,separation):
        print(">>in backbone_contacts", pdbfile, "cutoff =",cutoff, "scaling=",scaling )
        cutoff = float(cutoff);scaling=float(scaling)
        #opening pdbfile
        fin = open(pdbfile)
        backbone = dict()
        count_linres = 0
        #loading coordinates
        for i in fin:
            if count_linres == 0:
                #for adjusting residue number if res_num doesn't start with 1
                starting_residue = int(i[22:26])
            if i.startswith("ATOM"):
                #creating a backbone list
                #adjusting residue number 
                resi_num = int(i[22:26])-starting_residue
                if resi_num not in backbone:
                    backbone[resi_num] = dict()
                    count_linres = count_linres + 1
                atomname = i[12:16].strip()
                #loading backbone coordinates hashed to residue number and atom name
                backbone[resi_num][atomname] = np.float_([i[30:38],i[38:46],i[46:54]])
        #print (count_linres,len(backbone))
        assert(count_linres==len(backbone))
        fin.close()

        contacts = list()

        #atoms in backbone
        in_backbone = ["CA","C","O","N"]
        count = 0
        #total = len(backbone)-separation
        for r1 in range(0,len(backbone)-separation):
            #print ("Completed",int(100*float(count)/total),"%")
            count = count + 1
            #selecting residue 1
            for r2 in range(r1+separation,len(backbone)):
                #selecting residue 2
                dset = list()           #setting list of distances
                contact_flag = False    #no contact found
                for atomname1 in in_backbone:
                    #loading coordinates of atom in the residue
                    a1 = backbone[r1][atomname1]
                    for atomname2 in in_backbone:
                        #loading coordinates of atoms in residue2
                        a2 = backbone[r2][atomname2]
                        d = np.sum((a2-a1)**2)**0.5
                        dset.append(d)
                        if d<=cutoff:
                            #print (r1,atomname1,a1,r2,atomname2,a2,d)
                            contact_flag = True
                if contact_flag:
                    dmin = np.amin(np.array(dset))
                    CA1 = backbone[r1]["CA"]
                    CA2 = backbone[r2]["CA"]
                    ca_dist = (np.sum((CA2-CA1)**2)**0.5)
                    if dmin < scaling*ca_dist:
                        contacts.append([r1,r2,dmin,ca_dist])
        contacts.sort()
        contacts=np.asarray(contacts)
        np.savetxt('bb_new.t', contacts)
        return contacts
    def get_side_chain_contacts(self,pdbfile,cutoff,scaling,separation):
        print(">> in side_chain_contacts", pdbfile,"cutoff =",cutoff,"scaling=",scaling)
        cutoff = float(cutoff);scaling=float(scaling)
        fin = open(pdbfile)
        sidechain = dict()
        sidechain_CB = dict()
        count_linres = 0
        for i in fin:
            if count_linres == 0:
                starting_residue = int(i[22:26])
            if i.startswith("ATOM"):
                #creating a sidechain list
                resi_num = int(i[22:26])-starting_residue
                if resi_num not in sidechain:
                    sidechain[resi_num] = list()
                    count_linres = count_linres + 1
                atomname = i[12:16].strip()
                if atomname not in ["CA","C","O","N"]:
                    sidechain[resi_num].append(np.float_([i[30:38],i[38:46],i[46:54]]))
                    #if atomname == "CB":
                    #sidechain_CB[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
        assert(count_linres==len(sidechain))
        fin.close()

        fin = open("native_cb.pdb")
        for i in fin:
            resi_num = int(i[22:26])-starting_residue
            atomname = i[12:16].strip()
            if atomname == "CB":
                sidechain_CB[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
        fin.close()


        contacts = list()
        count = 0
        #total = len(sidechain)-separation
        for r1 in range(0,len(sidechain)-separation):
            #print ("Completed",int(100*float(count)/total),"%")
            count = count + 1
            #selecting residue 1
            for r2 in range(r1+separation,len(sidechain)):
                #selecting residue 2
                dset = list()           #setting list of distances
                contact_flag = False    #no contact found
                for a1 in sidechain[r1]:
                    for a2 in sidechain[r2]:
                        d =  (np.sum((a2-a1)**2)**0.5)
                        dset.append(d)
                        if d<=cutoff:
                            contact_flag = True
                
                if contact_flag:
                    dmin = np.amin(np.array(dset))
                    CB1 = sidechain_CB[r1]
                    CB2 = sidechain_CB[r2]
                    cb_dist = (np.sum((CB2-CB1)**2)**0.5)
                    if dmin < scaling*cb_dist:
                        contacts.append([r1,r2,dmin,cb_dist])
        contacts.sort()
        contacts=np.asarray(contacts)
        np.savetxt('sidechain_new.t', contacts)
        return contacts
    def get_bb_sc_contacts(self,pdbfile,cutoff,scaling,separation):
        print(">> in bb_sc_contacts", pdbfile, "cutoff =", cutoff, "scaling =",scaling)
        cutoff = float(cutoff);scaling=float(scaling);separation=int(separation)

        fin = open(pdbfile)
        #input pdbfile
        backbone = dict(); sidechain = dict()
        #for native Coarse grain file
        backbone_CA = dict(); sidechain_CB = dict()
        
        #loading all-arom pdbfile
        count_linres = 0
        for i in fin:
            if count_linres == 0:
                starting_residue = int(i[22:26])
            if i.startswith("ATOM"):
                #creating a backbone list
                resi_num = int(i[22:26])-starting_residue
                if resi_num not in backbone:
                    backbone[resi_num] = list()
                    sidechain[resi_num] = list()
                    count_linres = count_linres + 1
                atomname = i[12:16].strip()
                if atomname in ["CA","C","O","N"]:
                    backbone[resi_num].append(np.float_([i[30:38],i[38:46],i[46:54]]))
                    #if atomname == "CA":
                    #   backbone_CA[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
                else:
                    sidechain[resi_num].append(np.float_([i[30:38],i[38:46],i[46:54]]))
                    #if atomname == "CB":
                    #   sidechain_CB[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
        fin.close()
        #loading native co-ordinates
        fin = open("native_cb.pdb")
        for i in fin:
            resi_num = int(i[22:26])-starting_residue
            atomname = i[12:16].strip()
            if atomname == "CB":
                sidechain_CB[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
            elif atomname == "CA":
                backbone_CA[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
        
        assert(count_linres==len(backbone))
        assert(count_linres==len(sidechain))
        fin.close()
        
        contacts = list()
        count = 0
        #total = len(sidechain)
        for r1 in range(0,len(backbone)):
            count = count + 1
            #selecting residue 1
            for r2 in range(0,len(sidechain)):
                #selecting residue 2
                if r2-r1 >= separation or r1-r2 >= separation:
                    dset = list()           #setting list of distances
                    contact_flag = False    #no contact found
                    for a1 in backbone[r1]:
                        for a2 in sidechain[r2]:
                            d =  (np.sum((a2-a1)**2)**0.5)
                            dset.append(d)
                            if d<=cutoff:
                                contact_flag = True
                
                    if contact_flag:
                        dmin = np.amin(np.array(dset))
                        CA = backbone_CA[r1]
                        CB = sidechain_CB[r2]
                        ca_cb_dist = (np.sum((CB-CA)**2)**0.5)
                        if dmin < scaling*ca_cb_dist:
                            contacts.append([r1,r2,dmin,ca_cb_dist])
        contacts.sort()
        contacts=np.asarray(contacts)
        np.savetxt('bb_sc_new.t', contacts)
        return contacts
    def get_bb_contacts_SOPC(self,pdbfile,nativefile,cutoff,scaling,separation):
        #SOP-SC applies cut-off criterion to coarse-grained beads. Not to the all-atom PDB.
        #Slightly different from Cheung-Thirumalai handling.
        #NATIVE INTERACTIONS.
        print(">> in bb_contacts (SOP-SC)", nativefile, "cutoff =", cutoff)
        cutoff = float(cutoff); scaling = float(scaling)
        contacts = []
        #Matches paper results when only CA beads are taken into account. As done below!
        #Pairs of residues
        backbone_CA = list()
        fin = open(nativefile)
        for i in fin:
            if i.startswith("ATOM"):
                if i[12:16].strip() == "CA":
                    backbone_CA.append(np.float_([i[30:38],i[38:46],i[46:54]]))
        nca = len(backbone_CA)
        #Paira of atoms
        count=0
        for i in range(0,nca-separation):
            for j in range(i+separation,nca):
                    #CA-CA
                    x1 = backbone_CA[i]; y1 = backbone_CA[j]
                    distance = np.sum((y1-x1)**2)**0.5
                    # cutoff and scaling criterion.
                    if (distance <= cutoff):
                        #print (x1,y1,distance*10)
                        contacts.append([i, j,distance, distance])
                        count=count+1
        print ("Found",count,"bb-bb contacts");
        #np.savetxt ('bb_sopsc.dat',contacts,fmt='%-4d  %-4d  %10.6f  %10.6f')
        return contacts
    def get_ss_contacts_SOPC(self,pdbfile, nativefile, cutoff, scaling, separation):
        cutoff = float(cutoff);scaling = float(scaling)
        fin = open(nativefile)
        side_CB = dict()
        contacts = list()
        count = 0
        ncb = 0
        found_GLY_without_CB = False
        for i in fin:
            if i.startswith("ATOM"):
                if i[17:20] != "GLY":
                    if found_GLY_without_CB:
                        print ("Missing GLY CB!!!. Check your input")
                        found_GLY_without_CB = False
                    if i[12:16].strip() == "CB":
                        side_CB[ncb] = np.float_([i[30:38],i[38:46],i[46:54]])
                        ncb += 1
                if i[17:20] == "GLY":
                    if i[12:16].strip() == "CA":
                        side_CB[ncb] = np.float_([1000000,10000000,100000000])
                        found_GLY_without_CB = True
                        ncb += 1
                    if i[12:16].strip() == "CB":
                        side_CB[ncb-1] = np.float_([i[30:38],i[38:46],i[46:54]])
                        found_GLY_without_CB = False
                        
        for i in range(0,ncb-separation):
            for j in range(i+separation,ncb):
                x1 = side_CB[i] ;y1 = side_CB[j]
                distance = np.sum((y1-x1)**2)**0.5
                if distance < 10**(-10):
                    #if distance is 0
                    #using buffer gly co-ordinates skip
                    continue
                if (distance <= 8.0):
                    contacts.append([i, j, distance,  distance ])
                    count = count + 1
        print("Found", count, "sc-sc contacts")
        np.savetxt('scsc_sopsc.dat', contacts, fmt='%-4d  %-4d  %10.6f  %10.6f')
        return contacts
    def get_bs_contacts_SOPC(self,pdbfile, nativefile, cutoff, scaling, separation):
        cutoff = float(cutoff); scaling = float(scaling)
        fin = open(nativefile)
        bb = dict()
        sc = dict()
        nca = 0
        ncb = 0
        found_gly_without_cb = False
        for i in fin:
            if i.startswith("ATOM"):
                if i[12:16].strip() == "CA":
                    if found_gly_without_cb:
                        #Directly hit next CA after GLY CA
                        #Found gly in nativre without CB...Adding dummy co-ordinates
                        sc[ncb] = np.float_([1000000,100000,100000])
                        print ("Missing GLY CB!!!. Check your input")
                        ncb += 1
                        found_gly_without_cb = False
                    bb[nca] = np.float_([i[30:38],i[38:46],i[46:54]])
                    nca += 1
                    if i[17:20] == "GLY":
                        found_gly_without_cb = True
                elif i[12:16].strip() == "CB":
                    sc[ncb] = np.float_([i[30:38],i[38:46],i[46:54]])
                    ncb += 1
                    if i[17:20] == "GLY":
                        found_gly_without_cb = False
        count = 0
        contacts = list()
        #checking length of CA and CB list should be same
        assert(len(sc)==len(bb))
        assert(ncb==len(sc))
        assert(ncb==nca)
        #######_______TESTING________#
        for i in range(0,nca-separation):
            for j in range(i+separation,nca):
                #if np.absolute(j-i) < separation:
                #   continue
                #here nca =  ncb
                x1 = sc[i]; x2 = bb[i]
                y1 = bb[j]; y2 = sc[j]
                # CB-CA
                distance = np.sum((y1-x1)**2)**0.5
                if (distance <= 8.0):
                    contacts.append([i, j, distance, distance])
                    count += 1
                # x2 and y2:
                #
                ##CA-CB
                #distance = np.sum((y2-x2)**2)**0.5
                #cutoff and scaling criterion.
                #if (distance <= 8.0):
                #   contacts.append([i, j, distance * 10, distance * 10])
                #   count += 1
        print("Found", count, "bb-sc contacts. See SOP-SC.txt")
        np.savetxt('new_bb_sc_sopc.dat', contacts,fmt='%-4d  %-4d  %10.6f  %10.6f')
        contacts = np.asarray(contacts)
        return contacts
    def get_SOPC_non_native_angles(self,nativefile,CA_rad,epsl,CG_indices):
        print ('>> in get_sopc_non_native_angles')
        #SOP-SC model has non-native non-bonded interactions(angular repulsions!).
        #angles formed by triplets: Non-bonded part of triplet is written in pairs section.
        sigmabb=float(CA_rad)
        CA_indices = CG_indices[0]
        CB_indices = CG_indices[1]
        #select only backbone atoms
        contacts=[]
        count=0
        for i in CA_indices[:-2]:
            count=count+1
            contacts.append([i, i+4,sigmabb, sigmabb ])
        print ('Found',count,'non-native,non-bonded bb-bb angular interactions')
        #select CB atoms
        count=0
        #get side-chain radius.
        U=Utils(); sc_rad=U.amino_acid_radius_dict()

        CB_indices.pop(0) #remove first element from list.
        fin = open(nativefile)
        sc_residue = dict()
        for i in fin:
            if i.startswith("ATOM"):
                if int(i[6:11])-1 in CB_indices:
                    sc_residue[int(i[6:11])-1] = i[17:20].strip()
        fin.close()
        for i in CB_indices[:-1] :
            #assign radii to each side-chain
            name = sc_residue[i]
            sigmasc=float(sc_rad[name])
            sigmabs=0.4*(sigmasc+sigmabb)
            #print (sigmabs,sigmabb,sigmasc)
            count=count+2
            contacts.append([i,i+1,sigmabs,sigmabs]) #contacttype 1=LJ for GROMACS and SBM.#FTW
            contacts.append([i, i-3,sigmabs,sigmabs])

        print('Found', count, 'non-native,non-bonded bb-sc angular interactions')
        #np.savetxt('nn_angles_sopsc.dat', contacts, fmt='%-4d  %-4d  %10.6f  %10.6f')
        #must retturn residue numbers, and not atom numbers!
        return contacts

class protsbm(Select,object):
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
    # For help type python protSBM.py --help
        def __init__(self):
            return
        def globals(self,Ka1,Kb1,Kd1,CA_rad1,skip_glycine1,sopc1,dswap1,btparams1,CAcom1,hphobic1,hpstrength1,hpdist1,dsb1,mjmap1,btmap1,CBfar1,CBcharge1,):
            #Give default values to global variables.
            global Ka,Kb,Kd,CA_rad,hpdist,hpstrength
            Ka=Ka1;Kb=Kb1;Kd=Kd1;CA_rad=CA_rad1;hpdist=float(hpdist1);hpstrength=float(hpstrength1)
            global skip_glycine,sopc,dswap,btparams,CAcom,hphobic,dsb,btmap,mjmap,CBfar,CBcharge
            skip_glycine=skip_glycine1;dswap=dswap1;btparams=btparams1;sopc=sopc1;CAcom=CAcom1;hphobic=hphobic1;dsb=dsb1;CBfar=CBfar1;CBcharge=CBcharge1
            mjmap=mjmap1;btmap=btmap1
            #print ('Global logicals are as follows:\n--------------------------')
            #print (' skip_glycine=',skip_glycine,'\n','dswap=',dswap,'\n','bt=',btparams,'\n','sopc=',sopc,'\n','CAcom=',CAcom,'\n','hphobic=',hphobic)
            global radtodeg,kcaltokj
            radtodeg = 57.295779513;kcaltokj=4.184
            #print (dsb)
            return
############____________shift in separate tools_____________________#
        def two_lists_to_dict(self,keys,values):
            return dict(zip(keys, values))
        def get_dict_from_csv(self,csvfile):
            #read parameters as table and save as list
            import csv
            with open(csvfile, 'r') as f:
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
            inv_d = {v: k for k, v in d.items()}
            return inv_d
        def amino_acid_dict2(self):
            d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
            return d
        def accept_residue(self,residue):
            #Removing non amino acid residues, if any
            global aa_res 
            aa_res = self.amino_acid_dict2()
            if residue.get_resname().lstrip().rstrip() in aa_res:
                return 1
            else:
                return 0
        def amino_acid_radius_dict(self):
            #Radius of side-chains for two-bead model
            #Change here for custom radii.
            #http://www.pnas.org/cgi/doi/10.1073/pnas.1019500108. Table S2, SI.
             d = {'CYS': 2.74, 'ASP': 2.79, 'SER': 2.59, 'GLN': 3.01, 'LYS': 3.18,
            'ILE': 3.09, 'PRO': 2.78, 'THR': 2.81, 'PHE': 3.18, 'ASN': 2.84,
            'GLY': 2.25, 'HIS': 3.04, 'LEU': 3.09, 'ARG': 3.28, 'TRP': 3.39,
            'ALA': 2.52, 'VAL': 2.93, 'GLU': 2.96, 'TYR': 3.23, 'MET': 3.09}
             return d
        def amino_acid_charge_dict(self):
            d = {'CYS': 0, 'ASP': -1, 'SER': 0, 'GLN': 0, 'LYS': 1,
            'ILE': 0, 'PRO': 0, 'THR': 0, 'PHE': 0, 'ASN': 0,
            'GLY': 0, 'HIS': 1, 'LEU': 0, 'ARG': 1, 'TRP': 0,
            'ALA': 0, 'VAL': 0, 'GLU': -1, 'TYR': 0, 'MET': 0}
            return d
        def amino_acid_charge_dict2(self):
            d = {'C': 0, 'D': -1, 'S': 0, 'Q': 0, 'K': 1,
            'I': 0, 'P': 0, 'T': 0, 'F': 0, 'N': 0,
            'G': 0, 'H': 1, 'L': 0, 'R': 1, 'W': 0,
            'A': 0, 'V': 0, 'E': -1, 'Y': 0, 'M': 0}
            return d
        def get_aa_sequence(self,pdbfile):
            fin = open(pdbfile)
            code = self.amino_acid_dict2(); seq = ""
            prev_resi = 0
            for i in fin:
                if i.startswith("ATOM"):
                    new_resi = int(i[22:26])
                    if new_resi != prev_resi:
                        seq = seq + code[i[17:20].strip()]
                    prev_resi = new_resi
            return seq
        def split_chains(self,pdbfile):
            U=Utils()
            U.split_chain(pdbfile)
        def potential_dict(self,potential):
            #three potentials only implemented for now.
            d = {'lj612': 1, 'lj1012': 2,'dsb': 8}  
            return
        def get_back_bone_COM(self,pdbfile):
            print (">>> Calculating amino-acid backbone center of mass co-ordinates.")
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
                                #print atom._assign_atom_mass()
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
                        #print COM
                        COM1 += COM
                        count = count + 1
                COM1 = np.array(COM1).reshape(len(COM1) / 2, 2)
            return COM1
        def get_side_chain_COM(self,pdbfile,skip_glycine):
            print (">>> Calculating amino-acid side chain center of mass co-ordinates.")
            #returns COM of each side-chain as numppy array
            #SOPC puts CB on glycine hydrogen!!!!

            #assert (sopc!=skip_glycine)
            if not skip_glycine:
                glyname= 'GLY111'
                print ("Warning: Adding CB to Glycines as well!")
            else:
                glyname='GLY'
            # res numbering starts from 1
            
            #opeing structure
            p = PDBParser(PERMISSIVE=0)
            structure = p.get_structure('test', pdbfile)
            backbone=['N','CA','C','O']
            COM1=[]
            count = 0
            for model in structure:
                for chain in model:
                    for residue in chain:
                        coord = []
                        mass = []
                        #storing native residue number in count
                        count = residue.id[1]
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
                        if res==glyname:
                            print ("Excluding glycine, res = ",count+1)
                            #print atom_mass
                        else:
                            if len(coord) == 0:
                                print ("No co-ordinates found. Make sure to add Hydrogens if not skipping Glycines")
                                exit()
                            assert len(coord)==len(mass)
                            mass_sum=sum(mass)
                            #print mass
                            COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]),coord)],[count]]
                            #print COM
                            COM1+=COM
                        #count = count + 1
                        #print ("testing_count",count)
                COM1=np.array(COM1).reshape(len(COM1)/2,2)
                #coordinates,resnum
                # print (COM1)
            #print (len(COM1))
            return COM1
        def get_farthest_atom_in_residue_list(self,pdbfile,index_list):
            #returns coordinates of atom farthest from CA in residue.
            farthest_sc_coord = list()
            P = PDBParser()
            
            #loading list of residue in the input structure
            residues_list = P.get_structure('test',pdbfile).get_residues()
            for residue in residues_list:
                dist = 0
                ca = residue["CA"]      #loading CA atom object
                farthest_atom = ca      #intitializing as the farthest atom
                for atom in residue:    #looping ovar all sidechain atoms
                    if atom.get_name() not in ["C","O","N"] and not atom.get_name().startswith("H"):
                        if atom-ca > dist:
                            dist = atom-ca
                            farthest_atom = atom
                #if no atom found, CA will be saved as the farthest atom with dist = 0
                farthest_sc_coord.append((farthest_atom.get_name(),farthest_atom.get_coord()))
            return farthest_sc_coord
        def get_farthest_CB_coords(self,pdbfile,seq):
            print(">> Finding farthest amino-acid sidechain heavy atom co-ordinates")

            #number of residues
            n_res=len(seq)
            
            #reading residue number
            pdbstruc = PDBParser().get_structure('test',pdbfile)
            pdbresi = pdbstruc.get_residues()
            resi_num = [x.id[1] for x in pdbresi]

            #generating co-ordinate list
            CB_far_list = self.get_farthest_atom_in_residue_list(pdbfile,range(n_res))

            #filtering co-ordinates for GLY
            CB_list=[]
            for i in range(len(CB_far_list)):
                if CB_far_list[i][0] == "CA":
                    print ("Skipping GLY residue",resi_num[i])
                    continue
                else:
                    CB_list.append([ CB_far_list[i][1],resi_num[i]])
            
            #returning array
            return np.array(CB_list)
        def write_CB_to_native(self,pdbfile,sopc,atomtypes):
            print (">>> Writing Protein native files in native_ca.pdb and native_cb.pdb.\n")

            #protein sequence (chains appended)
            seq = self.get_aa_sequence(pdbfile)

            if CAcom:
                CA=self.get_back_bone_COM(pdbfile)      
                CA_coords = CA[:, 0].tolist()   #converting array to list
            
            #CB dictionary. Will remain empty if atomtyoes = 1
            CB_dict = dict()
            if atomtypes == 2:
                #checking if glycine is to be skipped
                if not skip_glycine:
                    G_list = list()
                    for i in range(len(seq)):
                        if seq[i] == "G":
                            G_list.append(i)
                    print ("Use Glycine set TRUE. Gly in residues: ",len(G_list))
                else:
                    sopc=False

                #if CB postion if extreme/farthest side chain atom
                if CBfar:
                    if not skip_glycine:
                        print ("CBFar option is not available with skip_glycine False ")
                        exit()
                    CB=self.get_farthest_CB_coords(pdbfile,seq)
                else:
                    #get side-chain center of mass
                    CB=self.get_side_chain_COM(pdbfile,skip_glycine)
                
                CB_coords = CB[:, 0].tolist()   #converting array to list

                if CAcom and not skip_glycine:
                    assert len(CB)==len(CA)
                res_num = CB[:, 1].tolist()     #converting array to list

                #if CA and CB are to be at the same position as the respective atom
                #in the all-atom file....they will be directly written from the all atom pdbfile itself.

                #converting list into dict()
                CB_dict = self.two_lists_to_dict(res_num, CB_coords)
                
                #opening CB native file
                f1 = open('native_cb.pdb', "w+")
            
            #opening CA native file
            f2 = open('native_ca.pdb', "w+")
            #Put each CB under the corresponding CA residue
            #Read CA from pdbfile

            #freading input all-atom pdbfile
            line = open(pdbfile)

            #counters
            count_CA=0
            count_CB=0
            count=0

            k = {}
            j = {}
            for i in line:
                if i.startswith("TER"):
                    #writing terminal
                    if atomtypes == 2:
                        f1.write("TER".ljust(6)+" "*(16)+k[5]+"\n")
                if i[12:16].strip()=="CA":
                    count=count+1
                    k[0] = "ATOM".ljust(6)  #                   atom#6s
                    k[1] = str(count).rjust(5)                  #atom num#5d
                    k[2] = 'CA'.center(4)                       # atomname$#4s
                    k[3] = i[17:20].strip().ljust(3)  # resname#1s
                    k[4] = i[21].rjust(1)                       # chain_id
                    k[5] = i[22:26].strip().rjust(4)  # resnum
                    if CAcom:   #CA COM co-ordinates
                        k[6] = str('%8.3f' % (float(CA_coords[count_CA][0]))).rjust(8)  # x
                        k[7] = str('%8.3f' % (float(CA_coords[count_CA][1]))).rjust(8)  # y
                        k[8] = str('%8.3f' % (float(CA_coords[count_CA][2]))).rjust(8)  # z
                    else:       #CA position in all-atom file
                        k[6] = str('%8.3f' % (float(i[30:38].rstrip().lstrip()))).rjust(8)  # x
                        k[7] = str('%8.3f' % (float(i[38:46].rstrip().lstrip()))).rjust(8)  # y
                        k[8] = str('%8.3f' % (float(i[46:54].rstrip().lstrip()))).rjust(8)  # z
                    k[9] = str('%6.2f' % (float(i[54:60].rstrip().lstrip()))).rjust(6)  # occ
                    k[10] = str('%6.2f' % (float(i[60:66].rstrip().lstrip()))).ljust(6)  # temp
                    k[11] = i[66:78].rjust(12)  # elname
                    #print("%s%s %s %s %s%s    %s%s%s%s%s%s" % (
                    #k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10], k[11]))
                    f2.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n" % (
                        k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10], k[11]))

                    #iterating CA count
                    count_CA = count_CA + 1
                    if atomtypes == 2:
                        #first riting CA in native CA+CB file (native_cb.pdb)
                        f1.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n" % (
                            k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10], k[11]))

                    #if the residue exits for CB
                    if int(i[22:26]) in CB_dict.keys():
                        #writing CB
                        #print (int(j[5]))
                        count=count+1
                        j[0] = "ATOM".ljust(6)#atom#6s
                        #j[1] = j[1].rjust(5)#aomnum#5d
                        j[1] = str(count).rjust(5)
                        j[2] = 'CB'.center(4)#atomname$#4s
                        j[3] = i[17:20].ljust(3)#resname#1s
                        j[4] = i[21].rjust(1) #Astring
                        j[5] = i[22:26].rstrip().lstrip().rjust(4) #resnum
                        #print (str('%8.3f' % (float(CB_coords[count_CB][0]))).rjust(8))
                        #print (int(i[22:26]),CB_dict[int(i[22:26])])
                        j[6] = str('%8.3f' % (float(CB_dict[int(i[22:26])][0]))).rjust(8)  # x                 
                        j[7] = str('%8.3f' % (float(CB_dict[int(i[22:26])][1]))).rjust(8)  # y
                        j[8] = str('%8.3f' % (float(CB_dict[int(i[22:26])][2]))).rjust(8)  # z
                        j[9] = str('%6.2f' % (float(i[54:60].rstrip().lstrip()))).rjust(6)  # occ
                        j[10] = str('%6.2f' % (float(i[60:66].rstrip().lstrip()))).ljust(6)  # temp
                        j[11] = i[66:78].rjust(12)  # elname
                        #uterating CB_count
                        count_CB=count_CB+1
                        f1.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n"%(j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10],j[11]))
                        #print("%s%s %s %s %s%s    %s%s%s%s%s%s"%(j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10],j[11]))
                        #print j[0]+"    "+j[1]+"  "+j[2]+"  "+j[3]+" "+j[4]+"  "+j[5]+"       "+j[6]+"  "+j[7]
                    else:
                        continue
            line.close()
            if atomtypes == 2:
                print ('See file: native_cb.pdb with ',count_CB+count_CA,'atoms')
                print ('CB=',count_CB,'CA=', count_CA)
                f1.close()
            print ('See file: native_ca.pdb with ',count_CA,'atoms')
            print ('CB=',count_CB,'CA=', count_CA)
            f2.close()
            return 1
        def get_CB_radius(self,atomtypes,nativefile,sopc,scaling,CBradii,seq):
            #excluded volume for each cb bead is distance from its C-alpha.
            #sopc gets side-chain radii from BT params!
            scaling=float(scaling)
            aa3_to_aa1 = self.amino_acid_dict()
            radius=[]
            if skip_glycine:
                sopc=False
            if CBradii:
                #using user defined redius
                print ('>>reading radii from file:','radii.dat')
                d=self.get_dict_from_csv('radii.dat')
                assert (len(d)==20)
                #while using the radius GLY will be skipped on it's own
                for i in seq:
                    radius.append(d[aa3_to_aa1[i]])
                return np.float_(radius)
            else:
                #using default radius
                d = self.amino_acid_radius_dict()
                #from greddy JPCB Supplementary.(https://pubs.acs.org/doi/suppl/10.1021/acs.jpcb.6b13100/suppl_file/jp6b13100_si_001.pdf)
                assert (atomtypes==2),'Atomtypes not set to two.'
                for i in seq:
                    radius.append(d[aa3_to_aa1[i]])
                radius=(0.8*(np.float_(radius)+3.8))/2
                return radius
        def get_atom_types(self,pdbfile,atomtypes,seq,CB_alltypes):
            d = self.amino_acid_dict()
            glyname='G'
            gly = [pos + 1 for pos, char in enumerate(seq) if char == 'G']
            ncb = len(seq) - len(gly)
            if not skip_glycine:
                glyname='G111'
                ncb=len(seq)
            natoms = 1 + ncb #one for CA

            count_cb=0
            atomname=[]
            atomname.append('CA')
            for i in seq:
                if i!=glyname:
                    if CB_alltypes:
                        x = 'CB' + str(count_cb)
                    else:
                        x = 'CB' + i
                    count_cb=count_cb+1
                    if x not in atomname:
                        atomname.append(x)

            #make dict
            if CB_alltypes:
                l1=np.arange(0,natoms)
            else:
                l1=np.arange(0,len(atomname))
            #print (len(l1),len(atomname))
            assert len(l1)==len(atomname)
            d=self.two_lists_to_dict(atomname,l1)
            return collections.OrderedDict(sorted(d.items()))
        def get_atom_names(self,pdbfile,atomtype,seq,CB_alltypes):
            d = self.amino_acid_dict()
            glyname='G'
            if atomtype==1:
                atomname=[]
                for i in seq:
                    #append CA for every residue
                    atomname.append('CA')
                return atomname
            elif atomtype==2:
                if not skip_glycine:
                    glyname='G111'      #G111 is just used as a tag
                #won't count glycine if the tag is not G
                gly = [pos + 1 for pos, char in enumerate(seq) if char == glyname]
                nca = len(seq)
                ncb = nca - len(gly)
                natoms = nca + ncb
                
                count_cb=0
                atomname=[]
                for i in seq:
                    atomname.append('CA')
                    if i!=glyname:  #glynam can be G or G111. G111 cnnot be there in the sequence hnce Gly is included
                        if CB_alltypes:
                            x = 'CB' + str(count_cb)
                        else:
                            x = 'CB' + i
                        count_cb=count_cb+1
                        atomname.append(x)
                assert len(atomname)==natoms
                return atomname
        def get_CG_index(self,nativefile,atomype):
                fin = open(nativefile)
                CA_indices = list(); CB_indices = list()
                AA_indeices = list()
                for i in fin:
                    if i.startswith("ATOM"):
                        if i[12:16].strip() == "CA":
                            CA_indices.append(int(i[6:11])-1)
                        elif i[12:16].strip() == "CB":
                            CB_indices.append(int(i[6:11])-1)
                        AA_indeices.append(int(i[6:11])-1)
                fin.close()
                assert (len(AA_indeices)==len(CA_indices)+len(CB_indices))
                return (AA_indeices,CA_indices,CB_indices)
        def write_header_SBM(self):
            f = open('SBM.INP', "w+")
            f.write('%s\n' % ('SBM topology file for use with OPTIM/PATHSAMPLE generated from Go-kit repo.'))
            f.write('%s\n' % ('Debye-Huckel Parameters: PREFACTOR, Dielectric Constant, Monovalent Ion Concentration, DH switching distance, DH cutoff distance'))
            f.write('%s\n' % (' 332.000   80.000    0.100   12.000   15.000'))

            f.write('%s\n' % ('nonbonded switching distance, truncate distance'))
            f.write('%s\n' % ('7.999999 8.999999'))
            f.close()
        def write_atomtypes_section(self, pdbfile, atomtype,sopc,CA_rad,CBcom,CBradii):
            CA_rad=float(CA_rad)        
            print ('>> in write_atomtypes_section',pdbfile,atomtype,sopc,CA_rad)
            assert atomtype<=2

            f = open('SBM.INP', "a")
            if skip_glycine:
                sopc=False
            
            seq = self.get_aa_sequence(pdbfile)
            d=self.amino_acid_dict()

            nca=len(seq)
            glyname='G'
            if atomtype==2:
                gly = [pos + 1 for pos, char in enumerate(seq) if char == 'G']
                ncb = nca - len(gly)
                if not skip_glycine:
                    glyname='G111'
                    ncb = nca 
                natoms = nca + ncb
                #print (nca,ncb)
            elif atomtype==1:
                natoms=len(seq)
            
            atomname=[]
            #getting atom mass
            atommass=np.ones(natoms,dtype=np.float32)
            
            #getting atom charge
            charge_in_atomtypes_section = False            
            if CBcharge and atomtype==2 and charge_in_atomtypes_section: 
                Q = self.amino_acid_charge_dict2()
                atomcharge = list()
                atomcharge.append(0.0)#for CA  
                for i in seq:
                    atomcharge.append(Q[i])
                atomcharge = np.float_(atomcharge)
            else:
                atomcharge=np.zeros(natoms,dtype=np.float32)
            
            #atomptype is set to A
            atomptype=['A']*natoms

            #CA_rad = np.ones(nca, dtype=np.float32) * CA_rad
            CA_rad = float(CA_rad)

            c6=0.000

            natomtypes = 21
            f.write('%d %s\n' % (natomtypes,' atomtypes'))
            f.write('%d\t%8.5f\t%s\n' % (1, CA_rad, '1.00'))
            

            #appending CA atomtype
            atomname.append(['CA', atommass[0], atomcharge[0], atomptype[0], c6, CA_rad])
            if atomtype==1:
                return atomname #only one atom type
            elif atomtype==2:
                count = 0
                count_cb = 0
                CB_rad = self.get_CB_radius(atomtype, 'native_cb.pdb', sopc, 1, CBradii,seq)
                count_gly = 0
                for i in seq:
                    if i==glyname:
                        count_gly = count_gly + 1
                CB_alltypes = False         ###_______TESTING___________#
                if CB_alltypes:
                    for i in seq:
                        #print (i)
                        if i != glyname:
                            x = 'CB' + str(count_cb) #str(residues[count].id[1]).lstrip().rstrip()
                            #print (atomcharge[count], atomptype[count], c6, CB_rad[0][count_cb], d[i])
                            #print (i,len(seq),count_cb,len(CB_rad[0]))
                            #print (x, atommass[count+1], atomcharge[count+1], atomptype[count], c6, CB_rad[0][count_cb], d[i])
                            atomname.append([x,atommass[count+1],atomcharge[count+1],atomptype[count],c6,CB_rad[count_cb],d[i]])
                            f.write('%d\t%8.5f\t%s\n' % (count_cb+2, CB_rad[count], '1.00'))
                        count = count + 1
                else:
                    cbtype_list = list()
                    for i in seq:
                        if i != glyname:
                            x = 'CB' + i        #CBK,CBA.CBP etc
                            if x not in cbtype_list:
                                cbtype_list.append(x)
                                atomname.append([x,atommass[count+1],atomcharge[count+1],atomptype[count],c6,CB_rad[count],d[i]])
                        count = count + 1
            return atomname
        def write_atoms_section(self,pdbfile,atomtype,sopc):
            print ('>> in write_atoms_section\t',atomtype,pdbfile,atomtype,sopc)
            assert atomtype<=2
            

            seq = self.get_aa_sequence(pdbfile)
            d=self.amino_acid_dict()
            
            f = open('SBM.INP', "a")
            nca=len(seq)

            gly= [pos+1 for pos,char in enumerate(seq) if char== 'G']
            print('Glycine in residues:', gly)
            ncb=len(seq) - len(gly)
            glyname='G'
            if not skip_glycine:
                ncb=len(seq)
                glyname='G111'
            
            str1= 'atoms (atomnum, atomtype, resnum, resname, atomname, charge, mass)'
            
            CB_alltypes = False #___________TESTING_____#
            g=self.get_atom_names(pdbfile,atomtype,seq,CB_alltypes)
            h=self.get_atom_types(pdbfile,atomtype,seq,CB_alltypes)
            
            atomtype_ca=1
            atomtype_cb=2
            
            #charges
            charge_in_atoms_section = True
            if atomtype==2 and CBcharge and charge_in_atoms_section:
                Q = self.amino_acid_charge_dict2()
            else:
                Q = dict()
                for i in seq:
                    Q[i]=0.0000
            charge=0.00000
            mass=1.00000
            if atomtype==2:
                natoms = nca + ncb
                f.write('\n%d %s\n' % (natoms, str1))
            if atomtype==1 and not dswap:
                f.write('\n%d %s\n' % (nca, str1))
                natoms=nca
            if atomtype==1 and dswap:
                f.write('\n%d %s\n' % (2*nca, str1))
                natoms=nca

            #opening native_cb.pdb for extracting residue numbers
            if atomtype==2:
                native = PDBParser().get_structure("test","native_cb.pdb")
            elif atomtype==1:
                native = PDBParser().get_structure("test","native_ca.pdb")
            else:
                print("Fatal error. Check attype input.")
                
            all = native.get_residues()
            residues = list()
            for i in all:
                residues.append(i)

            #write Calphas first
            rescount=0
            count=0
            all=[]
            for i in seq:
                res_num = residues[rescount].id[1]
                f.write(' %d\t%d\t%d\t%s\t%s\t%f\t%f\n' % (count+1, h[g[count].strip()]+1,res_num,d[i],g[count],charge,mass))
                all.append([count+1, atomtype_ca,res_num,d[i],'CA',charge,mass])
                count=count+1
                if atomtype==2:
                    #print count
                    if i!=glyname:
                        f.write(' %d\t%d\t%d\t%s\t%s\t%f\t%f\n' % (count + 1, h[g[count].strip()]+1, res_num, d[i], g[count], Q[i], mass))
                        all.append([count + 1, atomtype_cb, res_num, d[i], 'CB', Q[i], mass])               
                        count=count+1
                    else:
                        print ("skipped glycine",res_num)
                rescount=rescount+1     
            assert count==natoms
            return all,g,h
        def write_contacts_section(self,pdbfile,nativefile,btfile,cutoff,atomtype,sopc,btparams,dswap,CG_indices):
            #Cheung-Thirumalai protSBM. #Parameters from btfile are read.
            #should write a file called  CT_SBM.INP that can be used with OPTIM.
            ##################################native###############################
            assert atomtype <= 2
            U=Utils()
            print ('>>> writing contacts section,atomtype=',atomtype,'SOPC=',sopc)
            cutoff=float(cutoff)
            d = self.amino_acid_dict2()
            f = open('SBM.INP', "a")
            f1= open('contacts.txt', "w+")
            f2 = open('backbone.txt', "w+")

            contacttype=2
            strength_CA=float(1.0) #eps bb-bb
            strength_CB=float(1.0) #eps bb-sc

            # if U.file_exists(btfile) and btparams:
            #     btmap = self.get_dict_from_csv(btfile)
            #get CA-CA contacts.(CA atomtype=1,CB atomtype=2)
            #dict_seq = self.two_lists_to_dict(np.arange(0, len(seq)).tolist(), list(seq))
            #nc=native,sc=side-chain
            #--------------------------------
            if atomtype==1:
                Y = conmaps()
                L = pdb_kit()
                nc = Y.all_atom_contacts(pdbfile, cutoff,1.2)
                nc_2 =  L.get_backbone_contacts(pdbfile, cutoff,1.2,4)
                print (len(nc),len(nc_2))
                exit()
                dist_ca = nc[:, 3]
                #print (dist_ca)
                contacts_ca = nc[:, (0, 1)]
                l1 = CG_indices[0]
                num_atoms = len(l1)
                l2 = np.arange(0, len(l1))
                assert len(l1) == len(l2)
                dict1 = self.two_lists_to_dict(l2, l1)
                #print (dict1)
                num_contacts=len(contacts_ca)
                if hphobic:
                    hp_pairs=Y.write_hydrophobic_contacts('contacts.txt',pdbfile)
                    num_contacts=num_contacts+len(hp_pairs)
                f.write('\n%d\t%s\n' % (num_contacts, 'contacts'))
                contacts = []
                count=0
                print ("CA_contacts\n")
                if dsb:
                    contacttype=8;strength_CA=1.00000e+00
                for i in contacts_ca:
                    if dsb:
                        contacts.append([dict1[int(i[0])] + 1, dict1[int(i[1])] + 1, contacttype, dist_ca[count] * 10, strength_CA])
                    contacts.append([dict1[int(i[0])]+1,dict1[int(i[1])]+1,contacttype,dist_ca[count]*10,strength_CA])
                    f.write('%s\t%s\t%d\t%f\t%e\n'%(dict1[int(i[0])]+1,dict1[int(i[1])]+1,contacttype, dist_ca[count]*10,strength_CA))
                    f2.write(' %s\t%d\t%s\t%d\n' % ('1', dict1[int(i[0])] + 1, '1', dict1[int(i[1])] + 1))
                    count=count+1
                if hphobic:
                    print("hydrophobic contacts\n")
                    print ("Adding hydrophobic contacts: Distance =",hpdist, "Angstroms")
                    #hp_pairs=Y.write_hydrophobic_contacts('contacts.txt',pdbfile)
                    p=hp_pairs-1
                    #dist_hp = md.compute_distances(traj_hp, p, periodic=False, opt=True)[0][0]
                    counthp=0
                    for i in p:
                        print (i[0],i[1],counthp)
                        contacts.append([dict1[int(i[0])]+1,dict1[int(i[1])]+1,contacttype,hpdist,hpstrength])
                        f.write('%s\t%s\t%d\t%f\t%e\n' % (
                        dict1[int(i[0])] + 1, dict1[int(i[1])] + 1, contacttype, hpdist, hpstrength))
                        f2.write(' %s\t%d\t%s\t%d\n' % ('1', dict1[int(i[0])] + 1, '1', dict1[int(i[1])] + 1))
                        counthp = counthp + 1

                f2.close()
                f.close()
                return contacts
            #-------------------------------
            elif atomtype==2:
                contacts = []
                scaling = 1.2
                if skip_glycine:
                    sopc=False
                L = pdb_kit()
                if sopc:
                   nc = L.get_bb_contacts_SOPC(pdbfile,nativefile,4,1.0,4)                
                   sc = L.get_ss_contacts_SOPC(pdbfile,nativefile,4,1.0,2)
                   nc_sc=L.get_bs_contacts_SOPC(pdbfile,nativefile,4,1.0,2)
                   #nn_angles = Y.get_SOPC_non_native_angles(pdbfile, nativefile, 3.8, 1.0)
                   nn_angles = L.get_SOPC_non_native_angles(nativefile,3.8,1.0,CG_indices)
                else:
                    #separation 
                    cacasep=4;cacbsep=3;cbcbsep=3
                    
                    #Using pdb_kit
                    nc = L.get_backbone_contacts(pdbfile, cutoff,scaling,cacasep)
                    sc = L.get_side_chain_contacts(pdbfile, cutoff,scaling,cbcbsep)
                    nc_sc = L.get_bb_sc_contacts(pdbfile, cutoff,scaling,cacbsep)
                    
                    #print ("using conmaps")#__TESTING__
                    #Y = conmaps()
                    #nc = Y.get_backbone_contacts(pdbfile, cutoff,scaling,cacasep)
                    #sc = Y.get_side_chain_contacts(pdbfile, cutoff,scaling,cbcbsep)
                    #nc_sc = Y.get_bb_sc_contacts(pdbfile, cutoff,scaling,cacbsep)
                
                print ('Found', len(nc), 'backbone contacts')
                print ('Found', len(sc), 'side-chain contacts')
                print ('Found', len(nc_sc), 'bb-sc contacts')

                CA_indices = CG_indices[0]
                CB_indices = CG_indices[1]
                atoms_in_residue = dict()
                for i in range(0,len(CA_indices)):
                    atoms_in_residue[i] = []
                    ca = CA_indices[i]
                    atoms_in_residue[i].append(ca)
                    if ca+1 in CB_indices:
                        #CB for a residue is CA + 1, if CB exists
                        atoms_in_residue[i].append(ca+1)

                #residue_list=np.arange(0,num_ca)
                #write native CA-CA contacts first.
                pairs_bb = list()
                for i in range(0, len(nc)):
                    res1 = int(nc[i][0]); res2 = int(nc[i][1])
                    #Getting CA atoms in the resiude
                    res1 = int(CA_indices[res1])
                    res2 = int(CA_indices[res2])   
                    pairs_bb.append([res1, res2])
                #dist
                pairs_bb.sort()
                distances = Measure().distances(nativefile,pairs_bb)
                distances = distances[0]
                for i in range(0,len(pairs_bb)):                                
                    dist = distances[i]
                    (res1,res2) = pairs_bb[i]
                    f2.write(' %s\t%d\t%s\t%d\n' % ('1', int(res1)+1, '1', int(res2)+1))
                    contacts.append([res1+1, res2+1, int(contacttype), dist, strength_CA])

                # write side-chain contacts next
                for i in range(0, len(sc)):
                    ##map CB in residue to CB in two bead file.
                    res1 = int(sc[i][0]);res2=int(sc[i][1])

                    res1 = atoms_in_residue[res1]
                    res2 = atoms_in_residue[res2]
                    
                    #if len of res1 is not 2, this might be a Glycine

                    #print(sc[i],res1,res2,len(res1),len(res2))
                    assert len(res1)==2;assert len(res2)==2

                    #second atom is always CB. Put assertion.
                    res1=res1[1];res2=res2[1]
                    
                    dist = sc[i][3]#[0][0]

                    #add two-bead potentials here!!
                    #btparams=Flase
                    #------------------------------#
                    if btparams:
                        Y = conmaps()
                        print (" --interactions flag set to true. External interaction matrix will be read from file p.interaction.dat")
                        # boltzmann_kcal_mol=float(0.0019872041)
                        Kb=1
                        #print d[Y.get_residue_name(pdbfile,int(res1))[0]]
                        res1bt = d[Y.get_residue_name(nativefile, int(res1))[0]]
                        res2bt = d[Y.get_residue_name(nativefile, int(res2))[0]]
                        strength_CB=float(Y.get_btmap_val(res1bt,res2bt,'interaction.dat'))
                        strength_CB=0.5*(0.7- float(strength_CB))*300*Kb*0.001987 #kcal/mol
                        contacts.append([res1 + 1, res2 + 1, int(contacttype), dist,strength_CB])
                        print ([res1 + 1, res2 + 1, int(contacttype), dist, strength_CB, 'SCSCbt'])
                    #_---------------------------#
                    else:
                        contacts.append([res1+1, res2+1, int(contacttype), dist, strength_CB])
                        #print ([res1 + 1, res2 + 1, int(contacttype), dist, strength_CB,'SCSC'])

                #write bb-sc contacts next
                pair_bb_sc = list()
                for i in range(0,len(nc_sc)):
                    res1 = int(nc_sc[i][0])
                    res2 = int(nc_sc[i][1])
                    res1 = atoms_in_residue[res1] 
                    res2 = atoms_in_residue[res2]

                    assert len(res1)>=1 and len(res2)>=1
                    #res1 = CA, res2=CB
                    res1=res1[0];res2=res2[1]
                    #print (res1, res2)
                    #exit()
                    pair_bb_sc.append([res1, res2])
                #distances
                distances = Measure().distances(nativefile,pair_bb_sc)
                distances = distances[0]
                for i in range(0,len(pair_bb_sc)):           
                    (res1,res2) = pair_bb_sc[i]
                    dist = distances[i]
                    contacts.append([res1 + 1, res2 + 1, int(contacttype), dist, strength_CA])
                    #f4.write(' %s\t%d\t%s\t%d\n' % ('1', int(res1) + 1, '1', int(res2) + 1))
                
                #write angle-angle repulsion term for sop-sc model only!
                if sopc:
                    pair_sopc = list()
                    epsl=1 #kcal/mol
                    for i in range(0,len(nn_angles)):
                        contacttype=1
                        #Use the 6 term from LJ potential for repulsion!(supply negative to make repulsive)
                        # fudgeQQ,qi,qj,V,W for gromacs
                        res1=int(nn_angles[i][0]);res2=int(nn_angles[i][1])
                        #print res1,res2
                        pair_sopc.append([res1, res2])#;  pair = np.reshape(pair, (1, 2))
                    pair_sopc.sort()   
                    distances = Measure().distances(nativefile,pair_sopc)
                    for i in range(0,len(pair_sopc)):
                        (res1,res2) = pair_sopc[i]
                        dist = distances[i]
                        #dist = md.compute_distances(traj, pair, periodic=False, opt=True)[0][0] * 10
                        contacts.append([res1 + 1, res2 + 1, int(contacttype),dist,epsl])
            else:
                U.fatal_errors(7)

            # check if contacts are being repeated.
            contacts=np.asarray(contacts)
            contacts= contacts[np.argsort(contacts[:,0])]
            test=[]
            for i in contacts:
                j=tuple(i)
                #print (j)
                test.append(j)
            #remove repeat contacts.
            a=set(test)
            print ("Found",len(a)-len(contacts),"duplicates in list. Removing...")
            a=list(a)
            f.write('\n%d\t%s\n' % (len(contacts), 'contacts'))
            for i in a:
                #print i[0],i[1],i[2],i[3],i[4]
                #"%.10lf\n
                f.write(' %d\t%d\t%d\t%.12f\t%10.3f\n' % (i[0],i[1],i[2],i[3],i[4]))
                f1.write(' %s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('1',i[0],'1',i[1], i[2], i[3], i[4]))
            f.close()
            f1.close();f2.close()
            return contacts
        def write_exclusions_section(self,contacts):
            f = open('SBM.INP', "a")
            print (">>Exclusions are automatically accounted for in OPTIM.")
            f.write('\n%d\t%s\n' % (len(contacts), 'contacts'))
            for i in range(0,len(contacts)):
                f.write('\n%d\t%s\n' % (len(contacts), 'contacts'))
            f.close()
            return
        def get_chain_termianl_residue_index(self,nativefile):
            #returns the chain terminal residue index
            
            #reading chains
            P =PDBParser()
            chains = [x for x in P.get_structure("temp",nativefile).get_chains()]
            terminal_residue_index = list()
            
            #loading residue object from the chsin
            residues = [x.get_residues() for x in chains]
            prev_chain_length = 0
            for r in residues:
                #lenth of residues in each chain - 1 gives the index of terminal residue
                prev_chain_length = prev_chain_length + len([x for x in r]) 
                terminal_residue_index.append(prev_chain_length-1)
            #if 2 chains have same name, then the pdbparsar will consider as same chain
            return terminal_residue_index
        def write_bonds_section(self,nativefile,Kb,atomtype,CG_indices):
            print (">> writing bonds_section: Native= ",nativefile)
            
            f = open('SBM.INP', "a")
            #nativefile='native_cb.pdb'

            #if the model is CA only, CB_indices will be empty
            (CA_indices,CB_indices,terminal_residue) = CG_indices

            #getting CB indices
            pairs_CA=[]
            pairs_CAB=[]
            for residue_index in range(len(CA_indices)):
                if residue_index not in terminal_residue:
                    #Not including bonded term for CAs of two different chains
                    #print (terminal_residue)
                    print (residue_index,residue_index+1,[CA_indices[residue_index],CA_indices[residue_index+1]])
                    pairs_CA.append([CA_indices[residue_index],CA_indices[residue_index+1]])
                if CA_indices[residue_index] + 1 in CB_indices:
                    #As CB_indices is empty if atomype == 1
                    #if CA_indice + 1 is present in CB_indices
                    pairs_CAB.append([CA_indices[residue_index],CA_indices[residue_index] + 1])


            M = Measure()
            #values returned are in A
            pairs_total = pairs_CA + pairs_CAB
            dists = M.distances(nativefile,pairs_total)
            for i in range(0,len(pairs_total)):
                assert(pairs_total[i]==dists[i+1])
            
            num_bonds = len(dists[0])

            bondtype=1
            K_bond = np.array(Kb, dtype=float)
            btype = np.ones(num_bonds).T*bondtype

            assert len(pairs_total) == len(dists[0]) == len(dists)-1 == len(btype)
           
            all = []
            for i in range(0,len(dists[0])):
                all.append((dists[1+i][0],dists[1+i][1],dists[0][i]))
            all.sort()

            f.write('\n%d\t%s\n' % (num_bonds, 'bonds'))
            bonds=[]
            for i in all:
                bonds.append([i[0]+1, i[1]+1, bondtype, i[2], K_bond])
                f.write(' %d\t%d\t%d\t%8.5f\t%e\n' % (i[0]+1, i[1]+1, bondtype, i[2], K_bond))
            return bonds
        def write_angles_section(self,nativefile,Ka,atomtype,CG_indices):
            print ('>>writing angles_section',nativefile,Ka)


            #if the model is CA only, CB_indices will be empty
            #loading CA and CB indices
            (CA_in_native,CB_in_native,terminal_residue) = CG_indices


            triplets=[]
            if atomtype == 1:
                for i in range(1,len(CA_in_native)-1):
                    triplets.append([CA_in_native[i-1], CA_in_native[i], CA_in_native[i+1]])
            
            elif atomtype == 2:
                #terminal residues CA indices
                termianl_CA = list()
                
                #Adding CAi-1_CAi_CAi+1 triplets
                for i in range(1,len(CA_in_native)-1):
                    #i here is residue index
                    #adding a-1:a:a+1 CAs
                    if i-1 not in terminal_residue and i not in terminal_residue:
                        #avoiding angle between CAs of different chain
                        triplets.append([CA_in_native[i-1],CA_in_native[i],CA_in_native[i+1]])
                    elif i in terminal_residue:
                        termianl_CA.append(CA_in_native[i])
                    elif i-1 in terminal_residue:
                        termianl_CA.append(CA_in_native[i-1])
                
                #Adding CAi-1_CAi_CBi and CB_CAi_CAi+1
                a_prev = 0
                for b in CB_in_native:
                    #Adding a-1:a:b and b:a:a+1 angles
                    #avoiding if a-1 or a are terminal CAs
                    #b is CB atom number. Hence CA atom number a-1
                    a = b-1         #CA atom number of same residue
                    if a!=0 and a_prev not in termianl_CA:
                        #a-1:a:b #if a_prev is terminal CA, donot include the triplet    
                        triplets.append([a_prev,a,b])
                    if a not in termianl_CA and b!=CB_in_native[len(CB_in_native)-1]:                    
                        #b:a:a+1    #if a is in termial CA or b is last CB atom, donot include the triplet
                        a_next = b+1    #CA atom of next residue    
                        triplets.append([b,a,a_next])
                    a_prev = a      #CA atom number of prev residue for next cycle


            num_angles=len(triplets)
            f = open('SBM.INP', "a")
            f.write('\n%d\t%s\n' % (num_angles, 'angles'))

            M = Measure()
            CG_angles=M.angles(nativefile,triplets)
            K_angle=float(Ka)

            count=0
            all=[]
            count=0
            for i in range(len(triplets)):
                #print triplets[i],CG_angles[0][i]
                all.append([triplets[i][0]+1,triplets[i][1]+1,triplets[i][2]+1, CG_angles[0][i],K_angle])
                f.write(' %s\t%s\t%d\t%f\t%e\n' % (triplets[i][0]+1,triplets[i][1]+1,triplets[i][2]+1, CG_angles[0][i],K_angle))
                count=count+1
            if dswap:
                num_atoms = len(CA_in_native)+len(CB_in_native)
                count=0
                for i in range(len(triplets)):
                    all.append([triplets[i][0]+1+num_atoms,triplets[i][1]+1+num_atoms,triplets[i][2]+1+num_atoms, CG_angles[0][i], K_angle])
                    f.write(' %s\t%s\t%d\t%f\t%e\n' % (triplets[i][0]+1+num_atoms,triplets[i][1]+1+num_atoms,triplets[i][2]+1+num_atoms, CG_angles[0][i],K_angle))
                    count=count+1
            f.close()
            return  all
        def get_CB_chiral_dihedrals(self,nativefile,terminal_residue):
            print ('>>writing CB chiral dihedrals',nativefile)
            fnat = open(nativefile)
            lines = fnat.readlines()
            fnat.close()
            
            residues = dict()
            res_count = 0
            for i in lines:
                if res_count == 0:
                    start_residue = int(i[22:26]) 
                res_count = res_count + 1
                if i.startswith("ATOM"):
                    r = int(i[22:26]) - start_residue
                    if  not r in residues:
                        residues[r] = dict()
                    residues[r][i[12:16].strip()] = int(i[6:11]) - 1
            
            quadraplets = list()
            for r in range(1,len(residues)-1):
                if r-1 not in residues or "CB" not in residues[r]:
                    #skip the loop if there is no CB (GLY) ir r-1 is not a residue index
                    continue
                if r in terminal_residue or r-1 in terminal_residue:
                    #skip the loop if r or r-1 are the indices of termianl residues
                    continue
                #print (r)
                #avoid chiral shift 
                #adding Cai-1 : Cai+1 : Cai : Cbi improper dihedral term
                #ref: doi:10.1016/j.jmb.2009.08.010 Yaakov_Levy_2009
                ai = residues[r-1]["CA"]
                aj = residues[r+1]["CA"]
                ak = residues[r]["CA"]
                al = residues[r]["CB"]
                quadraplets.append([ai,aj,ak,al])

            #for i in quadraplets:
            #    print (i)
            return quadraplets 
        def write_dihedrals_section(self,nativefile,Kd,atomtypes,CG_indices):
            #Only over C-alpha atoms for Cheung-Thirumalai model.

            print ('>>> writing dihedrals section',nativefile,Kd)

            f = open('SBM.INP', "a")
            
            #Select only C-alpha atoms. #CA indices
            indices = CG_indices[0]
            num_atoms=len(indices)

            #get terminal residue index
            terminal_residue = CG_indices[2]

            quadruplets = []
            K_dihedral=float(Kd)
            for i in range(2, len(indices) - 1):
                #here i is amino-acid residue index
                if i-2 not in terminal_residue and i-1 not in terminal_residue and i not in terminal_residue:
                    #if i-2,i-1 or i are termianl CA, donot include in the quadraplet
                    quadruplets.append([indices[i-2],indices[i-1], indices[i], indices[i+1]])
            if atomtypes == 2:
            #improper dihedrals (updated version using virtual sites with harmonic bond)
                #_____TESTING_________#
                no_chiral = False
                #_____________________#
                if not no_chiral:
                    CB_quad=self.get_CB_chiral_dihedrals(nativefile,terminal_residue)
                    for i in CB_quad:
                        quadruplets.append(i)
            if dswap:
                num_dihedrals=2*len(quadruplets)
            else:
                num_dihedrals = len(quadruplets)

            f.write('\n%d\t%s\n' % (num_dihedrals, 'dihedrals'))

            count=0
            all=[]
            M = Measure()
            dihedrals = M.dihedrals(nativefile,quadruplets)
            #chirals   = M.dihedrals(nativefile,CB_quad)
            count=0
            for i in range(len(quadruplets)):
                #print triplets[i],CA_angles[0][i]
                all.append([quadruplets[i][0],quadruplets[i][1],quadruplets[i][2],quadruplets[i][3], '1 ', dihedrals[0][i],K_dihedral])
                f.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\n'% (quadruplets[i][0]+1,quadruplets[i][1]+1,quadruplets[i][2]+1,quadruplets[i][3]+1, '1 ', dihedrals[0][i],K_dihedral))
                count=count+1
            if dswap:
                for i in range(len(quadruplets)):
                    all.append([quadruplets[i][0]+num_atoms, quadruplets[i][1]+num_atoms, quadruplets[i][2]+num_atoms, quadruplets[i][3]+num_atoms, '1 ',
                                dihedrals[0][i], K_dihedral])
                    f.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\n' % (
                    quadruplets[i][0]+1+num_atoms, quadruplets[i][1]+1+num_atoms, quadruplets[i][2]+1+num_atoms, quadruplets[i][3]+1+num_atoms, '1 ', dihedrals[0][i], K_dihedral))
            all=np.asarray(all)
            f.close()
            #chi = list()
            #for i in range(len(CB_quad)):
            #    chi.append([CB_quad[i][0],CB_quad[i][1],CB_quad[i][2],CB_quad[i][3], '1 ', chirals[0][i],K_dihedral])
            #    count=count+1
            return all

#class protsbm(Select,object):
        def write_gro_moleculetype(self,topfilename):
            print ('>>writing GROMACS moleculetypesection', topfilename)
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ moleculetype ]'))
            f.write('%s\n' % ('; name            nrexcl'))
            f.write('%s\n' % ('  Macromolecule   3'))
            f.close()
            return 1
        def write_gro_header(self,topfilename,atomtypes):
            print ('>>> writing Protein GROMACS toptology header section', topfilename)
            f = open(topfilename, "w+")
            f.write('%s\n' % (';'))
            f.write('%s\n' % ('; Topology file generated from PYeSBM repo. '))
            f.write('%s\n' % ('; https://bitbucket.org/nsridhar/pyesbm/wiki/Home'))
            f.write('\n%s\n' % ('[ defaults  ]'))
            f.write('%s\n' % ('; nbfunc comb-rule gen-pairs'))
            f.write('%s\n' % ('  1      1         no   \n'))
            return
        def write_gro_nonbondparams(self,topfilename,atoms_in_top,CA_rad):
            f = open(topfilename, "a")
            #f.write('%s\n' % ('[ nonbond_params ]'))
            #add non-bonded r6 term in sop-sc model for all non-bonded non-native interactions.
            #f.write('%s\n' % ('; i j  func sigma(c10)    eps(c12)'))

            #CA-CA
            #r = ((CA_rad/10 + cA_rad/10)/2)**12
            r = ((CA_rad/10 + CA_rad/10))**12
            f.write('  %s\t%s\t1\t%e\t%e\n' % ('CA','CA',0,r))
            #CA-CB excl
            for i in (atoms_in_top):
                if i[0]!='CA':
                    #r = ((CA_rad/10 + i[5]/10)/2)**12
                    r = ((CA_rad/10 + i[5]/10))**12
                    f.write('  %s\t%s\t1\t%e\t%e\n' % ('CA',i[0].strip(),0,r))
            for x in range(0,len(atoms_in_top)):
                i = atoms_in_top[x]
                if i[0]!="CA":
                    for y in range(x,len(atoms_in_top)):
                        j = atoms_in_top[y]
                        if j[0]!='CA':
                            #r  = ((i[5]/10 + j[5]/10)/2)**12
                            r  = ((i[5]/10 + j[5]/10))**12 
                            f.write('  %s\t%s\t1\t%e\t%e\n' % (i[0].strip(),j[0].strip(),0,r))
            f.close()
        def write_gro_atomtypes(self,topfilename,atomtypes,pdbfile,sopc,CA_rad,CBcom,CBradii,excl_rule):
            print (">>> writing Protein GROMACS topology atomtypes section",topfilename,atomtypes,pdbfile,sopc,CA_rad)

            f = open(topfilename, "a")
            #1:CA model or 2:CA+CB model
            assert atomtypes<=2
            
            a=self.write_atomtypes_section(pdbfile,atomtypes,sopc,CA_rad,CBcom,CBradii)
            assert a[0][0]=='CA'    #make sure that tha fist element is CA
            ca=a[0]

            CA_rad12=(CA_rad/10)**12
            f.write('%s\n'%('[ atomtypes ]'))
            f.write('%s\n' % ('; name mass  charge ptype c6    c12'))
            f.write('   %s %8.3f %8.3f %s\t%e\t%e \n' % (ca[0], ca[1], ca[2], ca[3], ca[4], CA_rad12))
            for i in (a):
                if i[0]!='CA':
                    r=(i[5]/10)**12 #C12 term
                    f.write('  %s %8.3f %8.3f %s\t%e\t%e \n'%(i[0],i[1],i[2],i[3],i[4],r))
    
            #add non-bonded r6 term in sop-sc model for all non-bonded non-native interactions.
            f.write('\n%s\n' % ('[ nonbond_params ]'))
            f.write('%s\n' % ('; i j  func sigma(c10)    eps(c12)'))
            f.close()
            if excl_rule == 2:
                #The input files by default direct use of Combination rule 1 (Geometric mean) in Gromacs
                #To define Arith. mean using excl_rule as 2
                self.write_gro_nonbondparams(topfilename,a,CA_rad)
            return 0;
        def write_gro_atoms(self,topfilename,atomtypes,pdbfile,nativefile,sopc):
            print (">> writing Protein GROMACS tpology atom section",topfilename,atomtypes,pdbfile,nativefile,sopc)
            
            f = open(topfilename, "a")
            assert atomtypes<=2
            #d = atom names and e = atom types
            atoms,d,e = self.write_atoms_section(pdbfile,atomtypes,sopc)

            num_atoms=len(atoms)

            f.write('\n%s\n' % ('[ atoms ]'))
            f.write('%s\n' % (';nr  type  resnr residue atom  cgnr'))

            for i in range(0,len(atoms)):
                d[i]=d[i].rjust(4)
                f.write('\t%d\t%s\t%d\t%s\t%s\t%d\t%6.2f\t%6.2f \n' %(atoms[i][0],d[i],atoms[i][2],atoms[i][3],d[i].strip()[:2],atoms[i][0],atoms[i][5],atoms[i][6]))
            if dswap:
                for i in range(0, len(atoms)):
                    d[i] = d[i].rjust(4)
                    f.write('\t%d\t%s\t%d\t%s\t%s\t%d\t%6.2f\t%6.2f \n' % (
                    atoms[i][0]+num_atoms, d[i], atoms[i][2]+num_atoms,atoms[i][3], d[i].strip()[:2], atoms[i][0],atoms[i][5],atoms[i][6]))
            f.close()
            return atoms,d,e
        def write_gro_bonds(self,topfilename,atomtypes,nativefile,Kb,ptype,dsb,CG_indices):
            print (">> writing Protein GROMACS topology bonds section", topfilename, atomtypes,nativefile,Kb)
            #GROMACS IMPLEMENTS Ebonds = (Kx/2)*(r-r0)^2
            #Input units KJ mol-1 A-2
            #GROMACS units KJ mol-1 nm-1 (100 times the input value) 
            Kb = float(Kb*100)    

            #GROMACS 4.5.2 : FENE=7 AND HARMONIC=1
            allowed_pots=[1,7,9]
            ptype=1 #READ POTENTIAL FROM DICTIONARY HERE IF NEEDED.
            assert ((ptype in allowed_pots)),'Make sure potential is included in allowed_pots.'
            assert atomtypes <= 2

            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ bonds ]'))
            f.write('%s\t%s\t%s\t%s\t%s\n' % (';ai', 'aj', 'func', 'r0(nm)', 'Kb'))
            bonds =self.write_bonds_section(nativefile,Kb,atomtypes,CG_indices)
            
            if dsb: #global variable
                ptype=9
            
            if ptype==1:
                #harminoc 
                for i in range(0,len(bonds)):
                    f.write('\t%d\t%d\t%s %e %e\n' %(bonds[i][0],bonds[i][1], str(ptype),bonds[i][3]/10,Kb))
                f.close()
                return bonds
            
            elif ptype==7:
                #FENE potential.
                kcalAtokjA=418.4 #kcal/mol/A2 to kcal/mol/nm2
                Kb=float(20*kcalAtokjA) #Gromacs units
                R=float(0.2) #nm
                for i in range(0,len(bonds)):
                    #not sure of format for top file. (where does distance go?)
                    f.write('\t%d\t%d\t%s %e %e %e\n' %(bonds[i][0]+1,bonds[i][1]+1, str(ptype),bonds[i][2]/10,R,Kb))
                f.close()
            
            elif ptype==9:
                for i in range(0,len(bonds)):
                    f.write('\t%d\t%d\t%s %e %e\n' %(bonds[i][0],bonds[i][1], "1",bonds[i][3]/10,Kb))
                #contacts written in bonds section.
                #get_contacts from contacts.txt
                #loading external pairs
                dsb_contacts=np.int_(np.loadtxt('contacts.txt',unpack=True))
                dsb_contacts =(dsb_contacts[[1, 3], :])
                dsb_contacts=dsb_contacts.T - 1
                print (dsb_contacts[0][0])
                count_dsb=0
                for i in range(0,len(dsb_contacts)):
                    f.write('\t%d\t%d\t%s %s %s\n' % (dsb_contacts[i][0],dsb_contacts[i][1], str(ptype),str(count_dsb), '1'))
                    count_dsb=count_dsb+1
                #generate table_files here
                Z=tables()
                if hphobic:
                    Z.gen_many_table_file('contacts.txt', nativefile, 0.06, 350, 0.02,True)
                else:
                    Z.gen_many_table_file('contacts.txt',nativefile,0.06,350,0.02,False)
                return bonds
        def write_gro_angles(self,topfilename,atomtypes,nativefile,Ka,CG_indices):
            print (">>> writing Protein GROMACS ropology angles section", topfilename, atomtypes, nativefile, Ka)
            #Eangles = (Ktheta/2)*(r-r0)^2
            #Input units KJ mol-1
            ##GROMACS units KJ mol-1  
            Ka = float(Ka)

            assert atomtypes <= 2
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ angles ]'))
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (';ai', 'aj', 'ak','func', 'th0(deg)', 'Ka'))
            angles=self.write_angles_section(nativefile,Ka,atomtypes,CG_indices)
            for i in range(0, len(angles)):
                f.write('\t%d\t%d\t%d\t%s %e %e\n' % (angles[i][0],angles[i][1],angles[i][2],'1',angles[i][3]*radtodeg, Ka))
            f.close()
            return angles
        def write_gro_dihedrals(self,topfilename,atomtypes,nativefile,Kd,CG_indices):
            print (">>> writing Protein GROMACS topoplogy dihedrals section", topfilename, atomtypes, nativefile, Kd)
            #converting Kcal/mol to KJ/mol
            #K_dihedrals will be divided by these factprs
            #_____USER_CAN_EDIT____#
            #The Kd will be divided the factor below
            factor_1 = 1.0      #for term with multiplicity 1
            factor_2 = 1.0      #for term with multiplicity 3
            #______________________#


            #GROMACS IMPLEMENTATION: Edihedrals Kphi*(1 + cos(n(phi-phi0)))
            #Our implementaion: Edihedrals = Kphi*(1 - cos(n(phi-phi0)))
            #The negative sign is included by add phase = 180 to the phi0
            #Kphi*(1 + cos(n(phi-180-phi0))) = Kphi*(1 + cos(n180)*cos(n(phi-phi0)))
            #if n is odd i.e. n=1,3.... then cos(n180) = -1
            #hence Edihedrals = Kphi*(1 - cos(n(phi-phi0)))

            phase = 180
            radtodeg = 180/np.pi
            Kd=float(Kd)        #KJ/mol
            assert atomtypes <= 2

            #dihedrals section
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ dihedrals ]'))
            d=self.write_dihedrals_section(nativefile,Kd,atomtypes,CG_indices)
            f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (';ai','aj','ak','al','func','phi0(deg)','Kd','mult'))
            for i in range(0, len(d)):
                #GROMAC dihedrals in degrees. For dihedraltype 1 add 180.
                l1= float(d[i][5])*radtodeg+phase
                d1=int(d[i][0])         
                d2 = int(d[i][1])
                d3 = int(d[i][2])
                d4 = int(d[i][3])
                d5 = int(d[i][4])
                if d5!=1:
                    U.fatal_errors(8)
                #print l1,d1,d2,d3,d4
                f.write('\t%d\t%d\t%d\t%d %d %e %e %s\n' % (d1+1, d2+1, d3+1, d4+1,d5,l1,Kd/factor_1,'1'))
                f.write('\t%d\t%d\t%d\t%d %d %e %e %s\n' % (d1+1, d2+1, d3+1, d4+1,d5,l1*3,Kd/factor_2,'3'))
            #f.write(";;;;;;;;chirality_improper_dihed;;;;;;;;;\n")
            #for i in range(0, len(chi)):
            #    #GROMAC dihedrals in degrees. For dihedraltype 1 add 180.
            #    l1= float(d[i][5])*radtodeg
            #    d1 = int(chi[i][0])         
            #    d2 = int(chi[i][1])
            #    d3 = int(chi[i][2])
            #    d4 = int(chi[i][3])
            #   d5 = 2
            #    f.write('\t%d\t%d\t%d\t%d %d %e %e %s\n' % (d1+1, d2+1, d3+1, d4+1,d5,l1,Kd*40,' '))
            f.close()
        def write_gro_pairs(self,topfilename,atomtypes,nativefile,pdbfile,contacttype,cutoff,sopc,btparams,CG_indices):
            print (">> writing GROMACS pairs sections",topfilename,atomtypes,nativefile,pdbfile,contacttype,cutoff,sopc,btparams)
            contacts_allowed=[1,2]
            assert atomtypes <= 2
            if skip_glycine:
                sopc=False

            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ pairs ]'))

            f.write('\t%s\n' % ('; ai aj type A B'))
            btfile='interaction.dat'
            #cutoff=4.5
            contacts=self.write_contacts_section(pdbfile,nativefile,btfile,cutoff,atomtypes,sopc,btparams,dswap,CG_indices)
            assert (contacttype in contacts_allowed),'Only the 10-12 potential contacttype=2 is implemented'
            if dsb: #global variable
                #contacts will be added as bonds
                return contacts

            for i in range(0,len(contacts)):
                contacttype=int(contacts[i][2])
                #print contacttype
                if contacttype==2:
                    # gromacs contact type is 1 for table potentials.
                    # optim contact types differ.

                    epsilonij=float(contacts[i][4])
                    #print epsilonij
                    #12-19 potential.
                    B=((contacts[i][3]*0.1)**12)*5*epsilonij
                    A=((contacts[i][3]*0.1)**10)*6*epsilonij
                    #i-j
                    f.write('\t%d\t%d\t%s\t%e\t%e\n' % (contacts[i][0],contacts[i][1],'1',A,B))
                # elif contacttype==1 and sopc:
                #     count_angles=count_angles+1
                #     A=float(0.0)
                #     B=-(contacts[i][3]*0.1)**6 #sigmabb*6 and (sigma_bs)**6 (minus for repilsive)
                #     #print contacts[i][3]
                #     q1=0;q2=0;fudgeqq=1
                #     #Fudgeqq,qi,qj,V,W
                #     f.write('\t%d\t%d\t%s\t%d\t%e\t%e\t%e\t%e\n' % (contacts[i][0], contacts[i][1],'2',fudgeqq,q1,q2,B,A))
                #     print ("Writing", count_angles, "additional entries in pairs section. SOPC =",sopc)
            f.close()
            return contacts
        def write_gro_exclusions(self,topfilename,contacts):
            print (">> writing GROMACS exclusionss sections",topfilename)
            self.write_exclusions_section(contacts)
            #write exlusions
            f = open(topfilename, "a")
            f.write('\n\t%s\n' % ('[ exclusions ]'))
            f.write('\t%s\n' % ('; ai aj'))
            for i in range(0,len(contacts)):
                f.write('\t%d\t%d\n'%(contacts[i][0],contacts[i][1]))
            f.close()
            return 1
        def write_gro_tail(self,topfilename):
            print (">>> writing GROMACS tail section",topfilename)
            f = open(topfilename, "a")
            f.write('\n%s\n' % ('[ system ]'))
            f.write('%s\n' % (';name'))
            f.write('  %s\n' % ('Macromolecule'))
            f.write('\n%s\n' % ('[ molecules ]'))
            f.write('%s\n' % (';name    #molec'))
            f.write('%s\n' % ('Macromolecule     1'))
            f.close()
            return
        def write_gro_gro(self,grofilename,atomtypes,num_at):
            print ('>> write_gro_gro',grofilename)
            if atomtypes==1:
                native = open("native_ca.pdb")
            elif atomtypes==2:
                native = open("native_cb.pdb")
            else:
                print ('atomtypes must be 1 or 2')
                exit()

            f = open(grofilename, "w+")
            f.write('%s\n'%('Grofile generated from Lab16 repo. See https://bitbucket.org/nsridhar/lab16'))
            f.write('%d\n'%(num_at))
            k = {}
            for i in native:
                if i.startswith("ATOM"):
                    k[0] = i[22:26].strip().rjust(5)        #residue number grofile
                    k[1] = i[17:20].strip().ljust(5)        #residue name
                    k[2] = i[12:16].strip().rjust(5)        #atom name
                    k[3] = i[ 6:11].strip().rjust(5)        #atom number
                    k[4]= str('%8.3f' % (float(i[30:38])/10)).rjust(8)  #X 31-38, 8char 
                    k[5]= str('%8.3f' % (float(i[38:46])/10)).rjust(8)	#X 39-46, 8char
                    k[6]= str('%8.3f' % (float(i[46:54])/10)).rjust(8)	#X 47-54, 8char
                    f.write('%s%s%s%s%s%s%s\n' % (k[0],k[1],k[2],k[3],k[4],k[5],k[6]))
            f.write('%8.4f %8.4f %8.4f'%(10.000,10.000,10.0000))
            f.close()
            native.close()
            return
        #def get_gro_from_xyz(self, pdbfile,nativefile, xyzfile, grofilename):
            ## convert .xyz to .gro
            #Y=conmaps()
            #traj = Y.read_xyz_as_traj(xyzfile, nativefile)
            #xyz=traj.xyz
            #atomnames=self.get_atom_names(pdbfile,1)
            #residues=Y.get_sequence(pdbfile)
            ##print (residues)
            #l1=np.arange(0,len(residues))
            #d=self.two_lists_to_dict(l1,residues)
            #d1=self.amino_acid_dict()
            #f = open(grofilename, "w+")
            #f.write('%s\n' % ('Grofile generated from Go-kit'))
            #f.write('%d\n' % (len(atomnames)))
            #for i in range(0,len(residues)):
            #    #print i,d1[d[i]],atomnames[i],xyz[0][i][1]
            #    f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (
            #    i+1,d1[d[i]],atomnames[i],i+1,xyz[0][i][0],xyz[0][i][1],xyz[0][i][2]))

            #f.write('%8.4f %8.4f %8.4f' % (5.0000, 5.0000, 5.0000))
            #f.close()
            #return 1
        def write_gromacs_top(self,topfilename,atomtypes,pdbfile,nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii,excl_rule):
            #Order for SBM file
            assert atomtypes <= 2
            #2 = 10-12 potential
            # Note: Parameters in GROMACS units need tweaking.
            contacttype=2;bond_type=1
            assert contacttype==2

            #Generate important parameters
            #Generatin atom indices
            #if the model is CA only, CB_indices will be empty
            indices,CA_indices,CB_indices = self.get_CG_index(nativefile,atomtypes)
            num_ca=len(CA_indices)
            num_cb=len(CB_indices)        
            num_atoms=len(indices)
            assert (num_atoms==num_ca+num_cb)
            #List of terminal residues
            global terminal_residue 
            terminal_residue = self.get_chain_termianl_residue_index(nativefile)
            CG_indices = (CA_indices,CB_indices,terminal_residue)

            #write gromacs format files.
            self.write_gro_header(topfilename,atomtypes)
            self.write_header_SBM()
            self.write_gro_atomtypes(topfilename,atomtypes,pdbfile,sopc,CA_rad,CBcom,CBradii,excl_rule)
            self.write_gro_moleculetype(topfilename)
            atoms=self.write_gro_atoms(topfilename, atomtypes, pdbfile, nativefile, sopc)
            contacts=self.write_gro_pairs(topfilename, atomtypes, nativefile, pdbfile, contacttype, cutoff, sopc, btparams,CG_indices)
            self.write_gro_bonds(topfilename, atomtypes, nativefile, Kb, bond_type,dsb,CG_indices)
            self.write_gro_exclusions(topfilename,contacts)
            self.write_gro_angles(topfilename,atomtypes,nativefile,Ka,CG_indices)
            self.write_gro_dihedrals(topfilename,atomtypes,nativefile,Kd,CG_indices)
            #self.write_gro_pairs(topfilename, atomtypes, nativefile, pdbfile, contacttype, cutoff, sopc, btparams,CG_indices)

            #for dswap write original protein files again.
            #self.write_CB_to_native(pdbfile,sopc)
            #  bond_ptype = 1
            #self.write_gro_atoms(topfilename,atomtypes,pdbfile,nativefile)
            #  if not sopc:
            #      self.write_gro_angles(topfilename,atomtypes,nativefile,Ka)
            #      self.write_gro_dihedrals(topfilename,atomtypes,nativefile,Kd)

            self.write_gro_tail(topfilename)
            #self.write_gro_gro('gromacs.gro',atomtypes,len(atoms))
            print ('See file:',topfilename,'for GROMACS topology.')
            return atoms

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate GROMACS and OPTIM potential files for enhanced SBM models.")
    #CA_rad (float)
    parser.add_argument("--CA_rad","-CA_rad", help="Radius for C-alpha atom. Default=4.0")
    #CA at COm on/off: T/F
    parser.add_argument("--CAcom","-CAcom",action='store_true',help="Place C-alpha at COM of backbone")
    #CB radous
    parser.add_argument("--Cb_rad","-CB_rad", help="Statistically derived values by default for attype 2.")
    parser.add_argument("--Kb","-Kb", help="Kbond")
    parser.add_argument("--Ka","-Ka", help="Kangle")
    parser.add_argument("--Kd","-Kd", help="Kdihedral")
    parser.add_argument("--cutoff","-cutoff", help="Cut-off for contact-map generation")
    parser.add_argument("--scaling","-scaling", help="Scaling for mapping to all-atom contact-map.")
    parser.add_argument("--attype", "-attype",help="Number of atom types. E.g. 1 for CA, 2 for CA and CB")
    parser.add_argument("--interaction","-interaction",action='store_true', default=False, help='User defined interactions in file interaction.dat.')
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
    parser.add_argument('--CB_radii',"-CB_radii",action='store_true', help='External CB radius stored in file radii.dat')
    parser.add_argument("--CBcharge","-CBcharge", action='store_true', default=False, help='Put charges on CB for K,L,H,D,E')
    parser.add_argument("--excl_rule",help="Use 1: Geometric mean. 2: Arithmatic mean")
    parser.add_argument("--Kr", help="Krepulsion. Default=5.7A")

    args = parser.parse_args()
    X = protsbm()
    Y = conmaps()
    U = Utils()

    #Set default parameters
    sopc=True        
    skip_glycine=False
    btparams=False
    mjmap=False
    btmap=False
    dswap=False
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
    CBcharge=False
    #print (args.CB_radii)

    #Loading input parameters
    #Force constants
    if args.Kb:
        Kb=float(args.Kb)
    if args.Kd:
        Kd=float(args.Kd)
    if args.Ka:
        Ka=float(args.Ka)
    if args.cutoff:
        cutoff=float(args.cutoff)
    if args.excl_rule:
        excl_rule = int(args.excl_rule)
        if excl_rule not in (1,2):
            print ("Choose correct exclusion rule. Use 1: Geometric mean or 2: Arithmatic mean")
        else:
	        excl_rule = 1

    if args.interaction:
        U.file_exists('interaction.dat')
        btparams=True
    else:
        btparams=False
    
    #other parameters
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
    if args.CB_rad:
        CB_rad = float(args.CB_rad)
        CBradii = True
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
        import shutil
        shutil.copy2('btmap.dat','interaction.dat')
        btparams=True;btmap=True
    if args.mjmap:
        import shutil
        shutil.copy2('mjmap.dat', 'interaction.dat')
        btparams=True;mjmap=True
    if args.CBcharge:
        CBcharge = True
    #########
    X.globals(Ka,Kb,Kd,CA_rad,skip_glycine,sopc,dswap,btparams,CAcom,hphobic,hpstrength,hpdist,dsb,mjmap,btmap,CBfar,CBcharge)

    if not args.ext_conmap and args.aa_pdb:
        pdbfile=args.aa_pdb
        Y.all_atom_contacts(pdbfile,cutoff,1.2)

    if not args.attype:
        if args.w_native:
            pdbfile=args.w_native
            Y.check_hetatm(pdbfile)
            if not skip_glycine:
                Y.check_glycineh(pdbfile,False)
            X.write_CB_to_native(pdbfile,False,1)
        if args.aa_pdb:
            pdbfile=args.aa_pdb
            X.write_CB_to_native(pdbfile,sopc,1)

    if args.attype and args.aa_pdb:
        topfile = 'gromacs.top'
        atomtypes = int(args.attype)
        if int(args.attype)==2:
            nativefile='native_cb.pdb'
            X.write_CB_to_native(pdbfile,sopc,2)
            X.write_gromacs_top(topfile,int(atomtypes),pdbfile,nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii,excl_rule)
            # else:
            #     print ('need --CBcom argument.')
        if int(args.attype)==1:
            U=Utils()
            if U.file_exists('native_ca.pdb'):
                nativefile='native_ca.pdb'
                X.write_gromacs_top(topfile, int(atomtypes), pdbfile, nativefile,CA_rad,sopc,btparams,Ka,Kb,Kd,cutoff,CBcom,CBradii,excl_rule)
            else:
                print ('Need native_ca.pdb file. Just grep "CA" pdbfile.')
        if args.pdbgro:
            #nativefile='native_cb.pdb'
            grofile='gromacs.gro'
            #X.write_gro_gro(grofile,pdbfile,atomtypes,sopc)
        if args.pl_map:
            l=np.loadtxt('contacts.txt')
            xi=l[:,1]
            yi=l[:,3]
            Y.plot_map(xi,yi,'conmap','Res2','Res1')
        U.make_dir('PATH');U.make_dir('MD')
        U.make_dir_struc('PATH','MD')
if __name__ == '__main__':
    main()