from __future__ import print_function
from __future__ import print_function
from __future__ import division
#matplotlib inline!
from Bio.PDB import *
from util import Utils
#from plots import plot_1
import mdtraj as md
import numpy as np
import itertools
import matplotlib.pyplot as plt
from collections import Counter
from Bio.PDB.PDBParser import PDBParser
from scipy.spatial.distance import squareform
__author__ = "Sridhar N"
__email__ = "sridharn@ncbs.res.in"

class conmaps():
    #examples of extracting data from trajectories.
    #python gr.py --native native_ca.pdb --traj traj.xtc --
    #No. of atoms read =   92
    #No. of residues read = 92
    #No. of frames (in ns)= 3000
    #Writing xyz file to : native_ca.pdb.xyz

        def get_trajfromxtc(self,xtcfile,pdbfile):
            # mdtraj.load_xtc(filename, top=None, stride=None, atom_indices=None, frame=None)
            return md.load_xtc(xtcfile,top=pdbfile,stride=1)
        def read_xyz(self,xyzfile):
            return md.formats.XYZTrajectoryFile(xyzfile, mode='r').read()
        def read_xyz_as_traj(self,xyzfile,native):
            print ('>>read_xyz_as_traj:',xyzfile,native)
            topology=md.load(native).topology
            return md.formats.XYZTrajectoryFile(xyzfile, mode='r').read_as_traj(topology=topology)
        def download_pdb_from_id(self, pdbid):
            pdb1 = PDBList()
            pdb1.retrieve_pdb_file(pdbid, pdir='PDB')
        def get_angles(self,pdbfile,triplets):
            traj=md.load(pdbfile)
            angles=md.compute_angles(traj,triplets,periodic=False,opt=True)
            return angles
        def get_dihedrals(self,pdbfile,quadruplets):
            traj=md.load(pdbfile)
            return md.compute_dihedrals(traj,quadruplets,periodic=False,opt=True)
        def get_atom_name(self,topology,index):
            return ([atom.name for atom in topology.atoms if (atom.index==index)])
        def get_residue_number(self,topology,index):
            return ([atom.residue.index for atom in topology.atoms if (atom.index==index)])
        def get_atoms_in_residue(self,topology,index):
            return ([atom.index for atom in topology.atoms if atom.residue.index==index])
        def get_CA_atoms_in_residue(self,topology,index):
            return [atom.index for atom in topology.atoms if ((atom.name == 'CA') and (atom.residue.index == index))]
        def get_CB_atoms_in_residue(self,topology,index):
            return [atom.index for atom in topology.atoms if ((atom.name == 'CB') and (atom.residue.index == index))]
        def get_bb_atoms_in_residue(self,topology,index):
            return ([atom.index for atom in topology.atoms if atom.residue.index == index and atom.is_backbone])
        def get_sc_atoms_in_residue(self,topology,index):
            return ([atom.index for atom in topology.atoms if atom.residue.index == index and atom.is_sidechain])
        def get_sc_heavy_atoms_in_residue(self,topology,index):
#            topology=topology.select_atom_indices('heavy')
            return ([atom.index for atom in topology.atoms if atom.residue.index == index and atom.is_sidechain])
        def get_CB_pairs(self,pdbfile):
            #returns all possible pairs of beta carbons!
            index=self.get_CB_index(pdbfile)
            return list(itertools.combinations(index,2))
        def get_residue_name(self, pdbfile, atom_index):
            topology = md.load(pdbfile).topology
            name = ([atom.residue.name for atom in topology.atoms if (atom.index == atom_index)])
            return name
        def get_res_name_fr_idx(self,topology,index):
            name = ([residue.name for residue in topology.residues if (residue.index == index)])
            return name
        def get_total_number_of_atoms(self,pdbfile):
            topology=md.load(pdbfile).topology
            return topology.n_atoms
        def get_total_number_of_residues(self,pdbfile):
            topology=md.load(pdbfile)
            return topology.n_residues
        def get_all_atom_indices(self,pdbfile):
            traj = md.load(pdbfile)
            topology = traj.topology
            number = ([atom.index for atom in topology.atoms ])
            return number
        def get_backbone_atoms(self,pdbfile):
            #partition into backbone and side-chain.
            traj = md.load(pdbfile)
            topology=traj.topology
            backbone = ([atom.index for atom in topology.atoms if (atom.is_backbone)])
            return backbone
        def get_sidechain_atoms(self,pdbfile):
            traj = md.load(pdbfile)
            topology = traj.topology
            sidechain = ([atom.index for atom in topology.atoms if (atom.is_sidechain)])
            return sidechain
        def get_coordinates(self,pdbfile):
            traj=md.load_pdb(pdbfile)
            return traj.xyz[0]
        def get_CB_coordinates(self,pdbfile):
            traj=md.load(pdbfile)
            topology=traj.topology
            CB_index=([atom.index for atom in topology.atoms if (atom.name == 'CB')])
            CB_coords = []
            for i in CB_index:
                CB_coords.append(traj.xyz[0, i] * 10)
            return CB_coords
        def get_CA_coordinates(self,pdbfile):
            traj=md.load_pdb(pdbfile)
            topology = traj.topology
            #selection = topology.select_expression('name CA')
            CA_index=([atom.index for atom in topology.atoms if (atom.name == 'CA')])
            CA_coords=[]
            for i in CA_index:
                CA_coords.append(traj.xyz[0,i]*10)
            return CA_coords
        def get_pdb_coords(self,pdbfile):
            traj = md.load_pdb(pdbfile)
            topology = traj.topology
            index = ([atom.index for atom in topology.atoms])
            coords = []
            for i in index:
                coords.append(traj.xyz[0, i] * 10)
            return coords
        def get_CA_index(self,pdbfile):
            traj = md.load_pdb(pdbfile)
            topology = traj.topology  # selection = topology.select_expression('name CA')
            CA_index = ([atom.index for atom in topology.atoms if (atom.name == 'CA')])
            return CA_index
        def get_CB_index(self,pdbfile):
            traj = md.load_pdb(pdbfile)
            topology = traj.topology  # selection = topology.select_expression('name CA')
            CA_index = ([atom.index for atom in topology.atoms if (atom.name == 'CB')])
            return CA_index

        def check_glycineh(self,pdbfile,skiph):
            print ('>>in check_glycineh',skiph)
            traj = md.load_pdb(pdbfile); topology = traj.topology; seq=self.get_sequence(pdbfile)
            print (seq)
            gly = [pos + 1 for pos, char in enumerate(seq) if char == 'G'][0]
            print ('Glycine in residue:',gly)
            print(self.get_sc_atoms_in_residue(topology,gly))
            print (len(self.get_sc_atoms_in_residue(topology,gly)))

            if not skiph:
                 assert(len(self.get_sc_atoms_in_residue(topology,gly))>1),'No Hydrogens found on Glycine. Use --skip_glycine command to proceed.'
            return

        def get_btmap_val(self,res1,res2,mapfile):
            #returns btvalue for two residues (single alphabet notation)
            #print ('>>reading interaction parameters from file:',mapfile)
            U=Utils()
            btmap = U.get_dict_from_csv(mapfile)
            j = res1 + res2
            #print (j)
            if j not in btmap.keys():
                r = j[1] + j[0]
                assert r in btmap.keys(), "Residue pair not in file!"
                #print (r)
                return btmap[r]
            else:
                return btmap[j]

        def get_epsij_from_ext(self,contactfile,pdbfile):
            #read external contact file and generate epsij matrix.
            U=Utils()
            pairs=self.get_pairs_ext(contactfile)
            d1=U.amino_acid_3_1()
            eps=[]
            for i in pairs:
                res1=self.get_residue_name(pdbfile,i[0]);res2=self.get_residue_name(pdbfile,i[1])
                res1= (d1[res1[0]]);res2=(d1[res2[0]])
                eps+=[self.get_btmap_val(res1,res2,'mjmap.dat')]
            print (eps)
            return eps
        def get_distances(self,traj,pairs):
            return md.compute_distances(traj, pairs, periodic=False, opt=True)

        def get_bb_contacts_SOPC(self,pdbfile,nativefile,cutoff,scaling,separation):
            #SOP-SC applies cut-off criterion to coarse-grained beads. Not to the all-atom PDB.
            #Slightly different from Cheung-Thirumalai handling.
            #NATIVE INTERACTIONS.
            print(">> in bb_contacts (SOP-SC)", nativefile, "cutoff =", cutoff)
            cutoff = float(cutoff); scaling = float(scaling)
            topology = md.load('native_cb.pdb').topology
            traj = md.load('native_cb.pdb');natoms=traj.n_atoms
            seq=self.get_sequence(pdbfile)
            contacts = []
            #Matches paper results when only CA beads are taken into account. As done below!
            #Pairs of residues
            a = np.arange(0,len(seq)); b = list(itertools.combinations(a, 2))
            #Paira of atoms
            a1 = np.arange(0,natoms);b1=list(itertools.combinations(a1,2))
            #Loop thrice only for clarity. Should loop just once.
            count=0
            print ('Searching',len(b),'pairs')
            for i in b:
                if np.absolute(i[1]-i[0])>=separation:
                    #CA-CA
                    x1 = self.get_CA_atoms_in_residue(topology,i[0]);y1 = self.get_CA_atoms_in_residue(topology, i[1])
                    #print (x1,y1,i[0],i[1])
                    assert (len(x1)>=1) and (len(y1)>=1)
                    pair = ([x1, y1]);pair = np.reshape(pair, (1, 2))
                    distance = md.compute_distances(traj, pair)[0][0]
                    # cutoff and scaling criterion.
                    if (distance * 10 <= cutoff):
                        #print (x1,y1,distance*10)
                        contacts.append([i[0], i[1],distance*10, distance* 10])
                        count=count+1
            print ("Found",count,"bb-bb contacts")
            #np.savetxt ('bb_sopsc.dat',contacts,fmt='%-4d  %-4d  %10.6f  %10.6f')
            return contacts

        def get_ss_contacts_SOPC(self,pdbfile, nativefile, cutoff, scaling, separation):
            cutoff = float(cutoff);scaling = float(scaling)
            topology = md.load('native_cb.pdb').topology; traj = md.load('native_cb.pdb');
            natoms = traj.n_atoms; seq = self.get_sequence(pdbfile);  contacts = []
            a = np.arange(0, len(seq)); b = list(itertools.combinations(a, 2));count=0
            for i in b:
                x1 = self.get_CB_atoms_in_residue(topology, i[0]);y1 = self.get_CB_atoms_in_residue(topology, i[1])
                #if not x1 or not y1:
                 #   print ('Skipping',x1,y1,self.get_res_name_fr_idx(topology,i[0]),self.get_res_name_fr_idx(topology,i[1]))
                if np.absolute(i[1] - i[0]) >= separation and (x1 and y1):
                    # print (x1,y1,i[0],i[1])
                    #assert (len(x1) >= 1) and (len(y1) >= 1)
                    pair = ([x1, y1]);
                    pair = np.reshape(pair, (1, 2))
                    distance = md.compute_distances(traj, pair)[0][0]
                    # cutoff and scaling criterion.
                    if (distance * 10 <= 8.0):
                        # print (x1,y1,distance*10)
                        contacts.append([i[0], i[1], distance * 10, distance * 10])
                        count = count + 1
            print("Found", count, "sc-sc contacts")
            np.savetxt('scsc_sopsc.dat', contacts, fmt='%-4d  %-4d  %10.6f  %10.6f')
            return contacts

        def get_bs_contacts_SOPC(self,pdbfile, nativefile, cutoff, scaling, separation):
            cutoff = float(cutoff); scaling = float(scaling);topology = md.load('native_cb.pdb').topology; traj = md.load('native_cb.pdb');
            natoms = traj.n_atoms; seq = self.get_sequence(pdbfile);contacts = []; a = np.arange(0, len(seq)); b = list(itertools.combinations(a, 2))
            #pairs of residues
            count=0
            for i in b:
                x1 = self.get_CB_atoms_in_residue(topology, i[0]);x2=self.get_CA_atoms_in_residue(topology,i[0])
                y1 = self.get_CA_atoms_in_residue(topology, i[1]);y2=self.get_CB_atoms_in_residue(topology,i[1])
                assert len(y1)==1 and len(x2)==1;
                if x1 and y1:
                    if np.absolute(i[1]-i[0])>=separation:
                       # CB-CA
                        pair = ([i[0], i[1]]);
                        pair = np.reshape(pair, (1, 2))
                        distance = md.compute_distances(traj, pair)[0][0]
                        #cutoff and scaling criterion.
                        if (distance * 10 <= 8.0):
                            #print(pair, distance, x1, y1)
                            contacts.append([i[0], i[1], distance * 10, distance * 10])
                            count = count + 1
                # if x2 and y2:
                #      if np.absolute(i[1]-i[0])>=separation:
                #          #CA-CB
                #          pair = ([i[0], i[1]]);
                #          pair = np.reshape(pair, (1, 2))
                #          distance = md.compute_distances(traj, pair)[0][0]
                #          # cutoff and scaling criterion.
                #          if (distance * 10 <= 8.0):
                #              # print(pair, distance, x1, y1)
                #              contacts.append([i[0], i[1], distance * 10, distance * 10])
                #              count = count + 1
            print("Found", count, "bb-sc contacts. See SOP-SC.txt")
            np.savetxt('bb_sc_sopc.dat', contacts,fmt='%-4d  %-4d  %10.6f  %10.6f')
            contacts = np.asarray(contacts)
            return contacts


        def get_backbone_contacts(self,pdbfile,cutoff,scaling,separation):
            #partition all-atom pdb into backbone and side-chain components.
            print(">> in backbone_contacts", pdbfile, "cutoff =",cutoff, "scaling=",scaling )
            cutoff = float(cutoff);scaling=float(scaling)
            topology = md.load(pdbfile).topology; seq=self.get_sequence(pdbfile)
            traj=md.load(pdbfile)
            #all possible residue pairs
            a=np.arange(0,len(seq));b=list(itertools.combinations(a,2))
            #CA-CA-residues
            contacts=[]
            for i in b:
                if i[1]-i[0]>=separation:
                    #get all backbone atoms from all-atom pdb file
                    x1=self.get_bb_atoms_in_residue(topology,i[0]);y1=self.get_bb_atoms_in_residue(topology,i[1])
                    #compute distances between backbone atoms in all-atom pdb file.
                    pairs=list(itertools.product(x1,y1)); distance=md.compute_distances(traj,pairs)
                    #Compute distance between corresponding CA atoms.
                    CA1=self.get_CA_atoms_in_residue(topology, i[0]);CA2=self.get_CA_atoms_in_residue(topology, i[1])
                    # print (i[0],i[1],CA1,CA2)
                    # print (pairs)
                    pair = ([CA1,CA2]);pair = np.reshape(pair, (1, 2))
                    ca_dist = md.compute_distances(traj, pair, periodic=False, opt=True)
                    #cutoff and scaling criterion.
                    if (distance*10<cutoff).any():
                        if np.amin(distance)*10<scaling*ca_dist*10:
                            #print (distance,scaling,ca_dist*10)
                            #exit()
                            contacts.append([i[0],i[1],np.amin(distance)*10,ca_dist*10])
            contacts=np.asarray(contacts)
            np.savetxt('bb.t', contacts)
            return contacts

        def get_side_chain_contacts(self,pdbfile,cutoff,scaling,separation):
            print(">> in side_chain_contacts", pdbfile,"cutoff =",cutoff,"scaling=",scaling)
            cutoff = float(cutoff);scaling=float(scaling)
            #scaling=1.2
            traj = md.load(pdbfile)
            topology = traj.topology
            top1=md.load('native_cb.pdb').topology
            traj1=md.load('native_cb.pdb')
            seq = self.get_sequence(pdbfile)
            # all possible residue pairs
            a = np.arange(0, len(seq) - 1)
            b = list(itertools.combinations(a, 2))
            contacts = []
            #separation=3
            for i in b:
                if i[1]-i[0]>=separation and i:
                    x1=self.get_sc_atoms_in_residue(topology,i[0])
                    y1=self.get_sc_atoms_in_residue(topology,i[1])
                    if x1 and y1:
                        #get CB atoms from coarse-grained pdbfile
                        CB1 = self.get_CB_atoms_in_residue(top1, i[0])
                        CB2 = self.get_CB_atoms_in_residue(top1, i[1])
                        #assert (len(CB1)==len(CB2))
                        if CB1 and CB2:
                            pair = ([CB1, CB2]);
                            pair = np.reshape(pair, (1, 2))
                            cb_dist = md.compute_distances(traj1, pair, periodic=False, opt=True)
                            #print (cb_dist,pair,x1,y1)
                            pairs=list(itertools.product(x1,y1))
                            distance=md.compute_distances(traj,pairs)

                            #cutoff and scaling criterion
                            if (distance*10<cutoff).any():
                                if np.amin(distance)*10<scaling*cb_dist*10:
                                    contacts.append([i[0],i[1],np.amin(distance)*10,float(cb_dist)*10])
                                    #print (i[0],i[1],np.amin(distance)*10,float(cb_dist))
            contacts = np.asarray(contacts)

            np.savetxt('sidechain.t', contacts)
            return contacts

        def get_bb_sc_contacts(self,pdbfile,cutoff,scaling,separation):
            print(">> in bb_sc_contacts", pdbfile, "cutoff =", cutoff, "separation=",separation)
            cutoff = float(cutoff);scaling=float(scaling);separation=int(separation)
            traj = md.load(pdbfile)
            topology = traj.topology #all_atom pdb topology
            top_tb = md.load('native_cb.pdb').topology #native_all_topology
            traj_tb=md.load('native_cb.pdb')
            #top_ca = md.load('native_ca.pdb').topology #native_ca_topology
            seq = self.get_sequence(pdbfile)
            a = np.arange(0, len(seq) - 1)
            b = list(itertools.combinations(a, 2))
            contacts = []
            for i in b:
                if i[1] - i[0] >= separation:
                    x1 = self.get_bb_atoms_in_residue(topology,i[0]);y1=self.get_sc_atoms_in_residue(topology,i[1])
                    if x1 and y1:
                        pairs = list(itertools.product(x1, y1)); distance = md.compute_distances(traj, pairs)
                        #pairs from all atom pdb file!
                        CA1 = self.get_CA_atoms_in_residue(top_tb, i[0]);CB1 = self.get_CB_atoms_in_residue(top_tb, i[1])
                        if CA1 and CB1:
                            pair = ([CA1, CB1]); pair = np.reshape(pair, (1, 2))
                            ca_dist = md.compute_distances(traj_tb, pair, periodic=False, opt=True)
                            if (distance * 10 < cutoff).any():
                                if np.amin(distance)*10 < scaling * ca_dist * 10:
                                    contacts.append([i[0], i[1], np.amin(distance) * 10, ca_dist * 10])
            contacts = np.asarray(contacts)
            np.savetxt('bb_sc.t', contacts)
            return contacts


        def all_atom_contacts(self,pdbfile,cutoff,scaling):
            #Generate cut-off based contact map from pdb.:pdb_contacts
            #contacts.txt file as external pair to SMOG.
            print (">> in all_atom_contacts",pdbfile)
            cutoff=float(cutoff)*0.1
            #returns an Nx4 matrix containing resno.1 resno.2 dist_atom_12 dist_CA_12
            traj= md.load_pdb(pdbfile)
            topology = traj.topology
            #schemes can be ca closest closest-heavy
            print('Atoms in pdb    %s' % traj.n_atoms)
            print('Residues in pdb %s' % traj.n_residues)
            contacts = md.compute_contacts(traj, contacts='all', scheme='closest', ignore_nonprotein=False)
            #find contacts less than cutoff.
            indices=np.where((contacts[0]<cutoff))
            map=[]
            dist=[]
            ca_dist=[]
            #For contacts that match criterion, get residues and distances.
            for i in indices[1]:
#                map.append(contacts[1][i])
                resid1=contacts[1][i][0]
                resid2 = contacts[1][i][1]
                if (resid2 - resid1) > 3:
                    map.append(contacts[1][i])
                    #select alpha carbons in residue and compute CA-CA distances.
                    selection1=[atom.index for atom in topology.atoms if ((atom.name == 'CA') and (atom.residue.index==resid1) )]
                    selection2 = [atom.index for atom in topology.atoms if
                              ((atom.name == 'CA') and (atom.residue.index == resid2))]
                    pair=([selection1,selection2])
                    #print (pair)
                    pair=np.reshape(pair,(1,2))
                    #ca_dist is list of ca_distance between contacts that are within cut-off
                    ca_dist.append(md.compute_distances(traj,pair, periodic=False, opt=True))
                    #dist is list of atom-atom distance between contacts
                    dist.append(contacts[0][0][i])
                else:
                    continue
            map=np.array(map)
            #print (map.shape)
            #atom-atom distances between contacts
            dist=np.asarray(dist,dtype=np.float32)
            ca_dist=np.asarray(ca_dist,dtype=np.float32)
#           reshape distance array to stack  with residue numbers (and anything else one might want to pass on in contact map.
            ca_dist=np.reshape(ca_dist,(ca_dist.size,1))
            dist=np.reshape(dist,(dist.size,1))
            assert (len(dist)>0),'No contacts made. Check cutoff. '
            #print (map.shape, dist.shape, ca_dist.shape)
            conmap=np.hstack((map,dist,ca_dist))
            #conmap returns res_i res_j dist_ij
            print ('Found',len(map),'contacts in all-atom file',pdbfile,'with cut-off', cutoff, 'nm')

            save=True
            if save:
                np.savetxt('pdb_contacts',conmap,fmt='%f')
                map = np.array([x + 1 for x in map])
                ch1=np.ones(len(map))
                savecon=np.vstack((ch1,map[:,0],ch1,map[:,1]))
                savecon.reshape(len(map),4)
                np.savetxt( 'contacts.txt', savecon.T, fmt='%d')
                print('See file: contacts.txt for reading into other codes. \n See file: pdb_contacts for details.')
            #conmap contains the following:
            #resno1 resno2 atom-atom-dist ca-ca-dist
            return conmap

        def get_pairs_ext(self,filename):
            print ('>>in get_pairs_ext',filename)
            #read contact map generated from SMOG. Format : 1 i 1 j
            #check numbering of residues. Always start from 0 inside code and 1 outside code!!!!
            U=Utils()
            assert(U.file_exists(filename))
            pairs=np.loadtxt(filename,unpack=True)
            pairs=pairs.astype(int)
            #print (pairs)
            #choose column 2 and column 4. Change here for other columns.
            pairs=(pairs[[1, 3], :])
            if len(pairs[0]) != len(pairs[1]):
                U.fatal_errors(2)
                exit()
            print('No. of contacts read from file', filename, '=',len(pairs[0]))
            pairs=pairs.T - 1 #0 is first residue.

            return pairs

        def get_dihedral_energies(self,xyzfile):
            print ('>>in get_dihedral_energies')
            # calculate different energies for a given xyz trajectory file.
            # one single-point energy OPTIM run for each frame
            #xyzfile is file with lots of frames in xyz one after another
            U=Utils()
            import subprocess
            traj=self.read_xyz_as_traj(xyzfile,native='native_ca.pdb')
            print('Reading', len(traj), ' frames from xyz file')
            print('odata file is written for each frame')
            files='OPTIM SBM.INP finish'
            #check if all files needed are in folder
            for file in files.split():
                if U.file_exists(file):
                    print ('Found', file)
                else:
                    print ('Files missing. STOP')
                    exit()
            #remove useles files:
            #removefiles='points points.final odata.new'
            #os.remove(removefiles)
            #write odata file
            # loop over frames in xyz trajectory
            frame_count=0
            frame_count=0
            for j in xrange(0, len(traj)):
                frame_count = frame_count + 1
                xyz=traj[j].xyz
                #print("Writing odata file:")
                f = open('odata', 'w')
                f.write('STEPS 1\n')  # python will convert \n to os.linesep
                f.write('BFGSMIN 0.000001\n')  # python will convert \n to os.linesep
                f.write('POINTS\n')  # python will convert \n to os.linesep
                #write odata file
                for k in xrange(0,len(xyz[0])):
                    #NOTE: ,mdtraj works in nm and OPTIM in angstroms!
                    x=traj[j].xyz[0][k][0]*10
                    y = traj[j].xyz[0][k][1]*10
                    z = traj[j].xyz[0][k][2]*10
                    f.write('%s %10.6f %10.6f %10.6f\n' % ('SB',x,y,z))

                f.close()
                proc = subprocess.Popen(["OPTIM"], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
                print ("Frame", frame_count)
                f1=open('energies','a')
                for line in iter(proc.stdout.readline, ''):
                    #extract energy contributions. Remember each is summed up so subtract in end.
                    if 'Bond' in line.split():
                        e_bond=line.split()
                        e_bond=float(e_bond[1])
                    elif 'Angle' in line.split():
                        e_angle=float(line.split()[1]) - e_bond
                    elif 'Dihedral' in line.split():
                        e_dih=float(line.split()[1])-e_angle-e_bond
                    elif 'Contact' in line.split():
                        e_con=float(line.split()[1])-e_angle-e_bond-e_dih
                    elif 'force=' in line.split():
                        #print (line.split())
                        e_tot=float(line.split()[5])

                f1.write('%10.6f %10.6f %10.6f %10.6f %10.6f\n' % (e_bond,e_angle,e_dih,e_con,e_tot) )
            f1.close()
            return


        def get_contact_probability(self,traj,external_pair,start,finish,step):
            print (">>in get_contact_probability")
            # returns array containing residue pair, occurrence proability in trajectory
            # Recommended to alway
            # s give external_pair file.
            #print (traj)
            frames = np.arange(start, finish, step, dtype=int)
            #traj = self.get_trajfromxtc(xtcfile, pdbfile)

            #traj = traj.slice(frames, copy=True)
            # Give contacts externally with external_pair option.
            atom_pairs = self.get_pairs_ext(external_pair)
            #print ('Warning: First frame in file is assumed to be the native state.')
            #read xyz as traj
            #compute same pair distances in xyz trajectory
            dist=md.compute_distances(traj,atom_pairs,periodic=False,opt=True)
            dist_native=dist[0]
            print ('Found',len(dist[0]),' pair distances for',len(dist),'frames')
            frame_count=0
            scaling=1.2
            #loop over frames
            pair_count=[]
            for i in xrange(0,len(dist)):
                #loop over ALL contacts:
                #change here for specific contacts
                count=0
                for j in xrange(0,len(dist[0])):
                    if (dist[i][j] <= scaling* dist_native[j]):
                        pair_count.append([i, atom_pairs[j][0], atom_pairs[j][1]])
                        count=count+1
                    #else:
                        #save all missing pairs as frame, resid1,resid2
                frame_count=frame_count+1
#                print (frame_count, count)

            #pair count contains frame info too, Use it if necessary
            print (len(pair_count), 'pairs in total trajectory')
            pair_count=np.asarray(pair_count)
            resid1=pair_count[:,[1]].T[0].tolist()
            resid2 = pair_count[:, [2]].T[0].tolist()
            #count the occurrence of each residue pair in all formed contacts?
            #Counter returns dictionary of residue pair: Occurences.
            occurences=Counter(zip(resid1,resid2))
            #convert occurrences to np array and save file.
            occurences.keys()
            t1=[]
            for key,value in occurences.items():
                t1.append([key,value])
            t1=np.asarray(sorted(t1))
            scale=float(1/len(dist))
            #print (scale)
            t1[:,1]*=scale
            #print(t1)
            np.savetxt('contact_count_tr', t1, fmt='%s %s ')
            print ("See file: contact_count_tr")
            return t1


        def get_native_contacts(self,traj,native,aa_pdb,cutoff,scaling,Qfile,Qcut,pairs):
            #made ext_pair mandtory
            U=Utils()
            print (">>in get_native_contacts:", traj,native,aa_pdb,cutoff,scaling,Qfile,Qcut)
            scaling=float(scaling)
            traj_pdb= md.load_pdb(native)
            #atom pairs from external file:
            #assert(U.file_exists(pairfile))
            #remember counting starts from 0 in python!!!!
            #print ('Will not generate contact-map. Contacts in', pairfile)
            #print ('File name=', pairfile,'found. This will be used as contact-map. ')
            atom_pairs=pairs
            #compute distances in native-CA file.
            dist_pdb= md.compute_distances(traj_pdb,atom_pairs,periodic=False,opt=True)
            dist_pdb=dist_pdb[0]
            #n_contacts=len(atom_pairs)
            # get distances for all pairs in all frames from trajfile
            contacts_xtc = md.compute_distances(traj, atom_pairs, periodic=False, opt=True)
            frame_count = 0; contact_count = []
            # loop over frames in traj
            for j in xrange(0, len(traj)):
                frame_count = frame_count + 1
                # select contact distances in frame j
                dist_xyz = contacts_xtc[j]
                count = 0
                #loop over pairs
                for i in xrange(0, len(dist_xyz)):
                    if (dist_xyz[i] <= scaling * dist_pdb[i]):
                        #print (dist_xyz[i],dist_pdb[i],scaling,i)
                        count = count + 1
                contact_count.append(count)
            contact_count = np.asarray(contact_count)
            # max_contacts=np.max(contact_count)
            max = (np.amax(contact_count))
            fraction = np.divide(contact_count, max, dtype=float)
            print(fraction)
            np.savetxt(Qfile, (fraction), fmt='%10.6f')
            np.savetxt(Qfile+'count.test', (contact_count), fmt='%d')
            print('See file:',Qfile,Qfile+'count.test','cutoff=',cutoff,'scaling=',scaling)
            # returns count of contacts mapped per frame. contact_count[frame]=number of contacts formed.
            return contact_count



        def xyz_native_contacts(self,xyzfile,native,aa_pdb,cutoff,scaling):
            #DEPRECATED: Use get_native_contacts always.
            traj_xyz=self.read_xyz_as_traj(xyzfile,native)
            traj_native=md.load_pdb(native)
            print ('Reading',len(traj_xyz),' frames from xyz file')
            allatom_pdb=aa_pdb
            #topology = md.load(allatom_pdb).topology
            #get contact map for all atom file.
            #add here if external contact map is being supplied.
            pdb_conmap=self.all_atom_contacts(allatom_pdb, cutoff,scaling)
            print (pdb_conmap[:,3])
            #atom_pairs=pdb_conmap[:,[0,1]]
            #if not Utils.file_exists('external_pair'):
             #   print ('Fatal: external_pair file not found.')
             #   exit()
            ext_pair=True
            if ext_pair:
                atom_pairs=self.get_pairs_ext('external_pair')
                print (len(atom_pairs))
                #print (atom_pairs)
                contacts_native=md.compute_distances(traj_native,atom_pairs, periodic=False, opt=True)[0]
            else:
                U.fatal_errors(3)
                exit()
            # get contact distances for all frames in xyz trajectory
            contacts_xyz =  md.compute_distances(traj_xyz,atom_pairs, periodic=False, opt=True)

            #   print (contacts_xyz[0][5],pdb_conmap[5][3],pdb_conmap[5][0],pdb_conmap[5][1],atom_pairs[5])
            frame_count=0
            contact_count=[]
            #loop over frames in xyz trajector
            for j in xrange(0,len(traj_xyz)):
                frame_count=frame_count+1
                #select contact distances in frame j
                dist_xyz=contacts_xyz[j]
                #quit()
                count=0
                scaling = float(scaling)
                for i in xrange (0,len(dist_xyz)):
#                    print (dist_xyz[i], pdb_conmap[i][2])
                    if (dist_xyz[i]<=scaling*contacts_native[i]):
                        #print (dist_xyz[i],contacts_native[i],scaling)
                        count=count+1
                contact_count.append(count)

            contact_count=np.asarray(contact_count)
            #max_contacts=np.max(contact_count)
            max = (np.amax(contact_count))
            fraction=np.divide(contact_count,max,dtype=float)
#            print (fraction)
            print (contact_count)
            np.savetxt('Qmap', (fraction), fmt='%10.6f')
            print ('See file: Qmap')
            #returns count of contact mapped per frame. contact_count[frame]=number of contacts formed.
            return contact_count

        def write_as_xyz(self,xtcfile,pdbfile,start,finish,step,xyzfile):
            # nanometer, degree and picosecond.
            #write trajectory as xyz for pathsample readable file
            #traj=self.get_trajfromxtc(xtcfile,pdbfile)
            frames = np.arange(start, finish, step, dtype=int)
            traj = self.get_trajfromxtc(xtcfile, pdbfile)
            traj=traj.slice(frames, copy=True)
            print('No. of atoms read =   %s' % traj.n_atoms)
            print('No. of residues read = %s' % traj.n_residues)
            print('No. of frames (in ns)= %s' % traj.n_frames)
            print  ('Writing xyz file to : %s' % xyzfile)
            traj.center_coordinates().save_xyz((xyzfile),force_overwrite=True)
            return

        def get_dssp(self,traj):
            print ("Need dssp in folder.")
            return md.compute_dssp(traj)

        def get_SOPC_non_native(self,pdbfile,nativefile,epsl):
            natoms=self.get_total_number_of_atoms(nativefile)
            #list of all possible non-native contacts!
            print (type(natoms ))
            return 1

        def get_SOPC_non_native_angles(self,pdbfile,nativefile,CA_rad,epsl):
            print ('>> in get_sopc_non_native_angles')
            #SOP-SC model has weird non-native non-bonded interactions(angular repulsions!).
            #angles formed by triplets. Non-bonded part of triplet is written in pairs section.
            #All CA is odd number. All CB is even number.
            sigmabb=float(CA_rad)
            topology=md.load_pdb('native_cb.pdb').topology
            assert (self.get_atom_name(topology, 0)[0] == 'CA'),'First atom must be C-alpha in native file!'
            #select only backbone atoms
            bb=self.get_CA_index('native_cb.pdb'); contacts=[];count=0
            for i in bb[:-2]:
                count=count+1; assert(self.get_atom_name(topology,int(i))[0]=='CA')
                contacts.append([i, i+4,sigmabb, sigmabb ])
            print ('Found',count,'non-native,non-bonded bb-bb angular interactions')
            #select CB atoms
            count=0;sc=self.get_CB_index('native_cb.pdb')
            #get side-chain radius.
            U=Utils(); d=U.amino_acid_radius_dict()
            sc.pop(0) #remove first element from list.

            for i in sc[:-1] :
                atom_index=i
                #assign radii to each side-chain
                name = ([atom.residue.name for atom in topology.atoms if (atom.index == atom_index)])
                #print (d[name[0]],i,name)
                sigmasc=float(d[name[0]])
                sigmabs=0.4*(sigmasc+sigmabb)
                #print (sigmabs,sigmabb,sigmasc)
                count=count+2
                contacts.append([i,i+1,sigmabs,sigmabs]) #contacttype 1=LJ for GROMACS and SBM.#FTW
                contacts.append([i, i-3,sigmabs,sigmabs])

            print('Found', count, 'non-native,non-bonded bb-sc angular interactions')
            #np.savetxt('nn_angles_sopsc.dat', contacts, fmt='%-4d  %-4d  %10.6f  %10.6f')
            #must retturn residue numbers, and not atom numbers!

            return contacts

        def hierachical_cluster_rmsd(self,traj,pdbfile,start,finish,step):
            import scipy.cluster.hierarchy
            #get rmsds from xtc file and cluster hierarchically.
#            traj=self.get_trajfromxtc(trajectory,pdbfile)
            frames = np.arange(start, finish, step, dtype=int)
            traj = traj.slice(frames, copy=True)
            distances = np.empty((traj.n_frames, traj.n_frames))
            #distances=np.zeros((traj.n_frames, traj.n_frames))
            for i in xrange(traj.n_frames):
                distances[i] = md.rmsd(traj, traj, i)
                #print('Max pairwise rmsd: %f nm' % np.max(distances))
             #cluster with scipy linkage algorithm. takes only squareform.
            #assert np.all(distances-distances.T<1e-6)
            reduced_distances = squareform(distances, checks=False)
            print (reduced_distances)
            linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method='average')
            plt.title('RMSD Average linkage hierarchical clustering')
            _ = scipy.cluster.hierarchy.dendrogram(linkage, no_labels=True, count_sort='descendent')
            plt.show()

            return reduced_distances

        def check_pdb(self,pdbfile):
            traj = md.load_pdb(pdbfile)
            top=traj.topology
            traj.center_coordinates().save_pdb((pdbfile+'.1'), force_overwrite=True)
            seq=self.get_sequence(pdbfile)
            return pdbfile

        def check_hetatm(self,pdbfile):
            U=Utils()
            with open(pdbfile) as f1:
                a=[word for line in f1 for word in line.split()]
                if 'HETATM' in a:
                    U.fatal_errors(10)
            f1.close()
            return

        def pdb_read(self,pdbfile):
            #read pdb file.
            #returns X,Y,Z,res_name,res_nr. (add stuff if needed)
            x = [];y = [];z = []
            res_name = []
            res_nr=[]
            Nr = 0
            f = open(pdbfile, 'r')
            while 1:
                line = f.readline()
                if not line: break
                if line[0:6] == 'ATOM  ':
                    rx = float(line[30:38]);ry = float(line[38:46]);rz = float(line[46:54])
                    r1 = float(line[22:26])
                    if line[21] == 'A':
                        x.append(rx);
                        y.append(ry);
                        z.append(rz)
                        res_nr.append(r1)
                        Nr = Nr + 1
                        res_name.append(line[17:20])
            return (rx,ry,rz,res_name,res_nr)

        def pdb_renumber(self,pdbfile):
            print (">>in pdb_renumber")
            #returns value to be shifted by
            x=self.pdb_read(pdbfile)
            res_nr=np.array(x[4],dtype=int)
            min=np.amin(res_nr)
            if min !=1:
                c=min-1
                return c
            return c


        def get_hydrophobic(self, pdbfile):
            seq = self.get_sequence(pdbfile)
            #from Chan PNAS 2010 and wikipedia! .
            #alanine, valine, leucine, isoleucine, methionine, tryptophan, phenylalanine, tyrosine
            hydrophobic = {'A','V','L','I','M','F','W','Y'}
            positions = []
            for c in hydrophobic:
                positions += [pos for pos, char in enumerate(seq) if char == c]
            positions = np.sort(positions)
            # print positions
            print ("Found", len(positions), "hydrophobic residues")
            print (positions)
            return positions


        def get_polar(self, pdbfile):
            seq = self.get_sequence(pdbfile)
            polar = {'Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'M'}
            positions = []
            for c in polar:
                positions += [pos for pos, char in enumerate(seq) if char == c]
            positions = np.sort(positions)
            print   ("Found", len(positions), "polar residues")
            print (positions)
            return positions

        def get_charged(self, pdbfile):
            seq = self.get_sequence(pdbfile)
            charged = {'R','H', 'K', 'D', 'E'}
            positions = []
            for c in charged:
                positions += [pos for pos, char in enumerate(seq) if char == c]
            positions = np.sort(positions)
            print ("Found", len(positions), "charged residues", 'in', positions)
            print (positions)
            return positions

        def get_sequence(self, pdbfile):
            # Checks for duplicates, repeat atoms, general faults in pdbfile
            # get sequence as string from pdb
            p = PDBParser(PERMISSIVE=0)
            structure = p.get_structure('test', pdbfile)
            ppb = PPBuilder()
            seq = ''
            for pp in ppb.build_peptides(structure):
                seq += pp.get_sequence()
            #print (seq)
            #print (len(seq))
            return seq

        def plot_map(self,x,y,title,xaxis,yaxis):
            print ('>>in plot_map')
            #Simple x,y plot
            #X, Y are 1D numpy arrays
            maxx=max(x);maxy=max(y)
            plt.rcParams['backend'] = 'TkAgg'
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            plt.xlim(0,max(maxx,maxy))
            plt.ylim(0,max(maxx, maxy))#colors = ['k'] * len(x)
            ax.scatter(x, y ,alpha=0.5,s=10,linewidth=.05)
            plt.xlabel(xaxis)
            plt.ylabel(yaxis)
            #plt.savefig(title+'.png')
            plt.savefig(title + '.pdf')
            plt.rc('font', family='serif',size='20')
            #plt.show()
            print ('See: ',title+'.pdf')
            return 1

        def plot_scatter(self,x,y,title):
            #Simple x,y plot
            #X, Y are 1D numpy arrays
            import pylab
            plt.rcParams['backend'] = 'TkAgg'
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            #colors = ['k'] * len(x)
            ax.scatter(x, y,marker='o',markersize=1)
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.savefig(title+'.png')
            plt.rc('font', family='serif', size='20')
            #plt.show()
            print ('See: ',title+'.png')
            return 1


        def plot_heat_map(self,x,y,z):
            import matplotlib as mpl
            fig, ax = plt.subplots(figsize=(5,5))
            #ax.scatter(x, y, c=z, s=150, marker='<', edgecolor='none')
            colors = ['red', 'blue']
            levels = [0, 1]
            cmap, norm = mpl.colors.from_levels_and_colors(levels=levels, colors=colors, extend='max')
            ax.scatter(x, y, c='gray', s=150, marker='<', edgecolor='none', cmap=cmap, norm=norm)
            plt.show()

        def plot_histogram(self,x,title):
            import matplotlib.pyplot as plt
            #x=np.histogram(x, bins=np.arange(0,1,len(x)), density=True)
            noise = x
            num_bins = 100
            n, bins, _ = plt.hist(noise, num_bins, histtype='step')
            plt.hist(n,bins,normed=True,linestyle=('dashed'),color=('red'),label='pop')
            #!plt.legend(loc=1, ncol=2, borderaxespad=0.2, fontsize=15)
            plt.rc('font', family='serif', size='20')
            plt.style.use('ggplot')
            plt.savefig(title + '.pdf')
            print ("See file:", title+'.pdf')
            #plt.show()

        def plot_2d_hist(self,filename,title):
            from matplotlib.colors import LogNorm
            from matplotlib import cm
            import matplotlib.pyplot as plt
            #filename='rog_energy_top7.dat'
            font = {'family': 'serif',
                    'weight': 'normal',
                    'size': 10}
            plt.rc('font',**font)
            plt.rcParams["axes.linewidth"] = 2
            #plt.xlabel('Q',fontdict=font)

            plt.rc('ytick', labelsize=15)
            fig, ax = plt.subplots(figsize=(4,3))
            data = np.loadtxt(filename, dtype=float)
            x=data[:,0];y=data[:,1]
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.hist2d(x, y, bins=100, norm=LogNorm(),cmap=plt.cm.coolwarm)
            plt.xlabel('Q', fontdict=font)
            plt.colorbar()
            plt.show()

        def get_side_chain_COM(self,pdbfile):
            print (">> in get_side_chain_COM",pdbfile)
            #returns COM of each side-chain as numppy array
            #GLYCINES are skipped.
            # res numbering starts from 1
            p = PDBParser(PERMISSIVE=0)
            structure = p.get_structure('test', pdbfile)
            backbone=['N','CA','C','O']
            COM1=[]
            glyname='GLY111'
            for model in structure:
                for chain in model:
                    count=0
                    for residue in chain:
                        coord = []
                        mass = []
                        COM=[]
                        for atom in residue:
                            #get side chain only
                            if atom.name not in backbone:
                                #get coordinates and atomic mass.
                                coord.append(atom.get_coord())
                                #print len(coord)
                                mass.append(atom._assign_atom_mass())
                                print (atom._assign_atom_mass,atom.name)
                                #mass.append(float(0.00))
                                assert len(mass)!=0, "Zero length mass array!"
                                #print atom.get_parent()
                        #get COM for every
                        coord=np.array(coord)
                        mass=np.array(mass)
                        res=residue.get_resname().strip()
                        #Make sure GLYCINE is excluded.
                        if res=='GLY':
                            print ("Excluding glycine, res = ",count+1)
                        else:
                            assert len(coord)==len(mass)
                            mass_sum=sum(mass)
                            COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]),coord)],[count]]
                            #print (COM)
                            COM1+=COM
                        count = count + 1
                #print (np.array(COM1))

                COM1=np.array(COM1).reshape(len(COM1)/2,2)
                #coordinates,resnum
            return COM1


        def write_hydrophobic_contacts(self, pairfile, pdbfile):
            print ('\n-------------------\n>>in write_hydrophobic contacts',pairfile,pdbfile)
            print ('hp residues are A,V,L,I,M,F,W,Y')
            # Adding attractive interactions among hydrophobic residues to the SBM.INP file.
            # Three constraints are applied:1) hydrohpbic only, 2) not defined in native contacts.3) Distance > cut-off
            # Pick out positions of hydrophobic residues from pdbfile.
            # Todo: Add options for potential type(gaussian, db,10-5,6-12, etc.)
            hydrophobic = self.get_hydrophobic(pdbfile)

            # traj = md.load_pdb(pdbfile)
            nnc = []
            for a in itertools.combinations(hydrophobic, 2):
                # make list of non-native contacts i<j-3.
                if (a[1] - a[0] > 3):
                    nnc += a
            nnc = np.array(nnc)
            count_nnc=len(nnc)
            #print (count_nnc,int(len(nnc)/2))
            assert count_nnc==int(len(nnc))
            nnc = np.reshape(nnc, ((int(len(nnc)/2)), 2)).astype(int)
            # calculate distance between eacdist = md.compute_distances(traj, nnc)h contact.
            # dist = md.compute_distances(traj, nnc)
            #get pairs of native-contacts from external file 'contacts.txt'
            nc = self.get_pairs_ext(pairfile)+1
            # remove already existing native contact pairs.
            identical = np.empty((0, 2), int)
            #non-native contacts.
            nnc1 = np.array([x + 1 for x in nnc])
            count = 0
            count1 = 0
            indice = np.empty((0), int)
            for i in nnc1:
                count1 = count1 + 1
                for j in nc:
                    if np.array_equal(i, j):
                        # what is the index of nnc1
                        count = count + 1
                        identical = np.append(identical, np.array([i]), axis=0)
                        indice = np.append(indice, np.array([count1 - 1]), axis=0)
            print ('Found', len(nnc1), 'i<j-3 hydrophobic non-native contacts')
            print  ('Found', count, 'overlapping hydrophobic native and non-native contacts in file', pairfile,'from pdbfile:', pdbfile)
            nnc1 = np.delete(nnc1, indice, axis=0)
            # nnc1[nnc1[:,1],argsort()]
            #        np.sort(nnc1,order=['f1'],axis=0)
            print  ('After deleting overlapping contacts,', len(nnc1))

            # print in SBM.INP format to make it OPTIM readable.f
            np.savetxt('HYDROPHOBIC', nnc1, delimiter='\t', fmt='%d     %d 2   5.50000 1.00000e+00 ')
            np.savetxt('IDENTICAL', identical, delimiter='\t', fmt='%d')
            print ('See file: HYDROPHOBIC for hydrophobic contacts added.')
            #print (nnc1)
            return nnc1

        def get_rg(self,traj):
            print (">>get_rg")
            #traj=self.get_trajfromxtc(xtcfile,pdbfile)
            rg=md.compute_rg(traj)
            print (rg)
            np.savetxt('rg.dat',rg)
            return rg

        def get_energy(self,minfile):
            #get energies from min.data file
            print(">>get_energy")
            assert (minfile=='min.data')
            energy=np.loadtxt(minfile)
            energy=np.asarray(energy[:,0])
            # print (energy)
            return energy

        def get_energy_v_rg(self,xyztrajfile,minfile,nativefile):
            #make rog vs E curve
            print (">>get_energy_v_rg")
            traj=self.read_xyz_as_traj(xyztrajfile,nativefile)
            rg=np.asarray(self.get_rg(traj))
            e = np.asarray(self.get_energy(minfile))
            print (np.amax(e),np.amin(e))
            e = e-np.amin(e)
            #e = e-np.max(e)
            #assert (len(rg)==len(e)),'Are frames in trajectory equal to min.data indices?'
            rg.reshape(len(rg),1)
            e.reshape(len(rg), 1)
            print (rg.shape,e.shape)
            rge=np.column_stack((rg,e))
            np.savetxt('rge.dat', rge)
            self.plot_map(rg,e,'RGvE','X','Y')
            return True
        def get_dist_resids(self,xyzfile,nativefile,pair):
            #save distances between residue paurs
            traj = self.read_xyz_as_traj(xyzfile, nativefile)
            d1=md.compute_distances(traj,pair,periodic=False,opt=True)
            np.savetxt('pair_dist.dat',d1)
        def extract_below_energy(self,energyfile,trajfile,ecut):
            ecut=0
            a=np.loadtxt(energyfile,comments='@',skiprows=19)
            b=a[:,1]
        def get_gro_from_xyz(self,nativefile,xyzfile,grofile):
            #convert .xyz to .gro
            print ('>>get_gro_from_xyz')
            traj=self.read_xyz_as_traj(xyzfile,nativefile)
            print (traj.xyz)
        def xyz_to_pdb(self,xyzfile,nativefile):
            print ("Writing pdb for", xyzfile)
            traj=self.read_xyz_as_traj(xyzfile,nativefile)
            traj.save_pdb("xyz_to_pdb.pdb",force_overwrite=True, bfactors=None)
        def pdb_to_xyz(self,pdbfile):
            #convert pdb to xyz
            traj = md.load_pdb(pdbfile)
            traj.save_xyz("pdb_to_xyz.xyz", force_overwrite=True)

        def remove_hydrogens(self,pdbfile):
            import subprocess
            U=Utils()
            assert(U.file_exists('./reduce'));assert(U.file_exists(pdbfile))
            p = subprocess.Popen(["./reduce", pdbfile], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            out, err = p.communicate()
            print (out)
            return

        def cif_to_pdb(self,ciffile):
            p = MMCIFParser()
            structure = p.get_structure('test', ciffile)
            io=PDBIO()
            io.set_structure(structure)
            pdbfile=ciffile.replace('.cif','.pdb')
            io.save(pdbfile)
            print ('See file: ',pdbfile, 'and', ciffile)
            return 1

        # def count_transitions(self,file,framebin, Qthresh):
        #     #Read Qmap file to count transitions.
        #     print ('>>in count_transitions',file,framebin)
        #     d=np.loadtxt(file)
        #     d1=d[0:len(d):framebin]
        #     print("Framebin=", framebin, "Qthresh=", Qthresh, "\n... slicing", file, "with", len(d), "rows to",len(d1),"rows")
        #     diff = [y - x for x, y in zip(d1, d1[1:])];
        #     diff = np.asarray(diff)
        #     x = np.where(diff > Qthresh);
        #     print("FOUND", len(diff[x]),"TRANSITIONS")
        #
        # def split_chains(self,pdbfile):
        #     U=Utils()
        #     U.split_chain(pdbfile)


def main():
       #add argparse to read external i/o.
       import argparse
       parser = argparse.ArgumentParser(description=" ")
       parser.add_argument("--cutoff","-cutoff", help="Cut-off for calculating contact map.")
       parser.add_argument("--scaling","-scaling", help="Scaling for calculating contact map .")
       parser.add_argument('--get_seq',"-get_seq", help='Get sequence from pdb')
       parser.add_argument("--native","-native", help="Give CA-CA pdbfile as native file")
       parser.add_argument("--aa_pdb","-aa_pdb", help="Give all-atom pdbfile")
       parser.add_argument("--traj", "-traj", help="Give trajectory filename. Else traj.xtc assumed")
       parser.add_argument('--traj_list',"-traj_list", help='List of trajectories. Comma separated.No spaces.', type=str)
       parser.add_argument("--frames","-frames", help="Give three comma separated integers start, stop, step",type=str)
       parser.add_argument("--xyz_out","-xyz_out", help="Give XYZ output file name for extracting trajectory.")
       parser.add_argument("--xyz_traj","-xyz_traj", help="Give XYZ trajectory file containing lots of xyz structures.")
       parser.add_argument("--get_polar","-get_polar", help="Get positions of polar residues ")
       parser.add_argument("--get_charged","-get_charged", help="Get positions of charged residues ")
       parser.add_argument("--get_hp","-get_hp", help="Get positions of hydrophobic residues ")
       parser.add_argument("--ext_pair", "-ext_pair", help="Get contact pairs from any other program(say SMOG). Give smog.contacts(reads columns 2 and 4)")
       parser.add_argument("--ext_con_map","-ext_con_map", help="External contact map provided (0 or 1)")
       parser.add_argument("--con_prob","-con_prob", action='store_true', default=False,help="Get contact probabilities from xyz or xtc trajectory. ")
       parser.add_argument("--pl_fe","-pl_fe", help="Plot \" Free energy profile\". Input file containing only Q values. e.g. 1.map ")
       parser.add_argument("--cluster", "-cluster", action='store_true', default=False,help="Hierarchical clustering of structures in xtc file.")
       parser.add_argument('--get_energies',"-get_energies", action='store_true', default=False,help="Get energies from single point OPTIM runs for each frame in xyz trajectory ")
       parser.add_argument('--check_pdb',"-check_pdb",help="Give pdb  ID. Outputs smog readable pdb ")
       parser.add_argument('--get_pdb',"-get_pdb", help="Give pdb ID. Downloads pdb ")
       parser.add_argument('--gconmap', "-gconmap",help="Give pdb file. Generate contact map from all atom pdb ")
       parser.add_argument("--pl_map","-pl_map", action='store_true', default=False,
                           help="Plot contact map from gconmap.")
       parser.add_argument('--renumber_res', "-renumber_res", help="Give pdb file. Renumbers residues starting from 1 ")
       parser.add_argument('--hphobic',"-hphobic", help='Generate hydrophobic contacts. Give PDB file. Also need native_pairs file for eliminating duplicates.')
       parser.add_argument("--Qcut",
                           help="Extract as xyzfile structures from xtc file with Q>=value. Enter Q value. ")
       # parser.add_argument("--gen_table",
       #                     help="Give type of potential. E.g. 'db' for desolvation barrier.")
       # parser.add_argument("--gen_many_table",
       #                     help="Give type of potential. E.g. 'db' for desolvation barrier. And ext_pair option for pairs.")
       parser.add_argument("--rg","-rg", action='store_true', default=False,
                           help="Get radius of gyration for each frame in trajectory.")
       parser.add_argument("--rge","-rge", action='store_true', default=False,
                           help="Get E vs rog plot for entire trajectory.")
       parser.add_argument("--count_trans","-count_trans",
                           help="Count no. of transitions in MD trajectory. Input:filename,framebin,Qthresh")
       parser.add_argument("--mindata",help="Give minfile. Usually min.data.")
       parser.add_argument("--remove_H",
                           help="remove hydrogens from pdbdfile. Give pdbfile.")
       parser.add_argument("--cif2pdb","-cif2pdb", help="ciffile pdbfile")
       parser.add_argument("--pdb2xyz","-pdb2xyz", help="pdbfile xyzfile")
       parser.add_argument("--xyz2pdb","-xyz2pdb", help="xyzfile pdbfile;add --native native_ca.pdb")
       parser.add_argument("--dssp", "-dssp",help="Run dssp on xtc file. Needs dssp executable in folder.")


       #initialise classes .
       args = parser.parse_args()
       X = conmaps()
       U = Utils()
       #X.check_hetatm('1ris.pdb')
#           X.check_glycineh('1qys_H.pdb')
       #traj=X.read_xyz_as_traj('extract1.xyz','native_ca.pdb')
       #Y=X.get_dssp(traj)
       #print (Y)
       #X.split_chains('1.pdb')
       #get contact-map of chain 1
       #X.all_atom_contacts('chain_A.pdb',4.5)
       #get contact-map of chain 2
       #conmap=X.all_atom_contacts('chain_A.pdb', 4.5)
       #print (conmap[:,0]+1+59,conmap[:,0])
       #X.pdb_to_xyz('gb98_opt.pdb')
       #pairs=X.get_pairs_ext('ext_pair')
       #X.get_dist_resids('extract.xyz', 'native_ca.pdb',pairs)
       #X.get_energy('min.data')
       #dist=X.get_distances('/media/Data/sridhar/dswap/H2Z_Quench_f_146_0/OPTIM/extract1.pdb',[22,187])
       #np.savetxt('distances.dat')
       #X.xyz_to_pdb('extract_removed.xyz','native_all.pdb')
       #topology=md.load('1ubq.pdb').topology
       #X.get_SOPC_non_native_angles('1UBQ_PRO_WH.pdb','native_all.pdb',4.0,1)
       #X.get_SOPC_non_native('1UBQ_PRO_WH.pdb','native_all.pdb',1)
       #print (X.get_res_name_fr_idx(topology, 9))
       #exit()
       #X.get_native_contacts('80.xtc','native_ca.pdb','1qys.pdb',4.5,1.2,'Qfile',10000000)
       #X.gen_table_file('db',4.0,0.04,350,0.02,1)
       #traj=X.get_trajfromxtc('101_5.xtc','native_ca.pdb')
       #X.md_rmsd(traj)
       #X.gen_many_table_file('contacts.txt','native_ca.pdb')
       #X.get_bb_contacts_SOPC('1qys_H.pdb','native_all.pdb',8.0,1.0,4)
       #X.get_ss_contacts_SOPC('1qys_H.pdb','native_all.pdb',6.0,1.0,5)
       #X.get_epsij_from_ext('external_pair','native_all.pdb')

       #X.get_side_chain_COM('1qys_H.pdb')
       #X.count_md_transitions('traj.xtc','native_ca.pdb','1qys.pdb',4.5,1.2)
       #X.count_transitions('Qmap',100,0.2)

#for testing only
       #X.write_hydrophobic_contacts('contacts.txt', '1ris.pdb')
       #X.extract_below_energy('energy.xvg','traj.xtc',0)
       #bb=X.get_backbone_contacts('1qys.pdb',4.5)
       #sc=X.get_side_chain_contacts('1qys.pdb', 4.5)
       #print (len(bb)+len(sc),len(bb),len(sc))
       #X.check_pdb('2jvf.pdb')
       #X.check_pdb('2jvf.pdb')
       #X.write_as_xyz('traj.xtc','native_ca.pdb',1,50000,10,'sig5.xyz')
       #X.write_as_xyz(traj='traj.xtc', pdbfile='2jvf.pdb', st  art='1', finish='30000', step='10', xyzfile='M7.xyz')

       #X.all_atom_contacts('1qys.pdb',4.0)
       #X.get_pairs_ext('external_pair')

       # first frame in file is assumed as native state.
       #X.get_contact_probability(traj,'native_ca.pdb','external_pair')
       #x=np.loadtxt('80.map')
       #X.plot_2d_hist('rg_Q.dat','test')
       #X.get_rg('80.xtc','native_all.pdb')
       #X.cif_to_pdb('3spu.cif')

       #exit()
       if args.dssp:
           U.file_exists("native_ca.pdb")
           if not args.traj:
               U.fatal_errors(4)
           else:
               xtcfile=args.traj
               U.check_file_extension(xtcfile,'.xtc')
               traj=X.get_trajfromxtc(xtcfile,'native_ca.pdb')
               X.get_dssp(traj)

       if args.remove_H:
           pdbfile=str(args.remove_H)
           X.remove_hydrogens(pdbfile)
       if args.count_trans:
           arg1=[x for x in args.count_trans.split(',')]
           file=str(arg1[0]);framebin=int(arg1[1]);thresh=float(arg1[2])
           framebin=np.arange(1,10,1)
           for i in framebin:
            X.count_transitions(file,i,thresh)
           exit()
       if args.pdb2xyz:
           file=str(args.pdb2xyz)
           X.pdb_to_xyz(file)
       if args.xyz2pdb and not args.native:
           U.fatal_errors(6)
           print ("\nFatal: Give native pdb file with --native option. eg. --xyz2pdb 1.pdb --native native_ca.pdb\n")
       if args.xyz2pdb and args.native:
           xyzfile=str(args.xyz2pdb)
           nativefile=str(args.native)
           X.xyz_to_pdb(xyzfile,nativefile)
       if args.pl_fe:
           f=str(args.pl_fe)
           f = [item for item in args.pl_fe.split(',')]
           for i in f:
            print (i)
            if U.file_exists(i):
                x=np.loadtxt(i)
                X.plot_histogram(x,str(i))
           exit()

       if (args.cutoff):
          cutoff=args.cutoff
       else:
           cutoff=4.5
       if (args.scaling):
          scaling=args.scaling
       else:
           scaling=1.2

       if args.gconmap:
           pdbfile=args.gconmap
           #print (cutoff)
           X.check_hetatm(pdbfile)
           conmap=X.all_atom_contacts(pdbfile, cutoff,scaling)
           if args.pl_map:
               #plt([conmap[:,0].T.tolist()],[conmap[:,1].T.tolist()])
               x=conmap[:,0]
               y=conmap[:,1]
               assert len(x)==len(y)
               title='contact-map'
               X.plot_map(x,y,title,'Res2','Res1')
       if args.get_polar:
           pdbfile=args.get_polar
           X.get_polar(pdbfile)
       if args.get_charged:
           pdbfile=args.get_charged
           X.get_charged(pdbfile)
       if args.get_hp:
           pdbfile=args.get_hp
           X.get_hydrophobic(pdbfile)
       if args.get_seq:
           pdbfile=args.get_seq
           X.get_sequence(pdbfile)
           exit()
       if args.ext_pair:
           contactfile=args.ext_pair
           U.file_exists(contactfile)
           X.get_pairs_ext(contactfile)
       if args.get_pdb:
           pdbid=args.get_pdb
           pdb1=PDBList()
           pdb1.retrieve_pdb_file(pdbid, pdir='./', obsolete=False, overwrite=True)
           X.cif_to_pdb(pdbid+'.cif')
           exit()
       if args.aa_pdb:
        pdbfile = args.aa_pdb
       if args.cif2pdb:
           ciffile=str(args.cif2pdb)
           X.cif_to_pdb(ciffile)
       if args.frames:
           #frames read as string converted to list.
           list = [int(item) for item in args.frames.split(',')]
           start=list[0]
           finish= list[1]
           step= list[2]
           print ('start,finish,step',start,finish,step)
       else:
           #print ('Using default values for frames.')
           finish = 500000
           start = 1
           step=1
           #print('\n\n\n\n\n\nDefault: start,finish,step', start, finish, step)
       if args.rg:
           print ('Will calculate radius of gyration for trajectory')
           if args.traj:
            xtcfile=str(args.traj)
            pdbfile=args.native
            print ('Readin xtc file...',xtcfile)
            traj=X.get_trajfromxtc(xtcfile,pdbfile)
            X.get_rg(traj)
           elif args.xyz_traj:
            xyzfile=str(args.xyz_traj)
            nativefile=str(args.native)
            print('Reading xyz trajectory',xyzfile)
            traj=X.read_xyz_as_traj(xyzfile,nativefile)
            X.get_rg(traj)
           exit()
       if args.rge and args.mindata and args.native and args.xyz_traj:
           print ('Will calculate ROG vs Energy curves for trajectory')
           print('!!!beta: traj file should be in xyz format!!!')
           minfile=str(args.mindata)
           xyzfile=str(args.xyz_traj)
           nativefile=args.native
           U.file_exists(minfile);U.file_exists(xyzfile);U.file_exists(nativefile)
           X.get_energy_v_rg(xyzfile,minfile,nativefile)


       if args.get_energies and args.xyz_traj:
            xyzfile=args.xyz_traj
            X.get_dihedral_energies(xyzfile)
            exit()

       #write xyz from one trajectory
       if args.native and args.traj and args.xyz_out:
           pdbfile = args.native
           traj = args.traj
           U.file_exists(pdbfile)
           U.file_exists(traj)
           if args.xyz_out:
            xyzfile=args.xyz_out
            X.write_as_xyz(traj, pdbfile, start, finish, step,xyzfile)
            exit()
           else:
            print ("Not enough info. Give xyz_out")
            exit()

      #write many  xyzs from list of gromacs trajectories.
       if args.native and args.traj_list:
             if args.xyz_out:
                pdbfile=args.native
                traj = [item for item in args.traj_list.split(',')]
                U.file_exists(pdbfile)
                xyzfile = args.xyz_out
                count=0
                for i in traj:
                    count=count+1
                    count1=str(count)
                    xyzfile=str(count1+'.xyz')
                    U.file_exists(i)
                    X.write_as_xyz(i, pdbfile, start, finish, step, xyzfile)
                exit()
             elif args.aa_pdb:
                 traj = [item for item in args.traj_list.split(',')]
                 if (args.cutoff and args.scaling):
                     cutoff = args.cutoff
                     scaling = args.scaling
                 else:
                     #default values
                     cutoff=4.5
                     scaling=1.2
                     traj = [item for item in args.traj_list.split(',')]
                     print("Mapping contacts from xtc file", args.traj_list, "onto", args.aa_pdb)
                 count=0
                 for i in traj:
                     count=count+1
                     count1=str(count)
                     xtcfile = i
                     native = args.native
                     pdbfile = args.aa_pdb
                     U.file_exists(pdbfile)
                     U.file_exists(native)
                     U.file_exists(xtcfile)
                     Qfile=str(xtcfile +'.map')
                     if args.Qcut:
                         Qcut=args.Qcut
                     else:
                         Qcut=100000000
                         print ('Setting Qcut to', Qcut)
                     pairfile='external_pair'
                     traj=X.get_trajfromxtc(xtcfile,native)
                     pairs = X.get_pairs_ext(pairfile)
                     X.get_native_contacts(traj, native, pdbfile, cutoff, scaling,Qfile,Qcut,pairs)
                     xi=np.loadtxt(Qfile)
                     yi=np.arange(1,len(xi)+1)
                     title=xtcfile
                     X.plot_map(yi,xi,title,'X','Y')

                 exit()
             else:
                print("Not enough info. Give --xyz_out or --aa_pdb option.")
                exit()
       #trajectory=X.get_trajfromxtc(traj,pdbfile)
#      get Q for entire xyz trajectory
#        if args.xyz_traj and args.native and args.aa_pdb:
#            if (args.cutoff and args.scaling):
#              cutoff=args.cutoff
#              scaling=args.scaling
#            print ("Mapping contacts from xyz file", args.xyz_traj,"onto", args.aa_pdb)
#            xyzfile=args.xyz_traj
#            native=args.native
#            pdbfile=args.aa_pdb
#            contactfile = args.ext_pair
#            X.xyz_native_contacts(xyzfile,native,pdbfile,cutoff,scaling)

#      get Q from traj.xtc or extract.xyz
       if args.traj and args.native and args.aa_pdb:
           if (args.cutoff and args.scaling):
             cutoff=args.cutoff
             scaling=args.scaling
             print("Mapping contacts from xyz file", args.xyz_traj, "onto", args.aa_pdb)
             xtcfile = args.traj
             native = args.native
             pdbfile = args.aa_pdb
             pairs=X.get_pairs_ext(pairfile)
             X.get_native_contacts(traj,native,pdbfile,cutoff,scaling,'Qmap',100000,pairs)
       #pdbfile1='native.pdb'
       #xyzfile='native_ca.pdb.xyz'
       #X.get_native_contacts(pdbfile1,xyzfile)

       if args.check_pdb:
            pdbfile=args.check_pdb
            f=X.check_pdb(pdbfile)
            print ('See-file', f)

       if args.hphobic:
           pairfile = 'contacts.txt'
           pdbfile = args.hphobic
           if U.file_exists('contacts.txt'):
               X.write_hydrophobic_contacts(pairfile, pdbfile)
               exit()
           else:
               print ('Give blank file called native_pairs. See bitbucket examples.')
               exit()


#      Will take abot 2-3 minutes to cluster.

       if args.cluster and args.native:
          pdbfile=args.native
          if args.traj:
              xtcfile=args.traj
              traj=X.get_trajfromxtc(xtcfile,pdbfile)
          elif args.xyz_traj:
              xyzfile=args.xyz_traj
              traj=X.read_xyz_as_traj(xyzfile,pdbfile)

          else:
              print('Fatal:Need --xyz_traj or --traj argument')
              exit()
          X.hierachical_cluster_rmsd(traj,pdbfile,start,finish,step)

       if args.con_prob and args.native and args.ext_pair and args.traj:
           external_pair = str(args.ext_pair)
           xtcfile=str(args.traj)
           pdbfile = str(args.native)
           yesxtc = True if U.check_file_extension(xtcfile, '.xtc') else False
           yesxyz = True if U.check_file_extension(xtcfile, '.xyz') else False
           if yesxtc:
               traj=X.get_trajfromxtc(xtcfile,pdbfile)
           elif yesxyz:
               traj=X.read_xyz_as_traj(xtcfile,pdbfile)
           else:
               U.fatal_errors(1)
           a=X.get_contact_probability(traj,external_pair,start,finish,step)

           if args.pl_map:
               x = a[:, 0]
               x1=([s[0] for s in x])
               y1=([s[1] for s in x])
               z1= a[:,1].tolist()
               X.plot_heat_map(x1,y1,z1)

       # if args.gen_table and args.ext_pair and args.native:
       #     pottype=str(args.gen_table).trim()
       #     nativefile=args.native
       #     contactfile=args.ext_pair
       #     assert (pottype=='db')
       #     start=0.02
       #     stop=350
       #     step=0.02
       #     print ('start,stop,step=',start,stop,step)
       #     X.gen_many_table_file(contactfile,nativefile)


if __name__ == '__main__':
    main()

