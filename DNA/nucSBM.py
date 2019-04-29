#!/usr/bin/env python

"""
PyeSBM:
usage: python nucSBM.py --options

"""
from __future__ import print_function
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO,Select
import collections
import numpy as np
import math as m
import mdtraj as md
from eSBM import esbm
from gr import conmaps
from util import Utils

class PrePDB(Select):
	def __init__(self):
		return
	def accept_chain(self, chain):
		if chain in chain_list:
			return 1
		else:
			return 0
	def readPDB(self,file):
		fin = open(file)
		all = fin.read()
		pdblines = all.split("\n")
		pdblines = [x.rstrip().lstrip() for x in pdblines if x.startswith("ATOM")]
		return pdblines	
	def homoNmer2Monomer(self,file):
		print (">>converting homo-N-mer(if the pdb is) to a monomer.")
		P = PDBParser(PERMISSIVE=0)
		structure = P.get_structure('test', file)
		model = structure.get_chains()
		in_chain = dict()
		residues_in_chain = list()
		global chain_list; chain_list = list()		
		for chain in model:
			l = list()
			c = chain.get_id()
			for residue in chain:
				if [residue.get_resname(),residue.get_id()[1]] not in l:
					l.append([residue.get_resname(),residue.get_id()[1]])
			if l not in residues_in_chain:
				residues_in_chain.append(l)
				in_chain[c] = l		
				chain_list.append(chain)
		io = PDBIO()
		io.set_structure(structure)
		io.save("mono_"+file,PrePDB())
		return "mono_"+file
	def appendChain(self,file):
		print (">>>Appending the residues of the chains.\n>>>The chain_id remains same,", file)
		chain_terminal = list()
		fin = open(file)
		lines = [l for l in fin.readlines() if not l[12:16].strip().startswith('H')]
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		nuclines = list(); aalines = list()
		for l in lines:
			if l.startswith(("ATOM","TER")):
				if l.startswith("ATOM"):
					residue = l[17:21].strip()
				if residue in nuc_res:
					nuclines.append(l)
				else:
					aalines.append(l)
		if len(nuclines) == 0:
			print ("The PDB file lacks DNA/RNA. Prefer using the protein obly model......The code will continue")
		lines = nuclines + aalines
		rescount= 0; atcount = 0
		prev_res = 0
		outfile = "concat_"+file
		fout = open(outfile,"w+")
		for l in lines:
			if l.startswith("ATOM"):
				if l[22:26]!=prev_res:
					rescount = rescount + 1
				atcount = atcount + 1
				fout.write(l[0:6]+str(atcount).rjust(5)+l[11:22]+str(rescount).rjust(4)+l[26:])
				prev_res = l[22:26]
			elif l.startswith("TER"):
				fout.write(l[0:6]+str(atcount).rjust(5)+l[11:22]+str(rescount).rjust(4)+"\n")
				chain_terminal.append(rescount)
			else:
				fout.write(l)    
		return [outfile,chain_terminal]
	def sepNucPro(self,file):
		P = PDBParser(PERMISSIVE=0)
		structure = P.get_structure('test', file)
		print (">>writing nucleotide only file: ","nuc_"+file)
		io = PDBIO()
		io.set_structure(structure)
		io.save("nuc_"+file,nucsbm())
		print (">>writing protein only file: ","aa_"+file)
		io.set_structure(structure)
		io.save("aa_"+file,esbm())
		return ("nuc_"+file,"aa_"+file)
	def atomIn(self):
		atoms_in = dict()
		atoms_in["CA"] = ("N","C","CA","O")
		#atoms_in["CB"]  = ['CB', 'CD', 'CD1', 'CD2', 'CD3', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CG3', 'CH', 'CH1','CH2', 'CH', 'CZ', 'CZ1OG', 'CZ2', 'CZ3', 'ND1', 'ND2', 'NE', 'NE1NE2', 'NH1', 'NH2', 'NZ', 'OD1', 'OD2', 'OE1', 'OE2', 'OG1', 'OH', 'SD', 'SG']
		atoms_in["B"] = ("N1","C2","N2","N3","C4","C5","C6","O6","N6","N7","C8","N9","O2","O4","N4","C7")
		atoms_in["S"] = ("C1'","C2'","C3'","C4'","C5'","O2'","O3'","O4'","O5'")
		atoms_in["XP"] = ("P","OP1","OP2","O1P","O2P")
		return atoms_in
	def stackEps(self):
		#fin = open("stackEpsilon.dat")
		#file = [x.split() for x in fin.readlines() if not x.startswith("#") and len(x.rstrip().lstrip())!=0]
		#d=dict()
		#nuc = ["A","G","T","C","U"]
		#for l in file:
		#	if l[0][0] in nuc:
		#		a1 = "B"+l[0][0]
		#	else:
		#		a1 = l[0][0]
		#	if l[0][1] in nuc:
		#		a2 = "B"+l[0][1]
		#	else:
		#		a2 = l[0][1]
		#	if a1 not in d:
		#		d[a1]=dict()
		#	d[a1][a2] = l[1]
		#print("Stacking Epsilon")
		#for i in file:
		#	print(i)
		#for i in range(1,len(file)):
		#	d[file[i][0]] = dict()
		#	for j in range(1,len(file[i])):
		#		d[file[i][0]][file[0][j]] = file[i][j]
		d = {'BT': {'BG': '2', 'BA': '2', 'F': '2', 'H': '2', 'BC': '1', 'BT': '1', 'BU': '1', 'W': '3', 'Y': '2'}, 'BU': {'BG': '2', 'BA': '2', 'F': '2', 'H': '2', 'BC': '1', 'BT': '1', 'BU': '1', 'W': '3', 'Y': '2'}, 'BG': {'BG': '4', 'BA': '3', 'F': '2', 'H': '2', 'BC': '3', 'BT': '2', 'BU': '2', 'W': '3', 'Y': '2'}, 'BA': {'BG': '3', 'BA': '3', 'F': '2', 'H': '2', 'BC': '2', 'BT': '2', 'BU': '2', 'W': '3', 'Y': '2'}, 'BC': {'BG': '3', 'BA': '2', 'F': '2', 'H': '2', 'BC': '3', 'BT': '1', 'BU': '1', 'W': '3', 'Y': '2'}}
		return d

class nucsbm(esbm,Select):
	def __init__(self):
		return
	def accept_residue(self,residue):
		global nuc_res 
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		if residue.get_resname().lstrip().rstrip() in nuc_res:
			return 1
		else:
			return 0
	def getIndices(self,pdbfile):
		structure = PrePDB().readPDB(pdbfile)
		CG_indices = dict()
		for l in structure:
			atom_num = int(l[6:11].rstrip().lstrip())
			atom_name = l[12:16].rstrip().lstrip()
			res_id = int(l[22:26].rstrip().lstrip())
			X = float(l[30:38].rstrip().lstrip());Y = float(l[38:46].rstrip().lstrip());Z =	float(l[46:54].rstrip().lstrip())	
			xyz = np.array([X,Y,Z])
			CG_indices[(atom_num)] =[res_id,atom_name,xyz]
		return CG_indices
	def getCoordinates(self,pdbfile,CG_pos):
		print ("Getting Coarse-Grain Co-ordinates.")
		P = PDBParser(PERMISSIVE=0)
        
		structure = P.get_structure('test', pdbfile)
		atomname=dict()
		atomname["Bpu"] = ("N1","C2","H2-N2","N3","C4","C5","C6","O6-N6","N7","C8","N9")
		atomname["Bpy"] = ("N1","C2","O2","N3","C4","O4-N4","C5","C6")
		atomname["S"] = ("C1'","C2'","C3'","C4'","C5'","H2'-O2'","O3'","O4'","O5'")
		atomname["P"] = ("P","OP1","OP2","O5'")
		#Nucleotide residues in the PDB
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		bases = dict()
		bases["Bpy"] = ("T","C","U","DT","DC")
		bases["Bpu"] = ("A","G","DA","DG")
		
		count = 0
		CG_coord = dict()
		res_in_chain = list()
		
		for group,pos in CG_pos.items():
			CG_coord[group] = []
			count=0
			for chain in structure.get_chains():
				res_count = 0
				c=chain.get_id()
				#print(count)
				for residue in chain:
					r=residue.get_resname().rstrip().lstrip()
					coord = []
					mass = []
					#print (str(residue).split()[1].lstrip().rstrip())
					if r in nuc_res:
						if group[0]=="B" and r not in bases[group]:
							res_count = res_count + 1
							continue	#skipping this loop as group not in residue
						#print (r)
						if pos=="COM":
							#print(pos)
							for atom in residue:
								a=atom.get_name().rstrip().lstrip()
								if a in atomname[group]:
									#getting co-ordinates and mass
									coord.append(atom.get_coord())
									mass.append(atom._assign_atom_mass())
									assert len(mass)!=0, "Zero length mass array!"
							coord = np.array(coord)
							mass = np.array(mass)
							assert len(coord)==len(mass)
							mass_sum=sum(mass)
							COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]),coord)],[count+residue.id[1]],[(c,residue.id[1],r)]]
							#print(COM)
							#COM is a array of 2 elements in each row: 1)list of co-ordinates. 2) residue number
							CG_coord[group]+=COM
							#print ("COM",CG_coord[group])
						else:
							try:
								atom = residue[pos]
								#print(pos); print (atom)						
								coord.append(atom.get_coord())
								mass.append(atom._assign_atom_mass())
								coord = np.array(coord)
								mass = np.array(mass)
								mass_sum=sum(mass)
								COM = [[np.matmul(([float(atom_mass / mass_sum) for atom_mass in mass]),coord)],[count+residue.id[1]],[(c,residue.id[1],r)]]
								#in the below condition COM is not center of mass
								#COM = [[np.float_(coord)],[count+residue.id[1]],[(c,residue.id[1],r)]]
								#print (COM)
								CG_coord[group]+=COM
								#print ("ATOM",CG_coord[group])
							except:
								print ("Note: Atom "+pos+" missing in residue "+str(res_count+1)+" of chain "+c+". Skipping....")
						res_count = res_count + 1
				#print (res_count)
				count = count + res_count
				#print (count)
				#print(chain.get_id())
				if res_count > 0:
					cid = chain.get_id()
					if (cid,res_count) not in res_in_chain:
						res_in_chain.append((cid,res_count))
			CG_coord[group]=np.array(CG_coord[group]).reshape(len(CG_coord[group])/3,3)
		return CG_coord
	def writeNative(self,pdbfile,CG_pos,grofile):	
		if len(grofile)==0:
			grofile = pdbfile+".gro"
		print (">>Writing native and .gro files.", pdbfile, grofile)
		checksize = open(pdbfile)
		lines = checksize.readlines()
		if len(lines) <= 1:
			print ("As mentioned earlier, the ibpu file is PROTEIN ONLY.....Exiting now")
			return 0
		else:
			#if the file is not emplty
			terminal_residues = list()
			for i in lines:
				if i.startswith("ATOM"):
					resnum = int(i[22:26])
				elif i.startswith("TER"):
					terminal_residues.append(resnum)
			nchains = len(terminal_residues)
		#Defining Nitrogen Base types
		bases = dict()
		bases["Bpy"] = ("T","C","U","DT","DC")
		bases["Bpu"] = ("A","G","DA","DG")
		#Getting co-ordinates for beads
		CG_coord =dict()
		#function returns co-ordintates, residue, chain_id
		CG_coord = self.getCoordinates(pdbfile,CG_pos)
		CG_others =dict()
		natoms_CG = 0
		res_num  = list()
		for group,xyz in CG_coord.items():
			coord = dict(); others = dict()			
			for (x,y,z) in xyz:
				#separating co-rdinates from other parameters and storing at same index
				coord[int(y)] = x.tolist()     #[X,Y,Z] co-ordinates
				others[int(y)] = z				#[chain_id,res_number,res_name]			
				if int(y) not in res_num:
					res_num.append(int(y))
				#print (y); #print (coord[y])
				natoms_CG = natoms_CG + 1
			#segrigating based on group
			CG_coord[group] = coord
			CG_others[group] = others
		res_num.sort()
		#Writing PDB
		l = {}
		k = {}
		iter = 0
		#counting total atoms
		fgro = open(grofile,"w+")
		fgro.write('%s\n'%('Grofile generated from Lab16 repo. See https://bitbucket.org/nsridhar/lab16'))#.write('%s %s%s%s\n' % ("CG-MD for",pdbfile.split(".")[0],"\n",str(natoms_CG).rjust(5)))
		fgro.write('%s\n' % str(natoms_CG).rjust(5))
		fout = open('nuc_native_P-S-B.pdb', "w+")
		#for i in range(1,len(CG_coord["S"])+1):
		
		#getting terminal residues
		
		
		for i in res_num:
			print ("Writing residue ",100*i/len(CG_coord["S"]),"%")
			for group,coord in CG_coord.items():
				#print (group)
				try:
					#print(CG_others[group][i][2])
					r = CG_others[group][i][2]
					#print(r)
					if group[0]=="B" and r not in bases[group]:
						continue	#skipping this loop as group not in residue
					#print(CG_others[group][i])
					#print (group)
					#print(CG_coord[group][i])
					iter = iter + 1
					l[0] = "ATOM".ljust(6)						#ATOM 1-6,	6char
					l[1] = str(iter).rjust(5)					#atom number 7-11, 5char
					l[2] = "".ljust(1)							#blank at 12, 1char		
					if group[0] == "S":
						if r[0] == "D":
							l[3] = str("D"+group[0]).center(4)	#atom name (P,S or B) 13-16, 4char
						else:
							l[3] = str("R"+group[0]).center(4)
					elif group[0] == "P":
						l[3] = str("X"+group[0]).center(4)
					elif group[0] == "B":
						if r[0] == "D":
							l[3] = str(group[0]+r[1]).center(4)
						else:
							l[3] = str(group[0]+r[0]).center(4)
					l[4] = "".ljust(1)							#atom indicator/blank	17, 1char
					l[5] = r.ljust(3)							#residue name 18-20, 3char
					l[6] = "".ljust(1)							#blank 21, 1char
					l[7] = CG_others[group][i][0].ljust(1)		#chain 22, 1char
					l[8] = str(CG_others[group][i][1]).rjust(4)	#residue number 23-26, 4char
					l[9] = "".ljust(4)							# insertinal code 27, blank 28-30, 4char
					#print (CG_coord[group][i][0][0])
					l[10]= str('%8.3f' % float(CG_coord[group][i][0])).rjust(8)	#X 31-38, 8char 
					l[11]= str('%8.3f' % float(CG_coord[group][i][1])).rjust(8)	#X 39-46, 8char
					l[12]= str('%8.3f' % float(CG_coord[group][i][2])).rjust(8)	#X 47-54, 8char
					l[13]= "1.00".rjust(6)							#occuplancy factor 55-61, 6char
					l[14]= "00.00".rjust(6)							#temperature factor 62-68, 6char
					l[15]= "".rjust(12)							# element/labels/foontnotes name 69-80
					#print (CG_coord[group][i][0][0])
					#print("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" % (l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15]))
					fout.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" % (l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11],l[12],l[13],l[14],l[15]))
					k[0] = l[8].rjust(5)						#residue number grofile
					k[1] = l[5].ljust(5)						#residue name
					k[2] = l[3].rstrip().lstrip().rjust(5)		#atom name
					k[3] = l[1]									#atom number
					k[4]= str('%8.3f' % (float(CG_coord[group][i][0])/10)).rjust(8)	#X 31-38, 8char 
					k[5]= str('%8.3f' % (float(CG_coord[group][i][1])/10)).rjust(8)	#X 39-46, 8char
					k[6]= str('%8.3f' % (float(CG_coord[group][i][2])/10)).rjust(8)	#X 47-54, 8char
					fgro.write('%s%s%s%s%s%s%s\n' % (k[0],k[1],k[2],k[3],k[4],k[5],k[6]))
				except:
					buffer = 0
					#print ("\t"+group+" missing. Skipping")
			if int(l[8]) in terminal_residues:
				#writing terminal line
				fout.write("TER".ljust(6)+l[1]+l[2]+" "*(4)+l[4]+l[5]+l[6]+l[7]+l[8]+"\n")
		
		print ("Default box size is 10nm X 10nm X 10nm. Changes can be made directly in ", grofile," or using editconf.")
		fgro.write('%10.5f%10.5f%10.5f' % (10,10,10))
		print ("100% Done!!!!")
		fout.close(); fgro.close()
		return 1
	def nucAtomsType(self,pdbfile):
		#lading structure
		P = PDBParser(PERMISSIVE=0)
		structure = P.get_structure('test', pdbfile)
		#obtaining residues
		residues = structure.get_residues()
		r=[x.get_resname().lstrip().rstrip() for x in residues]
		#identifying unique nuclrotide residues
		temp = list()
		nuc_res = ["DA","DG","DC","DT","A","G","T","C","U"]
		for x in r:	
			if x in nuc_res and x not in temp:
				temp.append(x)
		r=temp
		#storing and counting deoxy nucleotides
		deoxy_count = [x for x in r if x.startswith("D")]
		#storing and counting oxy nucleotides
		oxy_count = [x for x in r if not x.startswith("D")]
		#identifying N-bases (x[1]) that are unique to oxy_count
		not_in_oxy = [x for x in deoxy_count if x[1] not in oxy_count]
		assert(len(r)-len(deoxy_count)==len(oxy_count))
		#types of coarsed grain atoms
		na_type=len(deoxy_count)+len(deoxy_count)-len(not_in_oxy) + 1 #1 for phosphate
		#(A union B)=A+B-(A intercept B)
		a_type = ["B"+x[1] for x in deoxy_count]		#B as base represntator
		a_type = a_type + ["B"+x for x in oxy_count if "B"+x not in a_type] + ["XP"]
		assert(len(a_type)==na_type)
		#Adding sugar bead
		if len(deoxy_count)!=0:
			na_type = na_type + 1 #deoxy sugar  
			a_type.append("DS")
		if len(oxy_count)!=0:
			na_type = na_type + 1 #deoxy sugar 
			a_type.append("RS")
		return (na_type,a_type)
	def getContacts(self,pdbfile,nativefile,cutoff,atomname,strength):
		print ('>>writing contacts section,atomtype')
		cutoff=float(cutoff)
		print (">>cutoff= ",cutoff," A")
		#loading all atom pdb file
		load_atom = PDBParser().get_structure("test",pdbfile).get_atoms()
		#storing atom objects in a list
		atom = [a for a in load_atom]
		#identifying all atom contacts
		#Note: unlike protein CG models, nucleotide models should incorporate interactions between adj. residues (resi2-resi1>=1)
		nuc_res = ("A","G","T","U","C","DA","DG","DT","DC")
		all_contacts = dict()
		nuc_contacts = dict()
		nucpro_contacts = dict()
		natoms = len(atom)
		print ("Determining contacts")
		for i in atom:
			for j in atom:
				d = float(j-i)#distance
				r1 = i.get_parent().id[1]
				r2 = j.get_parent().id[1]
				rname1 = i.get_parent().get_resname().lstrip().rstrip()
				rname2 = j.get_parent().get_resname().lstrip().rstrip()
				if d <= cutoff:
					if r2!=r1 and rname1 in nuc_res and rname2 in nuc_res:#
						if ((r2,r1),(j.get_id(),i.get_id())) not in nuc_contacts:
							nuc_contacts[(r1,r2),(i.get_id(),j.get_id())] = [d,(rname1,rname2)]
					elif (rname1 in nuc_res and rname2 not in nuc_res) or (rname1 not in nuc_res and rname2 in nuc_res):
						if ((r2,r1),(j.get_id(),i.get_id())) not in nucpro_contacts:
							nucpro_contacts[(r1,r2),(i.get_id(),j.get_id())] = [d,(rname1,rname2)]
					#for all in range < 4.5A
					if ((r2,r1),(j.get_id(),i.get_id())) not in all_contacts:
						all_contacts[(r1,r2),(i.get_id(),j.get_id())] = [d,(rname1,rname2)]
		print ("Total contacts in all atom pdb file = ",len(all_contacts))
		print ("Total contacts within nucleiotides = ",len(nuc_contacts))
		print ("Total nucleotide-protein interface contacs = ",len(nucpro_contacts))
		#loading nativefile (CG model file)
		CG_file = PrePDB().readPDB(nativefile)
		CG_coord = dict()
		for l in CG_file:
			key = (int(l[22:26].rstrip().lstrip()),l[12:16].rstrip().lstrip())
			X = float(l[30:38].rstrip().lstrip());Y = float(l[38:46].rstrip().lstrip());Z =	float(l[46:54].rstrip().lstrip())	
			A = int(l[6:11].lstrip().rstrip())
			value = np.array([X,Y,Z])
			CG_coord[key] = [value,A]
			#print (key,value)
		#defining atoms in CG bead
		atoms_in = dict()
		atoms_in = PrePDB().atomIn()
		stack_check = dict()	#bases are stacking if stleast 6*(6-1) = 30 in the interacting ditances
		nuc_cg = dict()
		for pair,distance in nuc_contacts.items():
			pair_id  = pair[0]
			pair_atom = pair[1]
			pair_resname = distance[1]
			pair_distance = distance[0]
			#identifyinf CG atoms for the representative all atoms
			for group,atoms in atoms_in.items():
				if pair_atom[0] in atoms:
					a1 = group
				if pair_atom[1] in atoms:
					a2 = group
			if a1 == "S":
				if pair_resname[0][0] == "D":
					a1 = "D" + a1	#adding D for deoxyribose
				else:
					a1 = "R" + a1	#adding R for ribose
			elif a1 == "B":
				if pair_resname[0][0] == "D":
					a1 = a1 + pair_resname[0][1]  #Adding A/T/G/C/
				else:
					a1 = a1 + pair_resname[0][0]  #Adding A/T/G/c/u
			if a2 == "S":
				if pair_resname[1][0] == "D":
					a2 = "D" + a2
				else:
					a2 = "R" + a2
			elif a2 == "B":
				if pair_resname[1][0] == "D":
					a2 = a2 + pair_resname[1][1] 
				else:
					a2 = a2 + pair_resname[1][0]
			pair_atom = (a1,a2)
			if a1[0]=="B" and a2[0]=="B":
				if (pair_id,pair_atom) not in stack_check:
					stack_check[(pair_id,pair_atom)] = 0
				stack_check[(pair_id,pair_atom)] = stack_check[(pair_id,pair_atom)] + 1
			#invoking CG coordinates from nativefile
			if (a1[1]=="P" and a2[1]=="S") or (a1[1]=="S" and a2[1]=="P"):
				if (int(pair_id[0])-int(pair_id[1]))**2 ==1:
					continue
					#Skipping the loop if P and S are adj,
			#storing CG coordinates from nativefile and calculating distances
			coord1 = CG_coord[pair_id[0],pair_atom[0]][0]
			coord2 = CG_coord[pair_id[1],pair_atom[1]][0]
			CG_dist = (np.sum((coord2-coord1)**2))**0.5 
			#storing atom id grom the nativefile
			atom1_id = CG_coord[pair_id[0],pair_atom[0]][1]
			atom2_id = CG_coord[pair_id[1],pair_atom[1]][1]
			#print (pair_id,pair_atom,CG_dist)	
			nuc_cg[(pair_id,pair_resname,pair_atom)] = (CG_dist,atom1_id,atom2_id)

		contacttype = 2
		#Writing aromatic contacts
		print (">>Aromatic contacts<<")
		for k,v in stack_check.items():
			print (k,v)
		#Storing a list of contacts for pair section
 		CG_contacts = list()
		for key,value in nuc_cg.items():
			distance = value[0]; atom1_id = value[1]; atom2_id = value[2]
			a1 = key[2][0]; a2 = key[2][1]
			#atom1_id,atom2_id,contacttype,distance,epsilon, residue 1, atomname1, residue 2, atomname2
			if a1[0]=="B" and a2[0]=="B" and stack_check[(key[0],key[2])] >= 30:
				CG_contacts.append((atom1_id,atom2_id,contacttype, distance, strength[key[2][0]][key[2][1]], key[0][0], a1, key[0][1], a2))
			CG_contacts.append((atom1_id,atom2_id,contacttype, distance, 1, key[0][0], a1, key[0][1], a2))			
		CG_contacts.sort()
		#writing pair
		
		#if len(nucpro_contacts)>0:
		#	fout = open("nucpro_inter_contacts.t","w+")
		#	for pair,distance in nucpro_contacts.items():
		#		pair_id  = pair[0]
		#		pair_atom = pair[1]
		#		pair_resname = distance[1]
		#		pair_distance = distance[0]
		#		fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(pair_atom[0],pair_atom[1],str(pair_id[0]),str(pair_id[1]),pair_resname[0],pair_resname[1],str(pair_distance)))
		#	fout.close()
		return CG_contacts
	def chargeContacts(self,pdbfile,nativefile,cutoff):
		print (">>Writig electrostatic pairs<<")
		native_charge_contacts = False	#will be user defined in future updates
		if native_charge_contacts:
			print ("Cutoff = ",cutoff)
			#loading all atom pdb file
			load_residues = PDBParser().get_structure("test",pdbfile).get_residues()
			#storing atom objects in a list
			residue = [a for a in load_residues]
			Q_pairs = list()
			count = 0
			for i in range(0,len(residue)-1):
				for atom1 in residue[i]:
					for j in range(i+1,len(residue)):
						for atom2 in residue[j]:
							if ((atom1 - atom2)**2)**0.5 <= 2*cutoff:
								r1 = atom1.get_parent().id[1]
								r2 = atom2.get_parent().id[1]
								if [r1,r2] not in Q_pairs: 
									Q_pairs.append([r1,r2])

			fnat = open(nativefile)
			#native atom number fo residues in pair
			native_XP = dict()
			for i in fnat:
				if i[12:16].strip() == "XP":
					native_XP[int(i[22:26])] = int(i[6:11])
			fnat.close()

			temp = list()
			for i in Q_pairs:
				try:
					temp.append([native_XP[i[0]],native_XP[i[1]]])
				except:
					print ("Seems like Residue",i[0],"or",i[1],"lack XP")
			Q_pairs = temp
		
			return Q_pairs
		else:
			all = PrePDB().readPDB(nativefile)
			q_list = list()
			q_coord = dict()
			q_contacts = list()
			for i in all:
				if i[12:16].strip() == "XP":
					q_list.append(int(i[6:11]))
					X = float(i[30:38].rstrip().lstrip());Y = float(i[38:46].rstrip().lstrip());Z =	float(i[46:54].rstrip().lstrip())	
					value = np.array([X,Y,Z])
					q_coord[int(i[6:11])] = value
			for i in range(0,len(q_list)-1):
				for j in range(i+1,len(q_list)):
					dist  = sum((q_coord[q_list[i]]-q_coord[q_list[j]])**2)**0.5
					if dist <= 9:
						q_contacts.append([q_list[i],q_list[j]])
			return q_contacts
	def writeGroAtomTypes(self,topfilename,atomtypes,pdbfile,rad):
		print (">> writing GROMACS atomtypes section",topfilename,atomtypes,pdbfile)
		#getting atom types
		(natomtypes,a_type) = self.nucAtomsType(pdbfile)
		#defining atomic mass
		atommass=np.ones(natomtypes,dtype=np.float32)
		#defining charge
		atomcharge = dict()
		charge_in_atomtypes = False #will be userdefined in future update
		for v in a_type:
			if charge_in_atomtypes:
				if v == "XP":
					atomcharge[v] = -1
				else:
					atomcharge[v] = 0
			else:
				atomcharge[v] = 0
		assert(len(atomcharge)==natomtypes)
		c6 = 0.000
		atomptype = ["A"]*natomtypes
		bases = dict()
		#Assiging radius N-bases based on pyrimidine and purine input radius
		bases["Bpy"] = ("T","C","U")
		bases["Bpu"] = ("A","G")
		for group,base in bases.items():
			for b in base:
				rad[b] = rad[group]
		atomname = list()
		for i in range(0,len(a_type)):
			atomname.append([a_type[i],atommass[1],atomcharge[a_type[i]],atomptype[i],c6,rad[a_type[i][1]]])
		f = open(topfilename, "a")
		assert atomtypes<=3
		a=atomname
		#writing section headrers
		f.write('%s\n'%('[ atomtypes ]'))
		f.write('%s\n' % ('; name mass  charge ptype c6    c12'))
		for i in (a):
			radius=(i[5]/10)**12 #C12 term
			i[0]=i[0].rjust(5)
			f.write('%s %8.3f %8.3f %s\t%e\t%e \n'%(i[0],i[1],i[2],i[3],i[4],radius))
		f.close()
		return a
	def writeGroAtoms(self,topfilename, pdbfile, nativefile, atomname,no_Pcharge):
		print ("writing atoms section")
		f1 = open(topfilename,"a")
		f1.write('\n%s\n'%('[ atoms ]'))
		f1.write('%s\n' % ("; nr  type  resnr  residue  atom  cgnr  charge  mass"))
		#P = PDBParser(PERMISSIVE=0)
		#structure = P.get_structure('test', nativefile)
		all = PrePDB().readPDB(nativefile)
		str1= 'atoms (atomnum, atomtype, resnum, resname, atomname, charge, mass)'
		natoms = len(all)
		l = {}
		iter = 0
		#atom_counter = dict()
		for i in all:
			#if i[12:16].lstrip().rstrip() not in atom_counter:
			#	atom_counter[i[12:16].lstrip().rstrip()] = 0
			#atom_counter[i[12:16].lstrip().rstrip()] = atom_counter[i[12:16].lstrip().rstrip()] + 1 
			iter = iter + 1
			l[0] = str(iter)
			l[1] = i[12:16].strip()
			l[2] = i[22:26].strip()
			l[3] = i[17:20].strip()
			l[4] = str(i[12:16]).strip()
			l[5] = str(iter)
			if i[12:16].lstrip().rstrip() == "XP":
				if not no_Pcharge:
					l[6] = str('%6.2f' % float(-1.000))
				elif no_Pcharge:
					l[6] = str('%6.2f' % float(0.000))
			else:
				l[6] = str('%6.2f' % float(0.000))
			l[7] = str('%6.2f' % float(1.000))
			f1.write("  ")
			for i in range(0,len(l)):
				f1.write(l[i]+"\t")
			f1.write("\n")
		f1.close()
		return all
	def writeGroPairs(self,topfilename, nativefile, pdbfile, contacttype, cutoff, atomname, btparams):
		print (">> writing GROMACS pairs sections",topfilename, nativefile, pdbfile,contacttype,cutoff,btparams)
    	#contacts_allowed=[1,2]
		#OPTIM works in kcal, GROMACS in kJ
		#Defining LJ-based stacking model Epsilon
		stackeps = PrePDB().stackEps()
		#definind Epsilon
		eps = 1.00
		print ("\t",end='')
		
		#Filling other eps values to stackpes
		for i in atomname:
			print (i[0].lstrip().rstrip()+"\t",end='')
		print ("\n",end='')
		for i in atomname:
			x = i[0].rstrip().lstrip()
			print (x,"\t",end='')
			for j in atomname:
				y = j[0].rstrip().lstrip()
				try:
					print (stackeps[x][y],"\t",end='')
				except:
					if x not in stackeps:
						stackeps[x] = dict()
					stackeps[x][y] = 1
					print (stackeps[x][y],"\t",end='')
			print ("\n",end='')
		#######################
		contacttype= 2	#10-12 potential used in this code
		contacts_allowed = [1,2] #1 for 6-10 and 2 for 10-12
		CG_contacts=self.getContacts(pdbfile,nativefile,cutoff,atomname,stackeps)
		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ pairs ]'))
		f1.write('%s\n' % ('; ai\taj\ttype\tA\tB'))
		for i in CG_contacts:
			epsilonij=float(i[4])
			B=((i[3]*0.1)**12)*5*epsilonij
			A=((i[3]*0.1)**10)*6*epsilonij
			f1.write('  %d\t%d\t%s\t%e\t%e\n' % (i[0],i[1],'1',A,B))
		#writing electrostatic pairs
		#q_contacts = self.chargeContacts(pdbfile,nativefile,cutoff)
		#for i in q_contacts:
		#	B=0; A=0 #no LJ_term needed
		#	f1.write('  %d\t%d\t%s\t%e\t%e\n' % (i[0],i[1],'1',A,B))
		#write exlusions
		f1.write('\n%s\n' % ('[ exclusions ]'))
		f1.write('%s\n' % ('; ai\taj'))
		for i in CG_contacts:
			f1.write('  %d\t%d\n'%(i[0],i[1]))
		f1.close()
		return
	def writeGroBonds(self,topfilename, nativefile, nKb, bond_ptype, CG_indices):
		print (">> writing GROMACS bonds section", topfilename, nativefile, nKb)
		#GROMACS 4.5.2 : FENE=7 AND HARMONIC=1
		allowed_pots = [1,7]
		bonds = list()
		for ai in range(1,len(CG_indices)+1):
			res1 = CG_indices[ai][0]
			atm1 = CG_indices[ai][1]
			xyz1 = CG_indices[ai][2]
			for aj in range(ai+1,len(CG_indices)+1):
				#aj>ai+1 avoids aj,ai as ai,aj allready exists
				res2 = CG_indices[aj][0]
				if res2-res1 in (0,1):
					atm2 = CG_indices[aj][1]
					xyz2 = CG_indices[aj][2]
					if (atm1[1]=="P" and atm2[1]=="S") or (atm1[1]=="S" and atm2[1]=="P"):
						bonds.append((ai,aj,(np.sum((xyz2-xyz1)**2))**0.5))
					elif (atm1[0]=="B" and atm2[1]=="S") or (atm1[1]=="S" and atm2[0]=="B"):
						bonds.append((ai,aj,(np.sum((xyz2-xyz1)**2))**0.5))
					else:
						continue
		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ bonds ]'))
		f1.write('%s\t%s\t%s\t%s\t%s\n' % ('; ai', 'aj', 'func', 'r0(nm)', 'Kb'))
		nKb = np.float_(nKb)*4.184	#1 Kcal/mol = 4.184 KJ/mol
		bonds.sort()
		for i in bonds:
			f1.write('  %d\t%d\t%d\t%e\t%e\n' % (i[0], i[1], bond_ptype, i[2]/10, nKb))
		f1.close()
		return bonds
	def writeGroAngles(self,topfilename, nativefile, nKa, CG_indices):
		print (">> writing GROMACS angle section", topfilename, nativefile, nKa)
		triplets = list()
		nKa = float(nKa)*4.184	#KCal/mol to KJ/mol
		#segregating indices based on atom
		resi_dict = dict()
		for key,value in CG_indices.items():
			resi = value[0]
			atom = value[1]
			xyz  = value[2]
			if resi not in resi_dict:
				resi_dict[resi] = dict()
			if atom[1] in ("P","S"):
				resi_dict[resi][atom[1]] = (key,xyz)
			elif atom[0] == "B":
				resi_dict[resi][atom[0]] = (key,xyz)

		triplets = list()
		for r in range(1,len(resi_dict)):
			try:#Set I P(i)S(i)B(i)
				triplets.append((resi_dict[r]["P"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r]["B"][0]-1))	#making indices start from 0
			except:
				buf = 0
			try:#Set II B(i)S(i)P(i+1)
				triplets.append((resi_dict[r]["B"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r+1]["P"][0]-1))	#making indices start from 0
			except:
				buf = 0
			try:#Set III P(i)S(i)P(i+1)
				triplets.append((resi_dict[r]["P"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r+1]["P"][0]-1))	#making indices start from 0
			except:
				buf = 0
			try:#Set IV S(i)P(i+1)S(i+1)
				triplets.append((resi_dict[r]["S"][0]-1,resi_dict[r+1]["P"][0]-1,resi_dict[r+1]["S"][0]-1))	#making indices start from 0
			except:
				buf=0

		Y = conmaps()
		CG_angles = Y.get_angles(nativefile,triplets)
		CG_angles = CG_angles[0]
		angles = list()
		for i in range(0,len(CG_angles)):
			 angles.append(((triplets[i]+np.array([1,1,1])).tolist(),CG_angles[i]))
		angles.sort()
		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ angles ]'))
		f1.write('%s\t%s\t%s\t%s\t%s\t%s\n' % ('; ai', 'aj', 'ak','func', 'th0(deg)', 'Ka'))
		for i in angles:
			f1.write('  %d\t%d\t%d\t%s\t%e\t%e\n' % (i[0][0],i[0][1],i[0][2],'1',i[1]*180/np.pi, nKa))
		f1.close()
		return angles
	def writeGroDihedrals(self,topfilename, nativefile, base_Kd, backbone_Kd, CG_indices):
		print (">> writing GROMACS dihedrals section", topfilename,nativefile,"Base Kd",base_Kd,"backbone Kd",backbone_Kd)
		#converting Kcal/mol to KJ/mol
		base_Kd = float(base_Kd)*4.184
		backbone_Kd = float(backbone_Kd)*4.184
		quadruplets_base = list()
		quadruplets_backbone = list()
		resi_dict = dict()
		for key,value in CG_indices.items():
			resi = value[0]
			atom = value[1]
			xyz  = value[2]
			if resi not in resi_dict:
				resi_dict[resi] = dict()
			if atom[1] in ("P","S"):
				resi_dict[resi][atom[1]] = (key,xyz)
			elif atom[0] == "B":
				resi_dict[resi][atom[0]] = (key,xyz)
		for r in range(1,len(resi_dict)):
			try: #Set I B(i)S(i)S(i+1)B(i+1)
				quadruplets_base.append((resi_dict[r]["B"][0]-1,resi_dict[r]["S"][0]-1,resi_dict[r+1]["S"][0]-1,resi_dict[r+1]["B"][0]-1))
			except:
				buf = 0
			try: #Set II P(i)P(i+1)P(i+2)P(i+3)				
				quadruplets_backbone.append((resi_dict[r]["P"][0]-1,resi_dict[r+1]["P"][0]-1,resi_dict[r+2]["P"][0]-1,resi_dict[r+3]["P"][0]-1))
			except:
				buf = 0
		
		Y = conmaps()
		CG_base_dihedrals = Y.get_dihedrals(nativefile,quadruplets_base)[0]
		#CG_backbone_dihedrals = Y.get_dihedrals(nativefile,quadruplets_backbone)[0]
		f1 = open(topfilename, "a")
		f1.write('\n%s\n' % ('[ dihedrals ]'))
		f1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('; ai','aj','ak','al','func','phi0(deg)','Kd','mult'))
		dihedrals = list()
		for i in range(0,len(quadruplets_base)):
			dihedrals.append(((quadruplets_base[i]+np.array([1,1,1,1])).tolist(),CG_base_dihedrals[i], base_Kd))
		for i in range(0,len(quadruplets_backbone)):
			dihedrals.append(((quadruplets_backbone[i]+np.array([1,1,1,1])).tolist(),np.pi, backbone_Kd))
		dihedrals.sort()
		for i in dihedrals:
			f1.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\t%s\n'% (i[0][0], i[0][1], i[0][2], i[0][3], '1 ', i[1]*180/np.pi, i[2], "1"))
			f1.write(' %d\t%d\t%d\t%d\t%s\t%e\t%e\t%s\n'% (i[0][0], i[0][1], i[0][2], i[0][3], '1 ', i[1]*180/np.pi, i[2], "3"))
		f1.close()
		return dihedrals
	def writeGromacsTop(self,topfilename,pdbfile,nativefile,btparams,rad,fconst,cutoff,no_Pcharge):
		#cheking if nucleotide exists
		checksize = open(pdbfile)
		if len(checksize.readlines()) <= 1:
			#if pdbfile size = 0, then there is no nucleotides im the sysytem
			#exiting with 0
			#print ("As mentioned earlier, the ibpu file is PROTEIN ONLY.....Exiting now")
			return 0
		checksize.close()

		print ("====>>>>> Generating Nucleotide CG files <<<<=========")
		#Order for SBM file
		#2 = 10-12 potential
		# Note: Parameters in GROMACS units need tweaking.
		#write gromacs format files.
		atomtypes = 3	#S,P,B
		contacttype = 2
		CG_indices = self.getIndices(nativefile)
		bond_ptype = 1
		self.write_gro_header(topfilename,atomtypes)  		#inheriting from eSBM
		self.write_header_SBM()						  		#inheriting from eSBM	
		atomname = self.writeGroAtomTypes(topfilename,atomtypes,pdbfile,rad)
		self.write_gro_moleculetype(topfilename)
		self.writeGroAtoms(topfilename, pdbfile, nativefile, atomname,no_Pcharge)
		self.writeGroBonds(topfilename, nativefile, fconst["nKb"], bond_ptype, CG_indices)
		self.writeGroAngles(topfilename, nativefile, fconst["nKa"], CG_indices)
		self.writeGroDihedrals(topfilename, nativefile, fconst["nKd"], fconst["P_nKd"], CG_indices)
		self.writeGroPairs(topfilename, nativefile, pdbfile, contacttype, cutoff, atomname, btparams)
		self.write_gro_tail(topfilename)
		print ('See file:',topfilename,'for GROMACS topology.')
		return 1
	def writeTablefile(self,debye,D,iconc,irad,T):
		print (">>Print table file",'table_nuc-aa.xvg')
		pi = np.pi
		irad = irad/10		#convertig A to nm
		if debye:
			#debye length
			#	 			e0*D*KB*T
			#dl**2 = -------------------------
			#		 2*el*el*NA*l_to_nm3*iconc

			e0 = 8.854e-21 			#C2 J-1 nm-1 #permitivity: The value have been converted into C2 J-1 nm-1
			D = D							#dielectric constant D
			KB = 1.3807*10**(-23) 	#J/K	#Boltzmann constant
			el = 1.6e-19			#C		#charge on electron
			NA = 6.022e+23          #n/mol  #Avogadro's number
			T = 298					#K		#Temperature
			l_to_nm3 = 1e-24		#nm3/l	#converting l to nm3

			dl_N = e0*D*KB*T	#numerator
			dl_D = 2*el*el*NA*l_to_nm3		#denom
			dl = (dl_N/dl_D)**0.5 #debye length nm M L-1
			#Note: still testing
			iconc = iconc #+ 0.5*(+12+7+3+7+12+6)/(6.022e+23*50e-24)
			dl = 1/dl			#nm-1 M-1 L
			dl = dl*(iconc**0.5) 		#nm-1
			Bk = m.exp(dl*irad)/(1+dl*irad)
		
		else:
			dl = 0
			Bk = 1
		Ke = 1/(4*pi*8.854e-21)                     #J nm C-2
		Ke = Ke*(2.56e-38)*(6.022e+23)*(1e-3)       #KJ nm C-2 e-2
		Ke=Ke*Bk/(D)
		fout = open('table_nuc-aa.xvg','w+')
		for i in range(0,29018/2+1):
			r = float(i)*0.002
			if r > 0.01:
				if r>=0.01:
					V 		= (Bk/D)*m.exp(-dl*r)/r
					V_1 	= (Bk/D)*(-m.exp(-dl*r)/r**2) + (Bk/D)*(-dl*m.exp(-dl*r)/r)
					V_1		= V_1*(-1)
				else:
					V = 0
					V_1 = 0
				C10 	= -1/r**10	
				C10_1	= -10/r**11 
				C12		= 1/r**12
				C12_1	= 12/r**13
				fout.write('%e %e %e %e %e %e %e\n' %(r,V,V_1,C10,C10_1,C12,C12_1))
			else:
				fout.write('%e %e %e %e %e %e %e\n' %(r,0,0,0,0,0,0))
		print ("testing",1/dl,"nm-1",iconc,"M")
		fout.close()
		return 	'table_nuc-aa.xvg'

def main():
	import argparse
	parser = argparse.ArgumentParser(description="Generate GROMACS and OPTIM potential files for enhanced SBM models.")
	parser.add_argument("--P_rad", help="Radius for Backbone Phosphate group bead. Default=3.7A")
	parser.add_argument("--S_rad", help="Radius for Backbone Sugar group bead. Default=3.7A")
	parser.add_argument("--Bpu_rad", help="Radius for N-Base Purine bead. Default=1.5A")
	parser.add_argument("--Bpy_rad", help="Radius for N-Base Pyrimidine bead. Default=1.5A")
	parser.add_argument("--nKb", help="Kbond for RNA/DNA")
	parser.add_argument("--nKa", help="Kangle for RNA/DNA. Default=20")
	parser.add_argument("--nKd", help="Kdihedral for Bi-Si-Si+1-Bi+1. Default=0.5")
	parser.add_argument("--P_nKd", help="Kdihedral for Backbone Pi-Pi+1-Pi+2-Pi+3. Default=0.7")
	parser.add_argument("--cutoff", help="Cut-off for contact-map generation")
	
	parser.add_argument("--Bpu_pos", help="Put input atom of Purine [N1,C2,H2-N2,N3,C4,C5,C6,O6-N6,N7,C8,N9,COM] as position of B. Default=COM(Center_of_Mass)")
	parser.add_argument("--Bpy_pos", help="Put input atom of Pyrimidine [N1,C2,O2,N3,C4,O4-N4,C5,C6,COM] as position of B. Default=COM(Center_of_Mass)")
	parser.add_argument("--S_pos", help="Put input atom of Sugar [C1',C2',C3',C4',C5',H2'-O2',O3',O4',O5',COM] as position of S   . Default=COM(Center_of_Mass)")
	parser.add_argument("--P_pos", help="Put input atom of Phosphate [P,OP1,OP2,O5',COM] group as position of P. Default=COM(Center_of_Mass)")
	
	parser.add_argument("--scaling", help="Scaling for mapping to all-atom contact-map.")
	parser.add_argument("--btmap", action='store_true', default=False, help='Use Betancourt-Thrumalai parameters for Bead interactions.')
	parser.add_argument("--w_native", help='Write native files, P-S-B from all atom PDB file.')
	parser.add_argument("--grotop",help='Gromacs topology file output name.')
	parser.add_argument("--pl_map", action='store_true', default=False, help='Plot contact map for two bead model')
	parser.add_argument("--aa_pdb", help='all-atom pdbfile e.g. 1qys.pdb')
	parser.add_argument("--all_chains", action='store_true', default=False, help='Will Not remove identical chains.')	
	parser.add_argument("--pdb2gro", help='Name for .gro file.')
	parser.add_argument("--pistacklen", help="pi-pi stacking length. Default=3.6A")

	parser.add_argument("--debye",action='store_true', help="Use debye electrostatic term")
	parser.add_argument("--T", help="System temperature. Default = 298K")
	parser.add_argument("--Kr", help="Krepulsion. Default=5.7A")
	parser.add_argument("--iconc", help="Solvant ion conc. Default=0.1M")  
	parser.add_argument("--irad", help="Solvant ion rad. Default=1.4A")  
	parser.add_argument("--dielec", help="Dielectric constant of solvant. Default=70")
	parser.add_argument("--no_Pcharge","-no_Pcharge", action='store_true', default=False, help='No negative charge on Phosphate bead')


	args = parser.parse_args()
	N = nucsbm()
	X = esbm()
	Y = conmaps()

	#Set default parameters
	no_Pcharge = False
	fconst = dict()					#force_csontant
	rad = dict()					#CG bead radii
	fconst["nKb"]= 200				#KCal/mol
	fconst["nKa"] = 20				#KCal/mol
	fconst["nKd"] = 0.5				#KCal/mol
	fconst["P_nKd"] = 0.7			#KCal/mol
	fconst["Kr"] = 5.7				#KCal/mol
	rad["P"] = 3.7					#A
	rad["S"] = 3.7					#A
	rad["Bpy"] = 1.5				#A
	rad["Bpu"] = 1.5				#A
	rad["stack"] = 3.6					#A
	cutoff = 4.5					#A
	irad = 1.4						#A (Na+ = 1.1A	Cl- = 1.7A)
	iconc = 0.01					#M
	D =	70							#dilectric constant (dimensionlsess)
	debye = False
	T = 120 #K
	pur_atom = ("N1","C2","H2-N2","N3","C4","C5","C6","O6-N6","N7","C8","N9","COM")
	pyr_atom = ("N1","C2","O2","N3","C4","O4-N4","C5","C6","H7-C7","COM")
	sug_atom = ("C1'","C2'","C3'","C4'","C5'","H2'-O2'","O3'","O4'","O5'","COM")
	phos_atom = ("P","OP1","OP2","O5'","COM")
	Bpu_pos = "COM"		#Center of Mass for purine
	Bpy_pos = "COM"		#Center of Mass for pyrimidine
	P_pos = "COM"			#Center of Mass for phosphate group
	S_pos = "COM"			#Center of Mass for sugar
	btparams = False	

	#Replacing default paramteres with input paramters
	if args.no_Pcharge:
		no_Pcharge = True
	if args.Bpu_pos:
		argsbuf = str(args.Bpu_pos)
		if argsbuf in pur_atom:
			Bpu_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Purine. The program will continue with default parameter")
	if args.Bpy_pos:
		argsbuf = str(args.Bpy_pos)
		if argsbuf in pyr_atom:
			Bpy_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Pyrimidine. The program will continue with default parameter")
	if args.P_pos:
		argsbuf = str(args.P_pos)
		if argsbuf in phos_atom:
			P_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Phosphate. The program will continue with default parameter")
	if args.S_pos:
		argsbuf = str(args.S_pos)
		if argsbuf in sug_atom:
			S_pos = argsbuf
		else:
			print("Warning!!! Wrong Atom name entered for Sugar. The program will continue with default parameter")
	
	CG_pos = {"Bpu":Bpu_pos,"Bpy":Bpy_pos,"S":S_pos,"P":P_pos}	
	print (CG_pos)
	#Force constants
	if args.nKb:
	    fconst["nKb"] = float(args.nKb)
	if args.nKa:
		fconst["nKa"] = float(args.nKa)
	if args.nKd:
		fconst["nKd"] = float(args.nKd)
	if args.P_nKd:
		fconst["P_nKd"] = float(args.P_nKd)

	if args.Kr:
		fconst["Kr"] = float(args.Kr)
	#CG radii
	if args.P_rad:
		rad["P"] = float(args.P_rad)
	if args.S_rad:
		rad["S"] = float(args.S_rad)
	if args.Bpu_rad:
		rad["Bpu"] = float(args.Bpu_rad)
	if args.Bpy_rad:
		rad["Bpy"] = float(args.Bpy_rad)
	if args.pistacklen:
		rad["stack"] = float(args.pistacklen)
	
	#solvant and ionic params
	if args.dielec:
		D = float(args.dielec)
	if args.iconc:
		iconc = float(args.iconc)
	if args.irad:
		irad = float(args.irad)
	sol_params = {"conc":iconc,"rad":irad,"D":D}
	if args.debye:
		debye = True
	if arg.T:
		T = float(args.T)
	if args.cutoff:
	    cutoff=float(args.cutoff)
	if args.btmap:
		btparams=True
	
	#cheking for read-write/input-output options
	grofile = ""
	if args.pdb2gro:
		grofile = str(args.pdb2gro)
	if args.w_native:
		pdbfile = args.w_native
		if not args.all_chains:
			pdbfile = PrePDB().homoNmer2Monomer(pdbfile)
		(pdbfile,terminal_residues) = PrePDB().appendChain(pdbfile)
		N.writeNative(pdbfile,CG_pos,grofile)
	if args.aa_pdb:
		pdbfile = args.aa_pdb
		if not args.all_chains:
			pdbfile = PrePDB().homoNmer2Monomer(pdbfile)
		(pdbfile,terminal_residues) = PrePDB().appendChain(pdbfile)
		N.writeNative(pdbfile,CG_pos,grofile)
	if args.grotop and args.aa_pdb :
		topfile = args.grotop
		nativefile="nuc_native_P-S-B.pdb"
		N.writeGromacsTop(topfile,pdbfile,nativefile,btparams,rad,fconst,cutoff,no_Pcharge)
	print ("done")
	


	#tobeCont
if __name__ == "__main__":
  main()


