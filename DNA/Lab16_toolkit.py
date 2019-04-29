from __future__ import print_function
import numpy as np 
import math as math
import sys


class PDB16_kit():
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
				#adjusting residue number if res_num doesn't start with 1
				starting_residue = int(i[22:26])
			if i.startswith("ATOM"):
				#creating a backbone list
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
		in_backbone = ["CA","C","O","N"]
		count = 0
		#total = len(backbone)-separation
		for r1 in range(0,len(backbone)-separation):
			#print ("Completed",int(100*float(count)/total),"%")
			count = count + 1
			#selecting residue 1
			for r2 in range(r1+separation,len(backbone)):
				#selecting residue 2
				dset = list()			#setting list of distances
				contact_flag = False	#no contact found
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
				dset = list()			#setting list of distances
				contact_flag = False	#no contact found
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
		backbone = dict(); sidechain = dict()
		backbone_CA = dict(); sidechain_CB = dict()
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
					#	backbone_CA[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
				else:
					sidechain[resi_num].append(np.float_([i[30:38],i[38:46],i[46:54]]))
					#if atomname == "CB":
					#	sidechain_CB[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
		fin = open("native_cb.pdb")
		for i in fin:
			resi_num = int(i[22:26])-starting_residue
			atomname = i[12:16].strip()
			if atomname == "CB":
				sidechain_CB[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
			elif atomname == "CA":
				backbone_CA[resi_num] = np.float_([i[30:38],i[38:46],i[46:54]])
		fin.close()

		assert(count_linres==len(backbone))
		assert(count_linres==len(sidechain))
		fin.close()
		contacts = list()
		count = 0
		#total = len(sidechain)
		for r1 in range(0,len(backbone)):
			#print ("Completed",int(100*float(count)/total),"%")
			count = count + 1
			#selecting residue 1
			for r2 in range(0,len(sidechain)):
				#selecting residue 2
				if r2-r1 >= separation or r1-r2 >= separation:
					dset = list()			#setting list of distances
					contact_flag = False	#no contact found
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
