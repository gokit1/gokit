from __future__ import print_function
#from __future__ import division
from eSBM import esbm
from Bio import PDB
from util import Utils
import numpy as np
import os
#from plots import plot_1
import itertools
from gr import conmaps
from pdb_reres import _renumber_pdb_residue

class dswap(object):
    #Generate input for domain swapped structures. Input must be a domain-swapped dimer
    #Order parameters for dswap: 1) Qinter:
    def rotation_matrix(self):
        return PDB.rotaxis2m(90, PDB.Vector(-2, 10, 10))
    def translation_matrix(self):
        return np.array((0.5, 0.5, 0), 'f')
    def transform_pdb(self,monomer_pdb,rotated_pdb):
        #make translated and rotated copy of monomer from pdb
        parser = PDB.PDBParser()
        io = PDB.PDBIO()
        #rotation=PDB.rotmat(pi, PDB.Vector(1, 0, 0))
        # translation=PDB.Vector((0,0,1),'f')
        # print translation
        struct = parser.get_structure('monomer',monomer_pdb)
        rotation_matrix = self.rotation_matrix()
        translation_matrix = self.translation_matrix()
        for atom in struct.get_atoms():
            atom.transform(rotation_matrix, translation_matrix)
        io=PDB.PDBIO()
        chains = list(struct.get_chains())
        #rename transformed structure with chain B
        chains[0].detach_parent()
        chains[0].id = 'B'
        io.set_structure(struct)
        print ('Rotated and translated monomer in file:',monomer_pdb)
        print ('rot_mat:\n',rotation_matrix)
        print ('tran_mat:\n', translation_matrix)
        io.save(rotated_pdb)
        print(' pdb file saved in:', rotated_pdb)
        return True
    def merge_pdb(self,monomer_pdb,transformed_pdb,merged_pdb):
            self.transform_pdb(monomer_pdb,transformed_pdb)
            parser = PDB.PDBParser(); io = PDB.PDBIO()
            s1 = parser.get_structure('monomer', monomer_pdb)
            s2 = parser.get_structure('transformed',transformed_pdb)
            chains = list(s2.get_chains());chains[0].detach_parent()
            s1[0].add(chains[0]);io.set_structure(s1);io.save("tmp")
            self.renumber_res("tmp",merged_pdb)
            print ("Added chains")
            print("Renumbered residues. See:",merged_pdb)
            os.remove('tmp')
            return True

    def renumber_res(self,pdbfile,pdbout):
        #picked code from pdb_reres.py from pdbtools.
        #original pdbfile is rewritten.
        user_opts = {'resid': 1, 'chain': None}  # defaults if no opt is called
        pdbfh = open(pdbfile)
        opt_dict = {'resid': 1, 'chain': None}
        new_pdb = _renumber_pdb_residue(pdbfh, user_opts)
        with open(pdbout, 'w') as f:
            for x in new_pdb:
                f.write(str(x))
        f.close()
        return

    def write_contact_map(self,monomerfile):
        #pdbfile is for monomer. Dimer is automatically generated.
        X=esbm()
        contacts=X.write_contacts_section(monomerfile,'native_ca.pdb','mjmap.dat',4.5,1,False,False,True)
        return contacts

    def split_chains(self,dimerpdb):
        U=Utils()
        U.split_chain(dimerpdb)
        return True

    def get_Q(self,traj,nativefile,interpairfile,intrapairfile):
        #no. of inter-molecular contacts formed per frame. Read inter-chain contacts externally for now.
        X=conmaps();U=Utils()
        #read trajectoryfile only once.
        U.file_exists(interpairfile);U.file_exists(intrapairfile)
        #Get Q_inter v Q_intra
        interpairs=X.get_pairs_ext(interpairfile)
        intrapairs=X.get_pairs_ext(intrapairfile)
        inter_contacts=X.get_native_contacts(traj, nativefile, 'test.pdb',4.5,1.2, 'Q_inter', 100000, interpairs)
        intra_contacts=X.get_native_contacts(traj, nativefile, 'test.pdb',4.5,1.2, 'Q_intra', 100000, intrapairs)
        #calculate probability contact maps.
        X.get_contact_probability(traj,interpairfile,0,1000000,1)
        X.get_contact_probability(traj, intrapairfile, 0, 1000000, 1)
        #X.plot_scatter(inter_contacts,intra_contacts,'Qintervintra')
        return (inter_contacts,intra_contacts)

    def get_chain_con_map(self,all_atompdb,cutoff,scaling):
        #get intra-pairfile and read first n/2 lines. (These are intra-chain contacts for chain-A)
        X=conmaps();Y=esbm()
        self.split_chains(all_atompdb)
        print ("Calculating monomer contacts:")
        monomerfile='chain_A.pdb'
        Y.write_CB_to_native(monomerfile,False)
        import mdtraj as md
        traj=md.load_pdb("chain_A.pdb")
        natoms=X.get_total_number_of_atoms("native_ca.pdb")
        contacts=X.all_atom_contacts(monomerfile,4.5)
        pairs=contacts[:,:2];dist=contacts[:,3]
        assert len(pairs)==len(dist)
        #All intra-chain contacts:
        pairs_i1j1=pairs+natoms
        #AB
        print (pairs_i1j1)



        # for i in xrange(0,len(pairs)):
        #     print (pairs[i]+len(pairs))







        #get contact_map for chain-A


        return









    #
    #
    # def get_atom_atom_distance(self,traj,pair):
    #     #atom numbering starts with 0 in mdtraj
    #     # h1m1=18 h4m2=200 distL1=0.9124 nm distL4=2.44761nm
    #     X=conmaps()
    #     dist=X.get_distances(traj,pair)
    #
    #     d0=[];d1=[];d2=[];d3=[];d4=[];d5=[];d6=[];d7=[]
    #
    #     for i in xrange(0,len(traj)):
    #         d0 += [dist[i][0]]
    #         d1 += [dist[i][1]]
    #         d2 += [dist[i][2]]
    #         d3 += [dist[i][3]]  # important one.
    #         d4 += [dist[i][4]]
    #         d5 += [dist[i][5]]
    #         d6 += [dist[i][6]]
    #         d7 += [dist[i][7]]
    #
    #     d0=np.asarray(d0)
    #     d1=np.asarray(d1)
    #     d2=np.asarray(d2)
    #     d3=np.asarray(d3)
    #     d4=np.asarray(d4)
    #     d5=np.asarray(d5)
    #     d6=np.asarray(d6)
    #     X.plot_histogram(d0, '0')
    #     X.plot_histogram(d1, '1')
    #     X.plot_histogram(d2, '2')
    #     X.plot_histogram(d3, '3')
    #     X.plot_histogram(d4, '4')
    #     X.plot_histogram(d5, '5')
    #     X.plot_histogram(d6, '6')
    #     X.plot_histogram(d7, '7')
    #     return dist
    #
    # def H2Z_helices(self):
    #     #track helix-helix distance in H2Z through atom no. only!
    #     #h1=helix1,m1=molecule1 #atomnumbering starts from 0
    #     #taken from  "L1_swap.pdb"
    #     #Chain1 has 4 helices, Chain2 has 5 helices in crystal structure.
    #     h1m1=18;h1m2=141
    #     h2m1=66;h2m2=164
    #     h3m1=81;h3m2=184
    #     h4m1=109;h4m2=200;h5m2=229
    #     C1=[h1m1,h2m1,h3m1,h4m1];C2=[h1m2,h2m2,h3m2,h4m2,h5m2]
    #     pairs=list(itertools.product(C1,C2))
    #     print (len(pairs))
    #     return pairs
    #
    def plot_2d_hist(self, filename, title):
        #2d histogram e,g rog vs energy
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
#        filename='rog_energy_top7.dat'
        #filename='plots/M7_rge.dat'
        font = {'family': 'sans-serif',
                'size': 16,
                'weight': 'normal'
                }

        plt.rc('font', **font)
        plt.rcParams["axes.linewidth"] = 2
        plt.rcParams["axes.linewidth"] = 2
        fig, ax = plt.subplots()
        data = np.loadtxt(filename, dtype=float)
        x = data[:, 0];
        y = data[:, 1]
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.hist2d(x, y, bins=150, norm=LogNorm(), cmap=plt.cm.coolwarm)

        #plt.ylim((0.5,6))
        plt.xlabel('$Q_{inter}$', fontdict=font)
        plt.ylabel('$Q_{intra}$')
        plt.xlim(100,458)
        xtics = np.arange(100, 500, 50)
        plt.xticks(xtics);plt.yticks(xtics)
        plt.ylim(100, 458)
        plt.colorbar(label='Count')
        plt.plot(267, 280, '-o',color='k')
        plt.savefig(title + '.png')
        plt.show()
        return True


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Domain-swapping set up and analysis.")
    parser.add_argument("--dimerise",nargs=2,help="Convert file monomerpdb to  dimerpdb.")
    parser.add_argument("--split_chain", nargs=1, help="split dimer into its chains. Output: Chain*.pdb")
    parser.add_argument("--gconmap", nargs=1, help="Get dimer contact-map from monomer.")
    parser.add_argument('--traj_list', help='List of trajectories. Comma separated.No spaces. (e.g.: 1.xtc,2.xtc...)', type=str)
    parser.add_argument('--pair',nargs=2,help='write atom-pair (eg.1 10). Numbering starts from 1.',type=int)
    parser.add_argument('--h2zpairs', action='store_true', default=False,
                        help="H2Z specific. Temp arg.")
    parser.add_argument('--Q_int', action='store_true', default=False,
                        help="Get Qinter and Qintra for trajectory.")
    parser.add_argument('--terfile',help="Inter-molecular contacts file")
    parser.add_argument('--trafile', help="Intra-molecular contacts file")
    parser.add_argument("--pl_map", action='store_true', default=False,
                        help="Plot stuff")
    args = parser.parse_args()
    #X.transform_pdb(a,"rotated.pdb")
    # X.merge_pdb(a,"rotated.pdb",b)
    # X.get_Q('extract1.xyz','native_ca.pdb',interchainpair,intrachainpair)
    # pair=([18,201])
    X=dswap()
    #X.get_chain_con_map("2H2Z_Cter_full.pdb",4.5,1.2)
    Y=conmaps()
    U=Utils()
    from timeit import default_timer as timer
    t1=timer()
    # traj=md.load_pdb('interesting.pdb')
    #traj=Y.read_xyz_as_traj('xaf.xyz','native_ca.pdb')
    # print('Reading', len(traj), ' frames from xyz file')
    # pairs=X.H2Z_helices()
    # print (pairs)
    # distance = md.compute_distances(traj, pairs)
#     #U.file_exists('native_ca.pdb')
    #X.plot_2d_hist('2dhist.dat','test')
    if args.pair:
        at1=int(args.pair[0])-1;at2= int(args.pair[1])-1
        pair=([at1, at2])
        print ("Atom pairs=",pair)
    elif args.h2zpairs:
        pair=X.H2Z_helices()
        print ("All pairs of possiblle distances between H2Z helices being tracked")
    else:
        print ("----")
    # if args.traj_list and args.h2zpairs:
    #     traj = [item for item in args.traj_list.split(',')]
    #     for i in traj:
    #         t1=Y.read_xyz_as_traj(i,'native_ca.pdb')
    #         X.get_atom_atom_distance(t1,pair)
    #         #11= alpha1-alpha4'
    #         #Y.plot_histogram(distance[0][11]*10, i+'.hist')

    if args.traj_list and args.Q_int:
        print ("---------")
        terfile=args.terfile;trafile=args.trafile
        U.file_exists('native_ca.pdb');U.file_exists(terfile);U.file_exists(trafile)
        file = [item for item in args.traj_list.split(',')]
        for i in file:
            if U.check_file_extension(i, '.xyz'):
                traj=Y.read_xyz_as_traj(i, 'native_ca.pdb')
                print (traj)
            elif U.check_file_extension(i,'.xtc'):
                traj=Y.get_trajfromxtc(i,'native_ca.pdb')
            else:
                U.fatal_errors(1)
                exit()

            x,y=X.get_Q(traj,'native_ca.pdb',terfile,trafile)

            # if args.pl_map:
            #     plot_1.plot_2d_hist('Q_intercount.test','Q_intracount.test','peace')





    #    X.split_chains('dimer.pdb')
    if args.dimerise:
        assert len(args.dimerise)==2
        a=str(args.dimerise[0]);b=str(args.dimerise[1])
        X.transform_pdb(a,'transformed.pdb')
        X.merge_pdb(a,'transformed.pdb',b)

    t2=timer()
    print ("Wall-time:", t2-t1)

#print (len(G98))
    #print list(X.get_three_letters(G98))

main()