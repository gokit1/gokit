from __future__ import print_function
from __future__ import division
import numpy as np
from util import Utils
import sys
sys.path.append('/usr/lib64/python2.7/site-packages')
import matplotlib.pyplot as plt
from gr import conmaps
import mdtraj as md


class tables(object):
    #class for generating tables as input for gromacs runs.
# Methods for energy and force calculations

    def get_10_12_lj(self,x,rcm):
        return  (5 * (rcm / x) ** 12) - (6 * (rcm / x) ** 10)
    def db_e1(self,x,rcm):
        return 5 * (rcm / x) ** 12 - 6 * (rcm / x) ** 10
    def db_e2(self,x,rcm,C,Y,rdb,epsdb):
        return (C * (Y ** 2) * (((Y ** 2) / 2) - (rdb - rcm) ** (2 * 2))) / (2 * 2) + epsdb
    def db_e3(self,x,rcm,B,Y,m,h1,h2):
        return -B * (Y - h1) / (Y ** m + h2)
    def db_f1(self,x,CO1,CO2):
        rm2=(1/x)**2;rm10=(1/x)**10
        return 60.0 * rm10 * (CO1 * rm2 / 5.0 - CO2 / 6.0)
    def db_f2(self,x,rcm,C,Y,k1,rdb):
        return (C) * (Y) * (Y ** 2 - k1) * (rdb - x)
    def db_f3(self,x,rcm,B,Y,m,h1,h2,rdb):
        return (B / ((Y ** m) + h2)) * (1 - (((Y - h1) * (m * Y ** (m - 1))) / ((Y ** m) + h2))) * (2 * (x - rdb))
    def hp_gaussian_e(self,x,eps,kappa,sigma_hp):
        k2=kappa*kappa
        rm12 = (sigma_hp/x) ** 12
        y = (x - sigma_hp - 1) ** 2
        #return eps*rm12 - k2*np.exp(-(1/2)*(y))
        return -k2*np.exp(-(1/2)*(y))
    def hp_gaussian_f(self,x,eps,kappa,sigma_hp):
        k2 = kappa * kappa
        y = (x - sigma_hp - 1) ** 2
        return 2 * k2 * np.sqrt(y) * np.exp(-(y ** 2) / 2) * y
    def gen_hp_gaussian(self,pottype,rcm,start,stop,step,suffix):
        x = np.arange(start, stop, step)
        sigma_hp=5.0 #Angstroms
        eps=1.0
        kappa=1.1*eps
        E1=self.hp_gaussian_e(x,eps,kappa,sigma_hp)
        E1=self.hp_gaussian_e(x,eps,kappa,sigma_hp)
        f1=self.hp_gaussian_f(x,eps,kappa,sigma_hp)
        E = E1
        f = f1
        table = np.array([x / 10, E, f])
        self.plot_table(x,E1,f1)
        return table

    def gen_db_table_file(self,pottype,hptrue,rcm,start, stop, step, suffix):
        # if hptrue:
        #     print ("Warning: Non-native hydrophobics will be added")
        assert (pottype=='db')
        if start == 0:
            print("set start to value slightly above zero.Say 0.02")
            exit()
        x = np.arange(start, stop, step)
        rssm = rcm + 3  # diameter of water molecule
        rdb = (rssm + rcm) / 2
        eps = 1
        epsdb = (0.1) * eps
        epsssm = (0.2) * eps
        k = 6
        m = 3
        n = 2
        # constants
        B = m * epsssm * (rssm - rdb) ** (2 * m - 2);
        h1 = (1 - (1 / m)) * ((rssm - rdb) ** 2) / ((epsssm / epsdb) + 1)
        h2 = (m - 1) * ((rssm - rdb) ** (2 * m)) / (1 + (epsdb / epsssm))
        C = (4 * n * (eps + epsdb)) / ((rdb - rcm) ** (4 * n))
        print('rcm,rssm,rdb,eps,epsssm,epsdb=', rcm, rssm, rdb, eps, epsssm, epsdb)
        #print (B,C,h1,h2)
        x1 = x[(x < rcm) & (x > 0)];x2 = x[(x >= rcm) & (x < rdb)];x3 = x[x >= rdb]
        CO1 = 5 * eps * (rcm ** 12);CO2 = 6 * eps * (rcm) ** 10
        E1 = self.db_e1(x1,rcm);f1 = self.db_f1(x1,CO1,CO2)
        Y2 = (x2 - rdb) ** 2
        k1 = (rdb - rcm) ** (4)
        E2 = self.db_e2(x2,rcm,C,Y2,rdb,epsdb);f2 = self.db_f2(x2,rcm,C,Y2,k1,rdb)
        Y3=(x3-rdb)**2
        E3 = self.db_e3(x3,rcm,B,Y3,m,h1,h2);f3 = self.db_f3(x3,rcm,B,Y3,m,h1,h2,rdb)
        E = np.hstack((E1, E2, E3))
        f=np.hstack((f1,f2,f3))
        if hptrue:
            sigma_hp = 5.0
            print("Adding gaussian hydrophobic term to potential, sigma_hp=", sigma_hp)
            kappa=1.1*eps
            ehp=self.hp_gaussian_e(x,eps,kappa,sigma_hp)
            fhp=self.hp_gaussian_f(x,eps,kappa,sigma_hp)
            E=E+ehp
            f=f+fhp

            #self.plot_table(x,Ehp,Fhp)
        assert len(E) == len(x)
        table = np.array([x / 10, E, f])
        #self.plot_table(x,E,f)
        #print('Number of entries in table file= ', len(x))
        # write table file.
        f1 = open('./MD/table_files/table_b' + suffix + '.xvg', 'w+')

        np.savetxt('./MD/table_files/table_b' + suffix + '.xvg', table.T, fmt='%10.3f %e %e')
        f1.close()
        return table

    def plot_table(self,x,E,f):
            # plt.style.use('ggplot')
            plt.rc('font', family='serif', size='20')
            zeroline = np.zeros(len(E))
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(1, 1, 1)
            plt.xlim([0, 4])
            plt.ylim([-2, 2])
            plt.xlim([1, 15])
            # plt.plot(x1,Eprime1)
            plt.plot(zeroline)
            ax.plot(x, E, color='k', ls='solid', label="E")
            ax.plot(x, f, color='b', ls='solid', label="Force")
            ax.plot(zeroline, color='g', ls='solid')
            plt.legend(loc=1, ncol=2, borderaxespad=0.2, fontsize=15)
            ax.set_xlabel('r')
            ax.set_ylabel('U(r)')
            plt.show()
            plt.close()
            return


    def gen_many_table_file(self, contactfile, nativefile, start, stop, step,hphobic):
        #print ('in gen_many_table')
        # typically something like contacts.txt generated from SMOG or --gconmap option.
        U=Utils()
        U.make_dir('MD');U.make_dir('MD/table_files')
        X = conmaps()
        pair = X.get_pairs_ext(contactfile)
        traj = md.load(nativefile)
        count = 0
        #print(len(pair))
        for i in pair:
            count = count + 1
            i = np.reshape(i, (1, 2))
            x = i[0][0];
            y = i[0][1]
            suffix = str(count)
            dist = md.compute_distances(traj, i)
            rcm = (dist[0][0] * 10)
            if hphobic:
                self.gen_db_table_file('db',True, rcm, start, stop, step, suffix)  # angstroms
            else:
                self.gen_db_table_file('db',False,rcm, start, stop, step, suffix)  # angstroms


#def main():
    #X=tables()
    #Y=conmaps()
    #X.gen_db_table_file('db', False,4, 0.01, 350,0.002,'peace')  # angstroms
    #X.gen_db_table_file('db', True, 4, 0.01, 350, 0.002, '')  # angstroms
    #X.gen_hp_gaussian('hpg',4,0.01,350,0.0002,'peace')
    #X.gen_fene('fene',20,2)

    #X.gen_many_table_file('smog.contacts.CG','native_ca.pdb',0.06,350,0.02)
#main()











