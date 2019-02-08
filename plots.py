import pylab
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm

class plot_1(object):
    #class for making publication quality plots with matplotlib

    def set_font(self):
        font = {'family': 'sans-serif',
                'size': 16,
                'weight': 'normal'
                }
        return font

    def read_twofl_columns(self,c1,c2,filename):
        #read two float columns from file.
        d1=np.loadtxt(filename, dtype=float);
        print (filename)
        x1=d1[:,c1];y1=d1[:,c2]
        return [x1,y1]




    def plot_roughness(self,f1,f2,f3,f4):
        #ebin 0.5,1,2,4 usually.
        print f1,f2,f3,f4
        d1 = np.loadtxt(f1, dtype=float);x1 = d1[:, 0];y1 = d1[:, 2]
        d2 = np.loadtxt(f2, dtype=float);x2 = d2[:, 0];y2 = d2[:, 2]
        d3 = np.loadtxt(f3, dtype=float);x3 = d3[:, 0];y3 = d3[:, 2]
        d4 = np.loadtxt(f4, dtype=float);x4 = d4[:, 0];y4 = d4[:, 2]
        font=self.set_font()
        plt.rc('font', **font)
        plt.rcParams["axes.linewidth"]=1
        n4_54_5=1360.30;n_tb=10;n_s6=70.09;n_m7=939.79
        #top7
        x1=(-x1+102)/2;x2=(-x2+51);x3=(-x3+26)*2;x4=(-x4+13)*4
        #tb
        # x1 = (-x1 + 290)/2;
        # x2 = (-x2 + 145);
        # x3 = (-x3 + 73)*2 ;
        # x4 = (-x4 + 37)*4
        # #S6
        # x1 = (-x1 + 112) / 2;
        # x2 = (-x2 + 56);
        # x3 = (-x3 + 28) * 2;
        # x4 = (-x4 + 14) * 4
        # M7
        # x1 = (-x1 + 160) / 2;
        # x2 = (-x2 + 80);
        # x3 = (-x3 + 40) * 2;
        # x4 = (-x4 + 20) * 4




        #y1=y1/n_tb;y2=y2/n_tb;y3=y3/n_tb;y4=y4/n_tb
        #y1=y1/n4_5;y2=y2/n4_5;y3=y3/n4_5;y4=y4/n4_5
        y1=y1/n_tb;y2=y2/n_tb;y3=y3/n_tb;y4=y4/n_tb
        y1 = y1 / n_m7;
        y2 = y2 / n_m7;
        y3 = y3 / n_m7;
        y4 = y4 / n_m7
        fig, ax = plt.subplots(figsize=(4, 2))
        xtics = np.arange(0,55,10)
        #xtics = np.arange(0,82,10)
        a1=ax.plot(x1,y1,'-',markersize=2,label='0.5')
        a2=ax.plot(x2, y2, '-', markersize=2, label='1')
        a3=ax.plot(x3, y3, '-', markersize=2, label='2')
        a4=ax.plot(x4, y4, '-', markersize=2, label='4')

        plt.xticks(xtics)
        plt.xlim((0,55 ))

        plt.legend(frameon=False,fontsize=10)
        plt.savefig('roughness' + '.eps', format='eps', dpi=1000)

        plt.show()

        return True


    def plot_prcontact_map(self,filename,file1):
        #plot probability contact-map
        font=self.set_font()
        plt.rc('font', **font)
        plt.rcParams["axes.linewidth"] = 1
        data = np.loadtxt(filename, dtype=float)
        data1 = np.loadtxt(file1,dtype=float)
        cm = matplotlib.cm.get_cmap('inferno_r')
        x = data[:, 0];y = data[:, 1];z = data[:,2]
        x1 = data1[:, 0];y1 = data1[:, 1]; z1 = data1[:, 2]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        sc = ax.scatter(x, y, c=z, vmin=0.20, vmax=1, marker='s', s=10, cmap=cm)
        ax.scatter(y1,x1, c=z1, vmin=0.20, vmax=1, marker='s', s=10, cmap=cm)
        #xtics=np.arange(0,92,20)
        #plt.xticks(xtics)
        plt.colorbar(sc)
        plt.grid()
        plt.savefig(filename+'.eps', format='eps', dpi=1000)
        plt.show()

    def plot_contact_map(self,f1,f2,f3,f4):
        #plot up to 5 files. #fig1 to show contact-maps.
        #files in SMOG.contacts format
        font=self.set_font()
        plt.rc('font', **font)
        plt.rcParams["axes.linewidth"] = 1
        #top half
        d1 = np.loadtxt(f1, dtype=float);x1 = d1[:, 1];y1 = d1[:, 3]
        d2 = np.loadtxt(f2, dtype=float);x2 = d2[:, 1];y2 = d2[:, 3]
        #bottom half
        d3 = np.loadtxt(f3, dtype=float);x3 = d3[:, 1];y3 = d3[:, 3]
        d4 = np.loadtxt(f4, dtype=float);x4 = d4[:, 1];y4 = d4[:, 3]
        #d5 = np.loadtxt(f5, dtype=float);x5 = d5[:, 0];y5 = d5[:, 1]
        # x3 = x3 + 1;
        # y3 = y3 + 1;
        # x4 = x4 + 1;
        # y4 = y4 + 1;
        #x5 = x5 + 1;
        #y5 = y5 + 1
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(1, 1, 1)

        #add diagonal
        diag1=np.arange(0,201,100)
        diag1=np.arange(0,100,10)

        a1 = ax.scatter(x1, y1, marker='s', s=10,color='#ff7f0e')
        a2 = ax.scatter(x2, y2, marker='s', s=20,color='#1f77b4')
        a3 = ax.scatter(y3, x3, marker='s', s=10,color='#ff7f0e')
        d1 = ax.plot(diag1,diag1, c="darkgray", marker=' ', linestyle='-',linewidth=2)

        a4 = ax.scatter(y4, x4, marker='s', s=20,color='#1f77b4')
        #a5 = ax.scatter(y5, x5, marker='s', s=25, color='black', facecolors='none')
        xtics=np.arange(0,101,40)
        ytics = np.arange(0, 101, 40)

        plt.xticks(xtics)
        plt.yticks(ytics)
        legend=plt.legend((a1,a2,a3,a4),('4.5','SCM'),fontsize=15,bbox_to_anchor=(-0.23,1.02,1,0.1), loc="lower left",ncol=5,frameon=False,handletextpad=0.0)
        plt.savefig('conmap'+'.eps', format='eps', dpi=1000)
        plt.show()

    def plot_scatter(self, x, y, title):
        # Simple x,y plot
        # X, Y are 1D numpy arrays
        plt.rcParams['backend'] = 'TkAgg'
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        # colors = ['k'] * len(x)
        ax.scatter(x, y, marker='o')
        ax.plot(x,y,'-o')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.savefig(title + '.png')
        plt.rc('font', family='serif', size='20')
        plt.show()
        print ('See: ', title + '.png')
        return True



    def plot_heat_map(self, x, y, z):
        #3d scatter plot.
        fig, ax = plt.subplots(figsize=(5, 5))
        # ax.scatter(x, y, c=z, s=150, marker='<', edgecolor='none')
        colors = ['red', 'blue']
        levels = [0, 1]
        cmap, norm = matplotlib.colors.from_levels_and_colors(levels=levels, colors=colors, extend='max')
        ax.scatter(x, y, c='gray', s=150, marker='<', edgecolor='none', cmap=cmap, norm=norm)
        plt.show()

    def plot_histogram(self,f1,f2,f3,f4,f5,f6,num_bins,title):
        #1d histogram for free-energy profile.
        self.set_font()
        plt.close('all')
        d1=np.loadtxt(f1);y1,x1=np.histogram((d1),bins=80,range=None,normed=True)
        d2=np.loadtxt(f2);y2,x2=np.histogram((d2),bins=80,range=None,normed=True)
        d3=np.loadtxt(f3);y3,x3=np.histogram((d3),bins=80,range=None,normed=True)
        d4=np.loadtxt(f4);y4,x4=np.histogram((d4),bins=80,range=None,normed=True)
        d5=np.loadtxt(f5);y5,x5=np.histogram((d5),bins=80,range=None,normed=True)
        d6=np.loadtxt(f6);y6,x6=np.histogram((d6),bins=80,range=None,normed=True)
        f, ((ax1, ax2,ax3,ax4,ax5,ax6)) = plt.subplots(6,1, sharex=True, sharey=False,figsize=(4,2))
        f.subplots_adjust(hspace=0.3)
        ax1.plot(x1[:-1],y1)
        ax2.plot(x2[:-1], y2)
        ax3.plot(x3[:-1], y3)
        ax4.plot(x4[:-1], y4)
        ax5.plot(x5[:-1], y5)
        ax6.plot(x6[:-1], y6)





        plt.savefig(title+'.eps', format='eps', dpi=1000)

        plt.show()

    def plot_2d_hist(self, filename, title):
        #2d histogram e,g rog vs energy
        from matplotlib.colors import LogNorm
#        filename='rog_energy_top7.dat'
        #filename='plots/M7_rge.dat'
        font=self.set_font()
        plt.rc('font', **font)
        plt.rcParams["axes.linewidth"] = 2
        plt.rcParams["axes.linewidth"] = 2
        fig, ax = plt.subplots()
        data = np.loadtxt(filename, dtype=float)
        x = data[:, 0];
        y = data[:, 2]
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.hist2d(x, y, bins=110, norm=LogNorm(), cmap=plt.cm.coolwarm)

        #plt.ylim((0.5,6))
        plt.xlabel('$Q_{inter}$', fontdict=font)
        plt.ylabel('$Q_{intra}$')
        plt.xlim(200,500);plt.ylim(350,500)
        xtics = np.arange(350, 500, 50)
        plt.xticks(xtics);plt.yticks(xtics)
        plt.colorbar(label='Count')
        #plt.plot(450, 450, '-o',color='k')
        #plt.plot(160, 308, '-o', color='k')
        plt.savefig(title + '.png')
        plt.show()
        return True

    def plot_lines(self,f1,f2,f3,f4,f5,f6,title):
        font = self.set_font()
        d1 = np.loadtxt(f1);x1 = d1[:, 0];y1 = d1[:, 1]
        d2 = np.loadtxt(f2);x2 = d2[:, 0];y2 = d2[:, 1]
        d3 = np.loadtxt(f3);x3 = d3[:, 0];y3 = d3[:, 1]
        d4 = np.loadtxt(f4);x4 = d4[:, 0];y4 = d4[:, 1]
        d5 = np.loadtxt(f5);x5 = d5[:, 0];y5 = d5[:, 1]
        d6 = np.loadtxt(f6);x6 = d6[:, 0];y6 = d6[:, 1]


        plt.rc('font', **font)
        plt.rcParams['backend'] = 'TkAgg'
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        # colors = ['k'] * len(x)
        #ax.scatter(x1, y1, marker='o')
        ax.plot(x1, y1+195.0197020000, '-',label='cut-off',color='#ff7f0e')
        ax.plot(x2, y2+199.0028432000, '--', label='cut-off+hp',color='#ff7f0e')
        ax.plot(x3, y3+254.9241883000, '-', label='SCM',color='#1f77b4')
        ax.plot(x4, y4+258.9067123000 ,'--', label='SCM+hp',color='#1f77b4')
        ax.plot(x5, y5+194.4277855000, '-', label='dsb',color='#2ca02c')
        ax.plot(x6, y6+198.4674764000, '--', label='dsb+hp',color='#2ca02c')
        plt.ylim((0, 200))
        plt.title(title)
        plt.savefig(title + '.eps', format='eps', dpi=1000)
        plt.rc('font', family='serif', size='20')
        plt.legend(frameon=False, fontsize=10)
        plt.show()
        return True

    def plot_3d(self):
        '''
        ======================
        3D surface (color map)
        ======================

        Demonstrates plotting a 3D surface colored with the coolwarm color map.
        The surface is made opaque by using antialiased=False.

        Also demonstrates using the LinearLocator and custom formatting for the
        z axis tick labels.
        '''

        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        import numpy as np

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make data.
        X = np.arange(-5, 5, 1)
        Y = np.arange(-5, 5, 1)
        X, Y = np.meshgrid(X, Y)
        R = np.sin(X ** 2 + Y ** 2)
        Z = np.sin(R)

        # Plot the surface.
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        # Customize the z axis.
        ax.set_zlim(-1.01, 1.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()

def main():
    X=plot_1()
    #X=plot_1()
    #X.plot_3d()
    #X.plot_3d()
    #X.plot_2d_hist('2dhist.dat', 'H2Z_2')
    #X.plot_2d_hist('test','test')
    X.plot_prcontact_map('MD_4_5_122', 'MD_dsb.map')

    #X.plot_prcontact_map('MD_4_5_122', 'MD_4_5_hp')
    #Top7 symmetric plot.
    #X.plot_contact_map('plots/4_5.map','plots/diff','plots/top7_SCM.map','plots/top7_SCM.map')

    #S6 and M7 plot.
    #X.plot_contact_map('plots/M7_4_5.map','plots/diff1M7','plots/S6_4.5','plots/diff_S6')

    #X.plot_contact_map('backbone.txt','bb_sc.txt','sidechain.txt','contacts.txt')
    folder=str('/home/sridhar/figures/histograms/')
    f1=folder+'S6_115.map'
    f2=folder+'S6_117_hp.map'
    f3=folder+'S6_135_6_5.map'
    f4=folder+'S6_139_6_5_hp.map'
    f5=folder+'S6_dsb.map'
    f6=folder+'tb_S6.map'

    #X.plot_histogram(f1,f2,f3,f4,f5,f6,100,'title')
    #X.plot_roughness('plots/4_5_hp/eps0.5', 'plots/4_5_hp/eps1', 'plots/4_5_hp/eps2', 'plots/4_5_hp/eps4')
    #X.plot_roughness('4_5/eps0.5','plots/4_5/eps1','plots/4_5/eps2','plots/4_5/eps4')
    #X.plot_roughness('plots/S6/eps0.5', 'plots/S6/eps1', 'plots/S6/eps2', 'plots/S6/eps4')
    #X.plot_roughness('plots/M7/eps0.5', 'plots/M7/eps1', 'plots/M7/eps2', 'plots/M7/eps4')
    #X.plot_roughness('plots/eps1', 'plots/tb/eps1', 'plots/tb/eps2', 'plots/tb/eps4')
    #x=np.loadtxt('plots/tb/eps1')[:,0];y=np.loadtxt('plots/tb/eps1')[:,2]

    #X.plot_scatter(x,y,'peace')
    #X.plot_roughness('plts/4_5_hp/eps0.5', 'plots/tb/eps1', 'plots/tb/eps2', 'plots/tb/eps4')
    #folder='/media/Data/sridhar/top7/new/trap_pdb/'
    #X.plot_lines(folder+'4_5/A2/EofS','4_5_hp/A2/EofS','6_5/A2/EofS','6_5_hp/A2/EofS','dsb/A2/EofS.9','dsb_hp/A2/EofS.8','A2 ')

    # X.plot_lines(folder + '4_5/A3/EofS', '4_5_hp/A3/EofS', '6_5/A3/EofS', '6_5_hp/A3/EofS', 'dsb/A3/EofS',
    #               'dsb_hp/A3/EofS', 'A3')
    # X.plot_lines(folder + '4_5/A4/EofS', '4_5_hp/A4/EofS', '6_5/A4/EofS', '6_5_hp/A4/EofS', 'dsb/A4/EofS',
    #               'dsb_hp/A4/EofS', 'A4')
    # X.plot_lines(folder + '4_5/A1/EofS', '4_5_hp/A1/EofS', '6_5/A1/EofS', '6_5_hp/A1/EofS', 'dsb/A1/EofS.8',
    #               'dsb_hp/A1/EofS.1', 'A1')
    # X.plot_lines(folder + '4_5/tb/EofS.14630', '4_5_hp/tb/EofS.14630', '6_5/tb/EofS.14630', '6_5_hp/tb/EofS.14630', 'dsb/A1/EofS.8',
    #               'dsb_hp/A1/EofS.8', 'tb14630')
    # X.plot_lines(folder + '4_5/tb/EofS.16908', '4_5_hp/tb/EofS.16908', '6_5/tb/EofS.16908', '6_5_hp/tb/EofS.16908',
    #               'dsb/tb/EofS.16908','dsb_hp/tb/EofS.5.16908', 'tb16908')
    # X.plot_lines(folder + '4_5/tb/EofS.35586', '4_5_hp/tb/EofS.35586', '6_5/tb/EofS.35586', '6_5_hp/tb/EofS.35586',
    #             'dsb/tb/EofS.2.35586','dsb_hp/tb/EofS.35586', 'tb35586')
    # #X.plot_lines(folder + '4_5/tb/EofS.16908', '4_5_hp/tb/EofS.16908', '6_5/tb/EofS.16908', '6_5_hp/tb/EofS.16908',
    #             #'dsb/A1/EofS.16908', 'dsb_hp/tb/EofS.16908', 'tb16908')


main()
