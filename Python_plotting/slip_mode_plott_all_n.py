# Slip mode plott for different 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob


A = [14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14]
path = '/Users/sebastiangrorud/matlab/Resultater/fiber_data/def90/MODE_data/'
file = glob.glob(path + '*.dat')
file = sorted(file)

x_1 = np.zeros(3)
x_2 = np.zeros(3)
x_3 = np.zeros(3)
x_4 = np.zeros(3)
x_5 = np.zeros(3)
xf_1 = np.zeros(48)
xf_2 = np.zeros(48)
xf_3 = np.zeros(48)
xf_4 = np.zeros(48)
xf_5 = np.zeros(48)

#files = [file[9], file[10],file[11]]
scale1 = [0.975, 0.99, 1.00]
scale2 = [0.99, 1.00, 1.01]
n_scale = [20, 50, 100, 1000]

ii = 0;
for n in range(0,1):
    for m in range (0,1):
        for i in range (0,5):
            tot = 0
            
            df = pd.read_fwf(file[ii+i], sep='  ',header=None, widths=A)
            df = df.apply(pd.to_numeric, errors='coerce')
            df.fillna(0,inplace=True)
            
            
            temp = np.zeros(96)
            temp = df.sum()
        
            x_res = np.zeros(48)
            for t in range (0,47):
                x_res[t] = temp[t] + temp[48+t]
        
            tot = sum(x_res[:])
            x_res[:] = x_res[:]/tot
            x_res
        
            x_plot = [sum(x_res[0:11]), sum(x_res[12:23]), sum(x_res[24:])]
            
            if i == 0:
                x_1 = x_plot
                xf_1 = x_res
                print(xf_1)
            if i == 1:
                x_2 = x_plot
                xf_2 = x_res
                print(xf_2)
            if i == 2:
                x_3 = x_plot
                xf_3 = x_res
                print(xf_3)
                
            if i == 3:
                x_4 = x_plot
                xf_4 = x_res
                print(xf_4)
                
            if i == 4:
                x_5 = x_plot
                xf_5 = x_res
                print(xf_5)
        ii = ii + 6
        print(ii)
        
        
        ## SINGEL PLOT
        # fig = plt.figure(1)
        # ax = fig.add_axes([0,0,1,1])
        # ax.bar(Y, x_plot)
        # ax.set_ylabel('Relative amount of shear by each slip mode')
        # ax.set_xlabel('Slip mode')
        # ax.set_title('Activated slip mode at 90% deformation')
        # ax.set_xticks(x)
        # ax.set_xticklabels(Y)
        # ax.legend()
        # plt.show()
        
        patterns = [ "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ]
        
        ##### MULTIPLOT
        Y = ['{110}<111>', '{211}<111>', '{321}<111>']
        x = np.arange(len(Y))  # the label locations
        width = 0.8  # the width of the bars
        
        fig, ax = plt.subplots()
         #fig.suptitle('Activated slip mode at $\epsilon_{vm}=20$ deformation', ha="center",rasterized=True)
        rects1 = ax.bar(x - width/2.5, x_1, width/5, label=r'$n=5$', edgecolor='grey', color= [(0.8980392156862745, 0.7686274509803922, 0.5803921568627451)], rasterized=True)
        rects2 = ax.bar(x - width/5, x_2, width/5, label=r'$n=20$', edgecolor='grey', color=[(0.4, 0.7607843137254902, 0.6470588235294118)], rasterized=True)
        rects3 = ax.bar(x, x_3, width/5, label=r'$n=50$', edgecolor='grey', color=[(0.9882352941176471, 0.5529411764705883, 0.3843137254901961)])
        rects4 = ax.bar(x + width/5, x_4, width/5, label=r'$n=100$', edgecolor='grey', color=[(0.5529411764705883, 0.6274509803921569, 0.796078431372549)],rasterized=True)
        rects5 = ax.bar(x + width/2.5, x_5, width/5, label=r'$n=1000$', edgecolor='grey', color= [(0.7019607843137254, 0.7019607843137254, 0.7019607843137254)], rasterized=True)
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Relative amount of shear by each slip mode')
        ax.set_xlabel('Slip mode', size=12)
        #ax.set_title(r'$\epsilon_{vm}=3.66$   $\zeta_1 = $'+str(scale1[n])+'   $\zeta_2 = $'+str(scale2[m]), ha="center")
        ax.set_xticks(x)
        ax.set_xticklabels(Y)
        ax.set_ylim([0.0, 1.0])
        ax.legend()
        
        #ax.bar_label(rects1, padding=3)
        #ax.bar_label(rects2, padding=3)
        name = 'slip_mode_full_90_'+str(ii)
        #plt.savefig(name+'.eps', format='eps')
        
