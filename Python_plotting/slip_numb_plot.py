# Chin plotting program 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#A = [14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
A = [14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
print(len(A))

#files = ['/Users/sebastiangrorud/GitHub_projects/ProjectWork_TaylorLin/Parameter/new_n_50_normal.dat']
files = ['/Users/sebastiangrorud/GitHub_projects/ProjectWork_TaylorLin/Parameter/Numb_slip_n100.dat']

#files = ['/Users/sebastiangrorud/matlab/Til rapport/CHIN_plotting/n_5.dat']
df = pd.read_fwf(files[0], sep='  ',header=None, widths=A)
df = df.apply(pd.to_numeric, errors='coerce')

df.fillna(0,inplace=True)

#x = np.linspace(0,100, 7)
#labels = ['0', '0.25', '0.50', '0.75', '1.00', '1.25', '1.50']
x = np.linspace(0,100, 9)
labels = ['0', '0.25', '0.50', '0.75', '1.00', '1.25','1.50','1.75', '2.00']





#=============================================================================
fig, ax = plt.subplots(figsize=(7.5,6), dpi=600)
plt.title(r'$n = 100$',size=14)
cmap = sns.color_palette("cubehelix", as_cmap=True)
ax = sns.heatmap(df,cmap=cmap,vmin=5.0, vmax=14.0, cbar=True)

ax.invert_yaxis()
plt.xlabel(r'$\zeta_2 = {\tau_c^{\{321\}\langle 111 \rangle}} / \tau_c^{\{110\}\langle 111 \rangle}$',size=16)
plt.ylabel(r'$\zeta_1 = {\tau_c^{\{211\}\langle 111 \rangle}} / \tau_c^{\{110\}\langle 111 \rangle}$',size=16)
plt.xticks(x, labels, rotation='horizontal', size=12)
plt.yticks(x, labels, rotation='horizontal', size=12)

for ax, spine in ax.spines.items():
    spine.set_visible(True)
#plt.savefig('numb_slip_n20.png', format='png')
plt.show()
