# Chin plotting program 
#import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
#matplotlib.rcParams['text.usetex'] = True

A = np.ones(100)
A = A*14
#A = [14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
A = [14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14,14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
print(len(A))

#files = ['/Users/sebastiangrorud/GitHub_projects/ProjectWork_TaylorLin/Parameter/new_n_50_normal.dat']
files = ['/Users/sebastiangrorud/matlab/Til rapport/CHIN_funker/new_n2.dat']

#files = ['/Users/sebastiangrorud/matlab/Til rapport/CHIN_plotting/n_5.dat']
df = pd.read_fwf(files[0], sep='  ',header=None, widths=A)
df = df.apply(pd.to_numeric, errors='coerce')

df.fillna(0,inplace=True)
df

#x = np.linspace(0,100, 7)
#labels = ['0', '0.25', '0.50', '0.75', '1.00', '1.25', '1.50']
x = np.linspace(0,100, 9)
labels = ['0', '0.25', '0.50', '0.75', '1.00', '1.25','1.50','1.75', '2.00']

# =============================================================================
# #plt.figure()
# fig, ax = plt.subplots(dpi=1200)
# plt.title(r'$n = 2$')
# ax.pcolormesh(df, shading='auto')
# plt.xlabel(r'$\zeta_1 = {\tau_c^{\{321\}\langle 111 \rangle}} / \tau_c^{\{110\}\langle 111 \rangle}$')
# plt.ylabel(r'$\zeta_2 {\tau_c^{\{211\}\langle 111 \rangle}} / \tau_c^{\{110\}\langle 111 \rangle}$')
# plt.xticks(x, labels, rotation='horizontal', size=8)
# plt.yticks(x, labels, rotation='horizontal', size=8)
# #plt.savefig('n_5_large_returnmap.png', dpi=2000)
# plt.show()
# 
# =============================================================================


#=============================================================================
fig, ax = plt.subplots(figsize=(5,5), dpi=1200)
#sns.color_palette("coolwarm", as_cmap=True)
plt.title(r'$n = 2$')
cmap = sns.diverging_palette(220, 20, s=60, as_cmap=True)
ax = sns.heatmap(df, vmin=4.0, vmax=15.0, cmap=cmap, cbar=False)
ax.invert_yaxis()
plt.xlabel(r'$\zeta_2 = {\tau_c^{\{321\}\langle 111 \rangle}} / \tau_c^{\{110\}\langle 111 \rangle}$')
plt.ylabel(r'$\zeta_1 = {\tau_c^{\{211\}\langle 111 \rangle}} / \tau_c^{\{110\}\langle 111 \rangle}$')
plt.xticks(x, labels, rotation='horizontal', size=12)
plt.yticks(x, labels, rotation='horizontal', size=12)

for ax, spine in ax.spines.items():
    spine.set_visible(True)
#plt.savefig('n_1000.png', format='png')
plt.show()

#=============================================================================
# $\zeta_1 = {\tau_c^{\{321\}\langle 111 \rangle}} / \tau_c^{\{110\}\langle 111 \rangle} $