import numpy as np
import readtreeHDF5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker

import illustris_python as il
"""
In this code be plot the stellar formation history for BCGs, where M_crir_200 > 5e14 Msun. 
To these BCGs I compute the median and percentile to 16 and 84, all this for different epocs. 

"""






#========================TNG300-1===========================#
ubication = 'IllustrisTNG'
simulation = 'L205n2500TNG'
simulation_name = 'IllustrisTNG300'
snapnum = 99; subfind_id = 1

basedir = '/n/hernquistfs3/%s/Runs/%s/output' % (ubication, simulation)

header = il.groupcat.loadHeader(basedir, snapnum)
nsubs = header['Nsubgroups_Total']
h = header['HubbleParam']

# Centrals only
group_first_sub = il.groupcat.loadHalos(basedir, snapnum, fields=['GroupFirstSub'])
sub_gr_nr = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloGrNr'])
nsubs = len(sub_gr_nr)
is_central = group_first_sub[sub_gr_nr] == np.arange(nsubs, dtype=np.uint32)

# Get M200 for the subhalos' parent halos
group_m200 = il.groupcat.loadHalos(basedir, snapnum, fields=['Group_M_Crit200'])
sub_m200 = group_m200[sub_gr_nr]

#====== only M200 > 5e14 ========#
locs = is_central & (sub_m200 > 5e14 / 1e10 * h)
subfind_ids = np.flatnonzero(locs)

#======================MERGER TREES=======================#
treedir = '/n/hernquistfs3/IllustrisTNG/Runs/L205n2500TNG/postprocessing/trees/SubLink_gal'
tree =  readtreeHDF5.TreeDB(treedir)

#========= Initialization from figure 
fig = plt.figure(figsize=(6,5))
ax = plt.subplot()

fila = 99  # 99 is the maximun of data of branch for mass
columna = len(subfind_ids)


mstellar_snap = np.zeros((fila , columna) )
mhalo_snap = np.zeros((fila , columna) )


#====== Charge main branch for the BGCs and add to the matrix vacuum
for k in range(len(subfind_ids)):

    branch = tree.get_main_branch(snapnum, subfind_ids[k], keysel = ['SubhaloMassType', 'SubhaloMass', 'SnapNum'])
        
    #=========================
    snapnums = branch.SnapNum[:]
        
    #==========================
    mhalo = branch.SubhaloMass[:]
    mhalo_Msun = mhalo* 1e10 / h
        
                
    #==========================
    mstar = branch.SubhaloMassType[:,4] 
    mstar_Msun = mstar * 1e10 / h
        
    #print(snapnums)
    
    for j in range(len(snapnums)): 
            
        mhalo_snap[:,k][j] = mhalo_Msun[j]
        mstellar_snap[:,k][j] = mstar_Msun[j]
    
    #============== PLOT ====================
    z = np.loadtxt('RedshiftsIllustrisTNG.txt')
    l = len(snapnums)
    print('size= %f '%l)
    # count k
    #k = snapnum -2
    # invert to z for tha snapnum = 99 --> z=0
    z = np.flip(z[-l:])

    #ax.xaxis.set_ticks_position('bottom')
    #plt.plot(snapnums, mstar)
    ax.plot(1.+z, mstar_Msun, c = 'b', linewidth = 0.2 )
    ax.plot(1.+z, mhalo_Msun,  c = 'g', linewidth = 0.2)


    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(0.9,10.5)
    
#======= create to array 
Mhalo_media = np.zeros(fila)
Mstellar_media = np.zeros(fila)

mhalo_percentile_16 = np.zeros(fila)
mhalo_percentile_84 = np.zeros(fila)

mstellar_percentile_16 = np.zeros(fila)
mstellar_percentile_84 = np.zeros(fila)


#========= Compute to median and percentiles
for i in range(fila):
    
    Mhalo_media[i] = np.median(mhalo_snap[i,:])
    Mstellar_media[i] = np.median(mstellar_snap[i,:])
    
    mhalo_percentile_16[i] = np.nanpercentile(mhalo_snap[i,:], 16)
    mhalo_percentile_84[i] = np.nanpercentile(mhalo_snap[i,:], 84)
    
    mstellar_percentile_16[i] = np.nanpercentile(mstellar_snap[i,:], 16)
    mstellar_percentile_84[i] = np.nanpercentile(mstellar_snap[i,:], 84)
    

#====== Plot 
ax.plot(1.+z, Mhalo_media, label = r'$M_{halo}$', c = 'g')
ax.plot(1.+z, Mstellar_media, label = r'$M_{stellar}$', c = 'b')


ax.fill_between(1.+z, mhalo_percentile_16, mhalo_percentile_84, alpha = 0.3, label = r'percentile 16-84')
ax.fill_between(1.+z, mstellar_percentile_16, mstellar_percentile_84, alpha = 0.3, label = r'percentile 16-84')

ax.axvline(3. , min(mstar_Msun), max(mhalo_Msun), linewidth = 0.5, c = 'black' )
plt.text(3.1, 5e14, r'$z = 2$', fontsize = 8)

ax.set_xlabel(r'$Redshift,\ 1 + z$')
ax.set_ylabel(r'$Mass\ [M_{\odot}]$')

plt.title(r'BCGs ($M_{200} > 5\times 10^{14} M_{\odot}$)')
plt.legend(loc=3, fontsize = 7, frameon = False , fancybox = True)
plt.subplots_adjust(wspace = 0.05, bottom = 0.19, right = 0.9, left = 0.2, top = 0.93)
plt.xticks([1,2,3,4,5,6,7,8,9,10], ["1","2","3","4","5","6","7","8","9","10"])#, '11', '12'] )

for ending in ['png', 'pdf']:
        #fig.savefig('/n/home10/dmontenegro/PythonModules/graficas/BCGS_more_massive_z%s.%s' %(int(z) , ending))
        fig.savefig('/n/home10/dmontenegro/PythonModules/graficas/merger_tree/formation_history_Mstar_median.%s' %( ending))
        


"""
for k in range(len(subfind_ids)):

    branch = tree.get_main_branch(snapnum, subfind_ids[k], keysel = ['SubhaloMassType', 'SubhaloMass', 'SnapNum'])

    #==========================
    mhalo = branch.SubhaloMass[:]
    mhalo_Msun = mhalo* 1e10 / h
    
 

    #==========================
    mstar = branch.SubhaloMassType[:,4] 
    mstar_Msun = mstar * 1e10 / h
    #median = np.median()
        

    #=========================
    snapnums = branch.SnapNum[:]


    #=========================
    z = np.loadtxt('RedshiftsIllustrisTNG.txt')
    l = len(snapnums)
    print('size= %f '%l)
    # count k
    #k = snapnum -2
    # invert to z for tha snapnum = 99 --> z=0
    z = np.flip(z[-l:])





    #ax.xaxis.set_ticks_position('bottom')
    #plt.plot(snapnums, mstar)
    ax.plot(1.+z, mstar_Msun, c = 'b', linewidth = 0.2 )
    ax.plot(1.+z, mhalo_Msun,  c = 'g', linewidth = 0.2)



    ax.set_xscale('log')
    ax.set_yscale('log')

#ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
#ax.set_xticks([1,2,3,4,5,6,7,8,9,10], ["1","2","3","4","5","6","7","8","9","10"] )
    ax.set_xlim(0.00001,10.5)
"""
