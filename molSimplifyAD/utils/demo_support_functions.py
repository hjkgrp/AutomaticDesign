import numpy as np
import scipy as sp
import pandas as pd
import sklearn as sk
from sklearn import decomposition
from matplotlib.pyplot import figure, imshow, axis
from matplotlib.image import imread
from molSimplifyAD.utils.report_tool.prepare_report import *
from molSimplifyAD.ga_init import *
from IPython.display import Image
from pymol.cgo import *
from pymol import cmd
from pymol import preset
from pymol import util
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from IPython import display
from scipy.spatial import ConvexHull
from scipy.interpolate import UnivariateSpline
from matplotlib import gridspec
def draw3d_big(load_path, save_path):
    cmd.delete('all')

    ## generate axes
    w = 0.06 # cylinder width 
    l = 5 # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter

    # start coords, end coords, w is radius, color1 with grad to color2 in rgb
    obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           CYLINDER, 0.0, 0.0, 0.0,   -l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   -l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   -l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
          ]
    # load axes into picture
    cmd.load_cgo(obj, 'axes')

    # load structure and define properties 
    cmd.load(load_path)
    cmd.bg_color(color="white")
    preset.ball_and_stick(selection='all', mode=1)

    util.cbaw
    cmd.color("gray","symbol c")
    cmd.select('hats', "symbol h")
    cmd.set("sphere_scale",0,'hats')
    cmd.set("stick_h_scale",1)
    cmd.set("stick_radius",0.70)

    #cmd.zoom('all', 10, 0, 0)
    #cmd.zoom('all',5,0,1)
    cmd.set('depth_cue', 0)

    # https://www.andre-gaschler.com/rotationconverter/
    view = [0.8535534, -0.3535534,  0.3826834,       0.5205981,  0.5497281, -0.6532815,       0.0205981,  0.7568349,  0.6532815,         0.000000000,    0.000000000,  -22.395387650,         0.000147820,    0.000066757,    0.000249624,        18.378971100,   26.411804199,  -25.000000000 ]

    cmd.set_view(view)
    cmd.zoom('all',0.5,0,0)

    cmd.png(save_path + '1.png', 350, 500, 500, 0, 0)
    
    #view2 from the other side
    for i in [2,3,4,5,6,7,8]:
        view[i] = view[i]*(-1)
    
    cmd.set_view(view)
    cmd.zoom('all',0.5,0,0)
    cmd.png(save_path + '2.png', 350, 500, 500, 0, 0)

def data_normalize(data, train_mean, train_var):
    data = data.astype(float)  # Make sure the data is always in float form
    d = np.shape(train_mean)[0]

    # print('normalizing with number of dimensions = ' +str(d))
    ### double check the variance in the training data
    delete_ind = list()
    # print(train_var)
    for idx, var in enumerate(np.squeeze(train_var)):
        if var < 1e-16:
            delete_ind.append(idx)
    if len(delete_ind)>0:
        print('Note: There are %d features with a variance smaller than 1e-16.' % len(delete_ind))
        print('Please double check your input data if this number is not what you expect...')
        data = np.delete(data, delete_ind, axis=1)
        train_mean = np.delete(train_mean, delete_ind, axis=0)
        train_var = np.delete(train_var, delete_ind, axis=0)
    # print(data.shape, train_mean.shape, train_var.shape)
    scaled_dat = np.divide((data.T - train_mean), np.sqrt(train_var), ).T
    return (scaled_dat)

def showImagesHorizontally(list_of_files):
    fig = figure(figsize = (20,20))
    number_of_files = len(list_of_files)
    for i in range(number_of_files):
        a=fig.add_subplot(1,number_of_files,i+1)
        image = imread(list_of_files[i])
        imshow(image,cmap='Greys_r')
        axis('off')
def find_gene_geo(new_tree,mad_path,target_gene):
    ## function to find the geo given a gene, no matter the gen it 
    ## was introduced
    for gen in range(0,new_tree.status_dictionary["gen"]+1):
        check_string = mad_path+'/initial_geo/gen_'+str(gen) + '/*'+target_gene +'*.xyz'
        matches = glob.glob(check_string)

        if len(matches) >0:
            break
    if len(matches)>0:
        return(matches[0])
    else:
        print('No match, error!  cannot find geo for gene ' + target_gene)
        return(False)
def get_info_from_ga(mad_path):
    ## function to load object from a specified mad path
    current_wd = os.getcwd()
    os.chdir(mad_path)
    GA_run = GA_run_defintion()
    GA_run.deserialize(os.getcwd() + '/.madconfig')
    new_tree = GA_generation('current_gen')
    ## read in info
    new_tree.read_state()
    os.chdir(current_wd)
    return(new_tree)
def translate_gene_names_into_chemistry(mad_path,genes):

    ## function to load object from a specified mad path
    current_wd = os.getcwd()
    os.chdir(mad_path)
    
    for g in genes:
        print('translating gene into words : ' + g)
        
        
    GA_run = GA_run_defintion()
    GA_run.deserialize(os.getcwd() + '/.madconfig')
    new_tree = GA_generation('current_gen')
    ## read in info
    new_tree.read_state()
    os.chdir(current_wd)
    return(new_tree)
def get_live_genes_sorted(new_tree):
    ## find genes that still live and have fitnessess 
    live_genes = dict()
    for gene in set(new_tree.genes.values()):
        if gene in new_tree.gene_fitness_dictionary.keys():
            live_genes.update({gene:new_tree.gene_fitness_dictionary[gene]})

    ## order them by fitness
    live_genes = sorted(live_genes.items(), key=operator.itemgetter(1),reverse=True)
    return(live_genes)

def visualize_best_genes(new_tree,mad_path,live_genes):
    ## make images of the top 3
    pics_to_show = [[],[]]

    if len(live_genes) == 1:
        lgs = [live_genes[0],live_genes[0],live_genes[0]]
    elif len(live_genes) == 2:
        lgs = [live_genes[0],live_genes[1],live_genes[1]]
    else:
        lgs = live_genes[0:3]
    ## loop over live geners
    for i, lg in enumerate(lgs):
        temp_path = os.getcwd()+ '/temp_' + str(i)+'_'
        orig_path = find_gene_geo(new_tree=new_tree,mad_path=mad_path,target_gene=lg[0])
        
        # make pictures
        draw3d_big(orig_path,temp_path)
        pics_to_show[0].append(temp_path+'1.png')
        pics_to_show[1].append(temp_path+'2.png')
    return(pics_to_show)
def update_data(new_tree,data,live_genes):
    ## update fitness
    mod_dat = copy.deepcopy(data)
    #print('have ' + str(len(new_tree.gene_fitness_dictionary.keys()))+' keys in tree')
    for gene in new_tree.gene_fitness_dictionary.keys():
        if gene in mod_dat["gene"].values:
            mod_dat.loc[mod_dat["gene"]==gene,"fitness"] = float(new_tree.gene_fitness_dictionary[gene])
            #print(gene,mod_dat.loc[mod_dat["gene"]==gene,"fitness"])
        else:
            print('gene ' + str(gene) + ' not in tree!')
    live_genes_to_match = [ live_genes[i][0] for i in range(0,len(live_genes))]
    
        
    live_db_inds = mod_dat.index[mod_dat["gene"].isin(live_genes_to_match)]
    all_inds = mod_dat.index[mod_dat["gene"].isin(new_tree.gene_fitness_dictionary.keys())]
    return(mod_dat,live_db_inds,all_inds)
def get_run_info(mad_path):
    gens = []
    mfs = []
    divs = []
    with open(mad_path+'statespace/_generation_meanFitness_diversity.csv') as f:
        for l in f.readlines():
            ll = l.strip().split(',')
            gens.append(float(ll[0]))
            mfs.append(float(ll[1]))
            divs.append(float(ll[2]))
    return gens, mfs, divs

def draw_double_plot(Xr,mod_dat,live_db_inds,all_inds,hull,gens,mfs,divs):
    font = {'weight' : 'bold',
            'size'   : 18}

    plt.rc('font', **font)
    #fig, axs = plt.subplots(nrows=1,ncols=2,figsize = (18,9))
    fig = figure(figsize=(20,10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.75, 1]) 
    ax1 =  plt.subplot(gs[0])
    ax2 =  plt.subplot(gs[1])
    font = {'weight' : 'bold',
            'size'   : 18}

    plt.rc('font', **font)
    #colors = [(0.3, 0.3, 0.3),(0.95,0.3,1)]  
    cmstr = 'winter'
    #cm = LinearSegmentedColormap.from_list(
    #    'atoms', colors,N=10)
    points  = Xr[:,0:2]
    for simplex in hull.simplices:
        ax1.plot(points[simplex, 0], Xr[simplex, 1], 'k-')

    o1 = ax1.tricontourf(Xr[all_inds,0], Xr[all_inds,1],mod_dat.loc[all_inds,'fitness'],
                   5,alpha=0.1,cmap=cmstr,linestyles='dashed', vmin=0, vmax=1)
    for c in o1.collections:
        c.set_edgecolor("face")

    o2 = ax1.scatter(Xr[all_inds,0], Xr[all_inds,1], c=mod_dat.loc[all_inds,"fitness"],
            cmap=cmstr,s=100,
                alpha=1, lw=1,marker=".",label='previous', vmin=0, vmax=1)
    o3 = ax1.scatter(Xr[live_db_inds,0], Xr[live_db_inds,1], c=mod_dat.loc[live_db_inds,"fitness"],
            cmap=cmstr,s=125,
                alpha=1.0, lw=3,marker="s",label='current', vmin=0, vmax=1,
                    edgecolor='black', linewidth=2)
    cbar = fig.colorbar(o1,ax=ax1)
    cbar.ax.set_ylabel('fitness',fontweight='bold')
    #cbar.add_lines(o1)
    ax1.set_xlabel('PC1',fontweight='bold')
    ax1.set_ylabel('PC2',fontweight='bold')
    ax1.legend()
    ax1.set_xlim(-19,21)
    ax1.set_ylim(-17,27)
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])

    

    ## second plot
    ax3 = ax2.twinx()
    ax2.plot(gens, mfs, 'go',markersize=6)
    ax3.plot(gens, divs, 'bs',markersize=6)

    if len(mfs) > 3:
        xnew = np.linspace(0,max(gens),500)
        mfs_spl = UnivariateSpline(gens, mfs,s=100)
        divs_spl = UnivariateSpline(gens, divs,s=100)
        ax2.plot(xnew, mfs_spl(xnew), 'g--',lw=4.5)
        ax3.plot(xnew, divs_spl(xnew), 'b--',lw=4.5)


    ax2.set_xlabel('generation',fontweight='bold')
    ax2.set_ylabel('mean fitness', color='g',fontweight='bold')
    ax3.set_ylabel('diversity', color='b',fontweight='bold')
    plt.setp(ax2.get_yticklabels(), color="g")
    plt.setp(ax3.get_yticklabels(), color="b")
    ax2.yaxis.set_tick_params(color='b')

    font = {'weight' : 'bold',
            'size'   : 18}

    plt.rc('font', **font)

    plt.show()
def save_double_plot(impath, Xr,mod_dat,live_db_inds,all_inds,hull,gens,mfs,divs):
    font = {'weight' : 'bold',
            'size'   : 12}

    plt.rc('font', **font)
    #fig, axs = plt.subplots(nrows=1,ncols=2,figsize = (18,9))
    fig = figure(figsize=(16,8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.75, 1]) 
    ax1 =  plt.subplot(gs[0])
    ax2 =  plt.subplot(gs[1])
    font = {'weight' : 'bold',
            'size'   : 12}

    plt.rc('font', **font)
    #colors = [(0.3, 0.3, 0.3),(0.95,0.3,1)]  
    cmstr = 'winter'
    #cm = LinearSegmentedColormap.from_list(
    #    'atoms', colors,N=10)
    points  = Xr[:,0:2]
    for simplex in hull.simplices:
        ax1.plot(points[simplex, 0], Xr[simplex, 1], 'k-')

    o1 = ax1.tricontourf(Xr[all_inds,0], Xr[all_inds,1],mod_dat.loc[all_inds,'fitness'],
                   5,alpha=0.1,cmap=cmstr,linestyles='dashed', vmin=0, vmax=1)
    for c in o1.collections:
        c.set_edgecolor("face")

    o2 = ax1.scatter(Xr[all_inds,0], Xr[all_inds,1], c=mod_dat.loc[all_inds,"fitness"],
            cmap=cmstr,s=100,
                alpha=1, lw=1,marker=".",label='previous', vmin=0, vmax=1)
    o3 = ax1.scatter(Xr[live_db_inds,0], Xr[live_db_inds,1], c=mod_dat.loc[live_db_inds,"fitness"],
            cmap=cmstr,s=125,
                alpha=1.0, lw=3,marker="s",label='current', vmin=0, vmax=1,
                    edgecolor='black', linewidth=2)
    cbar = fig.colorbar(o1,ax=ax1)
    cbar.ax.set_ylabel('fitness',fontweight='bold')
    #cbar.add_lines(o1)
    ax1.set_xlabel('PC1',fontweight='bold')
    ax1.set_ylabel('PC2',fontweight='bold')
    ax1.legend()
    ax1.set_xlim(-19,21)
    ax1.set_ylim(-17,27)
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])

    

    ## second plot
    ax3 = ax2.twinx()
    ax2.plot(gens, mfs, 'go',markersize=6)
    ax3.plot(gens, divs, 'bs',markersize=6)

    if len(mfs) > 3:
        xnew = np.linspace(0,max(gens),500)
        mfs_spl = UnivariateSpline(gens, mfs,s=100)
        divs_spl = UnivariateSpline(gens, divs,s=100)
        ax2.plot(xnew, mfs_spl(xnew), 'g--',lw=2.5)
        ax3.plot(xnew, divs_spl(xnew), 'b--',lw=2.5)


    ax2.set_xlabel('generation',fontweight='bold')
    ax2.set_ylabel('mean fitness', color='g',fontweight='bold')
    ax3.set_ylabel('diversity', color='b',fontweight='bold')
    plt.setp(ax2.get_yticklabels(), color="g")
    plt.setp(ax3.get_yticklabels(), color="b")
    ax2.yaxis.set_tick_params(color='b')

    font = {'weight' : 'bold',
            'size'   : 12}

    plt.rc('font', **font)
    
   
    fig.set_size_inches(16, 8)
    fig.savefig(impath+'-plot.png', dpi=600)
