from molSimplify.Classes import globalvars
from molSimplify.Classes import mol3D
from molSimplify.Scripts.io import *
from molSimplify.Informatics.autocorrelation import*
from molSimplify.Informatics.misc_descriptors import*
from molSimplify.Informatics.RACassemble import*
from molSimplify.Informatics.graph_analyze import*
from molSimplifyAD.ga_tools import*
from sklearn.model_selection import train_test_split, cross_val_score, LeaveOneOut
from sklearn.linear_model import LinearRegression
import scipy.stats as stats


###############################################################################################
# Below are the functions that are used to process the run results post dataframes from mAD.  #
# The run_results_processor just takes a runs_results_post file and matches things to get the #
# desired properties. A separate class (below) can be used to match complexes from two        #
# dataframes. Should work with external data as well, but needs to be adapted.                #
###############################################################################################

class run_results_processor:
    #### Class to take a run results post file and process it in various ways
    def __init__(self,frame,path_to_save=False, name = 'FirstRow', process_type='SSE',geo_dict = False):
        self.frame = frame.replace('undef',np.nan)
        self.frame = self.frame.replace('lig_mismatch',np.nan)
        #### This is done specifically for cases where we know the ligand mapping can be different. Makes it consistent.
        self.frame[['lig1','lig2','lig3','lig4','lig5','lig6']] = self.frame[['lig1','lig2','lig3','lig4','lig5','lig6']].replace(['co','cn'],['carbonyl','cyanide'])
        self.frame['name_without_HFX'] = (self.frame['metal']+'_'+self.frame['ox'].astype(str)+'_'+
                                          self.frame['lig1']+'_'+self.frame['lig2']+'_'+self.frame['lig3']+'_'+
                                          self.frame['lig4']+'_'+self.frame['lig5']+'_'+self.frame['lig6']+'_'+
                                          self.frame['spin'].astype(str))
        self.frame['chem_name'] = (self.frame['metal']+'_'+self.frame['ox'].astype(str)+'_'+
                                          self.frame['lig1']+'_'+self.frame['lig2']+'_'+self.frame['lig3']+'_'+
                                          self.frame['lig4']+'_'+self.frame['lig5']+'_'+self.frame['lig6']+'_ahf_'+
                                          self.frame['alpha'].astype(int).astype(str).str.zfill(2)+'_'+self.frame['spin'].astype(str))
        # sort_values puts most negative at top, so drop_duplicates will take first
        self.frame['energy'] = self.frame['energy'].apply(float)
        self.frame = self.frame.sort_values(['energy'],ascending=True)
        #### Path_to_save is the raw path such as the SI path, like this: /Dropbox/Paper/SI/.
        self.path_to_save = path_to_save
        self.process_type = process_type
        self.geo_dict = geo_dict
        #### This name should be descriptive.
        self.name = name
        
    def process_duplicates(self):
        # This function takes the duplicated DFT runs and takes the ones with the lowest energy.
        print('Initial dataframe shape: '+str(self.frame.shape))
        self.frame = self.frame.drop_duplicates(subset=['name_without_HFX','alpha'],keep='first')
        print('After dropping duplicates dataframe shape: '+str(self.frame.shape))
        
    def check_convergence(self):
        # This function eliminates the non converged jobs
        print('---- checking convergence now ----')
        print('Initial dataframe shape: '+str(self.frame.shape))
        dropped_frame = self.frame[self.frame['converged'] == False]
        self.frame = self.frame[self.frame['converged'] == True]
        print('After dropping convergence dataframe shape: '+str(self.frame.shape))
        return dropped_frame
    
    def geometry_check(self, custom_geo_dict = {}):
        print('---- checking geometries now ----')
        if len(custom_geo_dict) > 0:
            constraint_dictionary = custom_geo_dict
        else:
            if self.geo_dict == False:
                print('geo check dict not provided. Using defaults.')
                constraint_dictionary = {'num_coord_metal': 6,
                                         'rmsd_max': 0.3, 'atom_dist_max': 0.45,
                                         'oct_angle_devi_max': 12, 'max_del_sig_angle': 22.5,
                                         'dist_del_eq': 0.35, 'dist_del_all': 1,
                                         'devi_linear_avrg': 20, 'devi_linear_max': 28}
            else:
                constraint_dictionary = self.geo_dict
        print('Initial dataframe shape: '+str(self.frame.shape))
        dropped_frame = self.frame[~((self.frame['num_coord_metal'].astype(float)==constraint_dictionary['num_coord_metal'])&
                                   (self.frame['dist_del_eq'].astype(float)<=constraint_dictionary['dist_del_eq'])&
                                   (self.frame['dist_del_all'].astype(float)<= constraint_dictionary['dist_del_all'])&
                                   (self.frame['max_del_sig_angle'].astype(float)<= constraint_dictionary['max_del_sig_angle'])&
                                   (self.frame['oct_angle_devi_max'].astype(float)<= constraint_dictionary['oct_angle_devi_max'])&
                                   (self.frame['devi_linear_avrg'].astype(float)<= constraint_dictionary['oct_angle_devi_max'])&
                                   (self.frame['devi_linear_max'].astype(float)<= constraint_dictionary['devi_linear_max'])&
                                   (self.frame['rmsd_max'].astype(float)<= constraint_dictionary['rmsd_max']))
                                  ]
        self.frame = self.frame[((self.frame['num_coord_metal'].astype(float)==constraint_dictionary['num_coord_metal'])&
                                   (self.frame['dist_del_eq'].astype(float)<=constraint_dictionary['dist_del_eq'])&
                                   (self.frame['dist_del_all'].astype(float)<= constraint_dictionary['dist_del_all'])&
                                   (self.frame['max_del_sig_angle'].astype(float)<= constraint_dictionary['max_del_sig_angle'])&
                                   (self.frame['oct_angle_devi_max'].astype(float)<= constraint_dictionary['oct_angle_devi_max'])&
                                   (self.frame['devi_linear_avrg'].astype(float)<= constraint_dictionary['oct_angle_devi_max'])&
                                   (self.frame['devi_linear_max'].astype(float)<= constraint_dictionary['devi_linear_max'])&
                                   (self.frame['rmsd_max'].astype(float)<= constraint_dictionary['rmsd_max']))
                                  ]
        print('After dropping bad geo dataframe shape: '+str(self.frame.shape))
        print('----- Eliminated complexes due to geometric reason -----')
        return dropped_frame
    
    def homoleptic_geo_check(self,**kwargs): #### More stringent checks on homoleptic complexes
        print('---- checking homoleptic geometries now, making some checks more stringent ----')
        print('Initial dataframe shape: '+str(self.frame.shape))
        #### Made the dist_del_all check more stringent for homoleptic complexes
        dist_del_all_stringent = 0.5
        dist_del_eq_stringent = 0.2
        frame_copy = self.frame.copy().reset_index()
        remove_list = []
        full_row_list = range(0,len(frame_copy))
        for i, row in frame_copy.iterrows():
            if row['lig1'] == row['lig2'] == row['lig3'] == row['lig4'] == row['lig5'] == row['lig6']:
                if len(kwargs) > 0:
                    valid = True
                    for key, value in kwargs.items():
                        if row[key] > value:
                            valid = False
                            break
                    if not valid:
                        remove_list.append(i)
                else:
                    if (float(row['dist_del_all']) <= dist_del_all_stringent) and (float(row['dist_del_eq'])<=dist_del_eq_stringent):
                        continue
                    else:
                        remove_list.append(i)
            else:
                continue
        keep_list = list(set(full_row_list)-set(remove_list))
        dropped_frame = self.frame.iloc[remove_list]
        self.frame = self.frame.iloc[keep_list]
        print('After dropping homoleptic bad geo dataframe shape: '+str(self.frame.shape))
        print('----- Eliminated complexes due to more stringent checks on homoleptic complexes reason -----')
        return dropped_frame
        
    def spin_contamination_check(self, limit=1.0):
        print('---- checking spin contamination now ----')
        print('Initial dataframe shape: '+str(self.frame.shape))
        dropped_frame = self.frame[~(abs(self.frame['ss_act']-self.frame['ss_target'])<=limit)]
        self.frame = self.frame[(abs(self.frame['ss_act']-self.frame['ss_target'])<=limit)]
        print('After dropping spin contamination problems dataframe shape: '+str(self.frame.shape))
        print('----- Eliminated complexes due to spin contamination -----')
        return dropped_frame
    
    def empty_site_spin_contamination_check(self, limit=1.0): 
        print('---- checking spin contamination on empty site structures now ----')
        print('Initial dataframe shape: '+str(self.frame.shape))
        dropped_frame = self.frame[~((abs(self.frame['ss_act']-self.frame['ss_target'])<=limit)&
                                     (abs(self.frame['empty_ss_act']-self.frame['empty_ss_target'])<=limit))]
        self.frame = self.frame[(abs(self.frame['ss_act']-self.frame['ss_target'])<=limit)&
                                (abs(self.frame['empty_ss_act']-self.frame['empty_ss_target'])<=limit)]
        print('After dropping spin contamination problems dataframe shape: '+str(self.frame.shape))
        print('----- Eliminated complexes due to empty site spin contamination -----')
        return dropped_frame
    
    def spin_on_metal_check(self, limit=1.0): #### make sure to make this mode dependent (check empty ss)
        print('---- checking spin on metal now ----')
        print('Initial dataframe shape: '+str(self.frame.shape))
        dropped_frame = self.frame[~(abs(self.frame['net_metal_spin'].astype(float)-self.frame['spin'].astype(float)+1)<=limit)]
        self.frame = self.frame[(abs(self.frame['net_metal_spin'].astype(float)-self.frame['spin'].astype(float)+1)<=limit)]
        print('After dropping spin on metal problems dataframe shape: '+str(self.frame.shape))
        print('----- Eliminated complexes due to spin deviation from metal reason -----')
        return dropped_frame
    
    def get_ligand_energies(self):
        #### for LDE calculation
        ligand_energy_dict = {'acetonitrile':{0: -132.7384119779, 5: -132.7434703535, 10: -132.7487033151, 15: -132.7541064976, 20: -132.7596762509, 25: -132.7654090649, 30: -132.7712976021},
                         'ammonia':{0: -56.5383490619, 5: -56.5409735769, 10: -56.5436540385, 15: -56.5463883877, 20: -56.5491739809, 25: -56.5520095004, 30: -56.5548934166},
                         'carbonyl':{0: -113.2997886779 ,5: -113.3021967727,10:-113.3046900110 , 15: -113.3072685627, 20:-113.3099311791 , 25: -113.3126770301, 30: -113.3155056961},
                         'cyanide':{0: -92.8127841681,5: -92.8158354883,10: -92.8189530560, 15: -92.8221367255, 20: -92.8253845504, 25: -92.8286970857, 30: -92.8320722262},
                         'chloride':{0: -460.2048337746, 5: -460.2166452956,10: -460.2284826219, 15: -460.2403454252, 20: -460.2522334988, 25: -460.2641465881, 30: -460.2760844366},
                         'fluoride':{0: -99.7387107901, 5: -99.7425421104,10: -99.7463850300, 15: -99.7502395272, 20: -99.7541052427, 25: -99.7579818729, 30: -99.7618692971}, 
                         'misc': {0: -132.6995165613, 5: -132.7047649512,10: -132.7101997561, 15: -132.7158135271, 20: -132.7215990693, 25: -132.7275587014, 30: -132.7336798586},
                         'phosphine':{0: -343.1025218022, 5: -343.1124519769,10: -343.1224556638, 15: -343.1325324474,20: -343.1426794186, 25: -343.1528949569, 30: -343.1631777603},
                         'uthiol':{0: -399.3453112154, 5: -399.3556080928,10: -399.3659614254, 15: -399.3763700887,20: -399.3868329008, 25: -399.3973487455, 30: -399.4079166978},
                         'water': {0: -76.4009336829, 5: -76.4029944534,10: -76.4051137093, 15: -76.4072890001,20: -76.4095176518, 25: -76.4117985134, 30: -76.4141300174}}
        return ligand_energy_dict

    def limit_considered_symmetry(self, frame, symmetry=['trans','homoleptic']): # cis and five+1 limited due to ambiguity in RACs
        # This function takes in a frame, and only returns the allowed symmetries
        limited_frame = frame.copy()
        print('originally, frame shape: '+str(limited_frame.shape))
        symm_list = []
        for i, row in limited_frame.iterrows():
            if (row['lig1']==row['lig2'])&(row['lig1']==row['lig3'])&(row['lig1']==row['lig4']): #eq symm
                if (row['lig4']==row['lig5'])&(row['lig4']==row['lig6']):
                    symm_list.append('homoleptic')
                elif ((row['lig4']==row['lig5']) or (row['lig4']==row['lig6']))&(row['lig5']!=row['lig6']):
                    symm_list.append('5+1')
                else:
                    symm_list.append('trans')
            else:
                symm_list.append('cis')
        limited_frame['symm'] = symm_list
        dropped_frame = limited_frame[limited_frame['symm'].isin(symmetry)]
        dropped_frame.reset_index()
        print('after dropping specific symmetries, frame shape: '+str(dropped_frame.shape))
        return dropped_frame

    def limit_considered_ligands(self, frame, ligands=['ammonia','water','phosphine','uthiol','carbonyl','misc','acetonitrile']): #neutral ligands only
        # This function takes in a frame, and only returns the allowed ligands for dissociation (i.e. lig 6)
        limited_frame = frame.copy()
        print('originally, frame shape: '+str(limited_frame.shape))
        dropped_frame = limited_frame[limited_frame['lig6'].isin(ligands)]
        dropped_frame.reset_index()
        print('after dropping negatively charged ligands in lig6 position, frame shape: '+str(dropped_frame.shape))
        return dropped_frame
    
    def get_LDE(self, frame = False, symmetry = False, ligands = False):
        # if symmetry and ligands are provided (as lists), only those symmetries and ligands are considered for LDE
        ligand_energy_dict = self.get_ligand_energies()
        name_list = []
        name_list_no_HFX = []
        HFX_list = []
        LDE_list = []
        no_pair_list = []
        if not frame:
            pair_frame = self.frame.copy()
            pair_frame = pair_frame.reset_index()
        else:
            pair_frame = frame.reset_index()
        if symmetry:
            pair_frame = self.limit_considered_symmetry(pair_frame, symmetry=symmetry)
        if ligands:
            pair_frame = self.limit_considered_ligands(pair_frame, ligands=ligands)
        for i, row in pair_frame.iterrows():
            if i % 5000 == 0:
                print('counter',i)
            if ((row['lig1'] not in ligand_energy_dict.keys()) or (row['lig2'] not in ligand_energy_dict.keys()) or 
                (row['lig3'] not in ligand_energy_dict.keys()) or (row['lig4'] not in ligand_energy_dict.keys()) or 
                (row['lig5'] not in ligand_energy_dict.keys()) or (row['lig6'] not in ligand_energy_dict.keys())):
                continue
            complex_name = (row['metal']+'_'+str(row['ox'])+'_'+row['lig1']+'_'+row['lig2']+'_'+
                   row['lig3']+'_'+row['lig4']+'_'+row['lig5']+'_'+row['lig6']+'_s_'+str(row['spin'])+'_'+
                            str(int(row['alpha'])).zfill(2))
            complex_name_no_HFX = (row['metal']+'_'+str(row['ox'])+'_'+row['lig1']+'_'+row['lig2']+'_'+
                   row['lig3']+'_'+row['lig4']+'_'+row['lig5']+'_'+row['lig6']+'_s_'+str(row['spin']))
            if (not np.isnan(float(row['energy']))) and (not np.isnan(float(row['empty_sp_energy']))):
                name_list.append(complex_name)
                name_list_no_HFX.append(complex_name_no_HFX)
                HFX_list.append(int(row['alpha']))
                LDE_list.append((float(row['empty_sp_energy'])+float(ligand_energy_dict[row['lig6']][row['alpha']])-float(row['energy']))*627.509)
            else:
                no_pair_list.append(complex_name)
        LDE_frame = pd.DataFrame()
        LDE_frame['complex'] = name_list
        LDE_frame['complex_no_HFX'] = name_list_no_HFX
        LDE_frame['alpha'] = HFX_list
        LDE_frame['LDE'] = LDE_list
        no_pair_frame = pd.DataFrame()
        no_pair_frame['complex'] = no_pair_list
        if self.path_to_save:
            if not os.path.exists(self.path_to_save+'/MissingPairInfo/'):
                os.mkdir(self.path_to_save+'/MissingPairInfo/')
            LDE_frame.to_csv(self.path_to_save+'/'+str(self.name)+'_LDE.csv',index=False)
            no_pair_frame.to_csv(self.path_to_save+'/MissingPairInfo/'+str(self.name)+'_LDE_missing_pairs.csv',index=False)
            
        try:
            print('This LDE frame contains the following data types:')
            summarizeDataTypes(LDE_frame) #This function prints a data summary and does no other actions
        except:
            print('Summary for LDE frame could not be printed')
        return LDE_frame, no_pair_frame

    def pair_spin_states(self, frame = False, pair_type='LSHS'):
        # This function eliminates any rows that do not have the corresponding pair.
        if not frame:
            pair_frame = self.frame.copy()
            pair_frame = pair_frame.reset_index()
        else:
            pair_frame = frame.reset_index()
        if pair_type == 'LSHS':
            print('pairing LS HS, 4 e difference')
            # Gets SSEs that are 2 electrons different
            metal_spin_dictionary = {'cr': {2: [1, 5]},
                                     'mn': {2: [2, 6], 3: [1, 5]},
                                     'fe': {2: [1, 5], 3: [2, 6]},
                                     'co': {3: [1, 5]},
                                     'mo': {2: [1, 5]},
                                     'tc': {2: [2, 6], 3: [1, 5]},
                                     'ru': {2: [1, 5], 3: [2, 6]},
                                     'rh': {3: [1, 5]}}
            print('Using this pairing: ', metal_spin_dictionary)
        elif pair_type == 'ISHS':
            print('pairing IS HS, 2 e difference')
            metal_spin_dictionary = {'cr': {2: [3, 5]},
                                     'mn': {2: [4, 6], 3: [3, 5]},
                                     'fe': {2: [3, 5], 3: [4, 6]},
                                     'co': {3: [3, 5]},
                                     'mo': {2: [3, 5]},
                                     'tc': {2: [4, 6], 3: [3, 5]},
                                     'ru': {2: [3, 5], 3: [4, 6]},
                                     'rh': {3: [3, 5]}}
            print('Using this pairing: ', metal_spin_dictionary)
        elif pair_type == 'LSIS':
            print('pairing LS IS, 2 e difference')
            metal_spin_dictionary = {'cr': {2: [1, 3], 3: [2, 4]},
                                     'mn': {2: [2, 4], 3: [1, 3]},
                                     'fe': {2: [1, 3], 3: [2, 4]},
                                     'co': {2: [2, 4], 3: [1, 3]},
                                     'mo': {2: [1, 3], 3: [2, 4]},
                                     'tc': {2: [2, 4], 3: [1, 3]},
                                     'ru': {2: [1, 3], 3: [2, 4]},
                                     'rh': {2: [2, 4], 3: [1, 3]}}
            print('Using this pairing: ', metal_spin_dictionary)
        reference_frame = pair_frame.copy() # need this to find criteria that match
        name_list = []
        SSE_list = []
        no_pair_list = []
        name_list_no_HFX = []
        HFX_list = []
        for i, row in pair_frame.iterrows():
            if i % 5000 == 0:
                print('counter',i)
            ox_list = metal_spin_dictionary[row['metal']]
            if row['ox'] in ox_list.keys():
                spin_list = metal_spin_dictionary[row['metal']][row['ox']]
            else:
                continue
            if int(row['spin']) != spin_list[0]: #to avoid repeats...
                continue
            else:
                complex_name = (row['metal']+'_'+str(row['ox'])+'_'+row['lig1']+'_'+row['lig2']+'_'+
                       row['lig3']+'_'+row['lig4']+'_'+row['lig5']+'_'+row['lig6']+'_'+str(int(row['alpha'])).zfill(2))
                complex_name_no_HFX = (row['metal']+'_'+str(row['ox'])+'_'+row['lig1']+'_'+row['lig2']+'_'+
                       row['lig3']+'_'+row['lig4']+'_'+row['lig5']+'_'+row['lig6'])
                subframe = reference_frame[(reference_frame['metal']==row['metal'])&
                                           (reference_frame['ox']==row['ox'])&
                                           (reference_frame['lig1']==row['lig1'])&
                                           (reference_frame['lig2']==row['lig2'])&
                                           (reference_frame['lig3']==row['lig3'])&
                                           (reference_frame['lig4']==row['lig4'])&
                                           (reference_frame['lig5']==row['lig5'])&
                                           (reference_frame['lig6']==row['lig6'])&
                                           (reference_frame['alpha']==row['alpha'])&
                                           (reference_frame['spin']==spin_list[1])
                                          ]
                if subframe.shape[0] == 1:
                    name_list.append(complex_name)
                    name_list_no_HFX.append(complex_name_no_HFX)
                    HFX_list.append(int(row['alpha']))
                    SSE_list.append((float(subframe['energy'].values[0])-float(row['energy']))*627.509)
                else:
                    matching_complex_name = (row['metal']+'_'+str(row['ox'])+'_'+
                                             row['lig1']+'_'+row['lig2']+'_'+row['lig3']+'_'+
                                             row['lig4']+'_'+row['lig5']+'_'+row['lig6']+'_ahf_'+
                                             str(int(row['alpha'])).zfill(2)+'_'+str(spin_list[1]))
                    no_pair_list.append(matching_complex_name)
        SSE_frame = pd.DataFrame()
        SSE_frame['complex'] = name_list
        SSE_frame['complex_no_HFX'] = name_list_no_HFX
        SSE_frame['alpha'] = HFX_list
        SSE_frame['SSE'] = SSE_list
        no_pair_frame = pd.DataFrame()
        no_pair_frame['complex'] = no_pair_list
        if self.path_to_save:
            if not os.path.exists(self.path_to_save+'/MissingPairInfo/'):
                os.mkdir(self.path_to_save+'/MissingPairInfo/')
            SSE_frame.to_csv(self.path_to_save+'/'+str(self.name)+'_SSE_'+str(pair_type)+'.csv',index=False)
            no_pair_frame.to_csv(self.path_to_save+'/MissingPairInfo/'+str(self.name)+'_SSE_'+str(pair_type)+'_missing_pairs.csv',index=False)
        
        try:
            print('This '+pair_type+' SSE frame contains the following data types:')
            summarizeDataTypes(SSE_frame) #This function prints a data summary and does no other actions
        except:
            print('Summary for '+pair_type+' SSE frame could not be printed')
            
        return SSE_frame, no_pair_frame
    
    def filter_df_by(self, **kwargs):
        filtered_frame = self.frame.copy()
        if len(kwargs) > 0:
            for key, value in kwargs.items():
                filtered_frame = filtered_frame[filtered_frame[str(key)]==value]
        elif self.filter_dict:
            for key, value in self.filter_dict.items():
                print(key)
                print(value)
                filtered_frame = filtered_frame[filtered_frame[str(key)]==value]
        else:
            print('no dictionary provided for filter. Please revise')
            raise ValueError
        return filtered_frame
    
    def keep_specific_ligands(self,list_of_ligands):
        #### This function takes a frame and keeps only those with ligands in a certain list (list of ligands)
        filtered_frame = self.frame.copy()
        filtered_frame = filtered_frame[filtered_frame['lig1'].isin(list_of_ligands)]
        filtered_frame = filtered_frame[filtered_frame['lig2'].isin(list_of_ligands)]
        filtered_frame = filtered_frame[filtered_frame['lig3'].isin(list_of_ligands)]
        filtered_frame = filtered_frame[filtered_frame['lig4'].isin(list_of_ligands)]
        filtered_frame = filtered_frame[filtered_frame['lig5'].isin(list_of_ligands)]
        filtered_frame = filtered_frame[filtered_frame['lig6'].isin(list_of_ligands)]
        self.frame = filtered_frame #binds frame with dropped ligands first
        return filtered_frame
            
    def process_data(self, HFX=None, filter_ligands = True, metal_spin_limit = 1, spin_contam_limit = 1, homoleptic_geo_dict = {}, custom_geo_dict = {}):
        if filter_ligands:
            print('Currently filtering specific ligands... Initial shape: '+str(self.frame.shape))
            self.keep_specific_ligands(['ammonia','phosphine','water','uthiol','fluoride','chloride',
                                    'carbonyl','cyanide','misc','acetonitrile']) #if ligands are not these, will be eliminated
        if HFX != None:
            print('HFX '+str(HFX)+' provided. Filtering df...')
            self.frame = self.filter_df_by(alpha=HFX)
        print('Done filtering. Shape after filtering '+str(self.frame.shape))
        self.process_duplicates()
        convergence_drop = self.check_convergence()
        geo_drop = self.geometry_check(**custom_geo_dict)
        geo_drop['reason'] = 'bad geometry'
        geo_drop_homoleptic = self.homoleptic_geo_check(**homoleptic_geo_dict)
        geo_drop_homoleptic['reason'] = 'bad homoleptic geometry'
        spin_contam = self.spin_contamination_check(limit=spin_contam_limit)
        spin_contam['reason'] = 'spin contamination'
        metal_spin = self.spin_on_metal_check(limit=metal_spin_limit)
        metal_spin['reason'] = 'net metal spin'
        if self.process_type == 'LDE':
            spin_contam_empty = self.empty_site_spin_contamination_check()
            spin_contam_empty['reason'] = 'spin contamination on empty site'
            failed = pd.concat([convergence_drop,geo_drop,geo_drop_homoleptic,spin_contam,spin_contam_empty,metal_spin],axis=0,ignore_index=True)
        else:    
            failed = pd.concat([convergence_drop,geo_drop,geo_drop_homoleptic,spin_contam,metal_spin],axis=0,ignore_index=True)
        if self.path_to_save:
            self.frame.to_csv(self.path_to_save+'/'+str(self.name)+'_run_dataframe_after_elimination.csv',index=False)
            failed.to_csv(self.path_to_save+'/'+str(self.name)+'_failures_with_reasons.csv',index=False)
        self.failed = failed # bind the failure reason frames
        return failed ### RETURN PROCESSED FRAME HERE

###############################################################################################
# In the next class, we have methods that can take in one or two dataframes and match things  #
# across rows, down one across one, etc.                                                      #
###############################################################################################

class property_frame_processor:
    # This class takes in one or two dataframes and performs operations on them.
    # In particular, it makes matches between frames, calculates HFX sensitivities,
    # generates RACs. 
    def __init__(self, path_to_save=False, name = 'first_row'):
        self.path_to_save = path_to_save
        #### This name should be descriptive. For example, FirstRowISHS_+str(match_type) would suggest
        #### only first row data (so isoelectronics), and ISHS spin splitting energies. Be descriptive.
        self.name = name
        if path_to_save:
            if os.path.exists(self.path_to_save+'/'+name+'_'+match_type+'.csv'):
                decision = raw_input('this file already exists, press y to continue if overwriting.')
                if decision != 'y':
                    print('User decided not to overwrite present files')
                    raise ValueError
    def make_isovalent_match(self, frame1, frame2, prop='SSE'): #expects to operate with first row first
        metalset = set(frame1['complex'].str.split('_').str[0].tolist())
        if not any([val in ['cr','mn','fe','co'] for val in list(metalset)]):
            frame1, frame2 = frame2, frame1 #### swap so that frame1 has the first row.
        print('isovalent matching')
        match_dictionary = {'cr':'mo','mn':'tc','fe':'ru','co':'rh'}
        first_row_name = []
        second_row_name = []
        first_row_val = []
        second_row_val = []
        for i, row in frame1.iterrows():
            metal = row['complex'].split('_')[0]
            matching_metal = match_dictionary[metal]
            matching_complex = row['complex'].replace(metal,matching_metal)
            subframe = frame2[frame2['complex'] == matching_complex]
            if subframe.shape[0] == 1:
                first_row_name.append(row['complex'])
                second_row_name.append(matching_complex)
                first_row_val.append(float(row[str(prop)]))
                second_row_val.append(float(subframe[str(prop)].values[0]))
            else:
                continue
        isovalent_df = pd.DataFrame()
        isovalent_df['first_row'] = first_row_name
        isovalent_df['second_row'] = second_row_name
        isovalent_df['first_row_'+str(prop)] = first_row_val
        isovalent_df['second_row_'+str(prop)] = second_row_val
        return isovalent_df
      
    def make_down_one_across_one_match(self, frame1, frame2,prop='SSE'):
        metalset = set(frame1['complex'].str.split('_').str[0].tolist())
        if not any([val in ['cr','mn','fe','co'] for val in list(metalset)]):
            frame1, frame2 = frame2, frame1 #### swap so that frame1 has the first row.
        print('down and across matching')
        match_dictionary = {'cr':{2:['tc',3]},'mn':{2:['ru',3]},'fe':{2:['rh',3]}} #cannot make down one across one
        first_row_name = []
        second_row_name = []
        first_row_val = []
        second_row_val = []
        for i, row in frame1.iterrows():
            metal = row['complex'].split('_')[0]
            ox = int(row['complex'].split('_')[1])
            if ox == 3 or metal in ['co','rh']:
                continue
            matching_metal = match_dictionary[metal][ox][0]
            matching_ox = match_dictionary[metal][ox][1]
            matching_complex = row['complex'].replace(metal,matching_metal)
            matching_complex_split = matching_complex.split('_')
            matching_complex_split[1] = str(matching_ox)
            matching_complex = "_".join(matching_complex_split)
            subframe = frame2[frame2['complex'] == matching_complex]
            if subframe.shape[0] == 1:
                first_row_name.append(row['complex'])
                second_row_name.append(matching_complex)
                first_row_val.append(float(row[str(prop)]))
                second_row_val.append(float(subframe[str(prop)].values[0]))
            else:
                continue
        down_one_across_one_df = pd.DataFrame()
        down_one_across_one_df['first_row'] = first_row_name
        down_one_across_one_df['second_row'] = second_row_name
        down_one_across_one_df['first_row_'+str(prop)] = first_row_val
        down_one_across_one_df['second_row_'+str(prop)] = second_row_val
        return down_one_across_one_df
        
    def make_isoelectronic_match(self, frame, prop='SSE'): #assumes isoelectronic pairs are within same DF
        print('isoelectronic matching')
        match_dictionary = {'cr':{2:['mn',3]},'mn':{2:['fe',3]},'fe':{2:['co',3]},
                            'mo':{2:['tc',3]},'tc':{2:['ru',3]},'ru':{2:['rh',3]}} 
        compound_1_name = []
        compound_2_name = []
        compound_1_val = []
        compound_2_val = []
        frame2 = frame.copy()
        for i, row in frame.iterrows():
            metal = row['complex'].split('_')[0]
            ox = int(row['complex'].split('_')[1])
            if ox == 3 or metal in ['co','rh']:
                continue
            matching_metal = match_dictionary[metal][ox][0]
            matching_ox = match_dictionary[metal][ox][1]
            matching_complex = row['complex'].replace(metal,matching_metal)
            matching_complex_split = matching_complex.split('_')
            matching_complex_split[1] = str(matching_ox)
            matching_complex = "_".join(matching_complex_split)
            subframe = frame2[frame2['complex'] == matching_complex]
            if subframe.shape[0] == 1:
                compound_1_name.append(row['complex'])
                compound_2_name.append(matching_complex)
                compound_1_val.append(float(row[str(prop)]))
                compound_2_val.append(float(subframe[str(prop)].values[0]))
            else:
                continue
        isoelectronic_df = pd.DataFrame()
        isoelectronic_df['complex1'] = compound_1_name
        isoelectronic_df['complex2'] = compound_2_name
        isoelectronic_df['complex1_'+str(prop)] = compound_1_val
        isoelectronic_df['complex2_'+str(prop)] = compound_2_val
        return isoelectronic_df
    
    def do_LOOCV(self,sorted_x,sorted_y,off_line_tolerance=5):
        ### This function takes in an x and y list, and then does LOOCV  to see if points lie off line
        loo = LeaveOneOut()
        ytests = []
        ypreds = []
        loo_x = np.array([sorted_x]).T
        loo_y = np.array([sorted_y]).T
        X_val_corresponding_to_error = []
        for train_idx, test_idx in loo.split(loo_x):
            X_train, X_test = loo_x[train_idx], loo_x[test_idx]  # requires arrays
            y_train, y_test = loo_y[train_idx], loo_y[test_idx]
            X_val_corresponding_to_error.append(X_test)
            model = LinearRegression()
            model.fit(X=X_train, y=y_train)
            y_pred = model.predict(X_test)
            # there is only one y-test and y-pred per iteration over the loo.split, 
            # so we append them to respective lists.
            ytests += list(y_test)
            ypreds += list(y_pred)
        ytests = np.array(ytests)
        ypreds = np.array(ypreds)
        error_array = abs(ytests - ypreds)
        zipped_list = zip(X_val_corresponding_to_error,error_array)
        zipped_list.sort(key=lambda t:t[0])
        sorted_error_array = zip(*zipped_list)[1]
        temp_x = []
        temp_y = []
        for j, val in enumerate(sorted_error_array):
            if float(val) < off_line_tolerance:# These are outlying points
                temp_x.append(sorted_x[j])
                temp_y.append(sorted_y[j])
        return temp_x, temp_y
    
    def do_slope_sign_check(self,temp_x,temp_y):
        ### This function takes in and x and y list and checks the slope changes over the portions
        slope_signs = []
        for j2 in range(len(temp_y) - 1):
            slope, intercept, r_value, p_value, std_err = stats.linregress(temp_x[j2:j2 + 2], temp_y[j2:j2 + 2])
            slope_signs.append(np.sign(slope * 1000))
        signchange = ((np.roll(slope_signs, 1) - slope_signs) != 0).astype(int)
        signchange[0] = 0 #Do not want circular behavior
        idx_2_remove = np.where(signchange == 1)[0]
        if len(idx_2_remove) == 1: #one point on either end
            temp_x = list(temp_x)
            temp_y = list(temp_y)
            if idx_2_remove[0]>(len(temp_x)/2.0):
                temp_x = temp_x[idx_2_remove[0]:len(temp_x)]
                temp_y = temp_y[idx_2_remove[0]:len(temp_y)]
            else:
                temp_x = temp_x[0:idx_2_remove[0]]
                temp_y = temp_y[0:idx_2_remove[0]]
            if len(temp_x)<3:
                temp_x = []
                temp_y = []
        else:
            #### more than one discontinuity occurs in the line. Not worth calculating sensitivity.
            temp_x = []
            temp_y = []
        return temp_x, temp_y
        
        
    def calculate_property_HFX_sensitivity(self, frame, prop='SSE', off_line_tolerance = 5):
        R2_cutoff = 0.99
        grouped_frame = frame.groupby('complex_no_HFX')
        HFX_sensitivity_list = []
        complex_name_list = []
        failed_sensitivity_complex_list = []
        for name, group in grouped_frame:
            if group.shape[0]<3:
                #### less than 3 points. Not safe to compute sensitivity. Moving on.
                failed_sensitivity_complex_list.append(name)
            else:
                x = [int(alpha) for alpha in group['alpha'].tolist()]
                y = group[str(prop)].tolist()
                zipped_list = zip(x,y)
                zipped_list.sort(key = lambda t:t[0])
                [sorted_x,sorted_y] = zip(*zipped_list)
                #### sorted alpha and property, now doing LOOCV first 
                temp_x, temp_y = self.do_LOOCV(sorted_x,sorted_y,off_line_tolerance=off_line_tolerance)
                if len(temp_x) < 3:
                    failed_sensitivity_complex_list.append(name)
                    continue
                slope_signs = []
                slope, intercept, r_value, p_value, std_err = stats.linregress(temp_x, temp_y)
                if (r_value ** 2) < R2_cutoff:
                    # if the correlation of the line is lower than expected, do slope change check.
                    temp_x, temp_y = self.do_slope_sign_check(temp_x, temp_y)
                    if len(temp_x)<3:
                        failed_sensitivity_complex_list.append(name)
                    else:
                        slope, intercept, r_value, p_value, std_err = stats.linregress(temp_x,temp_y)
                        complex_name_list.append(name)
                        HFX_sensitivity_list.append(slope)
                else:
                    complex_name_list.append(name)
                    HFX_sensitivity_list.append(slope)
        HFX_sensitivity = pd.DataFrame()
        HFX_sensitivity['complex'] = complex_name_list
        HFX_sensitivity['sensitivity'] = HFX_sensitivity_list
        failed_sensitivity = pd.DataFrame()
        failed_sensitivity['complex'] = failed_sensitivity_complex_list
        if self.path_to_save:
            HFX_sensitivity.to_csv(self.path_to_save+'/'+str(self.name)+'_HFX_sensitivity.csv',index=False)
            failed_sensitivity.to_csv(self.path_to_save+'/'+str(self.name)+'_sensitivity_failure.csv',index=False)
        return HFX_sensitivity, failed_sensitivity

    def make_ligandclass_from_name(self, name):
        ligand_mol, emsg = lig_load(name)
        ligand_mol.convert2mol3D()
        lig_class = ligand(mol3D(), [],ligand_mol.denticity)
        lig_class.mol = ligand_mol
        conatoms = ligand_mol.cat
        return (lig_class,conatoms)
    
    def get_ligands_from_list(self,lignames):
        ### takes in a list of ligands and makes ligand classes.
        ligand_class_list  = []
        ligand_cons_list  = []
        mono_inds = []
        bi_inds = []
        tet_inds  = []
        for ii,name in enumerate(lignames):    
            lig, con = self.make_ligandclass_from_name(name)
            ligand_class_list.append(lig)
            ligand_cons_list.append(con)
        return ligand_class_list, ligand_cons_list
    
    def make_descriptors(self,complex_name, ligand_class_list, ligand_cons_list,lignames,prop_val=False,prop_type = 'SSE'):
        #### Assumes that it is passed in HFX containing name
        split_complex = complex_name.split('_')
        metal = split_complex[0].capitalize()
        ox = int(split_complex[1])
        metal_mol = mol3D()
        metal_mol.addAtom(atom3D(metal)) 
        ox_modifer = {metal: ox}
        if len(split_complex) == 11:
            if prop_type.lower() != 'sensitivity':
                prop_type = 'LDE'
            # This is the LDE case
            spin = int(split_complex[9])
            if prop_type.lower() != 'sensitivity':
                alpha = int(split_complex[10])
            else:
                alpha = False
        else:
            spin = False
            if prop_type.lower() != 'sensitivity':
                alpha = int(split_complex[8])
            else:
                alpha = False
        lig1, lig2, lig3, lig4, lig5, lig6 = (split_complex[2], split_complex[3], split_complex[4],
                                              split_complex[5],split_complex[6],split_complex[7])
        ind1, ind2, ind3, ind4, ind5, ind6 = (lignames.index(lig1), lignames.index(lig2), lignames.index(lig3),
                                              lignames.index(lig4), lignames.index(lig5), lignames.index(lig6))
        custom_ligand_dict = {"eq_ligand_list":[ligand_class_list[ind1], ligand_class_list[ind2],
                                                ligand_class_list[ind3], ligand_class_list[ind4]],
                              "ax_ligand_list":[ligand_class_list[ind5], ligand_class_list[ind6]],
                              "eq_con_int_list":[ligand_cons_list[ind1], ligand_cons_list[ind2],
                                                ligand_cons_list[ind3], ligand_cons_list[ind4]],
                              "ax_con_int_list":[ligand_cons_list[ind5], ligand_cons_list[ind6]]}
        # build the complex mol3D from parts...
        this_complex = assemble_connectivity_from_parts(metal_mol,custom_ligand_dict)
        descriptor_names, descriptors = get_descriptor_vector(this_complex,custom_ligand_dict,ox_modifer,NumB=True,Zeff=True)
        descriptor_names.append('complex')
        descriptors.append(complex_name)
        descriptor_names.append('ox')
        descriptors.append(ox)
        if alpha:
            descriptor_names.append('alpha')
            descriptors.append(alpha)
        if spin:
            descriptor_names.append('spin')
            descriptors.append(spin)
        if prop_val:
            descriptor_names.append(prop_type)
            descriptors.append(prop_val)
        RAC_dictionary = dict(zip(descriptor_names, descriptors))
        return RAC_dictionary
        
    def generate_RACs(self, frame, prop_type=False): ### Takes in a set of complex names and makes the geo free RACs
        lignames = ['acetonitrile','ammonia','carbonyl',
                       'cyanide','chloride','fluoride','misc',
                       'phosphine','uthiol','water']
        ligand_class_list, ligand_cons_list = self.get_ligands_from_list(lignames)
        RAC_list = []
        for i, row in frame.iterrows():
            if prop_type:
                RAC_dict = self.make_descriptors(row['complex'], ligand_class_list, ligand_cons_list,
                                  lignames,prop_val=row[str(prop_type)],prop_type=str(prop_type))
            else:
                # Do not append the property to the RAC frame
                RAC_dict = self.make_descriptors(row['complex'], ligand_class_list, ligand_cons_list, lignames)
            # print(RAC_dict)
            # sard
            RAC_list.append(RAC_dict)
        RACframe = pd.DataFrame(RAC_list)
        if self.path_to_save:
            if prop_type:
                RACframe.to_csv(self.path_to_save+'/'+str(self.name)+'_RACs.csv',index=False)
            else:
                RACframe.to_csv(self.path_to_save+'/'+str(self.name)+'_RACs_and_'+str(prop_type)+'.csv',index=False)
        return RACframe
        
                    
###############################################################################################
# Prints summaries of the type(s) of complexes included in a given dataframe                  #
###############################################################################################

def summarizeDataTypes(df):
    
    def identify_geometry(complex_name):
    #Takes the complex name, returns the perceived geometry of that complex
    
        lst = complex_name.split('_')[2:]
        
        #Clean the list so that it can accept complex names in a few different formats.
        #We want to reduce the list to just the ligand names
        lst.reverse()
        points_to_remove = 0
        for index in range(len(lst)):
            try:
                float(lst[index])
                points_to_remove = index+1
            except:
                if lst[index] == 's':
                    points_to_remove = index+1
        lst.reverse()
        lst = lst[:-1*points_to_remove]    
            
        reduced_list = list(set(lst))
        if len(reduced_list) > 2 or len(reduced_list) == 0:
            print lst
            print reduced_list
            raise Exception('Unexpected geometry with 3 unique ligands')
            
        if len(reduced_list) == 2:
            counted = count_and_combine_list(lst)
            entries = []
            for i in counted.keys():
                entries.append(counted[i])
            entries.sort()
            
            if entries[0]==1 and entries[1]==5:
                return '5+1'
            elif entries[0]==2 and entries[1]==4:
                if lst[0] != lst[2]:
                    return 'cis'
                elif lst[-1] == lst[-2]:
                    return 'trans'
                else:
                    raise exception('Could not identify bidentate geo')
            else:
                raise exception('Could not identify bidentate geo')
                
        if len(reduced_list) == 1:
            return 'homoleptic'
            
    def count_and_combine_list(lst):
    #Takes a list, counts the number of times each value occurs in the list
    #Returns a dictionary of the results
    
        results = {}
        for entry in lst:
            if entry not in results.keys():
                results[entry] = 1
            else:
                results[entry] += 1
        return results
        
        
    if 'complex' not in df.columns:
        print '---complex name not immediately found, attempting to name complex---'
        try:
            df['complex'] = (df['metal']+'_'+str(df['ox'])+'_'+df['lig1']+'_'+df['lig2']+'_'+
                       df['lig3']+'_'+df['lig4']+'_'+df['lig5']+'_'+df['lig6']+'_'+df['alpha'].apply(int).apply(str))
        except:
            raise Exception('Names could not be identified for complexes')
            
    complex_names = df['complex'].tolist()
    
    metals = [i.split('_')[0] for i in complex_names]
    geometries = [identify_geometry(i) for i in complex_names]
    
    print count_and_combine_list(metals)
    print count_and_combine_list(geometries)
            
    
