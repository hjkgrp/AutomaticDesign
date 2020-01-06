#NOTE: ORBKIT MUST BE INSTALLED FOR THIS SCRIPT TO BE USED
from orbkit import read, atomic_populations, options
from orbkit import analytical_integrals as ai
import matplotlib.pyplot as plt
import numpy as np 
options.quiet = True
class MoldenInfo:
    def __init__(self, moldenpath):
        self.moldenpath = moldenpath
        # print('Now Analyzing: ', moldenpath)
        self.qc = read.main_read(self.moldenpath, itype='molden')

    def find_homos(self):
        store_a = 0
        for idx, value in enumerate(self.qc.get_mo_labels()):
            number, letter = value.split()[1].split('.')
            if 'a' in letter:
                old = store_a
                new = int(number)
            else:
                old = 0
            if abs(new-old)>1:
                alphaidx = idx-1
                betaidx = len(self.qc.get_mo_labels())-1
                break
            store_a = new
        print(('alphaHOMO',self.qc.mo_spec[alphaidx]['energy']))
        print(('betaHOMO',self.qc.mo_spec[betaidx]['energy']))
        self.alphaidx = alphaidx
        self.betaidx = betaidx

    def get_MO_energies(self):
        self.list_of_alpha_energies = {i:self.qc.mo_spec[i]['energy'] for i, val in enumerate(self.qc.mo_spec) if self.qc.mo_spec[i]['spin']=='alpha'}
        self.list_of_beta_energies = {i:self.qc.mo_spec[i]['energy'] for i, val in enumerate(self.qc.mo_spec) if self.qc.mo_spec[i]['spin']=='beta'}
        self.list_of_alpha_energies.update((x, y*27.2114) for x, y in list(self.list_of_alpha_energies.items()))
        self.list_of_beta_energies.update((x, y*27.2114) for x, y in list(self.list_of_beta_energies.items()))

    def plot_MO_energies(self):
        plt.scatter(list(self.list_of_alpha_energies.keys()), list(self.list_of_alpha_energies.values()))
        plt.xlabel('MO Number')
        plt.ylabel('Total MO Energy [eV]')
        plt.ylim((-20,0)) #There are many more energy levels, just cutting plot off at -20 eV
        plt.show()

    def plot_dchar_vs_energies(self):
        energy_vals = [self.d_character_dict[i]['energy'] for i in list(self.d_character_dict.keys())]
        dchar_vals = [self.d_character_dict[i]['dchar'] for i in list(self.d_character_dict.keys())]
        maxval = max(energy_vals)
        limited_energy = [i for i in energy_vals if i >= maxval-5]
        limited_dchar = [x for x, y in zip(dchar_vals, energy_vals) if y >= maxval-5]
        print((np.mean(limited_dchar)))
        plt.scatter(limited_energy, limited_dchar)
        plt.xlabel('Total MO Energy [eV]')
        plt.ylabel('Percent d-character')
        plt.show()

    def get_overlap_matrix(self):
        self.ao_overlap_matrix = ai.get_ao_overlap(self.qc.geo_spec,self.qc.geo_spec,self.qc.ao_spec,ao_spherical=self.qc.ao_spherical)
        return self.ao_overlap_matrix

    def get_FMO_HOMO_coeffs(self):
        self.alphaHOMOcoeffs = self.qc.mo_spec[self.alphaidx]['coeffs']
        self.betaHOMOcoeffs = self.qc.mo_spec[self.betaidx]['coeffs']
        self.alphaHOMOenergy = self.qc.mo_spec[self.alphaidx]['energy']
        self.betaHOMOenergy = self.qc.mo_spec[self.betaidx]['energy']

    def get_d_orbital_coeffs(self):
        ################################################################################
        # This function assumes you are using an ECP with the metal listed first.      #
        # 3 lines for s, 3 for p, and 3 for d. Because the metal always comes first,   #
        # we can get the d-character of the molden file quite easily by taking the     #
        # coefficients that correspond to the d-orbitals, which will be the 12         #
        # coefficients that come after the first 12 (the first 3 correspond to 1s, 2s, #
        # and 3s, then the next 9 correspond to the 2p, 3p, and 4p (x, y, z) for each).#
        # The last 12 are the d-orbitals.                                              #
        ################################################################################
        dorb_alpha = self.alphaHOMOcoeffs.copy()
        dorb_beta = self.betaHOMOcoeffs.copy()
        idxlist = np.linspace(12,23,12) #This is basis set dependent. This is for LACVP* basis set.
        dorb_alpha = np.array([x if ind in idxlist else 0 for ind, x in enumerate(dorb_alpha)])
        dorb_beta = np.array([x if ind in idxlist else 0 for ind, x in enumerate(dorb_beta)])
        self.d_orbital_coeffs_alpha = dorb_alpha
        self.d_orbital_coeffs_beta = dorb_beta

    def get_all_orbital_d_character(self):
        self.get_overlap_matrix()
        idxlist = np.linspace(12,23,12)
        char_dict = dict()
        for i, val in enumerate(self.qc.mo_spec):
            full_coeffs = self.qc.mo_spec[i]['coeffs']
            d_orb = full_coeffs.copy()
            d_orb = np.array([x if ind in idxlist else 0 for ind, x in enumerate(d_orb)])
            MO_char = np.dot(d_orb, np.dot(self.ao_overlap_matrix, d_orb.T))
            if i in list(char_dict.keys()):
                print('Already in dictionary...')
            else:
                char_dict[i] = {'dchar': MO_char, 'energy': self.qc.mo_spec[i]['energy']*27.2114}
        self.d_character_dict = char_dict

    def check_normalization(self):
        self.find_homos()
        self.get_overlap_matrix()
        self.get_FMO_HOMO_coeffs()
        alpha_norm = np.dot(self.alphaHOMOcoeffs, np.dot(self.ao_overlap_matrix,self.alphaHOMOcoeffs.T))
        beta_norm = np.dot(self.betaHOMOcoeffs, np.dot(self.ao_overlap_matrix,self.betaHOMOcoeffs.T))
        # print('Norm of Alpha: ',alpha_norm)
        # print('Norm of Beta: ',beta_norm)

    def get_MO_character(self):
        self.find_homos()
        self.get_overlap_matrix()
        self.get_FMO_HOMO_coeffs()
        self.get_d_orbital_coeffs()
        MO_char_alpha = np.dot(self.d_orbital_coeffs_alpha, np.dot(self.ao_overlap_matrix,self.d_orbital_coeffs_alpha.T))
        MO_char_beta = np.dot(self.d_orbital_coeffs_beta, np.dot(self.ao_overlap_matrix,self.d_orbital_coeffs_beta.T))
        print(('Alpha d-character: ',MO_char_alpha))
        print(('Beta d-character: ', MO_char_beta))
        return MO_char_alpha, MO_char_beta

    def get_mullpop_and_charge(self):
        ################################################################################
        # This function assumes that the alpha spin is the majority spin and beta is   #
        # minority. Takes a while to calculate each separately.                        #
        ################################################################################
        self.qc_alpha = read.main_read(self.moldenpath, itype='molden',spin='alpha')
        self.qc_beta = read.main_read(self.moldenpath, itype='molden',spin='beta')
        self.pop = atomic_populations.mulliken(self.qc)
        self.pop_alpha = atomic_populations.mulliken(self.qc_alpha)
        self.pop_beta = atomic_populations.mulliken(self.qc_beta)
        self.spinpop = []
        self.alphapop = []
        self.betapop = []
        self.charge = []
        metallist = ['mn','cr','co','ni','fe']
        for idx, val in enumerate(self.qc.geo_info):
            if val[0].lower() in metallist:
                self.metalpopalpha = self.pop_alpha['population'][idx]
                self.metalpopbeta = self.pop_beta['population'][idx]
                self.metalcharge = self.pop['charge'][idx]
            self.spinpop.append((val[0].capitalize(), self.pop_alpha['population'][idx]-self.pop_beta['population'][idx]))
            self.alphapop.append((val[0].capitalize(), self.pop_alpha['population'][idx]))
            self.betapop.append((val[0].capitalize(), self.pop_beta['population'][idx]))
            self.charge.append((val[0].capitalize(), self.pop['charge'][idx]))
        self.oxygenpopalpha = self.pop_alpha['population'][-1]
        self.oxygenpopbeta = self.pop_beta['population'][-1]
        self.oxygencharge = self.pop['charge'][-1]
        self.totalspin = sum(self.pop_alpha['population']-self.pop_beta['population'])
        self.totalcharge = sum(self.pop['charge'])
        print(('Spin: ',self.totalspin))
        print(('Charge: ',self.totalcharge))
        print((self.pop['charge']))
        print(('Metal: ',self.metalcharge, 'Oxygen: ',self.oxygencharge))
        return self.totalspin, self.totalcharge, self.metalcharge, self.metalpopalpha, self.metalpopbeta, self.oxygencharge, self.oxygenpopalpha, self.oxygenpopbeta
            

#EXAMPLE
# moldenpath = '/Users/adityanandy/Desktop/catalysis.molden'
# my_df_final = MoldenInfo(moldenpath)
# my_df_final.get_MO_energies()
# my_df_final.get_all_orbital_d_character()
# my_df_final.plot_dchar_vs_energies()