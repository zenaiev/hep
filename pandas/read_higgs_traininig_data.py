# Based on https://www.kaggle.com/code/sugataghosh/higgs-boson-event-detection-part-1-eda
# Data to be retrieved from https://www.kaggle.com/competitions/higgs-boson/data?select=training.zip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math

data_train = pd.read_csv('higgs-boson/training/training.csv')
print(pd.Series({"Memory usage": "{:.2f} MB".format(data_train.memory_usage().sum()/(1024*1024)),
                 "Dataset shape": "{}".format(data_train.shape)}).to_string())
print(data_train.head())

# Splitting the training by target class
data_train_b = data_train[data_train['Label'] == 'b'] # Background events in the training set
data_train_s = data_train[data_train['Label'] == 's'] # System events in the training set

def hist_target(df, cols, target, bins, ncols = 3):
    nrows = math.ceil(len(cols) / ncols)
    fig, ax = plt.subplots(nrows, ncols, figsize = (5.5 * ncols, 2.5 * nrows), sharey = False)
    for i in range(len(cols)):
        sns.histplot(data = df, x = cols[i], bins = bins, hue = target, palette = ['red', 'grey'], ax = ax[i // ncols, i % ncols])
        ax[i // ncols, i % ncols].set_xlabel(cols[i])
        if i % ncols != 0:
            ax[i // ncols, i % ncols].set_ylabel(" ")
    plt.tight_layout()
    plt.show()

l1 = [
    'DER_mass_MMC', 'DER_mass_transverse_met_lep', 'DER_mass_vis',
    'DER_pt_h', 'DER_deltaeta_jet_jet', 'DER_mass_jet_jet',
]

l2 = [
    'DER_prodeta_jet_jet', 'DER_deltar_tau_lep', 'DER_pt_tot',
     'DER_sum_pt', 'DER_pt_ratio_lep_tau', 'DER_met_phi_centrality'
]

l3 = [
    'DER_lep_eta_centrality', 'PRI_tau_pt', 'PRI_tau_eta', 'PRI_tau_phi',
    'PRI_lep_pt', 'PRI_lep_eta', 'PRI_lep_phi', 'PRI_met', 'PRI_met_phi',
    'PRI_met_sumet', 'PRI_jet_leading_pt', 'PRI_jet_leading_eta'
]

l4 = [
    'PRI_jet_leading_phi', 'PRI_jet_subleading_pt',
    'PRI_jet_subleading_eta', 'PRI_jet_subleading_phi', 'PRI_jet_all_pt'
]

# Distributions of the float features in the training set by target class
hist_target(data_train.replace(-999, np.nan),
     l1,
     target = 'Label',
     bins = max(math.floor(len(data_train_b)**(1/3)), math.floor(len(data_train_s)**(1/3))),
     ncols = 3)

hist_target(data_train.replace(-999, np.nan),
     l2,
     target = 'Label',
     bins = max(math.floor(len(data_train_b)**(1/3)), math.floor(len(data_train_s)**(1/3))),
     ncols = 3)

hist_target(data_train.replace(-999, np.nan),
     l3,
     target = 'Label',
     bins = max(math.floor(len(data_train_b)**(1/3)), math.floor(len(data_train_s)**(1/3))),
     ncols = 3)

hist_target(data_train.replace(-999, np.nan),
     l4,
     target = 'Label',
     bins = max(math.floor(len(data_train_b)**(1/3)), math.floor(len(data_train_s)**(1/3))),
     ncols = 3)