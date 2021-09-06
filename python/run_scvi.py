## Init
import os
import numpy as np
import numpy.random as random
import pandas as pd

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.inference.autotune import auto_tune_scvi_model

from umap import UMAP

import torch
import scanpy as sc
import louvain

import logging
import pickle
from hyperopt import hp


# %matplotlib inline

use_cuda = True
n_epochs_all = None
save_path = ''
show_plot = True
os.chdir("/scratch/cs/csb/projects/lag3_scrnaseq/")


## Download data
f2059_0 = CsvDataset(filename='results/scvi/input_files/2059_0.csv', save_path='', sep=',', new_n_genes=False)
f2059_1 = CsvDataset(filename='results/scvi/input_files/2059_1.csv', save_path='', sep=',', new_n_genes=False)
f2059_3 = CsvDataset(filename='results/scvi/input_files/2059_3.csv', save_path='', sep=',', new_n_genes=False)

f2083_0 = CsvDataset(filename='results/scvi/input_files/2083_0.csv', save_path='', sep=',', new_n_genes=False)
f2083_1 = CsvDataset(filename='results/scvi/input_files/2083_1.csv', save_path='', sep=',', new_n_genes=False)
f2083_3 = CsvDataset(filename='results/scvi/input_files/2083_3.csv', save_path='', sep=',', new_n_genes=False)

f2342_0 = CsvDataset(filename='results/scvi/input_files/2342_0.csv', save_path='', sep=',', new_n_genes=False)
f2342_1 = CsvDataset(filename='results/scvi/input_files/2342_1.csv', save_path='', sep=',', new_n_genes=False)
f2342_3 = CsvDataset(filename='results/scvi/input_files/2342_3.csv', save_path='', sep=',', new_n_genes=False)

f2348_0 = CsvDataset(filename='results/scvi/input_files/2348_0.csv', save_path='', sep=',', new_n_genes=False)
f2348_1 = CsvDataset(filename='results/scvi/input_files/2348_1.csv', save_path='', sep=',', new_n_genes=False)
f2348_3 = CsvDataset(filename='results/scvi/input_files/2348_3.csv', save_path='', sep=',', new_n_genes=False)

f2350_0 = CsvDataset(filename='results/scvi/input_files/2350_0.csv', save_path='', sep=',', new_n_genes=False)
f2350_1 = CsvDataset(filename='results/scvi/input_files/2350_1.csv', save_path='', sep=',', new_n_genes=False)
f2350_3 = CsvDataset(filename='results/scvi/input_files/2350_3.csv', save_path='', sep=',', new_n_genes=False)

f2351_0 = CsvDataset(filename='results/scvi/input_files/2351_0.csv', save_path='', sep=',', new_n_genes=False)
f2351_1 = CsvDataset(filename='results/scvi/input_files/2351_1.csv', save_path='', sep=',', new_n_genes=False)
f2351_3 = CsvDataset(filename='results/scvi/input_files/2351_3.csv', save_path='', sep=',', new_n_genes=False)

## Combine
all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs = [f2059_0.X, f2059_1.X, f2059_3.X,
                                               f2083_0.X, f2083_1.X, f2083_3.X,
                                               f2342_0.X, f2342_1.X, f2342_3.X,

                                               f2348_0.X, f2348_1.X, f2348_3.X,
                                               f2350_0.X, f2350_1.X, f2350_3.X,
                                               f2351_0.X, f2351_1.X, f2351_3.X])

## Train, save and fin
# vae      = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels, n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')
# trainer  = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
# trainer.train(n_epochs=50)
# torch.save(trainer.model.state_dict(), 'results/scvi/output/lag3_all_genes.pkl')
# trainer.model.load_state_dict(torch.load('results/scvi/output/aml_all_genes.pkl'))



n_epochs = 50
max_evals = if_not_test_else(100, 1)
reserve_timeout = if_not_test_else(180, 5)
fmin_timeout = if_not_test_else(300, 10)

trainer, trials = auto_tune_scvi_model(
    gene_dataset=all_dataset,
    parallel=True,
    exp_key="all_dataset",
    train_func_specific_kwargs={"n_epochs": n_epochs},
    max_evals=max_evals,
    reserve_timeout=reserve_timeout,
    fmin_timeout=fmin_timeout,
)

torch.save(trainer.model.state_dict(), 'results/scvi/output/lag3_optim.pkl')



## Sample posterior to get latent representation and save those embeddings
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))

latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

np.savetxt("results/scvi/output/lag3_optim_batch_latent.csv", latent, delimiter=",")
np.savetxt("results/scvi/output/lag3_optim_batch_indices.csv", batch_indices, delimiter=",")
