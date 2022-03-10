%% Streptomycine FBA Project 
% BME 6315
clear all; close all;

addpath("../cobratoolbox")
initCobraToolbox;

%% Load data
load('Sco.mat')

%% Calculate Essential Genes
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution]...
    = singleGeneDeletion(model);
essential_genes_ind = find(grRateKO == 0);
essential_genes = model.genes(essential_genes_ind);