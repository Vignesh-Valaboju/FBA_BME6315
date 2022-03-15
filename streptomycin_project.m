%% Streptomycine FBA Project 
% BME 6315
clear all; close all;

addpath("../cobratoolbox")
initCobraToolbox;

clc
%% Load data
model = readCbModel('Sco.mat');

%% Alter optimization parameter for study of antibiotic production
%check current optimization parameter
printObjective(model);
%gather information about the reactions involved in actinorhodin
%surfNet(model, 'actinorhodin')
%surfNet(model, 'ACTS19')
%change optimization parameter to actinorhodin production reaction
model = changeObjective(model, 'ACTS19');

%% Calculate Essential Genes
[gene_grRatio, gene_grRateKO, gene_grRateWT, gene_hasEffect, gene_delRxns, gene_fluxSolution]...
    = singleGeneDeletion(model);
essential_genes_ind = find(gene_grRatio <= 0.01);
essential_genes = model.genes(essential_genes_ind);

%% Calculate Essential Reactions
[rxn_grRatio, rxn_grRateKO, rxn_grRateWT, rxn_hasEffect, rxn_delRxn, rxn_fluxSolution]...
    = singleRxnDeletion(model);
essential_rxn_ind = find(rxn_grRatio <= 0.01);
essential_rxns = model.rxnNames(essential_rxn_ind);