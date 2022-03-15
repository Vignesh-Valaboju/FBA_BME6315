%% Streptomycine FBA Project 
% BME 6315
clear all; close all;

addpath("../cobratoolbox")
initCobraToolbox;

%% Load data
% load('Sco.mat')
% load('ScoCombinedDraftModel.mat')
% model = ScoCombinedDraftModel;
% load('iMK1208.mat')
% model = iMK1208;

model = readCbModel('Sco.xml');
% writeCbModel(model, 'Sco-GEM.xml')

%% Calculate Essential Genes based on biomass formation
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution]...
    = singleGeneDeletion(model);
essential_genes_ind = find(grRatio < 1e-3);
essential_genes = model.genes(essential_genes_ind);

%% Compare Gene essentiality for table S9 from Wang et al. 
data = readtable('S9_data.xlsx');
for i=1:length(data.Genes)
    if data.BasedOnTn5MutangenesStudy(i)
        exp_true(i,1) = data.Genes(i);
    else 
        exp_true(i,1) = {''};
    end
end
for i=1:length(data.Genes)
    if data.BasedOnTn5MutangenesStudy(i)
        exp_false(i,1) = data.Genes(i);
    else 
        exp_true(i,1) = {''};
    end
end
TP = intersect(exp_true, essential_genes);
FP = setdiff(essential_genes, exp_true);

S9_neg = setdiff(model.genes, exp_true);
mod_neg = setdiff(model.genes, essential_genes);

TN = intersect(mod_neg, S9_neg);
for i=length(exp_true)+1:length(mod_neg)
    exp_true(i,1) = {''};
end
FN = intersect(exp_true, mod_neg);

accuracy = (length(TP)+length(TN))/(length(TP)+length(FP)+length(TN)+length(FN));
specificity = length(TN)/(length(TN)+length(FP));
sensitivity = length(TP)/(length(TP)+length(FN));



%% Alter optimization parameter for study of antibiotic production
% Check current optimization parameter
% printObjective(model);

% Gather information about the reactions involved in actinorhodin
%surfNet(model, 'actinorhodin')
%surfNet(model, 'ACTS19')

%change optimization parameter to actinorhodin production reaction
model = changeObjective(model, 'ACTS19');

%% Calculate Essential Genes based on ACTS19 rxn optimization 
[gene_grRatio, gene_grRateKO, gene_grRateWT, gene_hasEffect, gene_delRxns, gene_fluxSolution]...
    = singleGeneDeletion(model);
essential_genes_ind = find(gene_grRatio < 1e-3);
essential_genes = model.genes(essential_genes_ind);

%% Calculate Essential Reactions based on ACTS19 rx optimization
[rxn_grRatio, rxn_grRateKO, rxn_grRateWT, rxn_hasEffect, rxn_delRxn, rxn_fluxSolution]...
    = singleRxnDeletion(model);
essential_rxn_ind = find(rxn_grRatio < 1e-3);
essential_rxns = model.rxnNames(essential_rxn_ind);
essential_rxns_identifies = model.rxns(essential_rxn_ind);






