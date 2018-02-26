% This file explores the topology of structural connectivity in the 96 parcellation (structural matrix)
% thalamic nodes are 41, 42, 43, 89, 90, 91
% Load structural matrix
clear; clc; close all
Data=load('data_struct.mat')
net=Data.data_struct;
%imagesc(net)

%% degre measure 
[id,od,deg] = degrees_dir(net)

figure
info_plot=od
hist(info_plot)
thal.deg=info_plot([41, 42, 43, 89, 90, 91])
min(thal.deg)
max(thal.deg)
hold on

% Thalamic nodes do not seem to be that high degree
%% strength measure
[is,os,str] = strengths_dir(net)

figure
info_plot=str
hist(info_plot)
thal.deg=info_plot([41, 42, 43, 89, 90, 91])
min(thal.deg)
max(thal.deg)

%Thalamic nodes do not seem to be that high degree (is it in middle range)
%% Centrality on weighted, directed matrix
BC=betweenness_wei(net);
figure
info_plot=BC
hist(info_plot)
thal.deg=info_plot([41, 42, 43, 89, 90, 91])
min(thal.deg)
max(thal.deg)
% Two thalamic nuclei have a very high betweenness centrality 347

%% Centrality measure on the binarized matrix
W_nrm = weight_conversion(net, 'binarize')
BC=betweenness_bin(W_nrm);
figure
info_plot=BC
hist(info_plot)
thal.deg=info_plot([41, 42, 43, 89, 90, 91])
min(thal.deg)
max(thal.deg)

% thalamic nuclei show high values for clustering measure
%% Clustering measure
W_nrm = weight_conversion(net, 'normalize')
C = clustering_coef_wd(W_nrm)
figure
info_plot=C
hist(info_plot)
thal.deg=info_plot([41, 42, 43, 89, 90, 91])
min(thal.deg)
max(thal.deg)
% thalamic nuclei show high values for clustering measure


%% Community structure
nNodes=size(net,1);
ci = zeros(nNodes,1);
gamma = 1;
tau = 0.1;
nReps = 10;

for x = 1:500
    [ci_temp(:,x),q_temp(x,1)] = community_louvain(net,gamma,1:1:nNodes); 
end

%estimate a 'consensus' partition (tau and nReps can be altered to change threshold - see https://sites.google.com/site/bctnet/)
D = agreement(ci_temp);
ci = consensus_und(D,tau,nReps);
q = nanmean(q_temp);

%% Participation index (BA)  First do the community structure estimation before this step
BA = zeros(nNodes,1);
BA = participation_coef(net,ci);
BA([41, 42, 43, 89, 90, 91])

figure
info_plot=BA
hist(info_plot)
thal.deg=info_plot([41, 42, 43, 89, 90, 91])
min(thal.deg)
max(thal.deg)

%% within-module degree
Z=module_degree_zscore(net,ci)

figure
info_plot=Z
hist(info_plot)
thal.deg=info_plot([41, 42, 43, 89, 90, 91])
min(thal.deg)
max(thal.deg)

% sort of towards the higher end
% two thalamic nuclei have quite high within module degree

