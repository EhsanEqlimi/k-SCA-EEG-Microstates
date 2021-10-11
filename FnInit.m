function FnInit
Iter_Num=100; % The number of iteration
m_Set=[64]; % The number of sensors a.k.a the ambient dimension
n_Set=[4];% The number of sources (ith row correspond ith elemnet in m_set
% m_Set=[3]; % The number of sensors a.k.a the ambient dimension
% n_Set=[7];% The number of sources (ith row correspond ith elemnet in m_set

% SNR_Add=[30]; % SNR values of additive noise
Sigma_Off=[1e-4]; % Parameter controls the standard deviation of normal noise over the active lements

T= 2000;%c*30; %1000; % The Number of data poins a.k.a the number of samples (time points)
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A is satisfied
Iter_Num_A=1000;% The number of iteration to generate a good mixing matrix
Rank_Th_A=0.1; % to generate a good mixing matrix
Th_RansacSubspace=1e-5;% RANSAC Threshold in Distance Function for subspace estimation
ThBBC=1e-3; %BBC clustering for mixing idetification
Th_RansacMixing=1e-5;% RANSAC Threshold in Distance Function for mixing matrix identifiation

DualMode=0; %Twice subpspace clusteing
BAS_Set=[];NMSE_Set=[];Norm2Error_Set=[];Sub_Error_Set=[];Error_Set=[];
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='PermkSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
Orth=0; %if Orth=1 A is orth
p=0.95;
% Sigma=0;
CT=0;