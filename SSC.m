%--------------------------------------------------------------------------
% This is the main function for running SSC. 
% Load the DxN matrix X representing N data points in the D dim. space 
% living in a union of n low-dim. subspaces.
% The projection step onto the r-dimensional space is arbitrary and can 
% be skipped. In the case of using projection there are different types of 
% projections possible: 'NormalProj', 'BernoulliProj', 'PCA'. Please refer 
% to DataProjection.m for more information.
%--------------------------------------------------------------------------
% X: DxN matrix of N points in D-dim. space living in n low-dim. subspaces
% s: groundtruth for the segmentation
% n: number of subspaces
% r: dimension of the projection e.g. r = d*n (d: max subspace dim.)
% Cst: 1 if using the constraint sum(c)=1 in Lasso, else 0
% OptM: optimization method {'L1Perfect','L1Noise','Lasso','L1ED'}, see 
% SparseCoefRecovery.m for more information
% lambda: regularization parameter for 'Lasso' typically in [0.001,0.01] 
% or the noise level for 'L1Noise'. See SparseCoefRecovery.m for more 
% information.
% K: number of largest coefficients to pick in order to build the
% similarity graph, typically K = max{subspace dimensions} 
% Missrate: vector of misclassification rates
%--------------------------------------------------------------------------
% In order to run the code CVX package must be installed in Matlab. It can 
% be downlaoded from http://cvxr.com/cvx/download
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2010
%--------------------------------------------------------------------------

clc, clear all, close all
D = 64; %Dimension of ambient space
n = 4; %Number of subspaces
% d1 = 2; d2 = 1; %d1 and d2: dimension of subspace 1 and 2
d1=2; d2=2;d3=2;d4=3;d5=2;d6=3;d7=2;
N1 = 20; N2 = 20; %N1 and N2: number of points in subspace 1 and 2
N3=20;N4=20;N5=20;N6=20;N7=20;
X1 = randn(D,d1) * randn(d1,N1); %Generating N1 points in a d1 dim. subspace
X2 = randn(D,d2) * randn(d2,N2); %Generating N2 points in a d2 dim. subspace
X3 = randn(D,d3) * randn(d3,N3); %Generating N2 points in a d2 dim. subspace

X4 = randn(D,d4) * randn(d4,N4); %Generating N2 points in a d2 dim. subspace
X5 = randn(D,d5) * randn(d5,N5); %Generating N2 points in a d2 dim. subspace
X6 = randn(D,d6) * randn(d6,N6); %Generating N2 points in a d2 dim. subspace
X7 = randn(D,d7) * randn(d7,N7); %Generating N2 points in a d2 dim. subspace

% X = [X1 X2 X3 X4 X5 X6 X7];
% s = [1*ones(1,N1) 2*ones(1,N2) 3*ones(1,N3)  4*ones(1,N4) 5*ones(1,N5) 6*ones(1,N6) 7*ones(1,N7)]; %Generating the ground-truth for evaluating clustering results

X = [X1 X2 X3 ];
[coeff, score, latent, tsquared, explained, mu] =pca(X');
scatter3(score(1,:),score(2,:),score(3,:))
s = [1*ones(1,N1) 2*ones(1,N2) 3*ones(1,N3)  ]; %Generating the ground-truth for evaluating clustering results
%%
clear S
Fs=512;
L=128;
S=zeros(4,Fs);
Cosine= tukeywin(L,0)';
S(1,:)=[Cosine zeros(1,3*L)] ;%figure,plot(S);
S(2,:)=[zeros(1,60) Cosine zeros(1,Fs-(L+60))];
S(3,:)=[zeros(1,60*2) Cosine zeros(1,Fs-(L+2*60)) ];
S(4,:)=[zeros(1,60*3) Cosine zeros(1,Fs-(L+3*60))];
figure,
subplot(4,1,1)
plot(1:Fs,S(1,:))
subplot(4,1,2)

plot(1:Fs,S(2,:))
subplot(4,1,3)

plot(1:Fs,S(3,:))
subplot(4,1,4)

plot(1:Fs,S(4,:))
S=S(:,1:309);
figure,
subplot(4,1,1)
plot(1:309,S(1,:))
subplot(4,1,2)

plot(1:309,S(2,:))
subplot(4,1,3)

plot(1:309,S(3,:))
subplot(4,1,4)

plot(1:309,S(4,:))
X=randn(64, size(S,1))*S;
[coeff,score,latent,tsquared] =pca(X);
figure,scatter3(coeff(:,1),coeff(:,2),coeff(:,3));
%%
Iter_Num=100; % The number of iteration
m=[64]; % The number of sensors a.k.a the ambient dimension
n=[4];% The number of sources (ith row correspond ith elemnet in m_set
k=1;
% m_Set=[3]; % The number of sensors a.k.a the ambient dimension
% n_Set=[7];% The number of sources (ith row correspond ith elemnet in m_set

% SNR_Add=[30]; % SNR values of additive noise
Sigma_Off=[0]; % Parameter controls the standard deviation of normal noise over the active lements

T= 200;%c*30; %1000; % The Number of data poins a.k.a the number of samples (time points)
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A is satisfied
Iter_Num_A=1000;% The number of iteration to generate a good mixing matrix
Rank_Th_A=0.1; % to generate a good mixing matrix
Th_RansacSubspace=1e-5;% RANSAC Threshold in Distance Function for subspace estimation
ThBBC=1e-3; %BBC clustering for mixing idetification
Th_RansacMixing=1e-5;% RANSAC Threshold in Distance Function for mixing matrix identifiation

DualMode=0; %Twice subpspace clusteing
BAS_Set=[];NMSE_Set=[];Norm2Error_Set=[];Sub_Error_Set=[];Error_Set=[];
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='kSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
Orth=0; %if Orth=1 A is orth
p=0.95;
% Sigma=0;
CT=0;
c=nchoosek(n,k); % The number of all possible r-dim subspaces
                N=zeros(1,c);
                N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mode
                for j=1:k
                    for i=1:nchoosek(n,j)
                        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% N showes the number of each subspace in MSCA mode
                    end
                end
                e= 1-(1/c);%(c-1)/c;  % Probability that a point is an outlier
                %                     MaxIterEst=log(1-p)/log(1-(1-e)^k); %1e6;%log(1-p)/log(1-(1-e)^k);
                MaxIterEst=1e6;
                
                if MaxIterEst<1000% Sometimes if k=1, it happnes
                    MaxIterEst=1e4;
                end
                        Sigma=Sigma_Off(1);

[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Sigma,Iter_Num_A,Rank_Th_A,AMode,MixingMode);
T=size(S,2);
figure,
subplot(4,1,1)
plot(1:T,S(1,:))
subplot(4,1,2)

plot(1:T,S(2,:))
subplot(4,1,3)

plot(1:T,S(3,:))
subplot(4,1,4)

plot(1:T,S(4,:))

[coeff,score,latent,tsquared] =pca(X');
% figure,scatter3(X(1,:), X(2,:), X(3,:))

figure,scatter3(score(:,1), score(:,2), score(:,3))
%%
r = 0; %Enter the projection dimension e.g. r = d*n, enter r = 0 to not project
Cst = 0; %Enter 1 to use the additional affine constraint sum(c) == 1
OptM = 'Lasso'; %OptM can be {'L1Perfect','L1Noise','Lasso','L1ED'}
lambda = 0.001; %Regularization parameter in 'Lasso' or the noise level for 'L1Noise'
K = max(d1,d2); %Number of top coefficients to build the similarity graph, enter K=0 for using the whole coefficients
if Cst == 1
    K = max(d1,d2) + 1; %For affine subspaces, the number of coefficients to pick is dimension + 1 
end
r=c*k;
Xp = DataProjection(X,r,'NormalProj');
[A,B]=kmeans(Xp',2);
CMat = SparseCoefRecovery(Xp,Cst,OptM,lambda);
[CMatC,sc,OutlierIndx,Fail] = OutlierDetection(CMat,s);
if (Fail == 0)
    CKSym = BuildAdjacency(CMatC(1:T,1:T),k+1);
    [Grps , SingVals, LapKernel] = FnSpectralClustering(CKSym,c);
    figure,stem(Grps(:,1))
    Missrate = Misclassification(Grps,Labels)
   [ClusterinError1,ClusterinError2]=FnSubspaceClusteringErrorFinder(Grps(:,1),Labels)

    [idx,Clus]=kmeans(X',n);
        figure,stem(idx)
           [ClusterinError_km,ClusterinError_km]=FnSubspaceClusteringErrorFinder(idx,Labels)

    Missrate2 = Misclassification(idx',Labels)

%     save Lasso_001.mat CMat CKSym Missrate SingVals LapKernel Fail
else
    save Lasso_001.mat CMat Fail
end