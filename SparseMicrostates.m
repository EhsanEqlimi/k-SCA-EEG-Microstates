%--------------------------------------------------------------------------
% This is the main function for running EEG (simulated) microstate analysis
%based on sparse componnete analysis (SCA) 

% Load the mxT matrix X representing T data points in the m dim. space 
% living in a union of c low-dim. subspaces.
% The projection step onto the r-dimensional (r=c*k) space is arbitrary and  
% can be skipped. In the case of using projection there are different types of 
% projections possible: 'NormalProj', 'BernoulliProj', 'PCA'. Please refer 
% to DataProjection.m for more information.
%--------------------------------------------------------------------------
% Witten by Ehsan Eqlimi, @WAVES-UGent, Ghent, Belgium
%%
clc;clear;close all;
%% Initialization
m=64; %The dimension of ambient space/ the number of sensors/channels
n=4; %The number of sources
k=2; %Sparsity level or the dimension of subspaces
b=1; %If you input b ? 0, you obtain a rectwin window. If you input b ? 1,
% you obtain a hann window 
c=nchoosek(n,k); % the number of subspace if k is fixed for every single 
% time point
Fs=512; %Sampling Frequency
A=rand(m,n); %Normaly distrubuted Mixing vectors
L=Fs/4; %The length of each source/microstate (the number of samples) 
Cosine= tukeywin(L,b)';% The waveform/function for sources
GroundTruthLables=[];
OL=.5;% Overlap
RealL=Fs/2;
Sigma=0;%1e-4; %additive noise level
%% Generating Source Matrix
S=zeros(n,RealL);
for i=1:n
    Temp=[zeros(1,(i-1)*L*OL) Cosine zeros(1,4*L-L-(i-1)*L*OL)] ;%figure,plot(S);
    S(i,:)=Temp(:,1:RealL);
    subplot(4,1,i)
    plot(1:RealL,S(i,:))
end
currentFigure = gcf;
title(currentFigure.Children(end), 'Simulated Source (dipole) activation over time');
xlabel('Time [Sample]');
% xlim([0 RealL]);
%% Groundthruth labels for each time point
Actualc=4; 

for j=1:Actualc
GroundTruthLables=[GroundTruthLables j*ones(1,L/2)];
end
figure, plot(GroundTruthLables)
title('Groundtruth lables over time');
%% Mixing (X=AS+E)
X=A*S+Sigma*randn(m,RealL);
%% PCA projection
[LabelsKmenas,KmeansClusters]=kmeans(X',Actualc);
[~,StartIndClusters,~]=unique(LabelsKmenas, 'stable');
StartIndClusters(5)=RealL;
r=k*Actualc; %dimension of the space to project the data to
%%Enter the projection dimension e.g. r = d*n, enter r = 0 to not project
Colors={'.r', '.b', '.g', '.m','.y'}; % Colors for clusters
% PCA projection
[Coeff,Score,Latent,TSquared] =pca(zscore(X'));
% Data Projection based on PCA
Xp = FnDataProjection(X,r,'PCA');
figure,
for j=1:4   
scatter3(Xp(1,StartIndClusters(j):StartIndClusters(j+1)),Xp(2,StartIndClusters(j):StartIndClusters(j+1)),...
    Xp(3,StartIndClusters(j):StartIndClusters(j+1)),Colors{j});
hold on
end
title('PCA projection 1')

figure,
for j=1:4
scatter3(Score(StartIndClusters(j):StartIndClusters(j+1)-1,1),Score(StartIndClusters(j):StartIndClusters(j+1)-1,2),...
    Score(StartIndClusters(j):StartIndClusters(j+1)-1,3),Colors{j});
hold on
end
title('PCA projection 2')
%% Sparse Component Analysis
Cst=0; %Enter 1 to use the additional affine constraint sum(c) == 1
OptM='Lasso'; %OptM can be {'L1Perfect','L1Noisy','Lasso','L1ED'}
lambda=0.000001; %Sigma; %Regularization parameter in 'Lasso' or the noise level for 'L1Noise'
K = k; %Number of top coefficients to build the similarity graph, enter K=0 for using the whole coefficients
if Cst == 1
    K = k + 1; %For affine subspaces, the number of coefficients to pick is dimension + 1
end
% Xp = FnDataProjection(X,r,'NormalProj');
% CMat = SparseCoefRecovery(Xp,Cst,OptM,lambda);
% [CMatC,sc,OutlierIndx,Fail] = OutlierDetection(CMat,GroundTruthLables);
% if (Fail == 0)
%     CKSym = BuildAdjacency(CMatC,K);
%     [Grps , SingVals, LapKernel] = FnSpectralClustering(CKSym,Actualc);
%  
%     Missrate_Sparse = Misclassification(Grps,GroundTruthLables)
%     [ClusterinError1_SCA,ClusterinError2,EstLabels_SCA]=FnSubspaceClusteringErrorFinder(Grps(:,2),GroundTruthLables);
%     ClusterinError1_SCA
%        figure,stem(Grps(:,2));
%     title('SCA')
    %% ADMM
    
        r = 0; affine = false; outlier = true; rho = 1;
        [missrate,C] = SSC(X,r,affine,alpha,outlier,rho,GroundTruthLables);
    %%K-means
    [LabelsKmenas,KmeansClusters]=kmeans(X',Actualc);
   
    [ClusterinError1_Kmeans,ClusterinError2,EstLabels_Kmeans]=FnSubspaceClusteringErrorFinder(LabelsKmenas,GroundTruthLables);
    ClusterinError1_Kmeans
 figure,stem(EstLabels_Kmeans)
        title('k-Means')
    Missrate_Kmeans = Misclassification(LabelsKmenas,GroundTruthLables)
% end