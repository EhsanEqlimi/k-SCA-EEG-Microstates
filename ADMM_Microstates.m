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
A=FnColNormalizer(randn(m,n)); %Normaly distrubuted Mixing vectors
L=Fs/4; %The length of each source/microstate (the number of samples)
Cosine= tukeywin(L,b)';% The waveform/function for sources
GroundTruthLables=[];
OL=.5;% Overlap
RealL=256;%384;%Fs;
Sigma=0;%1e-4; %additive noise level
Actualc=4;
%% Generating Source Matrix
S=zeros(n,RealL);

for i=1:n
    Temp=[zeros(1,(i-1)*L*OL) Cosine zeros(1,n*L-L-(i-1)*L*OL)] ;%figure,plot(S);
    S(i,:)=Temp(:,1:RealL);
    subplot(n,1,i)
    plot(1:RealL,S(i,:))
end
S(2,:)=0;
currentFigure = gcf;
title(currentFigure.Children(end), 'Simulated Source (dipole) activation over time');
xlabel('Time [Sample]');
% xlim([0 RealL]);
%% Groundthruth labels for each time point
Actualc=4;
% S=S(:,1:384);
for j=1:Actualc
    GroundTruthLables=[GroundTruthLables j*ones(1,L/2)];
end
figure, plot(GroundTruthLables)
title('Groundtruth lables over time');
%% Mixing (X=AS+E)
X=A*S+Sigma*randn(m,RealL);
%% PCA projection
[LabelsKmenas,KmeansClusters]=kmeans(X',Actualc);
ChanLocs = readlocs('ChanLocsEE.loc'); % Read 64 EEG channels (Biosemi 10-10);

for mm=1:4
    figure,
    topoplot(KmeansClusters(mm,:), ChanLocs,'electrodes','on','style','map','plotrad',.54);
end

for mm=1:4
    figure,
    topoplot(A(:,mm), ChanLocs,'electrodes','on','style','map','plotrad',.54);
end
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
for j=1:n
    scatter3(Xp(1,StartIndClusters(j):StartIndClusters(j+1)),Xp(2,StartIndClusters(j):StartIndClusters(j+1)),...
        Xp(3,StartIndClusters(j):StartIndClusters(j+1)),Colors{j});
    hold on
end
title('PCA projection 1')

figure,
for j=1:n
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
figure,plot(X')
title('Mixture Signals');
 Mean=mean(X,1);
Xhat=X-repmat(mean(X,1),size(X,1),1);
figure,plot(Xhat')
title('Mixture Signals- Zero Mean');
alpha=20;
r = 0; affine = false; outlier = true; rho = 1;
[missrate,C,grps,Clus] = SSC(X,r,affine,alpha,outlier,rho,GroundTruthLables);
Sub=FnSubspaceCalc(X,grps',1,Actualc);
Subb=squeeze(Sub)
figure,stem(grps)
title('SCA-ADMM-Alpha-10')

[Sub,OrthSub]=FnSubspaceCalc(X,grps',2,4);
Subb=squeeze(Sub);
% Ahat=FnMixingCalc(OrthSub,2,6,4)
[Ahat,MinEV]=FnMixingCalc_Ehsan(OrthSub,2,4,4)
[MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error,Error] = FnMixingIdentificationError(A,Ahat)

%%K-means
 %%K-means
 [LabelsKmenas,KmeansClusters]=kmeans(Xhat',Actualc);
 figure,stem(LabelsKmenas)
 title('k-means')
 Missrate_Kmeans = Misclassification(LabelsKmenas,GroundTruthLables)
 %% Kmeans-microsates

 nSegments=1;
 ClustPars = struct('MinClasses',4,'MaxClasses',8,'GFPPeaks',false,'IgnorePolarity',false,'MaxMaps',500,'Restarts',5', 'UseAAHC',false,'Normalize',false);
 FitPars = struct('nClasses',4,'lambda',1,'b',20,'PeakFit',false, 'BControl',false,'Rectify',false,'Normalize',false);

  msinfo=FnMicrostateCalc(Xhat,nSegments,ClustPars);
  Maps = NormDimL2(msinfo.MSMaps(FitPars.nClasses).Maps,2);
  eegdata.data=Xhat;
  eegdata.event=0;
  eegdata.srate=512;
  [MSClass,gfp,ExpVar] = AssignMStates(eegdata,Maps,FitPars,ClustPars.IgnorePolarity);
figure,stem(MSClass)
title('Microstates-AAHC');
 Missrate_AAHC = Misclassification(MSClass,GroundTruthLables)
%% Test on rela EEG
 