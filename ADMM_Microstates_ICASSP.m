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
RealL=Fs;%256;%384;%Fs;
Sigma=0.01;%1e-4; %additive noise level
Actualc=4;
%% Generating Source Matrix
S=zeros(n,RealL);
for i=1:n
    Temp=[zeros(1,(i-1)*L*OL) Cosine zeros(1,n*L-L-(i-1)*L*OL)] ;%figure,plot(S);
    S(i,:)=Temp(:,1:RealL);
    subplot(n,1,i)
    plot(1:RealL,S(i,:))
end
S=[S zeros(n,512)];
S(2,256:256+127)=Cosine;
S(4,1:64)=Cosine(length(Cosine)/2+1:end);
S(1,320:320+length(Cosine)-1)=Cosine;
S(3,382:382+length(Cosine)-1)=Cosine;
S(4,447:447+length(Cosine)-1)=Cosine;
% S(4,447:447+length(Cosine)-1)=Cosine;
S(2,510:510+length(Cosine)-1)=Cosine;
S(1,638:638+length(Cosine)-1)=Cosine;
 S(4,697:697+length(Cosine)-1)=Cosine;
S(3,572:572+length(Cosine)-1)=Cosine;

ColorLab={'r','b','g','m'};
 S=S(:,1:757);
figure,
for i=1:n
%     Temp=[zeros(1,(i-1)*L*OL) Cosine zeros(1,n*L-L-(i-1)*L*OL)] ;%figure,plot(S);
%     S(i,:)=Temp(:,1:RealL);
    subplot(n,1,i)
    plot(1:size(S,2),S(i,:),ColorLab{i},'LineWidth',2)
    set(gca, 'fontsize', 10,'fontweight','bold','linewidth',1);
end

currentFigure = gcf;
title(currentFigure.Children(end), 'Simulated Overlapping Dipole Activation');
xlabel('Time [Sample]');

% xlim([0 637]);
%% Groundthruth labels for each time point
% SCat=[S(1,:) S(2,:) S(3,:)  S(4,:) ];
% figure, plot(SCat);
% bb=zeros(1,14);
% bb(1)=1;
% [aa,bb(2:end)]=findpeaks(SCat);

Labs=[1:5 2 6 4 5 3 6 1];
GroundTruthLables=[];
% for yy=1:length(bb)-1
% GroundTruthLables(bb(yy):bb(yy+1)-1)=Labs(yy);
% end


for tt=1:n
[aaa,bbc{tt}]=findpeaks(S(tt,:));
end
bbcMat=sort(cell2mat(bbc));
bbcMat=[1 bbcMat];
for yy=1:length(bbcMat)-1
GroundTruthLables(bbcMat(yy):bbcMat(yy+1)-1)=Labs(yy);
end
GroundTruthLables(bbcMat(end):size(S,2))=Labs(end);
% GroundTruthLables(1:bb(:,1))=1;
% GroundTruthLables(1:bb(:,1))=1;
% Actualc=6;
% % S=S(:,1:384);
% for j=1:5
%     GroundTruthLables=[GroundTruthLables j*ones(1,L/2)];
% end
% Labs=[1,6,4,5,2,6,3];
% for j=1:7
%     GroundTruthLables=[GroundTruthLables Labs(j)*ones(1,L/2)];
% end

figure, stem(GroundTruthLables)
title('Groundtruth subspace clustering');
ylim([0,6.3]);
%% Mixing (X=AS+E)
X=A*S+Sigma*randn(m,size(S,2));
%% PCA projection
[LabelsKmenas,KmeansClusters]=kmeans(X',4);
[~,StartIndClusters,~]=unique(LabelsKmenas, 'stable');
StartIndClusters(5)=size(S,2);
r=k*Actualc; %dimension of the space to project the data to
%%Enter the projection dimension e.g. r = d*n, enter r = 0 to not project
% Colors={'.r', '.b', '.g', '.m','-y','..r'}; % Colors for clusters
Colors=ColorLab;

% PCA projection
[Coeff,Score,Latent,TSquared] =pca(X');
% Data Projection based on PCA
% Xp = FnDataProjection(X,r,'PCA');
% figure,
% for j=1:Actualc
%     scatter3(Xp(1,StartIndClusters(j):StartIndClusters(j+1)),Xp(2,StartIndClusters(j):StartIndClusters(j+1)),...
%         Xp(3,StartIndClusters(j):StartIndClusters(j+1)),Colors{j});
%     hold on
% end
% title('PCA projection 1')
% 
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
alpha=10;
alpha=0.5*(1.25/std(X(:)))^2;
alpha=60;
r =n; affine = 0; outlier = true; rho = 1;
[missrate,C,grps,Clus] = SSC(X,r,affine,alpha,outlier,rho,GroundTruthLables);
Sub=FnSubspaceCalc(X,grps',1,Actualc);
Subb=squeeze(Sub)
figure,stem(grps)
title('SSC-ADMM')
ylim([0, 6.3]);
[Sub,OrthSub]=FnSubspaceCalc(X,grps',2,6);
Subb=squeeze(Sub);
% Ahat=FnMixingCalc(OrthSub,2,6,4)
[Ahat,MinEV]=FnMixingCalc_Ehsan(OrthSub,2,6,4)
[MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error,Error] = FnMixingIdentificationError(A,Ahat)

   for SubInd=1:size(OrthSub,3)
                    for DataInd=1:size(X,2)
                        DistADMM(SubInd,DataInd)=norm(OrthSub(:,:,SubInd)'*X(:,DataInd));
                    end
                end
                figure,stem(DistADMM')
%%K-means
 %%K-means
 [LabelsKmenas,KmeansClusters]=kmeans(Xhat',6);
 for SubInd=1:size(KmeansClusters,1)
                    for DataInd=1:size(X,2)
                        DistKmeans(SubInd,DataInd)=norm(Xhat(:,DataInd)-(Xhat(:,DataInd).*KmeansClusters(SubInd,:)').*KmeansClusters(SubInd,:)');
                    end
                end
                figure,stem(DistKmeans')
                
 [MixingIdentificatioError2, MixingVectorerror2,NMSE2,NMSSum2,AhatNew2,Norm2Error2,Error2] = FnMixingIdentificationError(A,KmeansClusters')

 figure,stem(LabelsKmenas)
 ClsusK=FnColNormalizer(KmeansClusters');
  [MixingIdentificatioError2, MixingVectorerror2,NMSE2,NMSSum2,AhatNew2,Norm2Error2,Error2] = FnMixingIdentificationError(A,ClsusK)

 title('k-Means')
 ylim([0 6.1]);
 Missrate_Kmeans = Misclassification(LabelsKmenas,GroundTruthLables)
 %% k-SCA RANSAC
 

Th_RansacSubspace=1e-5;% RANSAC Threshold in Distance Function for subspace estimation
ThBBC=1e-3; %BBC clustering for mixing idetification
Th_RansacMixing=1e-5;% RANSAC Threshold in Distance Function for mixing matrix identifiation

DualMode=0; %Twice subpspace clusteing
BAS_Set=[];NMSE_Set=[];Norm2Error_Set=[];Sub_Error_Set=[];Error_Set=[];
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='kSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
Orth=0; %if Orth=1 A is orth
p=0.95;
% Sigma=0;
CT=0;
SubspaceInds=cell(1,1);
MaxIterEst=1e5;
 [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFindV3(X,Th_RansacSubspace,k,c,SubspaceInds,MaxIterEst,p);
   [Ahat_MyKSCA,MinEV]=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,2,c,n)
 %                 [inliers,M]=FnDistanceClacBetweenSubspaceVector(ComplementOrthofSubSpaces(:,:,1),X(:,1),.1)
 for SubInd=1:size(ComplementOrthofSubSpaces,3)
     for DataInd=1:size(X,2)
         Dist(SubInd,DataInd)=norm(ComplementOrthofSubSpaces(:,:,SubInd)'*X(:,DataInd));
     end
 end
 figure, stem(Dist')
                %% Mixing Matrix Identification.
                % Option 1
                [Ahat,MinEV]=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,k,c,n); % Slow/High comp. cost/ needs all comb i.e. Order=Cx=nchoosek(1:c,g)
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
ChanLocs = readlocs('ChanLocsEE.loc'); % Read 64 EEG channels (Biosemi 10-10);

for mm=1:4
    figure,
    topoplot(AhatNew(:,mm), ChanLocs,'electrodes','on','style','map','plotrad',.54);
end

for mm=1:4
    figure,
    topoplot(A(:,mm), ChanLocs,'electrodes','on','style','map','plotrad',.54);
end

for mm=1:4
    figure,
    topoplot(AhatNew2(:,mm), ChanLocs,'electrodes','on','style','map','plotrad',.54);
end