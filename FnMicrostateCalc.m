function msinfo=FnMicrostateCalc(EEgMatrix,nSegments,ClustPar)

MapsToUse = [];
if ~isinf(ClustPar.MaxMaps)
    MapsPerSegment = hist(ceil(nSegments * rand(ClustPar.MaxMaps,1)),nSegments);
else
    MapsPerSegment = inf(nSegments,1);
end

for s = 1:nSegments
    if ClustPar.GFPPeaks == 1
        gfp = std(EEgMatrix(:,:,s),1,1);
        IsGFPPeak = find([false (gfp(1,1:end-2) < gfp(1,2:end-1) & gfp(1,2:end-1) > gfp(1,3:end)) false]);
        if numel(IsGFPPeak) > MapsPerSegment(s) && MapsPerSegment(s) > 0
            idx = randperm(numel(IsGFPPeak));
            IsGFPPeak = IsGFPPeak(idx(1:MapsPerSegment(s)));
        end
        MapsToUse = [MapsToUse EEgMatrix(:,IsGFPPeak,s)];
    else
        if (size(EEgMatrix,2) > ClustPar.MaxMaps) && MapsPerSegment(s) > 0
            idx = randperm(size(EEgMatrix,2));
            MapsToUse = [MapsToUse EEgMatrix(:,idx(1:MapsPerSegment(s)),s)];
        else
            MapsToUse = [MapsToUse EEgMatrix(:,:,s)];
        end
    end
end

flags = '';
if ClustPar.IgnorePolarity == false
    flags = [flags 'p'];
end
if ClustPar.Normalize == true
    flags = [flags 'n'];
end

if ClustPar.UseAAHC == false
    for nClusters = ClustPar.MinClasses:ClustPar.MaxClasses
        [b_model,~,~,exp_var] = eeg_kMeans(MapsToUse',nClusters,ClustPar.Restarts,[],flags);
        
        msinfo.MSMaps(nClusters).Maps = b_model;
        msinfo.MSMaps(nClusters).ExpVar = double(exp_var);
        msinfo.MSMaps(nClusters).ColorMap = lines(nClusters);
        msinfo.MSMaps(nClusters).SortMode = 'none';
        msinfo.MSMaps(nClusters).SortedBy = '';
        msinfo.MSMaps(nClusters).Communality= [];
    end
else
    [b_model,exp_var] = eeg_computeAAHC(double(MapsToUse'),ClustPar.MinClasses:ClustPar.MaxClasses,false, ClustPar.IgnorePolarity,ClustPar.Normalize);
    for nClusters = ClustPar.MinClasses:ClustPar.MaxClasses
        msinfo.MSMaps(nClusters).Maps = b_model{nClusters-ClustPar.MinClasses+1};
        msinfo.MSMaps(nClusters).ExpVar = exp_var(nClusters-ClustPar.MinClasses+1); %Ehsan
        msinfo.MSMaps(nClusters).ColorMap = lines(nClusters);
        msinfo.MSMaps(nClusters).SortMode = 'none';
        msinfo.MSMaps(nClusters).SortedBy = '';
        msinfo.MSMaps(nClusters).Communality= [];
    end
end

