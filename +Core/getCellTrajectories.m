function [intDistCurve] = getCellTrajectories(HIVMovie,scalingFactor)
    assert(isobject(HIVMovie),'Function input expected to be an object')
    assert(isa(HIVMovie,'Core.HIVCellMovie'),'Function input expected to be an object of class HIVCellMovie')
    
    filename = [HIVMovie.raw.movInfo.path2Cell filesep 'mask.mat'];
    mask = load(filename);
    field = fieldnames(mask);
    mask = mask.(field{1});
    
    traces = HIVMovie.traces3D;
    
    weight = [HIVMovie.raw.movInfo.pxSize/1e3 HIVMovie.raw.movInfo.pxSize/1e3 HIVMovie.info.dZ/1e3];
    
    %First, we need to have a weighted distance map from the cell and in the
    %cell
    distCell = zeros(size(mask));
    for i = 1:size(mask,4)
        currMask = mask(:,:,:,i);
        outCell = DistMap.calcWeightedDistMap(currMask,weight);
    
        inCell  = DistMap.calcWeightedDistMap(~currMask,weight);

        outCell(outCell==0) = -inCell(outCell==0);

        distCell(:,:,:,i) = outCell;
    end
    %clean up some memory
    clear mask
    %finally, for all trajectory we can build a trace as a distance map
    %from the cells and plot the intensity.
    
    %1 we interpolate the distance map to keep some resolution on
    %localization
    nDistCell = zeros([size(distCell,1)*scalingFactor,size(distCell,2)*scalingFactor,...
        size(distCell,3)*scalingFactor,size(distCell,4)]);
    for i = 1:size(distCell,4)
        nDistCell(:,:,:,i) = imresize3(distCell(:,:,:,i),[size(distCell,1)*scalingFactor,size(distCell,2)*scalingFactor,...
        size(distCell,3)*scalingFactor]);
    end
    distCell = nDistCell;
    clear nDistCell;
    
    intDistCurve = cell(size(traces));
    %loop through the traces
    for i =1:length(traces)
        currTrace = traces{i};
        
        currCol = round(currTrace.colM*scalingFactor);
        currRow = round(currTrace.rowM*scalingFactor);
        currZ   = round((currTrace.zM+1)*scalingFactor);
        t = currTrace.t;
        currInt = currTrace.intensity;
        currIntZ = currTrace.intZ;
        currDist = zeros(size(currCol));
        for j = 1: height(currTrace)
            
            currDist(j) = distCell(currRow(j),currCol(j),currZ(j),t(j));
            
        end
        
        
        intDistCurve{i} = table(currDist,currInt,currIntZ,'VariableNames',{'cellDist','intensity','intZ'});       
    end
end