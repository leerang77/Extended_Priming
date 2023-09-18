function [numActiveGCs, avgNumBcells] = getGCStat(gcnum, idx)
    numB = squeeze(gcnum{idx}(:,113,1))+squeeze(gcnum{idx}(:,226,1));
    avgNumBcells = mean(numB(numB>0));
    numActiveGCs = sum(sum(squeeze(gcnum{idx}(:,:,1)),2)>0);
end
