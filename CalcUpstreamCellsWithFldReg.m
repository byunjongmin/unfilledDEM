function upstreamCellsNo = CalcUpstreamCellsWithFldReg(DEM,chosenWatershed ...
    ,flood,SDSNbrY,SDSNbrX,floodedRegionIndex,floodedRegionCellsNo)
% @file CalcUpstreamCellsWithFldReg.m
% @brief Calculate upstream cells considering flooded region
%
% @version 0.1.1. / 2015-11-13
% @author Jongmin Byun
%==========================================================================

% constant
[mRows,nCols] = size(DEM);
Y = mRows - 2;
X = nCols - 2;
Y_INI = 2;
Y_MAX = Y+1;
X_INI = 2;
X_MAX = X+1;
FLOODED = 2;

% variables
[sArrayX,sArrayY] = meshgrid(X_INI:X_MAX,Y_INI:Y_MAX);
vectorY = reshape(sArrayY,[],1);
vectorX = reshape(sArrayX,[],1);
upstreamCellsNo = zeros(mRows,nCols);

% sort cells according to elevation
% note that it excludes flooded region and areas out of chosenWatershed
orgDEM = DEM;
DEM(flood == FLOODED | ~chosenWatershed) = -9999;
vectorDEM = reshape(DEM(Y_INI:Y_MAX,X_INI:X_MAX),[],1);

% sort cells according to whether it is outlet or not
fldRegOutlet = floodedRegionCellsNo > 0;
vectorFldRegOutlet = reshape(fldRegOutlet(Y_INI:Y_MAX,X_INI:X_MAX),[],1);

% for the cells with same elevation as well as of outlet, sort cells
% according to the averaged elevation of flooded region
outletIdx = find(floodedRegionIndex < 0);
meanElevFldReg = zeros(mRows,nCols);
for ithOutlet = 1:numel(outletIdx)
    
    ithFldRegID = -floodedRegionIndex(outletIdx(ithOutlet));
    ithFldRegIdx = floodedRegionIndex == ithFldRegID;
    meanElevFldReg(outletIdx(ithOutlet)) = mean(orgDEM(ithFldRegIdx));
    
end

vectorMeanElevFldReg ...
    = reshape(meanElevFldReg(Y_INI:Y_MAX,X_INI:X_MAX),[],1);
sortedDEMYX = [vectorY,vectorX,vectorDEM,vectorFldRegOutlet,vectorMeanElevFldReg];
sortedDEMYX = sortrows(sortedDEMYX,[-3,4,-5]);
consideringCellsNo = find(vectorDEM > -9999);
consideringCellsNo = size(consideringCellsNo,1);

% calculate upstream cells number in the order of elevation
for i = 1:consideringCellsNo

    ithCellY = sortedDEMYX(i,1);
    ithCellX = sortedDEMYX(i,2);
    
    % transfer the number of upstream cells to next cell
    
    % firstly, add one for itself and the number of upstream cells
    if floodedRegionIndex(ithCellY,ithCellX) >= 0
    
        upstreamCellsNo(ithCellY,ithCellX) ...
            = upstreamCellsNo(ithCellY,ithCellX) + 1;
        
    else
        
        % for the outlet of flooded region, plus the number of flooded
        % region cells
        upstreamCellsNo(ithCellY,ithCellX) ...
            = upstreamCellsNo(ithCellY,ithCellX) + 1 ...
            + floodedRegionCellsNo(ithCellY,ithCellX);
        
    end
    
    % transfer upstream cells number to the next cells
    downStreamNbrY = SDSNbrY(ithCellY,ithCellX);
    downStreamNbrX = SDSNbrX(ithCellY,ithCellX);
    
    % note that, if the downstream cell is flooded region, transfer the
    % number to the outelt of its outlet
    if flood(downStreamNbrY,downStreamNbrX) == FLOODED
        
       outletY = SDSNbrY(downStreamNbrY,downStreamNbrX);
       outletX = SDSNbrX(downStreamNbrY,downStreamNbrX);
       
       upstreamCellsNo(outletY,outletX) ...
        = upstreamCellsNo(outletY,outletX) ...
        + upstreamCellsNo(ithCellY,ithCellX);
        
    else
    
        upstreamCellsNo(downStreamNbrY,downStreamNbrX) ...
            = upstreamCellsNo(downStreamNbrY,downStreamNbrX) ...
            + upstreamCellsNo(ithCellY,ithCellX);
    
    end
    
end

% recognize the coordinate of the outlets of flooded regions
[tmpOutletY,tmpOutletX] = find(floodedRegionIndex < 0);
floodedRegionsNo = size(tmpOutletY,1);
for ithFloodedRegion = 1:floodedRegionsNo

    outletY = tmpOutletY(ithFloodedRegion,1);
    outletX = tmpOutletX(ithFloodedRegion,1);

    upstreamCellsNo(floodedRegionIndex ...
        == - floodedRegionIndex(outletY,outletX)) ...
        = upstreamCellsNo(outletY,outletX) - 1;

end % ithFloodedRegion = 1:
