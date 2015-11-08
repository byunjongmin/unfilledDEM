function [nUpstreamCells,doNotAccFlowCell] ...
    = mCalcUpstreamCells(nUpstreamCells,DEM,targetDrainage ...
    ,m2SDSNbrY,m2SDSNbrX,flatRegMap,nFlatRegCells,steepestDescentSlope ...
    ,iSubFldRegMap,doNotAccFlowCell)
% @file mCalcUpstreamCells.m
% @brief Calculate upstream cells of cells within a sub flooded region
%
% @param[in] DEM
% @param[in] targetDrainage
% @param[in] flood
% @param[in] m2SDSNbrY
% @param[in] m2SDSNbrX
% @param[in] flatRegMap
% @param[in) nFlatRegCells
% @param[in] steepestDescentSlope
% @param[in] iSubFldRegMap
% @param[in] doNotAccFlowCell
% @retval nUpstreamCells
% @retval doNotAccFlowCell
%
% @version 2.0
%==========================================================================

% constants
GONE_WSD = 4;

% sort cells according to sub-flooded region ID, elevation, flat outlet
targetDrainage ... % avoid flat without flow direction
    (abs(flatRegMap) > 0 & steepestDescentSlope <= 0) = false;

[vtY,vtX] = find(targetDrainage); % vector Y, vector X
vtDEM = DEM(targetDrainage); % vector DEM

flatOut = flatRegMap < 0; % outlet of flat
vtFlatOut = flatOut(targetDrainage); % vector of flat outlet

vtISubFldReg ... % vector of sub-flooded region ID
    = iSubFldRegMap(targetDrainage);

sortedDEMYX = [vtY,vtX,vtISubFldReg,vtDEM,vtFlatOut];
sortedDEMYX = sortrows(sortedDEMYX,[3,-4,5]);

nTotalCells = size(vtY,1);

% calculate upstream cells number
for iCell = 1:nTotalCells
    
    iCellY = sortedDEMYX(iCell,1); % Y and X coordinate of ith cell
    iCellX = sortedDEMYX(iCell,2);
    
    if iSubFldRegMap(iCellY,iCellX) == true
        % if ith cell is flat outlet, additionally add flat cells without
        % flow direction
        if nFlatRegCells(iCellY,iCellX) > 0
            nUpstreamCells(iCellY,iCellX) ...
                = nUpstreamCells(iCellY,iCellX) + 1 ...
                + nFlatRegCells(iCellY,iCellX);
        else
            nUpstreamCells(iCellY,iCellX) ...
                = nUpstreamCells(iCellY,iCellX) + 1;
        end
        
        % add the current sub-flooded region to the exception list
        doNotAccFlowCell(iCellY,iCellX) = GONE_WSD;
        
    end
    
    dStreamNbrY = m2SDSNbrY(iCellY,iCellX); % down stream neighbour Y
    dStreamNbrX = m2SDSNbrX(iCellY,iCellX); % down stream neighbour X

    if flatRegMap(dStreamNbrY,dStreamNbrX) > 0 ...
        && steepestDescentSlope(dStreamNbrY,dStreamNbrX) == 0
    
        flatOutY = m2SDSNbrY(dStreamNbrY,dStreamNbrX);
        flatOutX = m2SDSNbrX(dStreamNbrY,dStreamNbrX);
        
        nUpstreamCells(flatOutY,flatOutX) ...
            = nUpstreamCells(flatOutY,flatOutX) ...
            + nUpstreamCells(iCellY,iCellX);
        
    else
        
        if iSubFldRegMap(dStreamNbrY,dStreamNbrX) ~= 0
            nUpstreamCells(dStreamNbrY,dStreamNbrX) ...
                = nUpstreamCells(dStreamNbrY,dStreamNbrX) ...
                + nUpstreamCells(iCellY,iCellX);
        end
    end
end

% calculate the number of upstream cells of flat without flow direction
% in each flat
[tmpFlatOutY,tmpFlatOutX] ...
    = find(flatRegMap < 0 & nFlatRegCells > 0 ...
        & iSubFldRegMap);
nFlat = size(tmpFlatOutY,1):

for ithFlat = 1:nFlat

    iFlatOutY = tmpFlatOutY(ithFlat,1);
    iFlatOutX = tmpFlatOutX(ithFlat,1);

    iFlatID = abs(flatRegMap(iFlatOutY,iFlatOutX));
    iFlatWithoutFlwDir ...
        = flatRegMap == iFlatID & steepestDescentSlope <= 0;
    nUpstreamCells(iFlatWithoutFlwDir) ...
        = nUpstreamCells(iFlatOutY,iFlatoutX) - 1;
    
end
