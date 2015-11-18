function [nUpstreamCells,markForGoneCells] ...
    = mCalcUpstreamCells(nUpstreamCells,DEM,targetDrainage ...
    ,m2SDSNbrY,m2SDSNbrX,iSubFldRegMap,markForGoneCells)
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
% @param[in] markForGoneCells
% @retval nUpstreamCells
% @retval markForGoneCells
%
% @version 2.0.2 / 2015-11-14
% @author Jongmin Byun
%==========================================================================

% constants
GONE_WSD = 4;

% sort cells according to sub-flooded region ID, elevation
[vtY,vtX] = find(targetDrainage); % vector Y, vector X
vtDEM = DEM(targetDrainage); % vector DEM
vtISubFldReg ... % vector of sub-flooded region ID
    = iSubFldRegMap(targetDrainage);

sortedDEMYX = [vtY,vtX,vtISubFldReg,vtDEM];
sortedDEMYX = sortrows(sortedDEMYX,[3,-4]);

nTotalCells = size(vtY,1);

% calculate upstream cells number
for i = 1:nTotalCells
    
    iCellY = sortedDEMYX(i,1); % Y and X coordinate of ith cell
    iCellX = sortedDEMYX(i,2);
    
    % if iCell is involved in current sub-flooded region
    if iSubFldRegMap(iCellY,iCellX) == true

        nUpstreamCells(iCellY,iCellX) ...
            = nUpstreamCells(iCellY,iCellX) + 1;
     
    end
    
    dStreamNbrY = m2SDSNbrY(iCellY,iCellX); % down stream neighbour Y
    dStreamNbrX = m2SDSNbrX(iCellY,iCellX); % down stream neighbour X
        
    if iSubFldRegMap(dStreamNbrY,dStreamNbrX) == true
        
        nUpstreamCells(dStreamNbrY,dStreamNbrX) ...
            = nUpstreamCells(dStreamNbrY,dStreamNbrX) ...
            + nUpstreamCells(iCellY,iCellX);
    end
end