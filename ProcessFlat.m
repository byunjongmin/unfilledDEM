function flatRegMap = ProcessFlat(DEM,targetDrainage,slopeAllNbr)
% @file ProcessFlat.m
% @brief Identify flats
%
% @param[in] DEM
% @pararn[in] targetDrainage
% @param[in] slopeAllNbr
% @param[in] steepestDescentSlope
% @retval flatRegInfo
% @retval flatRegMap
%
% @version 3.0.2
% @author Jongmin Byun
%==========================================================================
% Constants
[mRows,nCols] = size(DEM);
Y = mRows - 2; X = nCols - 2;
Y_INI = 2; Y_MAX = Y+1;
X_INI = 2; X_MAX = X+1;
OUTER_BOUNDARY = true(mRows,nCols);
OUTER_BOUNDARY(Y_INI:Y_MAX,X_INI:X_MAX) = false;
flatRegMap = zeros(mRows,nCols);

for i=1:8
    iFlatRegMap = slopeAllNbr(:,:,i) == 0 & ~OUTER_BOUNDARY & targetDrainage;
    flatRegMap = iFlatRegMap | flatRegMap;        
end
