function [steepestDescentSlope,slopeAllNbr,SDSFlowDirection ...
    ,SDSNbrY,SDSNbrX] ...
    = CalcSDSFlow(DEM,dX,dY)
% @file CalcSDSFlow.m
% @brief Calculate flow direction at each cell based on the steepest
%   descent slope method
%
% @section INTRO CalcSDSFlow
%
% @version 0.1
%
% @retval steepestDescentSlope: steepest descent slope
% @retval slopeAllNbr: all slope values of 8 neighbors
% @retval SDSFlowDirection: flow direction
% @retval SDSNbrY: downstream neighbor Y coordinate
% @retval SDSNbrX: downstream neighbor X coordinate
% @param mRows
% @param nCols
% @param DEM
% @param dX
%
% =========================================================================

% Constants
[mRows,nCols] = size(DEM);
Y = mRows - 2;
X = nCols - 2;
Y_TOP_BND = 1;
Y_BOTTOM_BND = mRows;
Y_INI = 2;
Y_MAX = Y+1;
X_LEFT_BND = 1;
X_RIGHT_BND = nCols;
X_INI = 2;
X_MAX = X+1;

QUARTER_PI = 0.785398163397448;     % pi * 0.25
ROOT2 = 1.41421356237310;           % sqrt(2)
DISTANCE_RATIO_TO_NBR = [1 ROOT2 1 ROOT2 1 ROOT2 1 ROOT2];

[arrayX,arrayY] ...
    = meshgrid(X_LEFT_BND:X_RIGHT_BND,Y_TOP_BND:Y_BOTTOM_BND);

% Notet that calculation is performed within the boundary
[sArrayX,sArrayY] = meshgrid(X_INI:X_MAX,Y_INI:Y_MAX);

% matrix of local cell(e0) index
e0LinearIndicies = (arrayX-1) * mRows + arrayY;
% elevation of local cell(e0)
sE0LinearIndicies = (sArrayX-1) * mRows + sArrayY;

% (최대하부경사 유향을 구하는 과정에서) 8 방향 이웃 셀을 가리키는 3차원 색인 배열
s3IthNbrLinearIndicies = zeros(Y,X,8);

ithNbrYOffset = [0 -1 -1 -1  0  1  1  1];
ithNbrXOffset = [1  1  0 -1 -1 -1  0  1];

for ithDir = 1:8    
    s3IthNbrLinearIndicies(:,:,ithDir) ...
        = e0LinearIndicies(sE0LinearIndicies ...
        + (ithNbrXOffset(ithDir) * mRows + ithNbrYOffset(ithDir)));
end

SDSFlowDirection = nan(mRows,nCols);
sSDSFlowDirection = nan(Y, X);
steepestDescentSlope = nan(mRows,nCols);
sSteepestDescentSlope = -inf(Y, X);
slopeAllNbr = nan(mRows,nCols,8);

[SDSNbrX,SDSNbrY] = meshgrid(X_LEFT_BND:X_RIGHT_BND,Y_TOP_BND:Y_BOTTOM_BND);
[sSDSNbrX,sSDSNbrY] = meshgrid(X_INI:X_MAX,Y_INI:Y_MAX);

sE0Elevation = DEM(sE0LinearIndicies);

initialNbrX = sSDSNbrX;
initialNbrY = sSDSNbrY;

for ithNbr = 1:8

    sKthNbrElevation = DEM(s3IthNbrLinearIndicies(:,:,ithNbr));

    sKthNbrSlope = (sE0Elevation - sKthNbrElevation) ...
    / double((DISTANCE_RATIO_TO_NBR(ithNbr) * dX));

    slopeAllNbr(Y_INI:Y_MAX,X_INI:X_MAX, ithNbr) = sKthNbrSlope;

    biggerSlope = sKthNbrSlope > sSteepestDescentSlope;

    sSDSNbrY(biggerSlope) = initialNbrY(biggerSlope) ...
        + ithNbrYOffset(ithNbr);
    sSDSNbrX(biggerSlope) = initialNbrX(biggerSlope) ...
        + ithNbrXOffset(ithNbr);
    
    possitiveSlope = sKthNbrSlope > 0;

    sSDSFlowDirection(biggerSlope & possitiveSlope) ...
      = (ithNbr-1) * QUARTER_PI;

    sSteepestDescentSlope(biggerSlope) = sKthNbrSlope(biggerSlope);

end % for ithNbr = 1:8

SDSFlowDirection(Y_INI:Y_MAX,X_INI:X_MAX) = sSDSFlowDirection;
steepestDescentSlope(Y_INI:Y_MAX,X_INI:X_MAX) = sSteepestDescentSlope;
SDSNbrX(Y_INI:Y_MAX,X_INI:X_MAX) = sSDSNbrX;
SDSNbrY(Y_INI:Y_MAX,X_INI:X_MAX) = sSDSNbrY;