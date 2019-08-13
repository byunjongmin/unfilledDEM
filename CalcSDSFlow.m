function [steepestDescentSlope,slopeAllNbr,SDSFlowDirection ...
    ,SDSNbrY,SDSNbrX] = CalcSDSFlow(DEM,dX,dY)
% @file CalcSDSFlow.m
% @brief Calculate flow direction at each cell based on the steepest
%   descent slope method
%
% @section INTRO CalcSDSFlow
%
% @version 0.1.5 / 2015-11-26
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
Y = mRows - 2; X = nCols - 2;
Y_TOP_BND = 1; Y_BOTTOM_BND = mRows;
Y_INI = 2; Y_MAX = Y+1;
X_LEFT_BND = 1; X_RIGHT_BND = nCols;
X_INI = 2; X_MAX = X+1;

QUARTER_PI = pi * 0.25;
DIAGONAL_DIST = sqrt(dX^2 + dY^2);
DISTANCE_TO_NBR = [dX,DIAGONAL_DIST,dY,DIAGONAL_DIST ...
                    ,dX,DIAGONAL_DIST,dY,DIAGONAL_DIST];

[arrayX,arrayY] ...
    = meshgrid(X_LEFT_BND:X_RIGHT_BND,Y_TOP_BND:Y_BOTTOM_BND);

% Note that calculation is performed within the boundary
[sArrayX,sArrayY] = meshgrid(X_INI:X_MAX,Y_INI:Y_MAX);

% matrix of local cell(e0) index
e0LinearIndicies = (arrayX-1) * mRows + arrayY;
% elevation of local cell(e0)
sE0LinearIndicies = (sArrayX-1) * mRows + sArrayY;

% 3d array to indicate the index of the next downstream neighbor
s3IthNbrLinearIndicies = zeros(Y,X,8);

ithNbrYOffset = [0 -1 -1 -1  0  1  1  1];
ithNbrXOffset = [1  1  0 -1 -1 -1  0  1];

for i = 1:8    
    s3IthNbrLinearIndicies(:,:,i) ...
        = e0LinearIndicies(sE0LinearIndicies ...
        + (ithNbrXOffset(i) * mRows + ithNbrYOffset(i)));
end

SDSFlowDirection = nan(mRows,nCols);
sSDSFlowDirection = nan(Y, X); % without bounday
steepestDescentSlope = nan(mRows,nCols);
sSteepestDescentSlope = -inf(Y, X); % without boundary
slopeAllNbr = nan(mRows,nCols,8);

[SDSNbrX,SDSNbrY] = meshgrid(X_LEFT_BND:X_RIGHT_BND,Y_TOP_BND:Y_BOTTOM_BND);
[sSDSNbrX,sSDSNbrY] = meshgrid(X_INI:X_MAX,Y_INI:Y_MAX);

sE0Elevation = DEM(sE0LinearIndicies);

% main body to calculate the slope value for each neighbor
initialNbrX = sSDSNbrX;
initialNbrY = sSDSNbrY;

for k = 1:8

    sKthNbrElevation = DEM(s3IthNbrLinearIndicies(:,:,k));
    sKthNbrSlope = (sE0Elevation - sKthNbrElevation) / DISTANCE_TO_NBR(k);
    slopeAllNbr(Y_INI:Y_MAX,X_INI:X_MAX,k) = sKthNbrSlope;

    biggerSlope = sKthNbrSlope > sSteepestDescentSlope;
    sSDSNbrY(biggerSlope) = initialNbrY(biggerSlope) + ithNbrYOffset(k);
    sSDSNbrX(biggerSlope) = initialNbrX(biggerSlope) + ithNbrXOffset(k);
    
    possitiveSlope = sKthNbrSlope > 0;
    sSDSFlowDirection(biggerSlope & possitiveSlope) = (k-1) * QUARTER_PI;
    sSteepestDescentSlope(biggerSlope) = sKthNbrSlope(biggerSlope);

end % for i = 1:8

% output
SDSFlowDirection(Y_INI:Y_MAX,X_INI:X_MAX) = sSDSFlowDirection;
steepestDescentSlope(Y_INI:Y_MAX,X_INI:X_MAX) = sSteepestDescentSlope;
SDSNbrX(Y_INI:Y_MAX,X_INI:X_MAX) = sSDSNbrX;
SDSNbrY(Y_INI:Y_MAX,X_INI:X_MAX) = sSDSNbrY;