function [slopePerDist,slopePerDist_TF] ...
    = CalcSlope(considerNbr,distFromInit,streamProfElev,ithProf)
%
% function CalcSlope
% 
% 	Calculate slope according to the input window size and base DEM over
%   the defined stream path
%
% Input Variables
%
%	considerNbr: array consisting of numbers of considered neighbors
%	distFromInit: distance from initial point
%	streamProfElev: stream profile elevation smoothed according to window sizes
%   ithProf: row of streamProfElev
%	
% Output Variable
%
%	slopePerDist
%
% Version: 0.1.5 / 2016-01-05
%

if min(considerNbr) <= 1
    error('CalcSlope function requires a value for the condierNbr variable of more than 1.');
end

% define variables
nCNbr = numel(considerNbr);
nCells = numel(distFromInit);
slopePerDist = zeros(nCells,nCNbr);
slopePerDist_TF = false(nCells,nCNbr);

% calculate slope values
for ithDist = 1:nCNbr

    % boundary nodes
    % first node: forward difference
    firstElev = streamProfElev(1,ithProf);
	secondElev = streamProfElev(2,ithProf);
    
	distBtwCells = distFromInit(2) - distFromInit(1);
    slopePerDist(1,ithDist) = (firstElev - secondElev) / distBtwCells;
    % slopePerDist_TF(1,ithDist) = false;
    
    % upper boundary: central difference with variable ranges
    for ithCell = 2:considerNbr(ithDist)

        firstEdge = 1;
        secondEdge = ithCell + (ithCell - 1);

        firstElev = streamProfElev(firstEdge,ithProf);
        secondElev = streamProfElev(secondEdge,ithProf); 

        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
        slopePerDist(ithCell,ithDist) = (firstElev - secondElev) / distBtwCells;
        % slopePerDist_TF(ithCell,ithDist) = false;

    end
    
	% mid nodes: central difference	
    for ithCell = considerNbr(ithDist)+1:nCells-considerNbr(ithDist)

		firstEdge = ithCell - considerNbr(ithDist);
		secondEdge = ithCell + considerNbr(ithDist);

        firstElev = streamProfElev(firstEdge,ithProf);
        secondElev = streamProfElev(secondEdge,ithProf); 

        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
        slopePerDist(ithCell,ithDist) = (firstElev - secondElev) / distBtwCells;
        slopePerDist_TF(ithCell,ithDist) = true;
        
    end
    
    % lower boundary: central difference with variable ranges  
    for ithCell = nCells-considerNbr(ithDist)+1:nCells-1

        firstEdge = ithCell - (nCells - ithCell);
        secondEdge = nCells;

        firstElev = streamProfElev(firstEdge,ithProf);
        secondElev = streamProfElev(secondEdge,ithProf); 

        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
        slopePerDist(ithCell,ithDist) = (firstElev - secondElev) / distBtwCells;
        % slopePerDist_TF(ithCell,ithDist) = false;

    end

    % end node: backward difference
    firstElev = streamProfElev(nCells-1,ithProf);
    secondElev = streamProfElev(nCells,ithProf);
    
    distBtwCells = distFromInit(nCells) - distFromInit(nCells-1);
    slopePerDist(nCells,ithDist) = (firstElev - secondElev) / distBtwCells;
    % slopePerDist_TF(nCells,ithDist) = false;

end