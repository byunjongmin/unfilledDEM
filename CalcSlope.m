function slopePerDist = CalcSlope(considerNbr,distFromInit ...
    ,streamProfElev,ithProf)
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
%

nCNbr = numel(considerNbr);
nCells = numel(distFromInit);
slopePerDist = zeros(nCells,nCNbr);

for ithDist = 1:nCNbr

	% first node: forward difference
	firstElev = streamProfElev(1,ithProf);
    secondElev = streamProfElev(2,ithProf);
    distBtwCells = distFromInit(2) - distFromInit(1);

    slopePerDist(1,ithDist) = (firstElev - secondElev) / distBtwCells;

	% end node: backward difference
	firstElev = streamProfElev(nCells-1,ithProf);
    secondElev = streamProfElev(nCells,ithProf);
    distBtwCells = distFromInit(nCells) - distFromInit(nCells-1);

    slopePerDist(nCells,ithDist) = (firstElev - secondElev) / distBtwCells;

	% boundary nodes: central difference with variable ranges
	if considerNbr(ithDist) > 1
    
		% upper boundary
	    for ithCell = 2:considerNbr(ithDist)

	        firstEdge = ithCell - (ithCell - 1);
	        secondEdge = ithCell + (ithCell - 1);

	        firstElev = streamProfElev(firstEdge,ithProf);
	        secondElev = streamProfElev(secondEdge,ithProf); 

	        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
	        slopePerDist(ithCell,ithDist) = (firstElev - secondElev) / distBtwCells;

	    end

		% lower boundary    
	    for ithCell = nCells-considerNbr(ithDist)+1:nCells-1

	        firstEdge = ithCell - (nCells - ithCell);
	        secondEdge = ithCell + (nCells - ithCell);

	        firstElev = streamProfElev(firstEdge,ithProf);
	        secondElev = streamProfElev(secondEdge,ithProf); 

	        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
	        slopePerDist(ithCell,ithDist) = (firstElev - secondElev) / distBtwCells;

	    end
    
	end
    
	% mid nodes: central difference	
    for ithCell = considerNbr(ithDist)+1:nCells-considerNbr(ithDist)

		firstEdge = ithCell - considerNbr(ithDist);
		secondEdge = ithCell + considerNbr(ithDist);

        firstElev = streamProfElev(firstEdge,end);
        secondElev = streamProfElev(secondEdge,end); 

        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
        slopePerDist(ithCell,ithDist) = (firstElev - secondElev) / distBtwCells;
	        
    end

end