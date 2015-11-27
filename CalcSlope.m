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
% Version: 0.1.2 / 2015-11-27
%

% define variables
nCNbr = numel(considerNbr);
nCells = numel(distFromInit);
slopePerDist = zeros(nCells,nCNbr);

% calculate slope values

% % if you do not want negative slope values, remove depressions in an input
% % stram profile
% for i = nCells:-1:2
%     if streamProfElev(i-1) <= streamProfElev(i)
%         streamProfElev(i-1) = streamProfElev(i) + 0.000001;
%     end
% end

% calculate slope values
for ithDist = 1:nCNbr

    % first node: forward difference
    firstElev = streamProfElev(1,ithProf);
	secondElev = streamProfElev(2,ithProf);
	distBtwCells = distFromInit(2) - distFromInit(1);

    slopePerDist(1,ithDist) = (firstElev - secondElev) / distBtwCells;
    
    if slopePerDist(1,ithDist) <= 0.000001 * 2 / distBtwCells
        slopePerDist(1,ithDist) = nan;
    end

    % end node: backward difference
    firstElev = streamProfElev(nCells-1,ithProf);
    secondElev = streamProfElev(nCells,ithProf);
    distBtwCells = distFromInit(nCells) - distFromInit(nCells-1);

    slopePerDist(nCells,ithDist) = (firstElev - secondElev) / distBtwCells;
    
    if slopePerDist(nCells,ithDist) <= 0.000001 * 2 / distBtwCells
        slopePerDist(nCells,ithDist) = nan;
    end
    
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
            
            if slopePerDist(ithCell,ithDist) <= 0.000001 * (secondEdge - firstEdge) / distBtwCells
                slopePerDist(ithCell,ithDist) = nan;
            end

	    end

		% lower boundary    
	    for ithCell = nCells-considerNbr(ithDist)+1:nCells-1

	        firstEdge = ithCell - (nCells - ithCell);
	        secondEdge = ithCell + (nCells - ithCell);

	        firstElev = streamProfElev(firstEdge,ithProf);
	        secondElev = streamProfElev(secondEdge,ithProf); 

	        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
	        slopePerDist(ithCell,ithDist) = (firstElev - secondElev) / distBtwCells;
            
            if slopePerDist(ithCell,ithDist) <= 0.000001 * (secondEdge - firstEdge) / distBtwCells
                slopePerDist(ithCell,ithDist) = nan;
            end

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
        
        if slopePerDist(ithCell,ithDist) <= 0.000001 * (secondEdge - firstEdge) / distBtwCells
            slopePerDist(ithCell,ithDist) = nan;
        end
	        
    end

end