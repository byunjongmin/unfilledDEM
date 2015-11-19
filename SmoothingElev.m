function smoothedElev = SmoothingElev(considerNbr,distFromInit,baseDEM ...
    ,streamPath)
%
% function SmoothingElev
% 
% 	Provide smoothed elevation values according to the input window size
%
% Input Variables
%
%	considerNbr: array consisting of numbers for considering neighbour
%	distFromInit: distance from initial point
%	baseDEM: base DEM for calcuation
%	streamPath: extracted stream path
%	
% Output Variable
%
%	smoothedElev
%
% Version 0.1.1. / 2015-11-19

% output variable
nCNbr = numel(considerNbr); % number of different neighbors
nCells = numel(distFromInit);
smoothedElev = zeros(nCells,nCNbr);

for ithDist = 1:nCNbr

	% first node
	smoothedElev(1,ithDist) = double(baseDEM(streamPath(1)));
	% last node
    smoothedElev(nCells,ithDist) = double(baseDEM(streamPath(nCells)));

	% boundary nodes
	if considerNbr(ithDist) > 1
        
        for ithCell = 2:considerNbr(ithDist)

            firstEdge = ithCell - (ithCell - 1);
            secondEdge = ithCell + (ithCell - 1);

            smoothedElev(ithCell,ithDist) ...
                = mean(baseDEM(streamPath(firstEdge:secondEdge)));

        end
        
        for ithCell = nCells-considerNbr(ithDist)+1:nCells-1

            firstEdge = ithCell - (nCells - ithCell);
            secondEdge = ithCell + (nCells - ithCell);

            smoothedElev(ithCell,ithDist) ...
                = mean(baseDEM(streamPath(firstEdge:secondEdge)));

        end
        
    else
        
        error('Error. The number of considering neighbor should be larger than 1.');
        
    end % if considerNbr(ithDist) > 1
    
	% mid nodes
    for ithCell = considerNbr(ithDist)+1:nCells-considerNbr(ithDist)

		firstEdge = ithCell - considerNbr(ithDist);
		secondEdge = ithCell + considerNbr(ithDist);

		smoothedElev(ithCell,ithDist) ...
			= mean(baseDEM(streamPath(firstEdge:secondEdge)));
	        
    end  

end % for ithDist = 1:nCNbr