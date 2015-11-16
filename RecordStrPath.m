function [streamPath,distFromInit] = RecordStrPath(initY,initX ...
    ,endY,endX,mRows,nCols,SDSNbrY,SDSNbrX,dY,dX)
% @file RecordStrPath.m
% @brief Identify and record stream path using the input initial and end
% point coordinates. Meanwhile it calculates the distance from the
% initial point.
%
% @version 0.1.0. / 2015-11-12
% @author Jongmin Byun
%==========================================================================

initIdx = sub2ind([mRows,nCols],initY,initX);
endIdx = sub2ind([mRows,nCols],endY,endX);

% Identify the path determined by the new algorithm
streamPath = initIdx;
distFromInit = 0;
accDist = 0;
ithCellIdx = initIdx;
while ithCellIdx ~= endIdx
    
    nextCellIdx = sub2ind([mRows,nCols] ...
        ,SDSNbrY(ithCellIdx),SDSNbrX(ithCellIdx));
    
    if ismember(nextCellIdx,streamPath)
        % for debugging
        [tmpY2,tmpX2] = ind2sub([mRows,nCols],streamPath)
        [tmpY,tmpX] = ind2sub([mRows,nCols],nextCellIdx)
        error('Error: Infinite loop in RecordStrPath function.');
    end
    
    streamPath = [streamPath; nextCellIdx];
    
    % measure the distance between nodes  
    if abs(ithCellIdx - nextCellIdx) == mRows ...
            || abs(ithCellIdx - nextCellIdx) == 1
        
        accDist = accDist + dY;
        
    else
        
        accDist = accDist + sqrt(dY^2 + dX^2);
    
    end
    
    ithCellIdx = nextCellIdx;
    distFromInit = [distFromInit; accDist];
    
end