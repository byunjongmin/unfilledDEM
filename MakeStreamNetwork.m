function [streamNet,streamNetElement,streamNetNodeInfo] ...
    = MakeStreamNetwork(dX,dY,chosenWatershed,outletY,outletX ...
    ,SDSNbrY,SDSNbrX,floodedRegionIndex,floodedRegionCellsNo ...
    ,criticalRatio,upstreamCellsNo,flatRegionPathInfo,flatRegionMap)
%
% MakeStreamNetwork ver. 0.8
%     Make stream network on a DEM
%
% Outline
%    This function identifys every upstream cell begining with a main
%    outlet and record them in a Tree Data Structure (TDS) to make
%    a stream network using DEM.
%    During processing, in a junction cell in which has more than 2 upstream
%    cells which are equal to tributaries, it chooses one of them as first
%    upstream cell and find its more upstream cell. If it goes to the end
%    of tributary or terminal leaf node, it returns to the last visited
%    junction cell and chooses another upstream cell and find its more
%    upstream cell
%
% Note
%    2012/11/04: Includes all the cells within every flooded region
%    2012/09/19:I try to reduce the number of nodes in the DEM,
%    to choose the upstream cell linking main trunks.
%
%   Jong-Min Byun
%

% Constants
FIRST_MET = 1;
BOTH_FIRST_FINAL_MET = 3;
[mRows,nCols] = size(chosenWatershed);
% Offset for indicating 8 neighbours
ithOffset = [-1,mRows-1,mRows,mRows+1,1,-mRows+1, -mRows,-mRows-1];
% Tag describing the status of a node in a stream network
NOTDEFINED_NODE = 0;
MAINOUTLET_NODE = - 1; % root in a tree structure
JUNCTION_NODE = 2; % parent node in a tree structure
CHILD_NODE = 3; % child node in a tree structure
HEAD_NODE = 4; % terminal leaf node in a tree struectur

% Variables
% Linear indicies for downstream neighbour cell
% Note: Place NaN outside the chosen watershed
SDSNbrIdx = (SDSNbrX - 1) * mRows +  SDSNbrY;
SDSNbrIdx(~chosenWatershed) = nan;
% Variable describing the status of a cell in a stream network
streamNetElement = zeros(mRows,nCols);
% Additional property (ID of each node in a TDS, distance btw upstream node
% and current node) of a stream network
streamNetNodeInfo = zeros(mRows*nCols,3);
% Number of upstream cells at current cell 
nUpstreamNbr = nan(mRows,nCols);
% Number of times having gone upstream cell
% Note: In a junction cell, it identifies the number of upstream cells and
% chooses the cell one by one. So it need to remember the times having gone
% upstream cell.
nGoneUpstreamNbr = nan(mRows,nCols);
% List of indicies for returning to the last visited junction
% Note: This function record each junction index as a stack data structure
% (first in, last out). So if it goes to the terminal leaf node, it will
% return to the last visited junction node
idxListForReturnToJunction = zeros(mRows*nCols,1);

% Main body

% Initialization before iteration

% Insert the main outlet index to the current cell index 
outletIdx = (outletX - 1) * mRows + outletY;
currentCellIdx = outletIdx;
% TDS to make a stream network
streamNet = tree(currentCellIdx);
% Initialize the properties of the main outlet
streamNetNodeInfo(1,1) = 1; % current cell ID in a TSD (root node)
streamNetNodeInfo(1,2) = currentCellIdx; % current cell (main outlet) index

% Initialize the number of times having gone upstream cell
% Note: At the first time, it assumes that there is one upstream cell of
% the main outlet. When there are more than 2 upstream cells, it will be
% modified
nGoneUpstreamNbr(currentCellIdx) = 1;

% Initialize the list of indicies for returning to the last visited
% junction
ithParentNode = 1;
idxListForReturnToJunction(ithParentNode) =  currentCellIdx;

% Iterate until the cell for returning to the last visited junction will be
% the main outlet
do = true;
while do

    % for debug
    disp(currentCellIdx);        

    % A. Check whether the current cell is within the flat area
    
    % If it is not, check the flow direction of 8 neighbours and identify
    % the upstream cells and their number
    if flatRegionMap(currentCellIdx) == false
       
        % a. Check the flow direction of 8 neighbours at current cell and
        % identify the upstream cells and their number
        
        % Initialize the number of upstream cells
        nUpstreamNbr(currentCellIdx) = 0;
        % Initialize a variable for the information (idx,direction) of
        % upstream nbr cell
        upStreamNbrIdx = zeros(8,2);
        
        % 상류 셀들의 총 upstream 개수를 파악함
        totalUpstreamCellsNo = 0;
        for ithNbr = 1:8     

            % Check whether the flow direction of ith Nbr match with
            % current cell
            ithNbrIdx = currentCellIdx + ithOffset(ithNbr);
            if SDSNbrIdx(ithNbrIdx) == currentCellIdx
                
                    totalUpstreamCellsNo ...
                        = totalUpstreamCellsNo + upstreamCellsNo(ithNbrIdx);

            end
        end    
        
        for ithNbr = 1:8     

            % Check whether the flow direction of ith Nbr match with
            % current cell
            ithNbrIdx = currentCellIdx + ithOffset(ithNbr);
            if SDSNbrIdx(ithNbrIdx) == currentCellIdx

                if upstreamCellsNo(ithNbrIdx) ...
                        > (totalUpstreamCellsNo * criticalRatio)
                
                    % 상류 이웃 셀 개수를 1 증가 시키고, 이의 색인과 방향을
                    % 기록함
                    nUpstreamNbr(currentCellIdx) ...
                        = nUpstreamNbr(currentCellIdx) + 1;
                    upStreamNbrIdx(nUpstreamNbr(currentCellIdx),1) = ithNbrIdx;
                    upStreamNbrIdx(nUpstreamNbr(currentCellIdx),2) = ithNbr;

                end
            end
        end
    
    % flat region일 경우, flat region의 outlet인지 확인하고 만약 outlet이라면
    % flat region을 대표하는 juction 셀로 정의하고 flat region 주변 셀들 중에
    % upstream cell을 찾음. 하지만 outlet이 아닐 경우, flat region에 있는 outlet
    % 을 찾아 이를 upstream cell로 등록함. 만약 outlet마저 없을 경우는, 주변
    % 이웃 셀 중에 upstream cell을 찾음
    else 
        
        % flat region information
        flatRegID = flatRegionMap(currentCellIdx);
        flatRegIdx = find(flatRegionMap == flatRegID);
        
        if flatRegionPathInfo(currentCellIdx) == FIRST_MET ...
            || flatRegionPathInfo(currentCellIdx) == BOTH_FIRST_FINAL_MET
        
            % 이 경우 current cell은 flat region의 junction cell이 됨
            flatJunctionCellIdx = currentCellIdx;
            nFlatRegIdx = numel(flatRegIdx);

            % Initialize the number of upstream cells
            nUpstreamNbr(flatJunctionCellIdx) = 0;
            % Initialize a variable for the information (idx,direction) of
            % upstream nbr cell
            upStreamNbrIdx = [];
            
            % 상류 셀들의 총 upstream 개수를 파악함
            totalUpstreamCellsNo = 0;
            for ithFlatRegCell = 1:nFlatRegIdx
                
                tmpCurrentCellIdx = flatRegIdx(ithFlatRegCell);

                for ithNbr = 1:8     

                    % Check whether the flow direction of ith Nbr match with
                    % current cell
                    ithNbrIdx = tmpCurrentCellIdx + ithOffset(ithNbr);
                    if flatRegionMap(ithNbrIdx) ~= flatRegID
                        
                        if SDSNbrIdx(ithNbrIdx) == tmpCurrentCellIdx

                            totalUpstreamCellsNo ...
                                = totalUpstreamCellsNo + upstreamCellsNo(ithNbrIdx);

                        end
                    end
                end
            end
            
            for ithFlatRegCell = 1:nFlatRegIdx
                
                tmpCurrentCellIdx = flatRegIdx(ithFlatRegCell);

                for ithNbr = 1:8     

                    % Check whether the flow direction of ith Nbr match with
                    % current cell
                    ithNbrIdx = tmpCurrentCellIdx + ithOffset(ithNbr);
                    if flatRegionMap(ithNbrIdx) ~= flatRegID
                        
                        if SDSNbrIdx(ithNbrIdx) == tmpCurrentCellIdx

                            if upstreamCellsNo(ithNbrIdx) ...
                                    > (totalUpstreamCellsNo * criticalRatio)

                                % 상류 이웃 셀 개수를 1 증가 시키고, 이의 색인과
                                % 방향을 기록함
                                nUpstreamNbr(flatJunctionCellIdx) ...
                                    = nUpstreamNbr(flatJunctionCellIdx) + 1;
                                upStreamNbrIdx = [upStreamNbrIdx; ithNbrIdx,ithNbr];

                            end                            
                        end
                    end
                end       
            end
                        
        else
            
            flatRegOutletIdx ...
                = find( (flatRegionMap == flatRegID) ...
                            &  (flatRegionPathInfo == FIRST_MET ...
                                    | flatRegionPathInfo == BOTH_FIRST_FINAL_MET) );

            if numel(flatRegOutletIdx) > 0 ...
                    && SDSNbrIdx(currentCellIdx) == flatRegOutletIdx
                
                % Initialize the number of upstream cells
                nUpstreamNbr(currentCellIdx) = 0;
                % Initialize a variable for the information (idx,direction) of
                % upstream nbr cell
                upStreamNbrIdx = zeros(8,2);    

                % 상류 이웃 셀 개수를 1 증가 시키고, 이의 색인과 방향을
                % 기록함
                nUpstreamNbr(currentCellIdx) ...
                    = nUpstreamNbr(currentCellIdx) + 1;
                upStreamNbrIdx(nUpstreamNbr(currentCellIdx),1)...
                    = flatRegOutletIdx;
                
            else
                
                % Check the flow direction of 8 neighbours at current cell and
                % identify the upstream cells and their number

                % Initialize the number of upstream cells
                nUpstreamNbr(currentCellIdx) = 0;
                % Initialize a variable for the information (idx,direction) of
                % upstream nbr cell
                upStreamNbrIdx = zeros(8,2);
                
                % 상류 셀들의 총 upstream 개수를 파악함
                totalUpstreamCellsNo = 0;
                for ithNbr = 1:8     

                    % Check whether the flow direction of ith Nbr match with
                    % current cell
                    ithNbrIdx = currentCellIdx + ithOffset(ithNbr);
                    if SDSNbrIdx(ithNbrIdx) == currentCellIdx

                            totalUpstreamCellsNo ...
                                = totalUpstreamCellsNo + upstreamCellsNo(ithNbrIdx);

                    end
                end    

                for ithNbr = 1:8     

                    % Check whether the flow direction of ith Nbr match with
                    % current cell
                    if SDSNbrIdx(currentCellIdx + ithOffset(ithNbr)) ...
                            == currentCellIdx
                        
                        if upstreamCellsNo(currentCellIdx + ithOffset(ithNbr)) ...
                                > (totalUpstreamCellsNo * criticalRatio)

                            % 상류 이웃 셀 개수를 1 증가 시키고, 이의 색인과 방향을
                            % 기록함
                            nUpstreamNbr(currentCellIdx) ...
                                = nUpstreamNbr(currentCellIdx) + 1;
                            upStreamNbrIdx(nUpstreamNbr(currentCellIdx),1)...
                                = currentCellIdx + ithOffset(ithNbr);
                            upStreamNbrIdx(nUpstreamNbr(currentCellIdx),2) = ithNbr;

                        end
                    end
                end                
            end          
        end        
    end

    % b. Check whether there are upstream cells

    % 1) Check how many upstream cells there are?
    if nUpstreamNbr(currentCellIdx) > 0

        % If there is or are, do not hesitate to choose one (of them) as a
        % upstream cell
        % and find its more upstream cell

        % (1) Check whether there is one upstream cell
        if nUpstreamNbr(currentCellIdx) == 1

            % (a) If the number of upstream cell is one, find its more
            % upstream cell and record it as a child node

            % Change the current cell into the downstream cell which is
            % a parent node in a TDS. And change the newly founded upstream
            % cell into the upstream cell which is a child node in a TDS.
            downStreamCellIdx = currentCellIdx;
            upStreamCellIdx ...
                = upStreamNbrIdx(nUpstreamNbr(currentCellIdx),1);

            % 주의: 중복을 막기 위해
            if streamNetElement(upStreamCellIdx) == NOTDEFINED_NODE
                
                % Record the current cell as CHILD NODE
                streamNetElement(downStreamCellIdx) = CHILD_NODE;

                % Add the upstream cell as the child node under the downstream
                % cell in a TSD                
                [downStreamNodeID,tmp] ...
                    = find(streamNetNodeInfo(:,2) == downStreamCellIdx);
                [streamNet,nodeID] ...
                    = streamNet.addnode(downStreamNodeID,upStreamCellIdx);

                % Write additional property
                streamNetNodeInfo(nodeID,1) = nodeID;
                streamNetNodeInfo(nodeID,2) = upStreamCellIdx;

                % Calculate the distance btw downstream and upstream cell
                [upStreamY,upStreamX] = ind2sub([mRows,nCols],upStreamCellIdx);                
                [downStreamY,downStreamX] ...
                    = ind2sub([mRows,nCols],downStreamCellIdx);
                streamNetNodeInfo(nodeID,3) ...
                    = sqrt(double(((downStreamY-upStreamY)*dY)^2 ...
                    +((downStreamX-upStreamX)*dX)^2));

                % Change the current cell index into the upstreal cell index
                currentCellIdx = upStreamCellIdx;

            else
                
                % Record the current cell as the HEAD_NODE
                streamNetElement(currentCellIdx) = HEAD_NODE;
                
                % Chane the current cell index into the last visited
                % junction cell
                currentCellIdx = idxListForReturnToJunction(ithParentNode);

                % If the last visited junction cell is the main outlet,
                % end iteration
                % * 주의: outlet의 상류 이웃 셀이 두 개 이상일 경우,
                % 모든 상류 이웃 셀에 대한 탐색이 완료되어야 함.
                if currentCellIdx == outletIdx ...
                    && nGoneUpstreamNbr(currentCellIdx) ...
                        >= nUpstreamNbr(currentCellIdx)

                    do = false;

                end
                
            end
        % If the number of upstream cells is more than 2, current cell is a
        % junction node. So records the current cell as a junction node and
        % idenfys the number of upstream cells and puts them int the list
        % of the cells for returning to the last visited junction and
        % chooses one of them as upstream cell and find its more upstream
        % cell and record it as a child
        else

            % (a) Check whether this is first visit
            % Note: When the number of upstream cell is more than one, the
            % belows will be done repeatedly. So if this is first time, do
            % the followings
            if streamNetElement(currentCellIdx) == NOTDEFINED_NODE          

                % Record current cell as the JUNCTION_NODE
                streamNetElement(currentCellIdx) = JUNCTION_NODE;

                % Put the current cell index to the list of junction
                % Note: If the main outlet cell have more than 2
                % upstream cells, then add once more the index of the main
                % outlet to the list for juntion cells 
                if currentCellIdx == outletIdx

                    ithParentNode = ithParentNode + 1;
                    idxListForReturnToJunction(ithParentNode) ...
                        = currentCellIdx;

                end

                ithParentNode = ithParentNode + 1;
                idxListForReturnToJunction(ithParentNode) = currentCellIdx;

                % Initialize the number of time having gone upstream
                nGoneUpstreamNbr(currentCellIdx) = 0;

            end

            % (b) Check whether upstream cell to go remains
            if nGoneUpstreamNbr(currentCellIdx) ...
                    < nUpstreamNbr(currentCellIdx)

                % If upstream cells to go remains, puts it in the
                % stream network

                % Change the current cell into the downstream cell which is
                % a parent node in a TDS. And change the newly founded
                % upstream cell into the upstream cell which is a child
                % node in a TDS.
                % Note: MATLAB 색인 규약 때문에 '+1'
                downStreamCellIdx = currentCellIdx;
                upStreamCellIdx ...
                = upStreamNbrIdx(nGoneUpstreamNbr(currentCellIdx) + 1,1);

                % * 주의: 변경 후, 갔다 온 상류 셀 개수 1 증가 시킴
                nGoneUpstreamNbr(downStreamCellIdx) ...
                    = nGoneUpstreamNbr(downStreamCellIdx) + 1;

                % 현재 셀 색인을 하류 셀 아래 node로 추가함                
                [downStreamNodeID,tmp] ...
                    = find(streamNetNodeInfo(:,2) == downStreamCellIdx);
                [streamNet,nodeID] ...
                    = streamNet.addnode(downStreamNodeID,upStreamCellIdx);

                % Write additional property of upstream cell
                streamNetNodeInfo(nodeID,1) = nodeID;
                streamNetNodeInfo(nodeID,2) = upStreamCellIdx;

                % Calculate the distance btw downstream and upstream cell
                [upStreamY,upStreamX] ...
                    = ind2sub([mRows,nCols],upStreamCellIdx);                
                [downStreamY,downStreamX] ...
                    = ind2sub([mRows,nCols],downStreamCellIdx);
                streamNetNodeInfo(nodeID,3) ...
                    = sqrt(double(((downStreamY-upStreamY)*dY)^2 ...
                    +((downStreamX-upStreamX)*dX)^2));

                % Change the current cell index into the upstreal cell
                % index
                currentCellIdx = upStreamCellIdx;
                    
            % (c) If there is no left upstream cell, it returns to the last
            % visited junction cell
            else

                % Chane the current cell index into the last visited
                % junction cell
                ithParentNode = ithParentNode - 1;
                currentCellIdx = idxListForReturnToJunction(ithParentNode);

                % If the last visited junction cell is the main outlet,
                % end iteration
                % * 주의: outlet의 상류 이웃 셀이 두 개 이상일 경우,
                % 모든 상류 이웃 셀에 대한 탐색이 완료되어야 함.
                if currentCellIdx == outletIdx ...
                    && nGoneUpstreamNbr(currentCellIdx) ...
                        >= nUpstreamNbr(currentCellIdx)

                    do = false;

                end
            end %   if nGoneUpstreamNbr(currentCellIdx) 
        end % if nUpstreamNbr == 1

    % C. If there is no upstream cell, it is a terminal leaf node. So it
    % returns to the last visited junction cell
    else % elseif nUpstreamNbr <= 0

        % Record the current cell as the HEAD_NODE
        streamNetElement(currentCellIdx) = HEAD_NODE;

        % Change the current cell into last visited junction cell
        currentCellIdx = idxListForReturnToJunction(ithParentNode);

        % If the last visited junction cell is the main outlet,
        % end iteration
        if currentCellIdx == outletIdx ...
            && nGoneUpstreamNbr(currentCellIdx) ...
                == nUpstreamNbr(currentCellIdx)

            do = false;

        end
    end % if nUpstreamNbr > 0    
end % while do

% 유출구 색인을 MAINOUTLET_NODE 라고 기록함
streamNetElement(outletIdx) = MAINOUTLET_NODE;
