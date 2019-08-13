function [m2SDSNbrY,m2SDSNbrX,mFlowDir_SubFldReg,mFlowDir_Saddle ...
    ,root,fldRegInfo] ...
    = AssignFlowDirInFldReg(m2SDSNbrY,m2SDSNbrX,subFldRegOutInfo,DEM ...
    ,slopeAllNbr,regionalMin,fldRegID,subFldRegID,sharedOutlet)
% @file AssignFlowDirInFldReg.m
% @brief Assign flow direction to the cells along the maximum-depth route
% in flooded region to link to another regional minima or the next
% sub-flooded region
%
% @param[in] m2SDSNbrY: modified flow direction along the maximum-depth
%                       path to each regional minima
% @param[in] m2SDSNbrX
% @param[in] SDSFlowDirection
% @param[in] subFldRegOutInfo
% @param[in] DEM
% @param[in] slopeAllNbr
% @param[in] regionalMin
% @param[in] fldRegID
% @param[in] subFldRegID
% @param[in] fldRegID
% @param[in] subFldRegID
% @param[in] flatRegionMap
% @retval m2SDSNbrY
% @retval m2SDSNbrX
% @retval flatRegMap
% @retval mFlowDir_SubFldReg
% @retval root: tree database of sub-flooded regions
% @retval fldRegInfo
%
% @version 0.3.1 / 2019-08-13
% @author Jongmin Byun
%==========================================================================

% Constants
[mRows,nCols] = size(DEM);
ROOT_ID = 1;

% Variables
root = tree('Root Node'); % init a tree DB for sub-flooded regions
mFlowDir_SubFldReg ...  % flow direction modified cells within sub-flooded
    = nan(mRows,nCols); % regions
mFlowDir_Saddle = zeros(mRows,nCols); % flow direction modified cells on saddles

% Make flooded region info array: ID, outlet index,
% number of sub-flooded regions, outlet's elevation

% Calculate the number of flooded region
uniqueFldRegID = unique(fldRegID); % unique flooded region ID
uniqueFldRegID(uniqueFldRegID <= 0) = [];
nFldReg = numel(uniqueFldRegID); % number of flooded region
fldRegInfo = zeros(nFldReg,4); % info about the flooded region    

for i=1:nFldReg
    
    ithFldRegOutIdx = find(fldRegID == - uniqueFldRegID(i));

    if ~isnan(ithFldRegOutIdx)
        
        fldRegInfo(i,1) = uniqueFldRegID(i); % ith flooded region ID
        fldRegInfo(i,2) = ithFldRegOutIdx; % its outlet index
        % number of the cells on the ith sub-flooded region
        subFldRegInIthFldReg ...
            = unique(subFldRegID(fldRegID == uniqueFldRegID(i)));
        nSubFldRegInIthFldReg = numel(subFldRegInIthFldReg);
        fldRegInfo(i,3) = nSubFldRegInIthFldReg;
        fldRegInfo(i,4) = DEM(ithFldRegOutIdx);
        
    else % isnan(ithFldRegOutIdx)
        
        fprintf('Error in making the flooded region information.\n');
        fprintf('Flooded region ID does not match with the ID of its outlet.\n');
        error('Stop the AssignFlowDirInFldReg function');        

    end
    
end % for i=1:nFldReg

fldRegInfo = sortrows(fldRegInfo,4);

% Main body
for i = 1:nFldReg
    
    fprintf( '%d/%d is processing in the AssignFlowInFldReg function\n ' ...
        ,i,nFldReg) % for debug
    
    % Collect the info about the ith flooded region
    ithFldRegID = fldRegInfo(i,1); % ith flooded region ID
    ithFldRegOutIdx = fldRegInfo(i,2); % its outlet index
    [ithFldRegOutY,ithFldRegOutX] = ind2sub([mRows,nCols],ithFldRegOutIdx);
    nSubFldRegInIthFldReg = fldRegInfo(i,3);
    
    % Init the tree DB of sub-flooded regions in the ith flooded region
    [root,ithFldRegTreeID] = root.addnode(ROOT_ID,ithFldRegID);

    % Firstly, find a path from the true outlet to the minimum elevation of
    % the first visit sub-flooded region
    
    % Set variables for the first operation of the FindPathToMin function
    upStreamY = ithFldRegOutY;
    upStreamX = ithFldRegOutX;
    prevSubFldRegID = nan; % previously gone sub-flooded region ID
    targetSubFldRegID = nan; % target sub-flooded region ID

    goneSubFldRegID = []; % record for the already gone sub-flooded region ID
    outConnectedSubFldRegID = []; % true outlet connected sub-flooded region ID
    % or ID of the sub-flooded region adjacent to the true outlet
    nOutConnectedSubFldReg = sharedOutlet(upStreamY,upStreamX) + 1;
    
    for j = 1:nOutConnectedSubFldReg
        
        if j~= 1 % It means upStreamY,X is a shared true outlet.
            
            % Therefore, set the parameter the same again
            upStreamY = ithFldRegOutY;
            upStreamX = ithFldRegOutX;
            prevSubFldRegID = nan;
            targetSubFldRegID = outConnectedSubFldRegID(j);
        end
            
        [m2SDSNbrY,m2SDSNbrX,arrivedRegMinY,arrivedRegMinX ...
            ,mFlowDir_SubFldReg,mFlowDir_Saddle] ...
            = FindPathToMin(mRows,nCols,upStreamY,upStreamX,slopeAllNbr ...
            ,m2SDSNbrY,m2SDSNbrX,regionalMin,fldRegID,subFldRegID ...
            ,ithFldRegID,targetSubFldRegID ...
            ,prevSubFldRegID,mFlowDir_SubFldReg,mFlowDir_Saddle ...
            ,ithFldRegOutIdx,sharedOutlet,outConnectedSubFldRegID);

        % record the arrived sub-flooded region ID at the tree DB
        arrivedSubFldRegID = subFldRegID(arrivedRegMinY,arrivedRegMinX);
        [root,childrenID] ...
                        = root.addnode(ithFldRegTreeID,arrivedSubFldRegID);
                 
        % record the arrived sub-flooded region as gone sub-flooded region
        goneSubFldRegID = [goneSubFldRegID; arrivedSubFldRegID];
        % record the arrived sub-flooded region as true outlet connected
        % flooded region
        if ~ismember(arrivedSubFldRegID,outConnectedSubFldRegID)
            outConnectedSubFldRegID ...
                = [outConnectedSubFldRegID;arrivedSubFldRegID];
        end
        
        % if the number of sub-flooded regions is more than one,
        % identify adjacent sub-flooded regions for the arrivied one
        if nSubFldRegInIthFldReg > 1            
            
            % variables for indicating parent sub-flooded region
            indicatorForParent = 1; % reset the varialble

            pathNotDone = true;
            while pathNotDone
                
                % identify adjacent sub-flooded regions
                adjSubFldRegRow_2nd = find(subFldRegOutInfo(:,2) ...
                                        == goneSubFldRegID(indicatorForParent));
                adjSubFldRegID_2nd = subFldRegOutInfo(adjSubFldRegRow_2nd,3);
                adjSubFldRegRow_3rd = find(subFldRegOutInfo(:,3)  ...
                                        == goneSubFldRegID(indicatorForParent));
                adjSubFldRegID_3rd = subFldRegOutInfo(adjSubFldRegRow_3rd,2);
                
                adjSubFldRegRow = [adjSubFldRegRow_2nd; adjSubFldRegRow_3rd];
                adjSubFldRegID = [adjSubFldRegID_2nd; adjSubFldRegID_3rd];
                
                % if there is adjacent sub-flooded region, make lists for
                % true outlet connected sub-flooded region and saddle
                % connected sub-flooded region
                
                toGoSubFldRegID = []; % list of sub-flooded region
                nAdjSubFldReg = numel(adjSubFldRegRow);
                if nAdjSubFldReg > 0
                    
                    % sort the adjacent sub-flooded regions according to
                    % the elevation of each outlet
                    % note that choose the sub-flooded region with the highest
                    % saddle to prevent discontinuity of flow accumulation from
                    % a node sub-flooded region to a parent sub-flooded region
                    if nAdjSubFldReg > 1
                        
                        orderK = zeros(nAdjSubFldReg,3);
                        for z = 1:nAdjSubFldReg

                            orderK(z,1) = adjSubFldRegRow(z);
                            orderK(z,2) = adjSubFldRegID(z);
                            
                            outletY = subFldRegOutInfo(orderK(z,1),4);
                            outletX = subFldRegOutInfo(orderK(z,1),5);
                            orderK(z,3) = DEM(outletY,outletX);

                        end

                        orderK = sortrows(orderK,-3);
                        adjSubFldRegRow = orderK(:,1);
                        adjSubFldRegID = orderK(:,2);
                        
                    end
                    
                    for k = 1:nAdjSubFldReg
                    
                        row = adjSubFldRegRow(k);
                        outletY = subFldRegOutInfo(row,4);
                        outletX = subFldRegOutInfo(row,5);
                        
                        if sharedOutlet(ithFldRegOutIdx) == 0
                            
                            % not for the shared true outlet
                            if ~ismember(adjSubFldRegID(k),goneSubFldRegID)
                                
                                % list for saddle connected sub-flooded region
                                toGoSubFldRegID = [toGoSubFldRegID; ...
                                    outletY,outletX,adjSubFldRegID(k)];
                            end

                        else % sharedOutlet(ithFldRegOutIdx) ~= 0
                            
                            % for the shared true outlet
                            if outletY == ithFldRegOutY && outletX == ithFldRegOutX

                                if ~ismember(adjSubFldRegID(k),outConnectedSubFldRegID)
                                    % list for true outlet connected
                                    % sub-flooded region
                                    outConnectedSubFldRegID ...
                                        = [outConnectedSubFldRegID;adjSubFldRegID(k)];
                                    nOutConnectedSubFldReg = nOutConnectedSubFldReg + 1; 
                                end
                                    
                            else % ~(outletY == ithFldRegOutY && outletX == ithFldRegOutX)
                                    
                                if ~ismember(adjSubFldRegID(k),goneSubFldRegID)
                                    % list for saddle connected sub-flooded region
                                    toGoSubFldRegID = [toGoSubFldRegID; ...
                                        outletY,outletX,adjSubFldRegID(k)];
                                end
                            end                            
                        end % if sharedOutlet(ithFldRegOutIdx) == 0
                        
                    end % for k = 1:nAdjSubFldReg
                end % if nAdjSubFldReg > 0

                % if the list is made, assign flow direction on each
                % adjacent sub-flooded region
                if ~isempty(toGoSubFldRegID)

                    % find a parent sub-flooded region ID for tree DB
                    nodeIdx = depthfirstiterator(root,ithFldRegTreeID);
                    nodeIdx(1) = [];
                    for m = nodeIdx
                        tmpSubFldRegID = root.get(m);
                        if tmpSubFldRegID == goneSubFldRegID(indicatorForParent);
                            parentSubFldRegID = m;
                        end
                    end
                    
                    nToGoSubFldReg = numel(toGoSubFldRegID(:,1));
                    for l = 1:nToGoSubFldReg

                        upStreamY = toGoSubFldRegID(l,1);
                        upStreamX = toGoSubFldRegID(l,2);
                        prevSubFldRegID = goneSubFldRegID(indicatorForParent);
                        targetSubFldRegID = toGoSubFldRegID(l,3);

                        [m2SDSNbrY,m2SDSNbrX,arrivedRegMinY,arrivedRegMinX ...
                            ,mFlowDir_SubFldReg,mFlowDir_Saddle] ...
                            = FindPathToMin(mRows,nCols,upStreamY,upStreamX,slopeAllNbr ...
                            ,m2SDSNbrY,m2SDSNbrX,regionalMin,fldRegID,subFldRegID ...
                            ,ithFldRegID,targetSubFldRegID ...
                            ,prevSubFldRegID,mFlowDir_SubFldReg,mFlowDir_Saddle ...
                            ,ithFldRegOutIdx,sharedOutlet,outConnectedSubFldRegID);

                        % record the arrived sub-flooded region ID at the tree DB
                        arrivedSubFldRegID = subFldRegID(arrivedRegMinY,arrivedRegMinX);
                        [root,childrenID] ...
                                        = root.addnode(parentSubFldRegID,arrivedSubFldRegID);

                        % record the arrived sub-flooded region as gone sub-flooded region
                        goneSubFldRegID = [goneSubFldRegID; arrivedSubFldRegID];

                    end % for l = 1=nToGoSubFldReg
                end % if nToGoSubFldReg > 0

                % check the remaining nodes in the list for saddle
                % connected sub-flooded region
                indicatorForParent = indicatorForParent + 1;
                if indicatorForParent > numel(goneSubFldRegID)
                    pathNotDone = false;
                end
                
            end % while pathNotDone
        end % if nSubFldRegInIthFldReg > 1
    end % for j = 1:nOutConnectedSubFldReg
end % for i = 1:nFldReg