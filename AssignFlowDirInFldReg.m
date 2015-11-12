function [m2SDSNbrY,m2SDSNbrX,mFlowDir_SubFldReg,mFlowDir_Saddle ...
    ,root,fldRegInfo] ...
    = AssignFlowDirInFldReg(m2SDSNbrY,m2SDSNbrX,subFldRegOutInfo,DEM ...
    ,slopeAllNbr,regionalMin,fldRegID,subFldRegID,sharedOutlet)
% @file AssignFlowDirInFldReg.m
% @brief Assign flow direction to the cells along the deepest route in
% flooded region to link an anther regional minima in next sub flooded
% region
%
% @param(in] m2SDSNbrY
% @param[in) m2SDSNbrX
% @param[in] SDSFlowDirection
% @param[in] subFldRegOutInfo
% @param[in] DEM
% @param[in] slopeAllNbr
% @param[in] regionalMin
% @param[in] fldRegID
% @param[in] subFldRegID
% @param[in] fldRegID
% @param[in] subFldRegID
% @param[in) flatRegionMap
% @retval m2SDSNbrY
% @retval m2SDSNbrX
% @retval flatRegMap
% @retval mFlowDir_SubFldReg
% @retval root: tree database of sub-flooded regions
% @retval fldRegInfo
%
% @version 0.2.1 / 2015-11-12
% @author Jongmin Byun
%==========================================================================

% constants
[mRows,nCols] = size(DEM);
ROOT_ID = 1;

% variables
root = tree('Root Node'); % init the tree db of sub-flooded regions
mFlowDir_SubFldReg = nan(mRows,nCols); % flow direction modified cell within a sub-flooded region
mFlowDir_Saddle = zeros(mRows,nCols); % flow direction modified cell on a saddle

% make flooded region info array: ID, outlet index,
% number of sub-flooded regions, outlet's elevation
uniqueFldRegID = unique(fldRegID); % unique flooded region ID
uniqueFldRegID(uniqueFldRegID <= 0) = [];
nFldReg = numel(uniqueFldRegID); % number of flooded region
fldRegInfo ... % info about the flooded region
    = zeros(nFldReg,4);

for i=1:nFldReg
    
    ithFldRegOutIdx = find(fldRegID == - uniqueFldRegID(i));

    if ~isnan(ithFldRegOutIdx)
        
        subFldRegInIthFldReg ...
            = unique(subFldRegID(fldRegID == uniqueFldRegID(i)));
        nSubFldRegInIthFldReg = numel(subFldRegInIthFldReg);
        fldRegInfo(i,1) = uniqueFldRegID(i);
        fldRegInfo(i,2) = ithFldRegOutIdx;
        fldRegInfo(i,3) = nSubFldRegInIthFldReg;
        fldRegInfo(i,4) = DEM(ithFldRegOutIdx);
        
    else
        
        error('Error, \nFlooded region info.');

    end
end

fldRegInfo = sortrows(fldRegInfo,4);

% main body
for i = 1:nFldReg
    
    fprintf( '%d/%d in AssignFlowInFldReg function\n ' ...
        ,i,nFldReg) % for debug
    
    % collect the info about the ith flooded region
    ithFldRegID = fldRegInfo(i,1); % ith flooded region ID
    ithFldRegOutIdx = fldRegInfo(i,2); % its outlet index
    [ithFldRegOutY,ithFldRegOutX] = ind2sub([mRows,nCols],ithFldRegOutIdx);
    nSubFldRegInIthFldReg = fldRegInfo(i,3);
    
    % init the tree DB of sub-flooded regions in the ith flooded region
    [root,ithFldRegTreeID] = root.addnode(ROOT_ID,ithFldRegID);

    % firstly, find a path from the true outlet to the minimum of the first
    % visit sub-flooded region
    
    % set variables for FindPathToMin function
    upStreamY = ithFldRegOutY;
    upStreamX = ithFldRegOutX;
    prevSubFldRegID = nan; % previously gone sub-flooded region ID
    targetSubFldRegID = nan; % target sub-flooded region ID

    goneSubFldRegID = []; % record for already gone sub-flooded region ID
    outConnectedSubFldRegID = []; % true outlet connected sub-flooded region ID
    nOutConnectedSubFldReg = sharedOutlet(upStreamY,upStreamX) + 1;
    
    for j = 1:nOutConnectedSubFldReg
        
        if j~= 1 % it means upStreamY,X is a shared true outlet
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
        outConnectedSubFldRegID = [outConnectedSubFldRegID;arrivedSubFldRegID];
        
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
                
                % if there is adjacent sub-flooded region except for true
                % outlet connected sub-flooded region, make a list to visit
                
                toGoSubFldRegID = []; % list of sub-flooded region
                nAdjSubFldReg = numel(adjSubFldRegRow);
                if nAdjSubFldReg > 0
                    
                    for k = 1:nAdjSubFldReg
                    
                        row = adjSubFldRegRow(k);
                        outletY = subFldRegOutInfo(row,4);
                        outletX = subFldRegOutInfo(row,5);
                        
                        if sharedOutlet(ithFldRegOutIdx) == 0

                            if ~ismember(adjSubFldRegID(k),goneSubFldRegID)
                                toGoSubFldRegID = [toGoSubFldRegID; ...
                                    outletY,outletX,adjSubFldRegID(k)];
                            end

                        else % sharedOutlet(ithFldRegOutIdx) ~= 0

                            if outletY == ithFldRegOutY && outletX == ithFldRegOutX

                                if ~ismember(adjSubFldRegID(k),outConnectedSubFldRegID)
                                    outConnectedSubFldRegID = [outConnectedSubFldRegID; ...
                                                                adjSubFldRegID(k)];
                                    nOutConnectedSubFldReg = nOutConnectedSubFldReg + 1;
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

                indicatorForParent = indicatorForParent + 1;
                if indicatorForParent > numel(goneSubFldRegID)
                    pathNotDone = false;
                end
            end % while pathNotDone
        end % if nSubFldRegInIthFldReg > 1
    end % for j = 1:nOutConnectedSubFldReg
end % for i = 1:nFldReg