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
% @version 0.1.4 / 2015-11-10
% @author Jongmin Byun
%==========================================================================

% constants
[mRows,nCols] = size(DEM);

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
    
    ithFldRegID = fldRegInfo(i,1); % ith flooded region ID
    ithFldRegOutIdx = fldRegInfo(i,2); % its outlet index
    
    nSubFldRegInIthFldReg ... % number of sub-flooded region in ith flooded region
        = fldRegInfo(i,3);
    
    % init the tree of sub-flooded region in the ith flooded region
    [root,ithFldRegTreeID] = root.addnode(1,ithFldRegID);

    % find a path from the outlet to the first met regional minima
    [upStreamY,upStreamX] ... % initial upstream index
        = ind2sub([mRows,nCols],ithFldRegOutIdx);
        
    ithFldRegOutY = upStreamY;
    ithFldRegOutX = upStreamX;
    
    % for a shared true outlet
    nToGoSubFldRegID = sharedOutlet(upStreamY,upStreamX) + 1;
    toGoSubFldRegID = [];
    goneConnectedSubFldReg = []; % gone outlet connected sub-flooded region
    
    for j = 1:nToGoSubFldRegID
        
        % if shared outlet, reset upStreamY, X to revisit ithFldRegOutIdx
        if sharedOutlet(ithFldRegOutY,ithFldRegOutX) == true
            upStreamY = ithFldRegOutY;
            upStreamX = ithFldRegOutX;
        end
        
        % set the variable for determining a target flooded region as NULL
        arrivedRegMinY = 0; % regional minima at which a path will be ended
        arrivedRegMinX = 0;
        if j == 1
            targetSubFldRegID = nan(1,1); % target sub-flooded region ID
        else
            targetSubFldRegID = toGoSubFldRegID(j);
        end
        prevSubFldRegID = nan(1,1); % previously gone sub-flooded region ID
        
        [m2SDSNbrY,m2SDSNbrX,arrivedRegMinY,arrivedRegMinX ...
            ,mFlowDir_SubFldReg,mFlowDir_Saddle] ...
            = FindPathToMin(mRows,nCols,upStreamY,upStreamX,slopeAllNbr ...
            ,m2SDSNbrY,m2SDSNbrX,regionalMin,fldRegID,subFldRegID ...
            ,arrivedRegMinY,arrivedRegMinX,ithFldRegID,targetSubFldRegID ...
            ,prevSubFldRegID,mFlowDir_SubFldReg,mFlowDir_Saddle ...
            ,ithFldRegOutIdx,sharedOutlet,goneConnectedSubFldReg);
        
        arrivedSubFldRegID = subFldRegID(arrivedRegMinY,arrivedRegMinX);
        prevSubFldRegID ... % define previously gone sub-flooded region ID
            = arrivedSubFldRegID;
        goneConnectedSubFldReg = [goneConnectedSubFldReg; arrivedSubFldRegID];
        if j == 1
            % if first visit, fill toGoSubFldRegID
            toGoSubFldRegID = arrivedSubFldRegID;
        end

        % check if there is no adjacent sub flooded region in ith flooded region
        if nSubFldRegInIthFldReg == 1
            % if no, end the loop
            [root,childrenID] ...
                = root.addnode(ithFldRegTreeID,arrivedSubFldRegID);
        else
            % if there is, identify adjacent sub-flooded regions and find paths
            % from each outlet to another regional minima
            nHead = 1; % reset the number for indicating a parent one
            nSubHead = 1; % reset the number for indicating a child one
            goneSubFldRegID = toGoSubFldRegID(j); % gone sub-flooded region ID
            sortedSubFldRegOutYX = []; % outlets of sub-flooded regions
            
            % identify adjacent sub-flooded regions
            nSubFldRegOut = numel(subFldRegOutInfo(:,1));
            isOutletConnectedSubFldReg = true;
            pathNotDone = true;
            while pathNotDone
                
                % identify adjacent sub-flooded regions using sub-flooded
                % region outlet information
                nAdjSubFldReg = 0; % number of adjacent sub-flooded regions
                for k = 1:nSubFldRegOut
                    
                    tIdx = find(subFldRegOutInfo(k,2:3) ...
                                == goneSubFldRegID(nHead));
                    if numel(tIdx) > 0
                        if tIdx == 1 % 2nd column in subFldRegOutInfo
                            tAdjSubFldRegID = subFldRegOutInfo(k,3);
                        elseif tIdx == 2 % 3rd column in subFldRegOutInfo
                            tAdjSubFldRegID = subFldRegOutInfo(k,2);
                        end
                        
                        % remove repeated one except for the shared outlet
                        outletY = subFldRegOutInfo(k,4);
                        outletX = subFldRegOutInfo(k,5);
                        if sharedOutlet(ithFldRegOutIdx) == false
                            
                            if ~ismember(tAdjSubFldRegID,goneSubFldRegID) ...
                                    && tAdjSubFldRegID ~= 0
                                
                                goneSubFldRegID = [goneSubFldRegID;tAdjSubFldRegID];
                                nAdjSubFldReg = nAdjSubFldReg + 1;
                                sortedSubFldRegOutYX = [sortedSubFldRegOutYX ...
                                    ;subFldRegOutInfo(k,4),subFldRegOutInfo(k,5)];
                                
                            end
                        
                        else % sharedOutlet(ithFldRegOutIdx) == true
                                
                            if (outletY == ithFldRegOutY && outletX == ithFldRegOutX)
                            
                                if ~ismember(tAdjSubFldRegID,toGoSubFldRegID) ...
                                        && tAdjSubFldRegID ~= 0
                                    toGoSubFldRegID = [toGoSubFldRegID;tAdjSubFldRegID];
                                end
                            end
                            
                        end % if sharedOutlet(ithFldRegOutIdx) 
                    end % if numel(tIdx)
                end % for k = 1:
                
                % assign flow direction on each adjacent sub-flooded region
                for l = 1:nAdjSubFldReg
                    
                    % coordinate for the outlet of sub-flooded region
                    upStreamY = sortedSubFldRegOutYX(nSubHead,1);
                    upStreamX = sortedSubFldRegOutYX(nSubHead,2);
                    
                    % relocate the position for target sub-flooded region
                    nSubHead = nSubHead + 1;
                    targetSubFldRegID = goneSubFldRegID(nSubHead);
                    
                    arrivedRegMinY = 0;
                    arrivedRegMinX = 0;
                    
                    [m2SDSNbrY,m2SDSNbrX,arrivedRegMinY,arrivedRegMinX ...
                        ,mFlowDir_SubFldReg,mFlowDir_Saddle] ...
                        = FindPathToMin(mRows,nCols,upStreamY,upStreamX,slopeAllNbr ...
                        ,m2SDSNbrY,m2SDSNbrX,regionalMin,fldRegID,subFldRegID ...
                        ,arrivedRegMinY,arrivedRegMinX,ithFldRegID,targetSubFldRegID ...
                        ,prevSubFldRegID,mFlowDir_SubFldReg,mFlowDir_Saddle ...
                        ,ithFldRegOutIdx,sharedOutlet,goneConnectedSubFldReg);
                
                    arrivedSubFldRegID = subFldRegID(arrivedRegMinY,arrivedRegMinX);
    
                    % record the immediately gone sub-flooded region using tree
                    % data structure
                    
                    % find a parent tree ID
                    parentID = 0;
                    tmpindices = depthfirstiterator(root,ithFldRegTreeID);
                    tmpindices(1) = [];
                    for m = tmpindices
                        tmpSubFldRegID = root.get(m);
                        if tmpSubFldRegID == prevSubFldRegID
                            parentID = m;
                        end
                    end
                    
                    % check if parent node is the outlet of the ith
                    % flooded region
                    if parentID ~= 0 % parent is not the outlet!

                        [root,childrenID] ...
                            = root.addnode(parentID,arrivedSubFldRegID);

                    else % parentID == 0: parent is the outlet!

                        % check if the ouelet is an shared outlet
                        if sharedOutlet(upStreamY,upStreamX) == false

                            % record also the first met sub-flooded region
                            [root,childrenID] ...
                                = root.addnode(ithFldRegTreeID,prevSubFldRegID);
                            [root,grandChildrenID] ...
                                = root.addnode(childrenID,arrivedSubFldRegID);

                        elseif sharedOutlet(upStreamY,upStreamX) == true ...
                                && (ithFldRegOutY == upStreamY && ithFldRegOutX == upStreamX)

                            % shared outlet!
                            if isOutletConnectedSubFldReg == true
                                [root,childrenID] ...
                                    = root.addnode(ithFldRegTreeID,goneSubFldRegID(nHead));
                                isOutletConnectedSubFldReg = false;
                            end

                            [root,childrenID] ...
                                = root.addnode(ithFldRegTreeID,arrivedSubFldRegID);

                        end
                    end % if parentID -= 0

                end % for l = 1:nAdjSubFldReg
                
                if sharedOutlet(upStreamY,upStreamX) == true ...
                            && (ithFldRegOutY == upStreamY && ithFldRegOutX == upStreamX) ...
                            && nAdjSubFldReg == 0
                        
                        [root,childrenID] ...
                            = root.addnode(ithFldRegTreeID,goneSubFldRegID(nHead));
                        
                end

                % relocate the position for head
                nHead = nHead + 1;
                if nHead > numel(goneSubFldRegID)
                    
                    pathNotDone = false;
                    
                else
                    
                    prevSubFldRegID = goneSubFldRegID(nHead);
                    
                end
                
            end % while pathNotDone
        end % if nSubFldRegInIthFldReg == 1
    end % for j = 1:nToGoSubFldRegID
end % for i = 1:nFldReg