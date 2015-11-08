function [m2SDSNbrY,m2SDSNbrX,mFlowDir_SubFldReg,mFlowDir_Saddle ...
    ,subFldRegTree,fldRegInfo] ...
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
% @retval subFldRegTree
% @retval fldRegInfo
%
% @version 0.1.3 / 2015-11-08
% @author Jongmin Byun
%==========================================================================

% constants
[mRows,nCols] = size(DEM);

% variables
subFldRegTree = tree('Root Node'); % init the tree db of sub-flooded regions
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
    [subFldRegTree,ithFldRegTreeID] = subFldRegTree.addnode(1,ithFldRegID);

    % find a path from the outlet to the first met regional minima
    [upStreamY,upStreamX] ... % initial upstream index
        = ind2sub([mRows,nCols],fldRegInfo(i,2));
        
    initUpStreamY = upStreamY;
    initUpStreamX = upStreamX;
    
    nToGoSubFldRegID = sharedOutlet(upStreamY,upStreamX) + 1;
    toGoSubFldRegID = [];
    goneConSubFldReg = []; % gone outlet connected sub-flooded region
    
    for t4 = 1:nToGoSubFldRegID
        
        % if shared outlet, reset initUpStreamY, X to revisit
        % ithFldRegOutIdx
        if sharedOutlet(initUpStreamY,initUpStreamX)
            upStreamY = initUpStreamY;
            upStreamX = initUpStreamX;
        end
        
        % set the variable for determining a target flooded region as NULL
        arrivedRegMinY = 0; % regional minima at which a path will be ended
        arrivedRegMinX = 0;
        if t4 == 1
            targetSubFldRegID = nan(1,1); % target sub-flooded region ID
        else
            targetSubFldRegID = toGoSubFldRegID(t4);
        end
        prevSubFldRegID = nan(1,1); % previously gone sub-flooded region ID
        
        [m2SDSNbrY,m2SDSNbrX,arrivedRegMinY,arrivedRegMinX ...
            ,mFlowDir_SubFldReg,mFlowDir_Saddle] ...
            = FindPathToMin(mRows,nCols,upStreamY,upStreamX,slopeAllNbr ...
            ,m2SDSNbrY,m2SDSNbrX,regionalMin,fldRegID,subFldRegID ...
            ,arrivedRegMinY,arrivedRegMinX,ithFldRegID,targetSubFldRegID ...
            ,prevSubFldRegID,mFlowDir_SubFldReg,mFlowDir_Saddle ...
            ,ithFldRegOutIdx,sharedOutlet,goneConSubFldReg);
        
        arrivedSubFldRegID = subFldRegID(arrivedRegMinY,arrivedRegMinX);
        prevSubFldRegID ... % define previous sub-flooded region ID
            = arrivedSubFldRegID;
        goneConSubFldReg = [goneConSubFldReg; arrivedSubFldRegID];
        if t4 == 1
            toGoSubFldRegID = prevSubFldRegID;
        end

        % check if there is no adjacent sub flooded region in ith flooded region
        if nSubFldRegInIthFldReg == 1
            % if no, end the loop
            [subFldRegTree,childrenID] ...
                = subFldRegTree.addnode(ithFldRegTreeID,prevSubFldRegID);
        else
            % if there is, identify adjacent sub-flooded regions and find paths
            % from each outlet to another regional minima
            nHead = 1; % reset the number for indicating parent one
            nSubHead = 1; % reset the number for indicating child one
            goneSubFldRegID = toGoSubFldRegID(t4); % gone sub-flooded region ID
            sortedSubFldRegOutYX = []; % outlets of sub-flooded regions
            
            nSubFldRegOut = numel(subFldRegOutInfo(:,1));
            isOutletConnectedSubFldReg = true;
            pathNotDone = true;
            while pathNotDone
                
                % identify adjacent sub-flooded regions using sub-flooded
                % region outlet information
                nAdjSubFldReg = 0; % number of adjacent sub-flooded regions
                for t1 = 1:nSubFldRegOut
                    
                    tIdx = find(subFldRegOutInfo(t1,2:3) ...
                                == goneSubFldRegID(nHead));
                    if numel(tIdx) > 0
                        if tIdx == 1 % 2nd column in subFldRegOutInfo
                            tAdjSubFldRegID = subFldRegOutInfo(t1,3);
                        elseif tIdx == 2 % 3rd column in subFldRegOutInfo
                            tAdjSubFldRegID = subFldRegOutInfo(t1,2);
                        end
                        
                        % remove repeated one except for the shared outlet
                        outletY = subFldRegOutInfo(t1,4);
                        outletX = subFldRegOutInfo(t1,5);
                        if sharedOutlet(initUpStreamY,initUpStreamX) == false
                            
                            if ~ismember(tAdjSubFldRegID,goneSubFldRegID) ...
                                    && tAdjSubFldRegID ~= 0
                                
                                goneSubFldRegID = [goneSubFldRegID;tAdjSubFldRegID];
                                nAdjSubFldReg = nAdjSubFldReg + 1;
                                sortedSubFldRegOutYX = [sortedSubFldRegOutYX ...
                                    ;subFldRegOutInfo(t1,4),subFldRegOutInfo(t1,5)];
                                
                            end
                        
                        elseif sharedOutlet(initUpStreamY,initUpStreamX) == true ...
                                && (outletY == initUpStreamY && outletX == initUpStreamX)
                            
                            if ~ismember(tAdjSubFldRegID,toGoSubFldRegID) ...
                                    && tAdjSubFldRegID ~= 0
                                toGoSubFldRegID = [toGoSubFldRegID;tAdjSubFldRegID];
                            end

                        end
                    end
                end
                
                % assign flow direction on each adjacent sub-flooded region
                for t2 = 1:nAdjSubFldReg
                    
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
                        ,ithFldRegOutIdx,sharedOutlet,goneConSubFldReg);
                
                    arrivedSubFldRegID = subFldRegID(arrivedRegMinY,arrivedRegMinX);
    
                    % record the immediately gone sub-flooded region using tree
                    % data structure
                    
                    % find a parent tree ID
                    parentID = 0;
                    tmpindices = depthfirstiterator(subFldRegTree,ithFldRegTreeID);
                    tmpindices(1) = [];
                    for t1 = tmpindices
                        tmpSubFldRegID = subFldRegTree.get(t1);
                        if tmpSubFldRegID == prevSubFldRegID
                            parentID = t1;
                        end
                        % check if parent node is the outlet of the ith
                        % flooded region
                        if parentID ~= 0 % parent is not the outlet!
                            
                            [subFldRegTree,childrenID] ...
                                = subFldRegTree.addnode(parentID,arrivedSubFldRegID);
                            
                        else % parentID == 0: parent is the outlet!
                            
                            % check if the ouelet is an shared outlet
                            if sharedOutlet(upStreamY,upStreamX) == false
                            
                                % record also the first met sub-flooded region
                                [subFldRegTree,childrenID] ...
                                    = subFldRegTree.addnode(ithFldRegTreeID,prevSubFldRegID);
                                [subFldRegTree,grandChildrenID] ...
                                    = subFldRegTree.addnode(childrenID,arrivedSubFldRegID);
                                
                            elseif sharedOutlet(upStreamY,upStreamX) == true ...
                                    && (initUpStreamY == upStreamY && initUpStreamX == upStreamX)
                                
                                % shared outlet!
                                if isOutletConnectedSubFldReg == true
                                    [subFldRegTree,childrenID] ...
                                        = subFldRegTree.addnode(ithFldRegTreeID,goneSubFldRegID(nHead));
                                    isOutletConnectedSubFldReg = false;
                                end
                                
                                [subFldRegTree,childrenID] ...
                                    = subFldRegTree.addnode(ithFldRegTreeID,arrivedSubFldRegID);
                                
                            end
                        end % if parentID -= 0
                    end % for t1 = tmpindices
                end % for t2 = 1:nAdjSubFldReg
                
                if sharedOutlet(upStreamY,upStreamX) == true ...
                            && (initUpStreamY == upStreamY && initUpStreamX == upStreamX) ...
                            && nAdjSubFldReg == 0
                        
                        [subFldRegTree,childrenID] ...
                            = subFldRegTree.addnode(ithFldRegTreeID,goneSubFldRegID(nHead));
                        
                end

                % relocate the position for head
                nHead = nHead + 1;
                if nHead > numel(goneSubFldRegID)
                    
                    pathNotDone = false;
                    
                else
                    
                    prevSubFldRegID = goneSubFldRegID(nHead);
                    
                end
                
            end % while pathNotDone
        end % nSubFldRegInIthFldReg == 1
    end % for t4 = 1:nToGoSubFldRegID
end % for i = 1