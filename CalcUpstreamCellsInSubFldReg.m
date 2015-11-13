function nUpstreamCells ...
    = CalcUpstreamCellsInSubFldReg(root,mFlowDir_SubFldReg ...
    ,mFlowDir_Saddle,fldRegID,nUpstreamCells,fldRegInfo ...
    ,m2SDSNbrY,m2SDSNbrX,steepestDescentSlope ...
    ,subFldRegID,DEM,regionalMin)
% @file CalcUpstreamCellsInSubFldReg.m
% @brief Calculate the number of upstream cells of cells in each flooded
% region
%
% @param[in] root
% @param[in] mFlowDirMap
% @param[in] fldRegID
% @param[in] nUpstreamCells
% @param[in] fldRegInfo
% @param[in] m2SDSNbrY
% @param[in] m2SDSNbrX
% @param[in] steepestDescentSlope
% @param[in] subFldRegID
% @param[in] DEM
% @param[in] regionalMin
% @param[in] flatRegMap
% @param[in] nFlatRegCells
% @retval nUpstreamCells
%
% @version 2.0.0 / 2015-11-13
%
% @author Jongmin Byun
%==========================================================================
% constant
NOT_MODIFIED = 0;
FROM_REGMIN_TO_UP = 2;
FROM_TARWSD_TO_TARWSD = 3;
GONE_WSD = 4;
[mRows,nCols] = size(DEM);

% variables
s = strel('square',3); % structural element when dilating image
treeDepth = root.depthtree;
nChildNode = tree(root,0); % number of children nodes
nChildNode2 = tree(root,0); % number of children nodes to be processed
nUpstreamCells(fldRegID > 0) = 0; % reset upstream cells numbers for flooded regions

prevNUpstreamCells = nUpstreamCells;
doNotAccFlowMask = ... % cells not to be processed for the first time in mCalcUpstreamCells function
    mFlowDir_SubFldReg == FROM_TARWSD_TO_TARWSD ...
    | mFlowDir_SubFldReg == FROM_REGMIN_TO_UP;
doNotAccFlowCell = nan(mRows,nCols);
doNotAccFlowCell(doNotAccFlowMask) = mFlowDir_SubFldReg(doNotAccFlowMask);

nFldReg = numel(fldRegInfo(:,1));
for i = 1:nFldReg
% for debug
    fprintf(' %d/%d in CalcUpstreamCellsinSubFldReg function\n',i,nFldReg);
    
    % A. choose a sub-flooded region according to the distance from ith
    % flooded region outlet
    
    % a. identify ith flooded region ID in tree data structure
    ithFldRegID = fldRegInfo(i,1); % ith flooded region ID
    nodeIdx = find(treeDepth == 1);
    for j = nodeIdx
        if root.get(j) == ithFldRegID
            ithFldRegTreeID = j;
        end
    end
    
    % b. identify all leaf nodes in the ith flooded region and calculate
    % how many children nodes each parent node has
    leafNodeID = [];
    nodeIdx = depthfirstiterator(root,ithFldRegTreeID);
    nodeIdx(1) = [];
    for j = nodeIdx
        if root.isleaf(j)
            leafNodeID = [leafNodeID,j];
        end
        
        childrenID = root.getchildren(j);
        nChildren = numel(childrenID);
        
        if nChildren >= 1
            nChildNode = nChildNode.set(j,nChildren);
            nChildNode2 = nChildNode2.set(j,nChildren);
        end
    end
    % B. calculate the number of upstream cells for the cells in ith
    % flooded region starting from the leaf node
    
    % iterate until each leaf node goes up and meet its parent node
    nLeafNode = numel(leafNodeID);
    for j = 1:nLeafNode
        childID = leafNodeID(j); % ith leaf node
        
        pathANotDone = true;
        while pathANotDone
            
            childSubFldRegID ... % current sub-flooded region ID
                = root.get(childID);
            iSubFldRegMap ... % current sub-flooded region
                = subFldRegID == childSubFldRegID;
            dIthSubFldRegMap ... % dilated current sub flooded region
                = imdilate(iSubFldRegMap,s);
            
            % choose the area for calculating the number of upstream cells
            % over the current sub-flooded region
            % Note: avoid flow direction modified cells, but include the
            % previously calculated sub-flooded region outlet
            targetDrainage = ~isnan(DEM) ...
                & dIthSubFldRegMap ...
                & (doNotAccFlowCell ~= FROM_REGMIN_TO_UP ...
                    & doNotAccFlowCell ~= FROM_TARWSD_TO_TARWSD ...
                    & (mFlowDir_Saddle == NOT_MODIFIED ...
                        | mFlowDir_Saddle == childSubFldRegID));
                    
            nUpstreamCells = prevNUpstreamCells;
            % calculate the number of upstream cells over chosen area
            [nUpstreamCells,doNotAccFlowCell] ...
                = mCalcUpstreamCells(nUpstreamCells,DEM ...
                ,targetDrainage,m2SDSNbrY,m2SDSNbrX ...
                ,steepestDescentSlope,iSubFldRegMap,doNotAccFlowCell);
            
            % calculate the number of upstream cells of the remainder:
            % from the regional minima to its outlet
            
            % Identify regional minima in the sub-flooded region
            iRegMin = regionalMin & iSubFldRegMap;
            iSubFldRegMinIdx = find(iRegMin);            
            
            [upStreamY,upStreamX] = ind2sub([mRows,nCols],iSubFldRegMinIdx);
            pathBNotDone = true;
            while pathBNotDone
            
                if doNotAccFlowCell(upStreamY,upStreamX) ...
                        ~= GONE_WSD % avoid the flat already calculated
                    nUpstreamCells(upStreamY,upStreamX) ...
                        = nUpstreamCells(upStreamY,upStreamX) + 1;
                end
                
                downStreamNbrY = m2SDSNbrY(upStreamY,upStreamX);
                downStreamNbrX = m2SDSNbrX(upStreamY,upStreamX);

                if iSubFldRegMap(downStreamNbrY,downStreamNbrX) == false

                    pathBNotDone = false;
                    doNotAccFlowCell(upStreamY,upStreamX) = GONE_WSD;
                    
                else
                    
                    if doNotAccFlowCell(upStreamY,upStreamX) ...
                            ~= GONE_WSD % avoid the flat already calculated
                        
                        doNotAccFlowCell(upStreamY,upStreamX) = GONE_WSD;
                        % transfer upstream cells number to the downstream cell
                        
                        nUpstreamCells(downStreamNbrY,downStreamNbrX) ...
                            = nUpstreamCells(downStreamNbrY,downStreamNbrX) ...
                            + nUpstreamCells(upStreamY,upStreamX);
                        
                    end
                    
                    upStreamY = downStreamNbrY;
                    upStreamX = downStreamNbrX;
                    
                end
            end
                        
            % clean the upstream cells number outside ith sub-flooded
            % region
            prevNUpstreamCells(iSubFldRegMap) = nUpstreamCells(iSubFldRegMap);
            
            % remove one node from the remaining children nodes
            parentID = nChildNode.getparent(childID);
            nChildNode2 ...
                = nChildNode2.set(parentID,nChildNode2.get(parentID) - 1);
            % calculate the number of sibling nodes of the parent node.
            tmpNoChildren = nChildNode.get(parentID);
            % if there are more than 2 sibling nodes (that is, parent node)
            % and calculation of all sibling nodes are not finished
            if tmpNoChildren > 1 && nChildNode2.get(parentID) > 0
                % move on to another leaf node
                pathANotDone = false;
            else
                childID = parentID;
            end

            % if tree depth of the watershed in the flooded region is one
            if treeDepth.get(childID) == 1
                pathANotDone = false;
            end
            
        end
    end
end