function nUpstreamCells ...
    = CalcUpstreamCellsInSubFldReg(subFldRegTree,mFlowDir_SubFldReg ...
    ,mFlowDir_Saddle,fldRegID,nUpstreamCells,fldRegInfo ...
    ,m2SDSNbrY,m2SDSNbrX,steepestDescentSlope ...
    ,subFldRegID,DEM,regionalMin,flatRegMap,nFlatRegcells)
% @file CalcUpstreamCellsInSubFldReg.m
% @brief Calculate the number of upstream cells of cells in each flooded
% region
%
% @param[in] subFldRegTree
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
% @version 2.0
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
treeDepth = subFldRegTree.depthtree;
nChildNode = tree(subFldRegTree,0); % number of children nodes
nChildNode2 ... % number of children nodes to be processed
    = tree(subFldRegTree,0);

nUpstreamCells(fldRegID > 0) = 0; % clean upstream cells numbers of all flooded regions
prevNUpstreamCells = nUpstreamCells;

s = strel('square',3); % structural element when dilating image

doNotAccFlowMask = ... % cells not to be processed for the first time in mCalcUpstreamCells function
    mFlowDir_SubFldReg == FROM_TARWSD_TO_TARWSD ...
    | mFlowDir_SubFldReg == FROM_REGMIN_TO_UP;

doNotAccFlowCell = nan(mRows,nCols);
doNotAccFlowCell(doNotAccFlowMask) = mFlowDir_SubFldReg(doNotAccFlowMask);

flatIdx = {}; % indices for flat to be processed after all procedures
nFldReg = numel(fldRegInfo(:,1));
for iFldReg = 1:nFldReg
% for debug
    fprintf(' %d/%d in CalcUpstreamCellsinSubFldReg function\n',iFldReg,nFldReg);
    
    % A. choose a sub-flooded region according to the distance from ith
    % flooded region outlet
    
    % a. identify ith flooded region ID in tree data structure
    ithFldRegID = fldRegInfo(iFldReg,1); % ith flooded region ID
    tmpIndices = find(treeDepth == 1);
    for t1 = tmpIndices
        if subFldRegTree.get(t1) == ithFldRegID
            ithFldRegTreeID = t1;
        end
    end
    
    % b. idenity all leaf nodes in the ith flooded region and calculate how
    % many children nodes each parent node has
    leafNodeID = [];
    tmpIndices = depthfirstiterator(subFldRegTree,ithFldRegTreeID);
    tmpIndices(1) = [];
    for t1 = tmpIndices
        if subFldRegTree.isleaf(t1)
            leafNodeID = [leafNodeID,t1];
        end
        
        childrenID = subFldRegTree.getchildren(t1);
        nChildren = numel(childrenID);
        
        if nChildren >= 1
            nChildNode = nChildNode.set(t1,nChildren);
            nChildNode2 = nChildNode2.set(t1,nChildren);
        end
    end
    % B. calculate the number of upstream cells of cells in ith flooded
    % region starting from the leaf node
    
    % iterate until each leaf node goes up and meet its parent node
    nLeafNode = numel(leafNodeID);
    for t1 = 1:nLeafNode
        childID = leafNodeID(t1); % ith leaf node
        
        pathANotDone = true;
        while pathANotDone
            
            childSubFldRegID ... % current sub-flooded region ID
                = subFldRegTree.get(childID);
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
            [nUpstreamCells,doNotAccFlowCell] = mCalcUpstreamCells(nUpstreamCells,DEM ...
                ,targetDrainage,m2SDSNbrY,m2SDSNbrX ...
                ,flatRegMap,nFlatRegCells,steepestDescentSlope ...
                ,iSubFldRegMap,doNotAccFlowCell);
            
            % calculate the number of upstream cells of the remainder:
            % from the regional minima to its outlet
            
            % Identify regional minima in the sub-flooded region
            iRegMin = regionalMin & iSubFldRegMap;
            iSubFldRegMinIdx = find(iRegMin);
            if numel(iSubFldRegMinIdx) > 1
                flatIdx{end+1} ... % record the indices of flat for later processing
                    = iSubFldRegMinIdx;
                iSubFldRegMinIdx = find(iRegMin ...
                    & (steepestDescentSlope <= 0 & nFlatRegCells > 0)); % avoid flat without flow direction
            end
            
            
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
                            = nUpstreamCells(downStreamNbrY,downstreamNbrX) ...
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

nFlatRegMin = numel(flatIdx);
if nFlatRegMin > 0
    for t1 = 1:nFlatRegMin

        iFlatRegMin = cell2mat(flatIdx(1,t1));
        iFlatRegMinOutIdx ...
            = intersect(iFlatRegMin,find(nFlatRegCells > 0));
        iFlatRegMinNotOutIdx ...
            = intersect(iFlatRegMin,find(nFlatRegCells == 0));
        nUpstreamCells(iFlatRegMinNotOutIdx) = ...
            nupstreamCells(iFlatRegMinOutIdx) - 1;
        
    end
end