function nUpstreamCells ...
    = CalcUpstreamCellsInSubFldReg(root,mFlowDir_SubFldReg ...
    ,mFlowDir_Saddle,fldRegID,nUpstreamCells,fldRegInfo ...
    ,m2SDSNbrY,m2SDSNbrX,subFldRegID,DEM,regionalMin)
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
% @version 2.1.0 / 2015-11-18
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

treeDepth = root.depthtree; % coordinated tree where each node holds its depth
nChildNode = tree(root,0); % copy-constructor for the root tree. a tree for the number of children nodes
nChildNode2 = tree(root,0); % a tree for the number of children nodes to be processed

nUpstreamCells(fldRegID > 0) = 0; % mask for current sub-flooded region, reset upstream cells numbers for flooded regions
prevNUpstreamCells = nUpstreamCells; % for whole region

doNotAccFlowMask = ... % cells not to be processed for the first time in mCalcUpstreamCells function
    mFlowDir_SubFldReg == FROM_TARWSD_TO_TARWSD ...
    | mFlowDir_SubFldReg == FROM_REGMIN_TO_UP;
markForGoneCells = nan(mRows,nCols); % to mark trace
markForGoneCells(doNotAccFlowMask) = mFlowDir_SubFldReg(doNotAccFlowMask);

nFldReg = numel(fldRegInfo(:,1));
for i = 1:nFldReg
    
    % for debug
    fprintf(' %d/%d in CalcUpstreamCellsinSubFldReg function\n',i,nFldReg);
    
    % A. choose a sub-flooded region according to the distance from ith
    % flooded region outlet
    
    % a. identify the ith flooded region ID in tree data structure
    ithFldRegID = fldRegInfo(i,1); % ith flooded region ID
    nodeIdx = find(treeDepth == 1);
    for j = nodeIdx
        if root.get(j) == ithFldRegID
            ithFldRegTreeID = j;
        end
    end
    % ithFldRegID = find(strcmp(root,num2str(ithFldRegID));
    
    % b. identify all leaf nodes of the ith flooded region
    % and calculate how many children nodes each parent node has
    leafNodeID = [];
    % make a vector of nodes that traverse tree in a depth first manner
    nodeIdx = depthfirstiterator(root,ithFldRegTreeID); % for nodes of tree depth = 1
    nodeIdx(1) = []; % remove the head of ithFldRegTreeID
    for j = nodeIdx
        
        % make a list for leaf nodes
        if root.isleaf(j)
            leafNodeID = [leafNodeID,j];
        end
        
        % calculate the number of children nodes for each node
        childrenID = root.getchildren(j);
        nChildren = numel(childrenID);        
        if nChildren >= 1
            nChildNode = nChildNode.set(j,nChildren);
            nChildNode2 = nChildNode2.set(j,nChildren);
        end
        
    end
    % B. calculate the number of upstream cells for the cells in ith
    % flooded region starting from the leaf node
    
    % iterate until each leaf node goes up and meets its own parent node
    nLeafNode = numel(leafNodeID);
    for j = 1:nLeafNode
        
        childID = leafNodeID(j); % ith leaf node        
        pathANotDone = true;
        while pathANotDone
            
            currentSubFldRegID ... % current sub-flooded region ID
                = root.get(childID);
            iSubFldRegMap ... % current sub-flooded region
                = subFldRegID == currentSubFldRegID;
            dIthSubFldRegMap ... % dilated current sub flooded region map
                = imdilate(iSubFldRegMap,s);
            
            % define the area for calculating the number of upstream cells
            % over the current sub-flooded region
            % Note: avoid flow direction modified cells
            targetDrainage = ~isnan(DEM) ...
                & dIthSubFldRegMap ...
                & (markForGoneCells ~= FROM_REGMIN_TO_UP ... % flow direction modified cell
                    & markForGoneCells ~= FROM_TARWSD_TO_TARWSD);
                    
            nUpstreamCells = prevNUpstreamCells;
            % calculate the number of upstream cells over chosen area
            [nUpstreamCells,markForGoneCells] ...
                = mCalcUpstreamCells(nUpstreamCells,DEM ...
                ,targetDrainage,m2SDSNbrY,m2SDSNbrX ...
                ,iSubFldRegMap,markForGoneCells);
            
            % calculate upstream cells number for the remainder:
            % from the regional minima to its outlet
            
            % Identify the regional minima in the sub-flooded region
            iRegMin = regionalMin & iSubFldRegMap;
            iSubFldRegMinIdx = find(iRegMin);            
            
            % continue until the escape of the current sub-flooded region
            [downStreamY,downStreamX] = ind2sub([mRows,nCols],iSubFldRegMinIdx);
            pathBNotDone = true;
            while pathBNotDone
            
                % add one for itself
                if markForGoneCells(downStreamY,downStreamX) ...
                        ~= GONE_WSD % avoid the flat already calculated
                    nUpstreamCells(downStreamY,downStreamX) ...
                        = nUpstreamCells(downStreamY,downStreamX) + 1;
                end
                
                upStreamNbrY = m2SDSNbrY(downStreamY,downStreamX);
                upStreamNbrX = m2SDSNbrX(downStreamY,downStreamX);

                if dIthSubFldRegMap(upStreamNbrY,upStreamNbrX) == true
                    
                    if markForGoneCells(downStreamY,downStreamX) ...
                            ~= GONE_WSD % avoid the flat already calculated
                        
                        markForGoneCells(downStreamY,downStreamX) = GONE_WSD;
                        % transfer upstream cells number to the downstream cell
                        
                        nUpstreamCells(upStreamNbrY,upStreamNbrX) ...
                            = nUpstreamCells(upStreamNbrY,upStreamNbrX) ...
                            + nUpstreamCells(downStreamY,downStreamX);
                        
                    end
                    
                    downStreamY = upStreamNbrY;
                    downStreamX = upStreamNbrX;
                    
                else

                    pathBNotDone = false;
                    markForGoneCells(downStreamY,downStreamX) = GONE_WSD;
                     
                end % if iSubFldRegMap
            end % while pathBNotDeon
                        
            % update upstream cells number for the current flooded region
            prevNUpstreamCells(dIthSubFldRegMap) ...
                = nUpstreamCells(dIthSubFldRegMap);
            
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
            
        end % while pathANotDone
    end % for j = 1:nLeafNode
end