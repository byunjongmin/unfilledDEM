function [flood,m1SDSNbrY,m1SDSNbrX,m2SDSNbrY,m2SDSNbrX,fldRegID,nFldRegCells ...
    ,subFldRegOutlet,subFldRegOutInfo,subFldRegID,regionalMin,sharedOutlet] ...
    = ProcessSink(DEM,chosenWatershed,slopeAllNbr ...
    ,SDSNbrY,SDSNbrX,SDSFlowDirection)
% @file ProcessSink.m
% @brief Process sinks for flow to continue to move
% 
% @section Introduction - 1. Identify flooded regions originated from each
%   sink and their outlets. 2. Modify flow direction of each outlet for flow
%   to move out from each flooded region. During the process, identify
%   sub-flooded regions and their outlets where the flooded regions are
%   interconnected.
%
% @retval flood: flooded region
% @retval m1SDSNbrY,m1SDSNbrX: downstream neighbor coordinate modified for
%           both outlets and flooded regions
% @retval m2SDSNbrY,m2SDSNbrX: downstream neighbor coordinate modified only
%           fot outlets
% @retval fldRegID: flooded region's ID
% @retval subFldRegOutlet: sub-flooded region's outlet showing its status
% @retval subFldRegOutInfo: information on sub-flooded region outlets
% @retval sharedOutlet: number of shared depressions for an shared outlet
%
% @version: 2.2.0 / 2015-11-17
% @authoer: Jongmin Byun
%==========================================================================

% Constants
[mRows,nCols] = size(DEM);
Y = mRows - 2; X = nCols - 2;
Y_TOP_BND = 1; Y_BOTTOM_BND = mRows;
Y_INI = 2; Y_MAX = Y+1;
X_LEFT_BND = 1; X_RIGHT_BND = nCols;
X_INI = 2; X_MAX = X+1;
OUTER_BOUNDARY = true(mRows,nCols);
OUTER_BOUNDARY(Y_INI:Y_MAX,X_INI:X_MAX) = false;
ithNbrYOffset = [0 -1 -1 -1  0  1  1  1];
ithNbrXOffset = [1  1  0 -1 -1 -1  0  1];
VERY_HIGH = inf;
UNFLOODED = 0; CURRENT_FLOODED = 1; OLD_FLOODED = 2; SINK = 3;
% Note that OLD_FLOODED don't recognize the goneSubFldRegID.
OUTLET_AT_BOUNDARY = 1; NEW_OUTLET = 4;
% for subFldRegOutlet
OUTLET_DRY_NBR = 1; % outlet to a dry neighbor
OUTLET_WET_NBR = 2; % outlet to a wet neighbor
OUTLET_WET_NBR_FALSE = 3; % outlet linked to another flooded region with the same elevation outlet. so not true oulet
SHARED_OUTLET_DRY_NBR = 4; % shared outlet to a dry neighbor
SHARED_OUTLET_WET_NBR = 5; % shared outlet to a wet neighbor
SHARED_OUTLET_WET_NBR_FALSE = 6; % shared outlet linked to another flooded region with the same elevation outlet
SHARED_OUTLET_SURROUNDED_HIGH_ELEV = 7; % shared outlet surrounded by the cells with higher elevation. so not true outlet. so not true outlet.

% Variables
% for flooded regions and their outlets
m1SDSNbrY = SDSNbrY; m1SDSNbrX = SDSNbrX;
% for flooded region's outlets and flat
m2SDSNbrY = SDSNbrY; m2SDSNbrX = SDSNbrX;

flood = zeros(mRows,nCols); % status of flooding
ithFldReg = 0;              % ith flooded region
fldRegID = zeros(mRows,nCols);      % flooded region ID
nFldRegCells = zeros(mRows,nCols);  % each flooded region cells number
LISTMAX = mRows * nCols;
crtFldRegCell ...           % current fooded region cell
    = struct('Y',zeros(LISTMAX,1),'X',zeros(LISTMAX,1));
subFldRegOutlet = zeros(mRows,nCols); % sub-flooded region outlet location
subFldRegOutInfo = []; % information on sub-flooded region outlets
subFldRegID = zeros(mRows,nCols); % sub-flooded region ID
% number of shared depressions for a shared outlet possessing totally
% different sub-flooded regions
sharedOutlet = zeros(mRows,nCols); 

% Main body ---------------------------------------------------------------

% Define sink
noFlowDirection = isnan(SDSFlowDirection) & ~OUTER_BOUNDARY & chosenWatershed;
regionalMin = imregionalmin(DEM,8);
flood(noFlowDirection & regionalMin) = SINK; % only sink without flat

tmpElev = DEM(noFlowDirection & regionalMin);
[sinkCellsY,sinkCellsX] = find(noFlowDirection & regionalMin);
sortedSinks = [tmpElev,sinkCellsY,sinkCellsX];
sortedSinks = sortrows(sortedSinks,1);

allSinksNo = size(sinkCellsY,1);
for ithSink = 1:allSinksNo
    
    fprintf('%d/%d in ProcessSink function\n',ithSink,allSinksNo) % for debug
    
    OUTLET_FOUNDED = false;
    
    crtSinkY = sortedSinks(ithSink,2); % current sink y coordinate
    crtSinkX = sortedSinks(ithSink,3); % current sink x coordinate
    crtSinkIdx = sub2ind([mRows,nCols],crtSinkY,crtSinkX);
    
    ithFldReg = ithFldReg + 1; % flooded region ID
    flood(crtSinkIdx) = CURRENT_FLOODED; % mark the flood status
    subFldRegID(crtSinkIdx) = ithFldReg; % sub-flooded region ID
    mergedGoneSubFldRegID = ithFldReg; % reset the record of gone sub-flooded regions
    
    % initial list for the cells belonging to the ith flooded region
    crtFldRegCell.Y(:) = 0; crtFldRegCell.X(:) = 0;
    % total number of cells belonging to the current flooded region
    nCrtFldRegCells = 1;
    crtFldRegCell.Y(nCrtFldRegCells) = crtSinkY;
    crtFldRegCell.X(nCrtFldRegCells) = crtSinkX;
    
    while (OUTLET_FOUNDED == false)
        
        lowerElev = VERY_HIGH; % creiterion for the outlet candidate
        i = 1;
        
        % A. Find an outlet candidate of which neighbor is dry and
        % elevation is slightly higher. Note that the output variables for
        % the following while loop are flood, nCrtFldRegCells,
        % crtFldRegCellY, -X, outletCandY, -X, crtCellY, -X.
        while (i <= nCrtFldRegCells)
            
            % index of current cell
            crtCellY = crtFldRegCell.Y(i);
            crtCellX = crtFldRegCell.X(i);
            
            for j = 1:8
                
                % index of neighbor
                nbrY = crtCellY + ithNbrYOffset(j);
                nbrX = crtCellX + ithNbrXOffset(j);
                nbrIdx = sub2ind([mRows,nCols],nbrY,nbrX);
                
                % if the neighbor is dry, ... Note that when meeting a SINK,
                % it should be recognised as UNFLOODED region
                if (flood(nbrIdx) == UNFLOODED || flood(nbrIdx) == SINK)
                    
                    if DEM(nbrIdx) <= lowerElev
                        % outlet candidate
                        outletCandY = nbrY;
                        outletCandX = nbrX;
                        lowerElev = DEM(nbrIdx);
                        % record for the cell immediately before the outlet
                        % candidate
                        beforeOutletCandY = crtCellY;
                        beforeOutletCandX = crtCellX;
                    end
                
                % if the neighbor is in the OLDFLOODED,
                elseif flood(nbrIdx) == OLD_FLOODED
                    
                    % if the neighbor is included in the mergedGoneSubFldRegID,
                    % it is not the outlet candidate
                    if ismember(subFldRegID(nbrIdx),mergedGoneSubFldRegID)
                        
                        % add a new neighbor to the list
                        flood(nbrIdx) = CURRENT_FLOODED;
                        nCrtFldRegCells = nCrtFldRegCells + 1;
                        crtFldRegCell.Y(nCrtFldRegCells) = nbrY;
                        crtFldRegCell.X(nCrtFldRegCells) = nbrX;

                    % if the neighbor is a newly encountered old flooded
                    % region, it should be regarded as an outlet candidate
                    % to find a connected sub-flooded region because the
                    % outlet candidate can be a saddle rather than outlet.
                    else
                        
                        if DEM(nbrIdx) <= lowerElev
                            % outlet candidate
                            outletCandY = nbrY;
                            outletCandX = nbrX;
                            lowerElev = DEM(nbrIdx);
                            % record for the cell immediately before the outlet
                            % candidate
                            beforeOutletCandY = crtCellY;
                            beforeOutletCandX = crtCellX;
                        end
                        
                    end  
                % elseif flood(nbrIdx) == CURRENT_FLOODED
                    % skip!
                end
            end % for j = 1:8
            i = i + 1;     
        end % while (i <= nCrtFldRegCells)
        
        % B. Identify the status of the outlet candidate
        % TO DO: Remove the IsBoundary function
        outletCandIdx = sub2ind([mRows,nCols],outletCandY,outletCandX);
        if IsBoundary(outletCandY,outletCandX ...
                ,Y_TOP_BND,Y_BOTTOM_BND,X_LEFT_BND,X_RIGHT_BND)
            OUTLET_CAND_STATUS = OUTLET_AT_BOUNDARY;
        else
            OUTLET_CAND_STATUS = NEW_OUTLET;
        end
        
        % C. Process according to the status of the outlet candidate
        if OUTLET_CAND_STATUS == OUTLET_AT_BOUNDARY
            OUTLET_FOUNDED = true;
        else % OUTLET_CAND_STATUS == NEW_OUTLET
            
            % types for outlet candidate
            OUTLET_DRY_NBR_idx = 0; % outlet to dry neighbor
            OUTLET_WET_NBR_idx = 0; % outlet to wet neighbor
            OUTLET_WET_NBR_FALSE_tf = false; % outlet to wet neighbor false
            OUTLET_SURROUNDED_HIGH_ELEV_tf = false; % outlet surrounded by the cells with higher elevation
            SHARED_OUTLET_DRY_NBR_idx = 0; % shared outlet to dry neighbor
            SHARED_OUTLET_WET_NBR_idx = 0; % shared outlet to wet neighbor
            SHARED_OUTLET_WET_NBR_FALSE_tf = false; % shared outlet to wet neighbor
            SHARED_OUTLET_SURROUNDED_HIGH_ELEV_tf = false; % shared outlet surrouned by the cells with higher elevation

            
            % if the candiate is not an already defined outlet,
            if fldRegID(outletCandIdx) >= 0
                
                % searching for the neighbor enabling to drain
                steeperSlope = 0;
                outConnectedSubFldRegID = []; % outlet connected sub-flooded region ID
                for j = 1:8
                    
                    outletCandNbrY = outletCandY + ithNbrYOffset(j);
                    outletCandNbrX = outletCandX + ithNbrXOffset(j);
                    outletCandNbrIdx = sub2ind([mRows,nCols],outletCandNbrY,outletCandNbrX);
                    
                    % if a neighbor of the candidate is dry,
                    if flood(outletCandNbrIdx) == UNFLOODED ...
                            || flood(outletCandNbrIdx) == SINK
                        
                        if slopeAllNbr(outletCandY,outletCandX,j) ...
                                >= steeperSlope
                            
                            steeperSlope = slopeAllNbr(outletCandY,outletCandX,j);
                            maxSDSNbrY = outletCandNbrY;
                            maxSDSNbrX = outletCandNbrX;
                            OUTLET_DRY_NBR_idx ...
                                = sub2ind([mRows,nCols],maxSDSNbrY,maxSDSNbrX);
                        end
                        
                    % if it is wet, it is a cell involved in another linked
                    % flooded region    
                    elseif flood(outletCandNbrIdx) == OLD_FLOODED
                        
                        oldOutletY = m1SDSNbrY(outletCandNbrIdx);
                        oldOutletX = m1SDSNbrX(outletCandNbrIdx);
                        oldOutletIdx = sub2ind([mRows,nCols],oldOutletY,oldOutletX);
                        
                        if DEM(oldOutletIdx) ~= DEM(outletCandNbrIdx)
                        
                        
                            if slopeAllNbr(outletCandY,outletCandX,j) ...
                                    >= steeperSlope

                                steeperSlope ...
                                    = slopeAllNbr(outletCandY,outletCandX,j);
                                maxSDSNbrY = outletCandNbrY;
                                maxSDSNbrX = outletCandNbrX;
                                OUTLET_WET_NBR_idx ...
                                    = sub2ind([mRows,nCols],maxSDSNbrY,maxSDSNbrX);

                            end
                            
                        else % DEM(oldOutletIdx) == DEM(outletCandNbrIdx)
                            
                            % for a flooded region with the same elevtion outlet,
                            % gather the information of the flooded region
                            if ~ismember(subFldRegID(outletCandNbrIdx),outConnectedSubFldRegID)
                                outConnectedSubFldRegID = [outConnectedSubFldRegID ...
                                    ;subFldRegID(outletCandNbrIdx)];
                            end
                            
                        end
                        
                    % else flood(outletCandNbrIdx) == CURRENT_FLOODED
                        % skip!
                    
                    end
                end % for j = 1:8
                
                if steeperSlope == 0
                    
                    if numel(outConnectedSubFldRegID) == 0
                        
                        OUTLET_SURROUNDED_HIGH_ELEV_tf = true;
                        
                    else % numel(outConnectedSubFldRegID) ~= 0
                        
                        OUTLET_WET_NBR_FALSE_tf = true;
                        
                    end
                    
                    isTrueOutlet = false;
                    
                else
                    
                    % determine whether the outlet candidate is really true
                    maxSDSNbrIdx = sub2ind([mRows,nCols],maxSDSNbrY,maxSDSNbrX);
                    if maxSDSNbrIdx == OUTLET_DRY_NBR_idx ...
                            || maxSDSNbrIdx == OUTLET_WET_NBR_idx
                    
                        isTrueOutlet = true;
                        
                    end
                end
                    
            % if the candidate is an already defined outlet,
            % it is a shared outlet
            else % fldRegID(outletCandIdx) < 0

                steeperSlope = 0;
                outConnectedSubFldRegID = []; % outlet connected sub-flooded region ID
                isSameElev = false;
                for j = 1:8
                    
                    outletCandNbrY = outletCandY + ithNbrYOffset(j);
                    outletCandNbrX = outletCandX + ithNbrXOffset(j);
                    outletCandNbrIdx = sub2ind([mRows,nCols],outletCandNbrY,outletCandNbrX);
                    
                    % if the steepest neighbor of the candidate is
                    % UNFLOODED, the shared outlet is linked to dry neighbor
                    if flood(outletCandNbrIdx) == UNFLOODED ...
                            || flood(outletCandNbrIdx) == SINK
                        
                        if slopeAllNbr(outletCandY,outletCandX,j) ...
                                >= steeperSlope
                            
                            steeperSlope = slopeAllNbr(outletCandY,outletCandX,j);
                            maxSDSNbrY = outletCandNbrY;
                            maxSDSNbrX = outletCandNbrX;
                            SHARED_OUTLET_DRY_NBR_idx ...
                                = sub2ind([mRows,nCols],maxSDSNbrY,maxSDSNbrX);
                        end

                    % if the steepest neighbor of the candidate is in old
                    % flooded region, firstly check whether this flooded
                    % region is the flooded region sharing the outlet or
                    % another linked flooded region
                    elseif flood(outletCandNbrIdx) == OLD_FLOODED
                        
                        oldOutletY = m1SDSNbrY(outletCandNbrIdx);
                        oldOutletX = m1SDSNbrX(outletCandNbrIdx);
                        oldOutletIdx = sub2ind([mRows,nCols],oldOutletY,oldOutletX);

                        if oldOutletIdx == outletCandIdx
                            
                            % for a flooded region sharing the outlet,
                            % gather the information of the flooded region
                            if ~ismember(subFldRegID(outletCandNbrIdx),outConnectedSubFldRegID)
                                outConnectedSubFldRegID = [outConnectedSubFldRegID ...
                                    ;subFldRegID(outletCandNbrIdx)];
                            end
                            
                        else % oldOutletIdx ~= outletCandIdx
                            
                            if DEM(oldOutletIdx) ~= DEM(outletCandIdx)
                                
                                if slopeAllNbr(outletCandY,outletCandX,j) ...
                                    >= steeperSlope

                                    steeperSlope = slopeAllNbr(outletCandY,outletCandX,j);
                                    maxSDSNbrY = outletCandNbrY;
                                    maxSDSNbrX = outletCandNbrX;
                                    SHARED_OUTLET_WET_NBR_idx ...
                                        = sub2ind([mRows,nCols],maxSDSNbrY,maxSDSNbrX);
                                end
                                
                            else % DEM(oldOutletIdx) ~= DEM(outletCandIdx) 
                                
                                % for a shared outlet linked to another flooded region with the same elevation outlet
                                % gather the information of the flooded region
                                if ~ismember(subFldRegID(outletCandNbrIdx),outConnectedSubFldRegID)
                                    outConnectedSubFldRegID = [outConnectedSubFldRegID ...
                                        ;subFldRegID(outletCandNbrIdx)];
                                end
                                isSameElev = true;
                            
                            end
                            
                        end
                        
                    % elseif flood(outletCandNbrIdx) ~= CURRENT_FLOODED
                        % skip!
                        
                    end % if flood(outletCandNbrIdx)
                end % for j = 1:8
                
                if steeperSlope == 0
                    
                    if isSameElev == false

                        SHARED_OUTLET_SURROUNDED_HIGH_ELEV_tf = true;
                    
                    else % isSameElev == true
                        
                        SHARED_OUTLET_WET_NBR_FALSE_tf = true;
                        
                    end 
                    
                    isTrueOutlet = false;
                    
                else
                
                    % determine whether the outlet candidate is really true
                    maxSDSNbrIdx = sub2ind([mRows,nCols],maxSDSNbrY,maxSDSNbrX);
                    if maxSDSNbrIdx == SHARED_OUTLET_DRY_NBR_idx ...
                            || maxSDSNbrIdx == SHARED_OUTLET_WET_NBR_idx

                        isTrueOutlet = true;
                    
                    end
                end
            end % if fldRegID(outletCandIdx) >= 0

            % Perform a following procedure according to whether the outlet
            % candidate is true or not
            if isTrueOutlet == true

                OUTLET_FOUNDED = true;

                % a. mark the type of flooded region outlet
                if maxSDSNbrIdx == OUTLET_DRY_NBR_idx
                    
                    % just for recording
                    subFldRegOutlet(outletCandY,outletCandX) ...
                        = OUTLET_DRY_NBR;
                    
                elseif maxSDSNbrIdx == OUTLET_WET_NBR_idx
                    
                    % just for recording
                    subFldRegOutlet(outletCandY,outletCandX) ...
                        = OUTLET_WET_NBR;
                
                % note that if the true outlet is a shared outlet, find out
                % the sub-flooded regions sharing the outlet.
                elseif maxSDSNbrIdx == SHARED_OUTLET_DRY_NBR_idx ...
                        || maxSDSNbrIdx == SHARED_OUTLET_WET_NBR_idx
                        
                    % a. mark the type of the shared outlet
                    if maxSDSNbrIdx == SHARED_OUTLET_DRY_NBR_idx
                        
                        subFldRegOutlet(outletCandY,outletCandX) ...
                            = SHARED_OUTLET_DRY_NBR;
                        
                    elseif maxSDSNbrIdx == SHARED_OUTLET_WET_NBR_idx
                        
                        subFldRegOutlet(outletCandY,outletCandX) ...
                            = SHARED_OUTLET_WET_NBR;
                        
                    end

                    nOutConnectedSubFldReg = numel(outConnectedSubFldRegID);
                    for i=1:nOutConnectedSubFldReg
                        
                        % b. count the number of connected sub-flooded
                        % regions
                        sharedOutlet(outletCandY,outletCandX) ...
                            = sharedOutlet(outletCandY,outletCandX) + 1;
                        
                        % c. make subFldRegOutInfo for each outlet connected
                        % sub-flooded region
                        subFldRegOutInfo = [subFldRegOutInfo ...
                            ;ithFldReg ...
                            ,subFldRegID(beforeOutletCandY,beforeOutletCandX) ... % adjacent sub-flooded region ID
                            ,outConnectedSubFldRegID(i) ... % outlet connected sub-flooded region ID
                            ,outletCandY,outletCandX];

                        % d. add outlet connected flooded region cells
                        tFldRegID = unique(fldRegID(subFldRegID == outConnectedSubFldRegID(i)));
                        tmpIdx = find(fldRegID == tFldRegID(tFldRegID > 0));
                        [tmpY,tmpX] = ind2sub([mRows,nCols],tmpIdx);
                        crtFldRegCell.Y(nCrtFldRegCells+1:nCrtFldRegCells+numel(tmpIdx)) = tmpY;
                        crtFldRegCell.X(nCrtFldRegCells+1:nCrtFldRegCells+numel(tmpIdx)) = tmpX;
                        nCrtFldRegCells = nCrtFldRegCells + numel(tmpIdx);
                        
                    end
                end % if maxSDSNbrIdx

            else % if isTrueOutlet == false
          
                OUTLET_FOUNDED = false;
                
                % a. mark the type of outlet cand
                if OUTLET_WET_NBR_FALSE_tf == true

                    subFldRegOutlet(outletCandY,outletCandX) ...
                        = OUTLET_WET_NBR_FALSE;
                    % reset old flooded region outlet
                    fldRegID(oldOutletIdx) = 0;


                elseif SHARED_OUTLET_WET_NBR_FALSE_tf == true

                    subFldRegOutlet(outletCandY,outletCandX) ...
                        = SHARED_OUTLET_WET_NBR_FALSE;
                    % reset old flooded region outlet
                    fldRegID(oldOutletIdx) = 0;

                elseif SHARED_OUTLET_SURROUNDED_HIGH_ELEV_tf == true

                    subFldRegOutlet(outletCandY,outletCandX) ...
                        = SHARED_OUTLET_SURROUNDED_HIGH_ELEV;

                end

                nOutConnectedSubFldReg = numel(outConnectedSubFldRegID);
                for i=1:nOutConnectedSubFldReg
                    
                    % b. count the number of connected sub-flooded regions
                    sharedOutlet(outletCandY,outletCandX) ...
                        = sharedOutlet(outletCandY,outletCandX) + 1;
                    
                    % c. make sub-flooded region outlet info
                    subFldRegOutInfo = [subFldRegOutInfo ...
                        ;ithFldReg ...
                        ,subFldRegID(beforeOutletCandY,beforeOutletCandX) ... % adjacent sub-flooded region ID
                        ,outConnectedSubFldRegID(i) ... % shared sub-flooded region ID
                        ,outletCandY,outletCandX];
                    
                    % d. add outlet connected flooded region cells
                    tFldRegID = unique(fldRegID(subFldRegID == outConnectedSubFldRegID(i)));
                    tmpIdx = find(fldRegID == tFldRegID(tFldRegID > 0));
                    [tmpY,tmpX] = ind2sub([mRows,nCols],tmpIdx);
                    crtFldRegCell.Y(nCrtFldRegCells+1:nCrtFldRegCells+numel(tmpIdx)) = tmpY;
                    crtFldRegCell.X(nCrtFldRegCells+1:nCrtFldRegCells+numel(tmpIdx)) = tmpX;
                    nCrtFldRegCells = nCrtFldRegCells + numel(tmpIdx);

                    % e. reset the ID of connected sub-flooded region
                    flood(tmpIdx) = CURRENT_FLOODED;

                    % f. add the shared flooded region into the gone
                    % sub-flooded region ID
                    goneSubFldRegID = outConnectedSubFldRegID(i);
                    mergedGoneSubFldRegID ...
                        = [mergedGoneSubFldRegID;goneSubFldRegID];                    

                end % for i = 1:nOutConnectedSubFldReg
                
                % common tasks for the outlet cand     
                
                % mark its sub-flooded region ID
                % note that subFldRegID of the outlet candidate can be
                % already defined when the candidate is involed in a newly
                % encountered flooded region, and so check the status of
                % subFldRegID.
                if subFldRegID(outletCandY,outletCandX) == 0
                    subFldRegID(outletCandY,outletCandX) ...
                        = subFldRegID ...
                        (SDSNbrY(outletCandY,outletCandX) ...
                        ,SDSNbrX(outletCandY,outletCandX));
                end
                
                flood(outletCandY,outletCandX) = CURRENT_FLOODED;

                nCrtFldRegCells = nCrtFldRegCells + 1;
                crtFldRegCell.Y(nCrtFldRegCells) = outletCandY;
                crtFldRegCell.X(nCrtFldRegCells) = outletCandX;
                
            end % if isTrueOutlet
        end % if OUTLET_CAND_STATUS
    end % while (OUTLET_FOUNDED == false)
    
    % Tasks for true outlet
    
    % a. modify the flow direction of the outlet candidate
    m1SDSNbrY(outletCandY,outletCandX) = maxSDSNbrY;
    m1SDSNbrX(outletCandY,outletCandX) = maxSDSNbrX;
    m2SDSNbrY(outletCandY,outletCandX) = maxSDSNbrY;
    m2SDSNbrX(outletCandY,outletCandX) = maxSDSNbrX;
    
    % b. update flow directions for flooded regions in m1SDSNBRY, -X.
    floodedCellY = crtFldRegCell.Y(1:nCrtFldRegCells);
    floodedCellX = crtFldRegCell.X(1:nCrtFldRegCells);
    floodedCellIndex = sub2ind([mRows,nCols],floodedCellY,floodedCellX);
    m1SDSNbrY(floodedCellIndex) = outletCandY;
    m1SDSNbrX(floodedCellIndex) = outletCandX;
    
    % c. update flood as a OLD_FLOODED
    flood(floodedCellIndex) = OLD_FLOODED;
    
    % d. update fldRegID, nFldRegCells
    fldRegID(floodedCellIndex) = ithFldReg;
    fldRegID(outletCandY,outletCandX) = -ithFldReg;
    nFldRegCells(outletCandY,outletCandX) = nCrtFldRegCells;
    
end % for ithSink = 1

% to differ sub-flooded region ID with flooded region ID
fldRegID = fldRegID.*10^(numel(num2str(max(fldRegID(:)))));