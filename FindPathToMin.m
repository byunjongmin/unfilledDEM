function [m2SDSNbrY,m2SDSNbrX,arrivedRegMinY,arrivedRegMinX ...
    ,mFlowDir_SubFldReg,mFlowDir_Saddle] ...
    = FindPathToMin(mRows,nCols,upStreamY,upStreamX,slopeAllNbr ...
    ,m2SDSNbrY,m2SDSNbrX,regionalMin,fldRegID,subFldRegID,ithFldRegID ...
    ,targetSubFldRegID,prevSubFldRegID,mFlowDir_SubFldReg,mFlowDir_Saddle ...
    ,ithFldRegOutIdx,sharedOutlet,outConnectedSubFldRegID)
% @file FindPathToMin.m
% @brief Function to find a path to a regional minima
%
% @param[in] mRows
% @param[in] nCols
% @param[in] upStreamIdx
% @param[in] slopeAllNbr
% @param[in] flatRegMap
% @param(in] flatRegBndInfo
% @param[in] SDSFlowDirection
% @param[in] SDSNbrY
% @param[in] SDSNbrX
% @param[in] regionalMinima
% @param[in] flatRegionPathInfo
% @param[in] DEM
% @param(in] fldRegWithNoRegMin
% @param[in] fldRegID
% @param[in] subFldRegID
% @param(in] arrivedRegMinIdx
% @param(in] ithFldRegID
% @param[in] targetSubFldRegID
% @param[in] prevSubFldRegID
% @param[in] mFlowDir_SubFldReg
% @param[in] mFlowDir_Saddle
% @param[in] ithFldRegOutIdx
% @param(in] flatRegInfo
% @retval SDSNbrY
% @retval SDSNbrX
% @retval flatRegionPathInfo flow direction information in the flat regions
% @retval arrivedRegMinIdx finally arrived regional minima index
% @retval mFlowDir_SubFldReg
% @retval mFlowDir_Saddle
%
% @version 0.2.0. / 2015-11-12
% @author Jongmin Byun
%==========================================================================

% constant
Y = mRows - 2; X = nCols - 2;
NOT_MODIFIED = 0; % for mFlowDir_Saddle
FROM_REGMIN_TO_UP = 2; % for mFlowDir_SubFldReg, from regional minima to up
FROM_TARWSD_TO_TARWSD = 3; % from target to target sub flooded region
ithOffset ... % offset for ith neighbour
    = [mRows,mRows-1,-1,-mRows-1,-mRows,-mRows+1,1,mRows+1];

% identify the type of an initial upstream cell
IS_SHARED_TRUE_OUTLET = false;
IS_TRUE_OUTLET = false;
IS_SADDLE = false;
upStreamIdx = sub2ind([mRows,nCols],upStreamY,upStreamX);
if upStreamIdx == ithFldRegOutIdx
    
    if sharedOutlet(upStreamIdx) > 0
        IS_SHARED_TRUE_OUTLET = true;
    else
        IS_TRUE_OUTLET = true;    
    end
    
else
    
    if sharedOutlet(upStreamIdx) > 0
        IS_SADDLE = true;
    end
    
end

% find the Steepest Downstream Cell (SDC) at the initUpstreamIdx according
% to the type of the initial upstream cell
steeperSlope = 0; % for the target sub-flooded region
iSteeperSlope = 0; % for the previous gone sub-flooded region
firstMet_tf = true; % for a shared true outlet, the variable to indicate
                    % that it has already visited a sub-flooded region
for i = 1:8

    ithNbrIdx = upStreamIdx + ithOffset(i);
    [ithNbrY,ithNbrX] = ind2sub([mRows,nCols],ithNbrIdx);

    % for the cell within the domain boundary            
    if (1 < ithNbrY && ithNbrY < Y+2) ...
           && (1 < ithNbrX && ithNbrX < X+2)
       
        % for the ith flooded region
        if fldRegID(ithNbrIdx) == ithFldRegID

            % if the init upstream cell is true outlet
            if IS_TRUE_OUTLET == true

                if slopeAllNbr(upStreamY,upStreamX,i) > steeperSlope
                            
                    steeperSlope = slopeAllNbr(upStreamY,upStreamX,i);
                    steepestNbrY = ithNbrY;
                    steepestNbrX = ithNbrX;

                end
                
            elseif IS_SHARED_TRUE_OUTLET == true
                                   
                if isnan(targetSubFldRegID)
                    
                    % if targetSubFldRegID is not defined,
                    % firstly, figure out a first-met sub-flooded region,
                    % and then, define it as a target
                    if firstMet_tf == true
                        
                        if ~ismember(subFldRegID(ithNbrIdx),outConnectedSubFldRegID)

                            if slopeAllNbr(upStreamY,upStreamX,i) ...
                                > steeperSlope

                                steeperSlope ...
                                    = slopeAllNbr(upStreamY,upStreamX,i);
                                steepestNbrY = ithNbrY;
                                steepestNbrX = ithNbrX;

                                if firstMet_tf == true
                                    firstMetSubFldRegID = subFldRegID(ithNbrIdx);
                                    firstMet_tf = false;
                                end
                            end
                        end % if ~ismember(subFldRegID(ithNbrIdx)

                    else % firstMet_tf == false

                        if subFldRegID(ithNbrIdx) == firstMetSubFldRegID

                            if slopeAllNbr(upStreamY,upStreamX,i) ...
                                    > steeperSlope

                                steeperSlope ...
                                    = slopeAllNbr(upStreamY,upStreamX,i);
                                steepestNbrY = ithNbrY;
                                steepestNbrX = ithNbrX;

                            end
                        end % if subFldRegID(ithNbrIdx)                         
                    end % if firstMet_tf

                else % ~isnan(targetSubFldRegID)

                    % check if ith neighbor is within the target sub flooded region
                    if subFldRegID(ithNbrIdx) == targetSubFldRegID

                        if slopeAllNbr(upStreamY,upStreamX,i) ...
                                > steeperSlope

                            steeperSlope ...
                                = slopeAllNbr(upStreamY,upStreamX,i);
                            steepestNbrY = ithNbrY;
                            steepestNbrX = ithNbrX;

                        end
                    end % if subFldRegID(ithNbrIdx)
                end % if isnan(targetSubFldRegID)
                
            elseif IS_SADDLE == true
                        
                % check if ith Nbr is within the previous or target
                % sub-flooded region
                if subFldRegID(ithNbrIdx) == prevSubFldRegID

                    if slopeAllNbr(upStreamY,upStreamX,i) > iSteeperSlope

                        iSteeperSlope = slopeAllNbr(upStreamY,upStreamX,i);
                        iSteepestNbrY = ithNbrY;
                        iSteepestNbrX = ithNbrX;

                    end

                elseif subFldRegID(ithNbrIdx) == targetSubFldRegID

                    if slopeAllNbr(upStreamY,upStreamX,i) > steeperSlope

                        steeperSlope = slopeAllNbr(upStreamY,upStreamX,i);
                        steepestNbrY = ithNbrY;
                        steepestNbrX = ithNbrX;

                    end
                end % if subFldRegId(ithNbrIdx)
            end % if IS_TRUE_OUTLET
        end % if fldRegID(ithNbrIdx)
    end % if (1 < ithNbrY
end % for i = 1:8

% if upstream cell is a (shared) saddle, change flow
% direction of the cell to go in the gone watershed
if IS_SADDLE == true

    if mFlowDir_Saddle(upStreamIdx) == NOT_MODIFIED ...
            && isnan(mFlowDir_SubFldReg(upStreamIdx)) % except for the already corrected col
        m2SDSNbrY(upStreamIdx) = iSteepestNbrY;
        m2SDSNbrX(upStreamIdx) = iSteepestNbrX;
    end
    mFlowDir_Saddle(upStreamIdx) = subFldRegID(iSteepestNbrY,iSteepestNbrX);

end

% continue to find a path to regional minima

% if targetSubFldRegID is nan, upstream cell would be the outlet of
% flooded region. therefore define the targetSubFldReg using the sub
% flooded region ID of steepestNbrIdx
if isnan(targetSubFldRegID)
    targetSubFldRegID = subFldRegID(steepestNbrY,steepestNbrX);
end

pathNotDone = true;
while pathNotDone

    % C. check if the SDC is the regional minima
    % a. if it is, assign the flow direction and end the loop
    if regionalMin(steepestNbrY,steepestNbrX) == true

        % a) assign flow direction to the downstream cell
        m2SDSNbrY(steepestNbrY,steepestNbrX) = upStreamY;
        m2SDSNbrX(steepestNbrY,steepestNbrX) = upStreamX;
        mFlowDir_SubFldReg(steepestNbrY,steepestNbrX) = FROM_REGMIN_TO_UP;

        % b) identify the regional minima
        arrivedRegMinY = steepestNbrY;
        arrivedRegMinX = steepestNbrX;

        % c) end the loop
        pathNotDone = false;

    % b. If the DOC is not a regional minima, assign the flow direction to
    % upstream cell and change upstream cell to the downstream cell
    else

        % (a) assign flow direction to the downstream cell
        m2SDSNbrY(steepestNbrY,steepestNbrX) = upStreamY;
        m2SDSNbrX(steepestNbrY,steepestNbrX) = upStreamX;
        mFlowDir_SubFldReg(steepestNbrY,steepestNbrX) = FROM_TARWSD_TO_TARWSD;

        % (b) Change upstream cell with downstream cell
        upStreamY = steepestNbrY;
        upStreamX = steepestNbrX;

        % continue to find the SDC
        steeperSlope = 0; % for the target sub-flooded region
        for i = 1:8

            ithNbrIdx = upStreamIdx + ithOffset(i);
            [ithNbrY,ithNbrX] = ind2sub([mRows,nCols],ithNbrIdx);

            if subFldRegID(ithNbrIdx) == targetSubFldRegID

                if slopeAllNbr(upStreamY,upStreamX,i) > steeperSlope

                    steeperSlope = slopeAllNbr(upStreamY,upStreamX,i);
                    steepestNbrY = ithNbrY;
                    steepestNbrX = ithNbrX;

                end
            end
        end % for i = 1
    end % if regionalMin(steepestNbrY,
end % while pathNotDone
