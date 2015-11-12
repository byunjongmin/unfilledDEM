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
% @version 0.1.4. / 2015-11-12
% @author Jongmin Byun
%==========================================================================

% constant
Y = mRows - 2; X = nCols - 2;
NOT_MODIFIED = 0; % for mFlowDir_Saddle
% for mFlowDir_SubFldReg
FROM_REGMIN_TO_UP = 2; % from regional minima to upward
FROM_TARWSD_TO_TARWSD = 3; % from target to target sub flooded region

ithOffset ... % offset for ith neighbour
    = [mRows,mRows-1,-1,-mRows-1,-mRows,-mRows+1,1,mRows+1];

% main body ---------------------------------------------------------------
% initial upstream coordinate
initUpStreamIdx = sub2ind([mRows,nCols],upStreamY,upStreamX);
firstMet_tf = true; % for a shared true outlet, the variable to indicate
                    % that it has already visited a sub-flooded region
pathNotDone = true;
while pathNotDone
    
    % A. find the Steepest Downstream Cell (SDC) from the initUpstreamIdx
    upStreamIdx = sub2ind([mRows,nCols],upStreamY,upStreamX);
    steeperSlope = 0; % for the target sub-flooded region
    iSteeperSlope = 0; % for the previous gone sub-flooded region
    for ithNbr = 1:8
        
        ithNbrIdx = upStreamIdx + ithOffset(ithNbr);
        [ithNbrY,ithNbrX] = ind2sub([mRows,nCols],ithNbrIdx);
        
        % for the cell within the domain boundary            
        if (1 < ithNbrY && ithNbrY < Y+2) ...
               && (1 < ithNbrX && ithNbrX < X+2)
           
            % if the init upstream cell is not a shared outelt
            if sharedOutlet(initUpStreamIdx) == false
                
                % check if ith neighbor is within the ith flooded region
                if fldRegID(ithNbrIdx) == ithFldRegID
                    
                    % check if the init upstream cell is the outlet of ith
                    % flooded region
                    if upStreamIdx == ithFldRegOutIdx
                    
                        if slopeAllNbr(upStreamY,upStreamX,ithNbr) > steeperSlope
                            
                            steeperSlope = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                            steepestNbrY = ithNbrY;
                            steepestNbrX = ithNbrX;
                            
                        end
                    
                    % check whether this condition exist
                    else % upStreamIdx ~= ithFldRegOutIdx
                        
                        % check if ith Nbr is within the previous or target
                        % sub-flooded region
                        if subFldRegID(ithNbrIdx) == prevSubFldRegID
                        
                            if slopeAllNbr(upStreamY,upStreamX,ithNbr) > iSteeperSlope
                            
                                iSteeperSlope = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                                iSteepestNbrY = ithNbrY;
                                iSteepestNbrX = ithNbrX;
                                
                            end
                                
                        elseif subFldRegID(ithNbrIdx) == targetSubFldRegID
                                
                            if slopeAllNbr(upStreamY,upStreamX,ithNbr) > steeperSlope

                                steeperSlope = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                                steepestNbrY = ithNbrY;
                                steepestNbrX = ithNbrX;
                                    
                            end
                        end % if subFldRegId(ithNbrIdx)
                    end % if upStreamIdx == ithFldRegOutIdx
                end % if fldRegID(ithNbrIdx)
                
            else % sharedOutlet(initUpStreamIdx) == true
                
                if upStreamIdx == ithFldRegOutIdx
                   
                    if isnan(targetSubFldRegID)
                        % if targetSubFldRegID is not defined,
                        % firstly, figure out a first-met sub-flooded region,
                        % and then, define it as a target
                        if fldRegID(ithNbrIdx) == ithFldRegID

                            if firstMet_tf == true

                                if subFldRegID(ithNbrIdx) > 0

                                    if ~ismember(subFldRegID(ithNbrIdx),outConnectedSubFldRegID)

                                        if slopeAllNbr(upStreamY,upStreamX,ithNbr) ...
                                                > steeperSlope

                                            steeperSlope ...
                                                = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                                            steepestNbrY = ithNbrY;
                                            steepestNbrX = ithNbrX;

                                            if firstMet_tf == true
                                                firstMetSubFldRegID = subFldRegID(ithNbrIdx);
                                                firstMet_tf = false;
                                            end
                                        end
                                    end % if ~ismember(subFldRegID(ithNbrIdx)
                                end % if subFldRegID(ithNbrIdx)

                            else % firstMet_tf == false

                                if subFldRegID(ithNbrIdx) == firstMetSubFldRegID

                                    if slopeAllNbr(upStreamY,upStreamX,ithNbr) ...
                                            > steeperSlope

                                        steeperSlope ...
                                            = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                                        steepestNbrY = ithNbrY;
                                        steepestNbrX = ithNbrX;

                                    end
                                end % if subFldRegID(ithNbrIdx)                         
                            end % if firstMet_tf
                        end % if fldRegID(ithNbrIdx)

                    else % ~isnan(targetSubFldRegID)

                        % check if ith neighbor is within the target sub flooded region
                        if subFldRegID(ithNbrIdx) == targetSubFldRegID

                            if slopeAllNbr(upStreamY,upStreamX,ithNbr) ...
                                    > steeperSlope

                                steeperSlope ...
                                    = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                                steepestNbrY = ithNbrY;
                                steepestNbrX = ithNbrX;

                            end
                        end % if subFldRegID(ithNbrIdx)
                    end % if isnan(targetSubFldRegID)
                    
                else
                    
                    % check if ith neighbor is within the ith flooded region
                    if fldRegID(ithNbrIdx) == ithFldRegID
                       
                        % check if ith Nbr is within the previous or target
                        % sub-flooded region
                        if subFldRegID(ithNbrIdx) == prevSubFldRegID

                            if slopeAllNbr(upStreamY,upStreamX,ithNbr) > iSteeperSlope

                                iSteeperSlope = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                                iSteepestNbrY = ithNbrY;
                                iSteepestNbrX = ithNbrX;

                            end

                        elseif subFldRegID(ithNbrIdx) == targetSubFldRegID

                            if slopeAllNbr(upStreamY,upStreamX,ithNbr) > steeperSlope

                                steeperSlope = slopeAllNbr(upStreamY,upStreamX,ithNbr);
                                steepestNbrY = ithNbrY;
                                steepestNbrX = ithNbrX;

                            end
                        end % if subFldRegId(ithNbrIdx)
                    end % if fldRegID(ithNbrIdx)
                    
                end % if upStreamIdx
            end
        end
    end
        
    % if targetSubFldRegID is nan, upstream cell would be the outlet of
    % flooded region. therefore define the targetSubFldReg using the sub
    % flooded region ID of steepestNbrIdx
    if isnan(targetSubFldRegID)

        targetSubFldRegID = subFldRegID(steepestNbrY,steepestNbrX);

    end

    % B. if upstream cell is a (shared) saddle, change flow
    % direction of the cell to go in the gone watershed
    if ~isnan(prevSubFldRegID) ... % in the case of sub flooded region
            && (initUpStreamIdx == upStreamIdx) % initial upstream cell is the outlet of sub flooded region

        if mFlowDir_Saddle(upStreamIdx) == NOT_MODIFIED ...
                && isnan(mFlowDir_SubFldReg(upStreamIdx)) % except for the already corrected col
            m2SDSNbrY(upStreamIdx) = iSteepestNbrY;
            m2SDSNbrX(upStreamIdx) = iSteepestNbrX;
        end
        mFlowDir_Saddle(upStreamIdx) = subFldRegID(iSteepestNbrY,iSteepestNbrX);

    end

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
        m2SDSNbrY(steepestNbrY,steepestNbrX) = upStreamY; % นÝบน
        m2SDSNbrX(steepestNbrY,steepestNbrX) = upStreamX;
        mFlowDir_SubFldReg(steepestNbrY,steepestNbrX) = FROM_TARWSD_TO_TARWSD;

        % (b) Change upstream cell with downstream cell
        upStreamY = steepestNbrY;
        upStreamX = steepestNbrX;

    end % if regionalMin(steepestNbrY,
end  % while pathNotDone