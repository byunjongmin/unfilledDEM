function [rawDEMwithBnd,rndDEMwithBnd,nUpstreamCells,m2SDSNbrY,m2SDSNbrX ...
    ,mFlowDir_SubFldReg,mFlowDir_Saddle,subFldRegTree,fldRegInfo ...
    ,outletX,outletY] ...
    = CalcFlowPathFromUnfilledDEM(INPUT_DIR,DEMFileName)
% @file main.m
%
% @brief Function to derive flow path and calculate upstream cells number
%        from unfilled DEM
%
% @param[in] INPUT_DIR
% @param[in] DEMFileName
%
% @retval rawDEMwithBnd: raw DEM with nan boundary
% @retval rndDEMwithBnd: ramdon values added DEM with nan boundaries. It is
%         required to regenerate outputs
% @retval nUpstreamCells: flow accumulation including flooded regions
% @retval m2SDSNbrY,m2SDSNbrX: flow direction modified cells along the path
%         to each regional minima
% @retval mFlowDir_SubFldReg: flow direction modified cell within a
%         sub-flooded region 
% @retval mFlowDir_Saddle: flow direction modified cell on a saddle
% @retval subFldRegTree
% @retval fldRegInfo
%
% @version 1.0.0 / 2019-08-13
% @author Jongmin Byun
%==========================================================================

%% Load DEM

% Constants

DEMFilePath = fullfile(INPUT_DIR,DEMFileName);
[DEM,R] = geotiffread(DEMFilePath);
% Note that R is MapCellsReference or GeographicPostingsReference
% depending on the CoordinateSysteamType

% Basic properties of the imported DEM
mRows = R.RasterSize(1,1);
nCols = R.RasterSize(1,2);

% Check whether the DEM is projected to get the cell size: dX and dY
% 개선: 투영한 것만 사용할 것. geographic 인 경우에 셀 크기가 일정하지 않음
if strcmp(R.CoordinateSystemType,'planar')
    
    dX = R.CellExtentInWorldX;
    dY = R.CellExtentInWorldY;
    
else % strcmp(R,'geographic')
    
    CIR_EQUATOR = 40075; % circumference at equator [km]
    CIR_MERIDIAN = 40008; % meridional circumference [km]
    rLat = (R.LatitudeLimits(2) - R.LatitudeLimits(1)) * 0.5;
    
    dY = CIR_MERIDIAN * 1000 / 360 * R.SampleSpacingInLatitude;
    dX = CIR_EQUATOR * 1000 / 360 * cos(rLat * pi / 180) * R.SampleSpacingInLatitude;

end

% Save the imported raw DEM
rawDEM = DEM;

% Add boundaries for DEM to avoid boundary-related errors
mRows = mRows + 2;
nCols = nCols + 2;
DEM = nan(mRows,nCols);
DEM(2:end-1,2:end-1) = rawDEM;

rawDEMwithBnd = DEM; % raw DEM with boundaries

% Make a mask for nan
nanMask = isnan(DEM);
nanMask(DEM == 0) = true; % include 0 as null value
DEM(nanMask) = nan;

DEM = double(DEM); % for integer based DEM

% Draw DEM
figure(1); clf;
imagesc(DEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Digital Elevation Model');

%% Remove flat cells in DEM by adding random elevation values and smoothing
% using a specific filter
% To be improved: Smoothing should be avoided because of modification of
% raw DEM and generation of isolated area near the boundary of DEM

% Check whether flat cells exhibit

% Assign flow directions to the DEM using D8 algorithm
% Note that, before calculation of flow direction of the target drainage,
% fill NaN areas with higher elevation values
% For future renference, filling with lower elevation is used for GPSS
% because all flow should be directed for model outside.
DEM(nanMask) = inf;
[~,slopeAllNbr,~,~,~] = CalcSDSFlow(DEM,dX,dY);
DEM(nanMask) = nan;

% Make a map of flat cells
flatRegMap = ProcessFlat(DEM,~nanMask,slopeAllNbr);

nFlatCells = sum(flatRegMap(:));
while nFlatCells > 0
    
    fprintf('Flat areas are found in the study DEM.\n');
    fprintf('Random values (<10^6) are added to DEM to exclude the flat area.\n');
    
    % Add very little randomly generated values to DEM
    randElev = rand(mRows,nCols);
    DEM = DEM + randElev * 0.00001;
    
    % Or smooth DEM
    % 
    % % make a bell shaped weight window, h
    % fSize = 2; % radius of window
    % transRatio = 1; % ratio of transfering values to next neighbors
    % h = ones(fSize*2-1,fSize*2-1);
    % for i = 1:fSize
    %     h(i:end-(i-1),i:end-(i-1)) = transRatio^(fSize-i);
    % end
    % h = 1/sum(h(:)) * h;
    % 
    % % smooth using the filter h
    % nSmooth = 1; % number of times
    % for i=1:nSmooth
    %     DEM = filter2(h,DEM);
    % end
    
    % display difference between the original and smoothed dEM
    diffDEM = rawDEMwithBnd - DEM;

    figure(2); clf;

    subplot(1,2,1);
    imagesc(flatRegMap);
    set(gca,'DataAspectRatio',[1 1 1]);
    colorbar;
    title('Flat cells in DEM');

    subplot(1,2,2);
    diffDEM(nanMask) = nan;
    imagesc(diffDEM);
    set(gca,'DataAspectRatio',[1 1 1]);
    colorbar;
    title('Difference between the original and smoothed DEM');   
    
    % Assign flow directions to the DEM using D8 algorithm
    % Note that, before calculation of flow direction of the target drainage,
    % fill NaN areas with higher elevation values
    % For future renference, filling with lower elevation is used for GPSS
    % because all flow should be directed for model outside.
    DEM(nanMask) = inf;
    [~,slopeAllNbr,~,~,~] = CalcSDSFlow(DEM,dX,dY);
    DEM(nanMask) = nan;

    % Make a map of flat cells
    flatRegMap = ProcessFlat(DEM,~nanMask,slopeAllNbr);
    nFlatCells = sum(flatRegMap(:));
    
    rndDEMwithBnd = DEM; % ramdon values added DEM with boundaries
end

%% Make a test domain for debugging

% IS_IT_PART = false;
% if IS_IT_PART == true
% 
%     % coordinates of the test domain
%     tYMin = 558; tYMax = 658;
%     tXMin = 117; tXMax = 168;
%     
%     % cut DEM
%     rawDEMwithBnd = DEM(tYMin:tYMax,tXMin:tXMax); % original DEM
%     DEM = DEM(tYMin:tYMax,tXMin:tXMax); % smoothed DEM
%     
%     [mRows,nCols] = size(DEM);
%     
%     % mark a mask of nan
%     nanMask = (DEM == 32767 | DEM == -9999); % avoid null values of DEM
%     nanMask(DEM == 0) = true; % include 0 as null value
%     
%     % for debug
%     figure(3); clf;
%     set(gcf, 'Color',[1,1,1]);
% 
%     DEM(nanMask) = nan;
%     imagesc(DEM);
%     DEM(nanMask) = rawDEMwithBnd(nanMask);
%     
%     set(gca,'DataAspectRatio',[1 1 1]);
%     colorbar;
%     title('DEM of the test domain');
%     
% end

%% Pick the outlet of the study drainage and replace the outlet with nan
% Note that, to use ProcessSink function, a mask should be made for the
% outer region of the target drainage including its outlet

targetDrainage = ~nanMask; % study areas

% Pick the outlet of the target drainage
DEM(nanMask) = inf;
% Extract the boundary of the target drainage
s = strel('square',3); % structural element when eroding image
dilatedTarget = imerode(targetDrainage,s); 
targetBnd = targetDrainage & ~dilatedTarget;
targetBndIdx = find(targetBnd);

% Identify the coordinate of the outlet
targetBndElev = DEM(targetBndIdx);
[~,minElevIdx] = min(targetBndElev);
outletIdx = targetBndIdx(minElevIdx); % outlet on the boundary of drainage
[outletY,outletX] = ind2sub([mRows,nCols],outletIdx);

% Note that the elevation of the main outlet should be the lowest.
DEM(outletY,outletX) = min(DEM(:)) - 0.1;

% Locate the outlet of target drainage
targetDrainage(outletY,outletX) = false;

%% Assign flow direction to the target area in DEM using D8 algorithm

[~,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
    = CalcSDSFlow(DEM,dX,dY);

%% Process sinks
% Identify depressions and their outlets, then modify each depression
% outlet's flow direction to go downstream if flows are overspilled
% over the depression

[flood ... % flooded region map
,m1SDSNbrY,m1SDSNbrX ... % modified for flooed region and its outlet
,m2SDSNbrY,m2SDSNbrX ... % modified only for flooded region's outlets 
,fldRegID ... % flooded region ID
,nFldRegCells ... % each flooded region cells number
,subFldRegOutlet ... % sub-flooded region outlet location
,subFldRegOutInfo ... % information on sub-flooded region outlets
,subFldRegID ... % sub-flooded region ID
,regionalMin ... % regional minimum map
,sharedOutlet] ... % shared outlet map
    = ProcessSink(DEM,targetDrainage,slopeAllNbr ...
        ,SDSNbrY,SDSNbrX,SDSFlowDirection);

% Frequency distribution of the types of sub-flooded region's outlet
figure(3); clf;
set(gcf, 'Color',[1,1,1]);
histogram(subFldRegOutlet(subFldRegOutlet > 0));
xlabel('Type of Sub-flooded Region Outlet');
ylabel('Frequency');
grid on

% For reference
% 1: outlet to a dry neighbor
% 2: outlet to a wet neighbor
% 3: outlet linked to another flooded region with the same elevation outlet
%    so not true oulet
% 4: shared outlet to a dry neighbor
% 5: shared outlet to a wet neighbor
% 6: shared outlet linked to another flooded region with the same elevation outlet
% 7: shared outlet surrounded by the cells with higher elevation. so not
%    true outlet.

% Set a temporary boundary of DEM for debugging
fXMin = 1; fXMax = nCols;
fYMin = 1; fYMax = mRows;

figure(4); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(fldRegID(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Flooded Region ID');

figure(5); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(subFldRegID(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Sub-flooded Region ID');

figure(6); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(subFldRegOutlet(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Sub-flooded Region Outlet');

%% Modify flow direction of the cells within each depression

% Assign flow direction to the cells in each depression region
[m2SDSNbrY,m2SDSNbrX ... % flow direction modified cells along the path to each regional minima
,mFlowDir_SubFldReg ... % flow direction modified cell within a sub-flooded region 
,mFlowDir_Saddle ... % flow direction modified cell on a saddle
,subFldRegTree ...
,fldRegInfo] ... % informaiton on flooded regions: ID, outlet index, number of sub-flooded region, outlet's elevation
    = AssignFlowDirInFldReg(m2SDSNbrY,m2SDSNbrX,subFldRegOutInfo,DEM ...
    ,slopeAllNbr,regionalMin,fldRegID,subFldRegID,sharedOutlet);

% display the tree DB of flooded region
disp(subFldRegTree.tostring);

%% Calculate upstream cells number of all cells within every depression

% A. Calculate the number of upstream cells out of flooded region
nUpstreamCellsWithFldReg = CalcUpstreamCellsWithFldReg(DEM,targetDrainage ...
    ,flood,m1SDSNbrY,m1SDSNbrX,fldRegID,nFldRegCells);

% B. Calculate flow accumulation within each flooded region
nUpstreamCells ...
    = CalcUpstreamCellsInSubFldReg(subFldRegTree,mFlowDir_SubFldReg ...
    ,fldRegID,nUpstreamCellsWithFldReg,fldRegInfo ...
    ,m2SDSNbrY,m2SDSNbrX,subFldRegID,DEM,regionalMin);   

% Visualization for the number of upstream cells
figure(7);
set(gcf,'Color',[1 1 1])

subplot(1,3,1)
imagesc(nUpstreamCellsWithFldReg);
title('Upstream Cells No. with Depressions');
axis image
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,2)
imagesc(nUpstreamCells);
title('Upstream Cells No.');
axis image
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,3)
imagesc(log(nUpstreamCells));
title('Upstream Cells No. [Log]');
axis image
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar