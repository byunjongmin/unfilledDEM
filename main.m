function main
% @file main.m
% @brief Function to extract stream longitudinal profiles from unfilled
% DEMs
%
% @version 0.1.1. / 2015-11-16
% @author Jongmin Byun
%==========================================================================

%% Load DEM

% Constants
INPUT_DIR = '../data/input';
OUTPUT_DIR = '../data/output';
DEMFileName = 'tareaR50m.tif';
DEMFilePath = fullfile(INPUT_DIR,DEMFileName);
[DEM,R] = geotiffread(DEMFilePath);

% DEM basic properties
mRows = R.RasterSize(1,1);
nCols = R.RasterSize(1,2);
dX = R.DeltaX;
dY = -R.DeltaY;
DEM = double(DEM);

%% For debug for a part of DEM
tYMin = 311; tYMax = 370;
tXMin = 146; tXMax = 195;

DEM = DEM(tYMin:tYMax,tXMin:tXMax);
[mRows,nCols] = size(DEM);

%% Define target drainage
% Note that you should check the location of the outlet of the test domain.
% It should be located on the boundary of a drainage. If it is within the
% drainage, you should modify the DEM

% % for the domain surrounded by null
% nanMask = (DEM == 32767);
% DEM(nanMask) = inf;
% DEMArea = ~nanMask;

% for the domain filled with only elevations
nanMask = true(mRows,nCols);
nanMask(2:mRows-1,2:nCols-1) = false;
DEMArea = ~nanMask;
DEM(nanMask) = inf;

% extract the boundary of the target drainage
s = strel('square',3); % structural element when eroding image
eDEMArea = imerode(DEMArea,s); 
DEMAreaBnd = DEMArea & ~eDEMArea;
DEMAreaBndIdx = find(DEMAreaBnd);

% identify the coordinate of the outlet
DEMAreaBndElev = DEM(DEMAreaBndIdx);
[~,minElevIdx] = min(DEMAreaBndElev);
outletIdx = DEMAreaBndIdx(minElevIdx); % outlet on the boundary of drainage
[outletY,outletX] = ind2sub([mRows,nCols],outletIdx);

% Note that the elevation of the main outlet should be the lowest.
DEM(outletY,outletX) = min(DEM(:)) - 0.1;

% define target drainage
targetDrainage = (~nanMask);
targetDrainage(outletY,outletX) = false;

%% Main body

%% Smooth the imported DEM, until the flat is removed
orgDEM = DEM; % original DEM

% Assign flow directions to the DEM using D8 algorithm
% Note that, to use CalcSDSFlow function, elevation of outer region of the
% target drainage should be inf.
[steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
    = CalcSDSFlow(DEM,dX,dY);

% Process flat and sinks for flow to continue to move
flatRegMap = ProcessFlat(DEM,targetDrainage,slopeAllNbr);

rT = 0;
afterNFlat = inf;
while afterNFlat > 0

    fprintf('The remaining flat cell is %4.2f\n',afterNFlat);
    % smooth only flat. To prevent an infinite loop, change the size of
    % moving window when flat region does not reduce 
    if rT > 5
        nNbr = 5;
    else
        nNbr = 3;
    end
    h = ones(nNbr,nNbr)/nNbr^2;
    smtDEM = filter2(h,DEM);
    DEM(flatRegMap == true) = smtDEM(flatRegMap == true);
    [steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
        = CalcSDSFlow(DEM,dX,dY);
    flatRegMap = ProcessFlat(DEM,targetDrainage,slopeAllNbr);
    
    oldNFlat = afterNFlat;
    afterNFlat = numel(find(flatRegMap == true));
    if oldNFlat == afterNFlat
        rT = rT + 1;
    else
        rT = 0;
    end
    
end

% for debug
figure(1);
subplot(1,2,1);
imagesc(flatRegMap);
set(gca,'DataAspectRatio',[1 1 1]);
title('Distribution of Flat Region');

subplot(1,2,2);
diffDEM = orgDEM - DEM; % for debug
imagesc(diffDEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Difference in Elevation after Smoothing');

% %% For debug for a whole DEM
% INPUT_DIR = '../data/input';
% dataFileName = 'b_PSink_2015-09-15.mat';
% dataFilePath = fullfile(INPUT_DIR,dataFileName);
% load(dataFilePath);

%% Remove isolated area

DEM_BW = false(mRows,nCols); % binary of DEM
DEM_BW(~isinf(DEM)) = true;
CC = bwconncomp(DEM_BW); % identify connected components of valid cells
if CC.NumObjects > 1
    DEM(CC.PixelIdxList{2:end}) = inf; % remove the connected components except for the target area
    targetDrainage(CC.PixelIdxList{2:end}) = false;
end

%% Identify depressions and their outlets, then modify each depression
% outlet's flow direction to go downstream when flows are overspilled
% over the depression

% Note : To use ProcessSink function, make a mask excepting the main outlet
% and outer region of the target drainage
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
    
% % For debug
% INPUT_DIR = '../data/input';
% dataFileName = 'a_PSink_2015-11-09.mat';
% dataFilePath = fullfile(INPUT_DIR,dataFileName);
% load(dataFilePath);
    
% frequency distribution of the types of sub-flooded region's outlet
figure(2); clf;
set(gcf, 'Color',[1,1,1]);
h = histogram(subFldRegOutlet(subFldRegOutlet > 0));
xlabel('Type of Sub-flooded Region Outlet');
ylabel('Frequency');
grid on

% set boundary
fXMin = 1; fXMax = nCols;
fYMin = 1; fYMax = mRows;
figure(3)

subplot(2,2,1)
imagesc(DEM(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('DEM');

subplot(2,2,2)
imagesc(fldRegID(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Flooded Region ID');

subplot(2,2,3)
imagesc(subFldRegID(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Sub-flooded Region ID');

subplot(2,2,4)
imagesc(subFldRegOutlet(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Sub-flooded Region Outlet');

%% Modify flow direction of the cells within each depression

% Assign flow direction to the cells in each depression region
[m2SDSNbrY,m2SDSNbrX ... % modified for flooded region's outlet
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
    ,mFlowDir_Saddle,fldRegID,nUpstreamCellsWithFldReg,fldRegInfo ...
    ,m2SDSNbrY,m2SDSNbrX,subFldRegID,DEM,regionalMin);   

% Visualization for the number of upstream cells
figure(5);
set(gcf,'Color',[1 1 1])

subplot(1,3,1)
imagesc(nUpstreamCellsWithFldReg(2:end-1,2:end-1));
title('Upstream Cells No. with Depressions');
axis image
set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,2)
imagesc(nUpstreamCells(2:end-1,2:end-1));
title('Upstream Cells No.');
axis image
set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,3)
imagesc(log(nUpstreamCells(2:end-1,2:end-1)));
title('Upstream Cells No. [Log]');
axis image
set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

%% draw a stream longitudinal profile on the interesting stream path

% record stream path using the input initial and end point coordinates
initY = 22; initX = 39;
endY = 29; endX = 46;
[streamPath,distFromInit] = RecordStrPath(initY,initX,endY,endX ...
    ,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX);

% stream longitudinal profile using raw DEM
elev = DEM(streamPath(:));

% draw figure
figure(4); clf;
set(gcf, 'Color',[1,1,1]);
stairs(distFromInit,elev);

% grid on
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(elev),max(elev)])

% % c. export the number of upstream cells for ArcGIS
% nUpstreamCellsMap = nUpstreamCells;
% tFilePath = fullfile(OUTPUT_DIR,'nUpstreamCellsMap');
% ExportRasterAsArcGrid(tFilePath,nUpstreamCellsMap,dY,xllcorner(1),yllcorner(1));

%% Make stream network
% 유역면적이 일정면적 이상인 셀을 하천으로 정의함
% 유역 면적
upstreamArea = upstreamCellsNo * double(dX * dY);

% 하천 정의
chanInitArea = 6000; % 하천 시작 임계점
channel = upstreamArea >= chanInitArea;

% 하천인 셀에 한해서 stream network을 구함
% * 주의: 유출구가 하나 이상일 경우도 있음

% Make the database with tree structure for the stream network

% Note: Eliminate the flow direction of the outlet cell
SDSNbrY(outletY,outletX) = 0;
SDSNbrX(outletY,outletX) = 0;

% Note: contributing upstream area 가 다음 비율을 넘을 때만 stream network로
% 인정함
% criticalRatio

[streamNet,streamNetElement,streamNetNodeInfo] ...
    = MakeStreamNetwork(dX,dY,channel,outletY,outletX,SDSNbrY,SDSNbrX ...
    ,floodedRegionIndex,floodedRegionCellsNo,criticalRatio ...
    ,upstreamCellsNo,flatRegionPathInfo,flatRegionMap);

% display the database with tree structure
% disp(streamNet.tostring);

%% Make stream longitudinal profiles and analyze them

% Make a stream network with elevation
streamNetWithElev = tree(streamNet,'clear');

iterator = streamNetWithElev.breadthfirstiterator;		

for i = iterator
    
    ithIdx = streamNet.get(i);
    streamNetWithElev = streamNetWithElev.set(i,DEM(ithIdx));

end

% Make a stream network with every distances btw upstream and downstream
streamNetWithInterval ...
    = CalcDistBtwUpDownStream(streamNet,streamNetElement ...
    ,streamNetNodeInfo,dY,dX);

% Make a stream network with flooded region index
streamNetWithFldReg = tree(streamNet,'clear');

iterator = streamNetWithFldReg.breadthfirstiterator;		
for i = iterator
    
    ithIdx = streamNet.get(i);
    streamNetWithFldReg ...
        = streamNetWithFldReg.set(i,floodedRegionIndex(ithIdx));

end

% Identify major stream paths and calculate gradients of them and make a
% stream gradient map

% Input Variables
% Settting the number of nodes from upstream cell to downstream cell
noNodeBtwUpDown = 1;

[streamPathInfo,streamNetDepth,streamGradMap,streamPathMap] ...
    = IdentifyStreamPathInfo(noNodeBtwUpDown,streamNet,streamNetElement ...
    ,streamNetNodeInfo,streamNetWithInterval,DEM,upstreamCellsNo);


%% Draw stream longitudinal profile
[nPath,tmp] = size(streamPathInfo);

largestStreamLength = max(max(cell2mat(streamPathInfo(:,3))));

figure;
cc = jet(nPath);
legendNames = cell(nPath,1);
for ithPath = 1:nPath

    elev = streamPathInfo{ithPath,5};
    distanceFromDivide = streamPathInfo{ithPath,6};
    
    % 가장 긴 하천종단곡선의 거리에서 현 하천종단곡선의 거리를 뺌.
    distBtwEndOfProfileToMainOutlet ...
        = largestStreamLength - distanceFromDivide(end);

    plot(distanceFromDivide + distBtwEndOfProfileToMainOutlet ...
        ,elev,'color',cc(ithPath,:));
    
    legendNames{ithPath,1} = num2str(streamPathInfo{ithPath,1});
    
    hold on  
        
end

% legend(legendNames);
    
grid on
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')


%% Draw slope - area graph
[nPath,tmp] = size(streamPathInfo);

figure;
cc = jet(nPath);
legendNames = cell(nPath,1);
for ithPath = 1:nPath

    upstreamCellsNoP = streamPathInfo{ithPath,8};
    upstreamAreaP = upstreamCellsNoP * double(dX * dY);
    streamGradient = streamPathInfo{ithPath,7};

    
    % 가장 긴 하천종단곡선의 거리에서 현 하천종단곡선의 거리를 뺌.

    scatter(upstreamAreaP,streamGradient,'Marker','*','MarkerFaceColor',cc(ithPath,:));
        
    legendNames{ithPath,1} = num2str(streamPathInfo{ithPath,1});
    
    hold on  
        
end

% legend(legendNames);
set(gca,'XScale','log','YScale','log');
grid on
title('Area - Slope')
xlabel('Upstream Area')
ylabel('Slope')