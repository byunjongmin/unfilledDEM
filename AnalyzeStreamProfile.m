function [upstreamAreaProf,slopePerDist,Rd,streamPath,distFromInit ...
    ,area_bin,slope_bin] ...
    = AnalyzeStreamProfile(initY,initX,endY,endX,chosenProfile ...
    ,nNbrForSmooth,nNbrForCalcSlope,initXforRd,endXforRd ...
    ,chosenSlopeForAreaSlope,contInterval,logBinSize...
    ,rawDEMwithBnd,nUpstreamCells,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX)
%
% function AnalyzeStreamProfile
% 
% 	Draw longitudinal profile, calculate relative slope according to given
%   window sizes, draw area-slope relationship and log-binned relationship
%
%	
% Output Variable
%
%	upstreamAreaProf
%   slopePerDist
%
% Version: 0.1.1 / 2019-08-13
%

%% Identify stream path based on the given initial and end point coordinates
[streamPath,distFromInit] = RecordStrPath(initY,initX,endY,endX ...
    ,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX);

% Stream longitudinal profile using raw DEM
profElev = rawDEMwithBnd(streamPath(:));
nCells = numel(streamPath);

% Draw figure
figure(21); clf;

subplot(3,1,1);
set(gcf, 'Color',[1,1,1]);
stairs(distFromInit,profElev);
grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(profElev),max(profElev)])

%% Generate stream profiles smoothed using different moving windows
% Note that, if you input multiple values for the number of considering
% neighbor cells, you can notice their differences 

nCNbr = numel(nNbrForSmooth); % number of cases
smoothedProfElev ...
    = SmoothingElev(nNbrForSmooth,distFromInit,rawDEMwithBnd,streamPath);

subplot(3,1,2)

cc = jet(nCNbr);
for ithLine = 1:nCNbr
    plot(distFromInit,smoothedProfElev(:,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
set(gca,'FontSize',13,'fontWeight','bold' ...
    ,'XMinorTick','on','YMinorTick','on')
title('Smoothed Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(profElev),max(profElev)])

%% Draw stream gradient

% Chose base elevation of profile
inputProfElev = smoothedProfElev;
nCNbrForSlope = numel(nNbrForCalcSlope);
% 너무 크게 구간을 잡으면 연산이 안됨.
[slopePerDist,slopePerDist_TF] ...
    = CalcSlope(nNbrForCalcSlope,distFromInit,inputProfElev,chosenProfile);

subplot(3,1,3)

cc = jet(nCNbrForSlope);
for ithLine = 1:nCNbrForSlope
    % exclude boundary nodes for slope calculation
    RdIdx = slopePerDist_TF(:,ithLine) == true;
    plot(distFromInit(RdIdx),slopePerDist(RdIdx,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
set(gca,'FontSize',13,'fontWeight','bold' ...
    ,'XMinorTick','on','YMinorTick','on')
title('Stream Gradient')
xlabel('Distance From Initiaion [m]')
ylabel('Slope')
xlim([0 distFromInit(end)])
ylim([min(slopePerDist(:,1)),max(slopePerDist(:,1))])

%% Identify Knickpoints: Relative Slopes (Rd)

Rd = nan(nCells,1); % relative slope
nValuesForRd = nan(nCells,1); % number of values for curve fitting
Rsquare = nan(nCells,1); % coefficient of determination

rangeX = initXforRd:endXforRd;
for ithCell = 1:nCells
    
    rXIdx = find(slopePerDist(ithCell,slopePerDist_TF(ithCell,rangeX) == true));
    if numel(rXIdx) > 1
        
        y = slopePerDist(ithCell,rangeX(rXIdx));
        x = nNbrForCalcSlope(rangeX(rXIdx)) .* 2 .* double(dX);
        [p,s] = polyfit(x,y,1); % curve fitting
        
        c = polyfit(x,y,1);
        if c(1) ~= p(1)
            error('error');
        end
            
        Rd(ithCell,1) = -p(1);
        nValuesForRd(ithCell,1) = numel(y);
        Rsquare(ithCell,1) = 1 - (s.normr / norm(y - mean(y)))^2;
    end

end

figure(22);

subplot(3,1,1)
plot(distFromInit,Rd);
grid on
set(gca,'FontSize',13,'fontWeight','bold' ...
    ,'XMinorTick','on','YMinorTick','on')
title('Relative Slope')
xlabel('Distance From Divide [m]')
ylabel('Relative Slope [m^-1]')
xlim([0 distFromInit(end)])
ylim([-1E-5,1E-5])

subplot(3,1,2)
plot(distFromInit,nValuesForRd);
grid on
set(gca,'FontSize',13,'fontWeight','bold' ...
    ,'XMinorTick','on','YMinorTick','on')
title('Number of values for curve fitting')
xlabel('Distance From Divide [m]')
ylabel('Number of values')
xlim([0 distFromInit(end)])

subplot(3,1,3)
plot(distFromInit,Rsquare);
grid on
set(gca,'FontSize',13,'fontWeight','bold' ...
    ,'XMinorTick','on','YMinorTick','on')
title('R^2')
xlabel('Distance From Divide [m]')
ylabel('R^2')
xlim([0 distFromInit(end)])

%% Simple statistics of slopePerDist for determining the number of cells
% for calculation of the local and trend slope

slopePerDistMin = zeros(1,nCNbrForSlope);
slopePerDistMax = zeros(1,nCNbrForSlope);
slopePerDistMean = zeros(1,nCNbrForSlope);
slopePerDistStd = zeros(1,nCNbrForSlope);

for i=1:nCNbrForSlope
    
    RdIdx = slopePerDist_TF(:,i) == true;
    slopePerDistMin(1,i) = min(slopePerDist(RdIdx,i));
    slopePerDistMax(1,i) = max(slopePerDist(RdIdx,i));
    slopePerDistMean(1,i) = mean(slopePerDist(RdIdx,i));
    slopePerDistStd(1,i) = std(slopePerDist(RdIdx,i));
    
end

figure(23); clf;

subplot(2,2,1)
scatter(nNbrForCalcSlope .* dX,slopePerDistMin);
grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Minimum Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,2)
scatter(nNbrForCalcSlope .* dX,slopePerDistMax);
grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Maximum Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,3)
scatter(nNbrForCalcSlope .* dX,slopePerDistMean);
grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Mean Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,4)
scatter(nNbrForCalcSlope .* dX,slopePerDistStd);
grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Standard Deviation Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Standard Deviation')


%% Draw corrected upstream area profiles on the stream paths of interest

% calculate upstream area
upstreamAreaProf = nUpstreamCells(streamPath(:)) .* dX .* dY;

figure(24); clf;
set(gcf, 'Color',[1,1,1]);

subplot(2,1,1)
plot(distFromInit,upstreamAreaProf);
set(gca,'YScale','log');
grid off
title('Upstream Contributing Area - Distance')
xlabel('Distance From Divide [m]')
ylabel('Area [m^2]')
xlim([0 distFromInit(end)])
ylim([min(upstreamAreaProf),max(upstreamAreaProf)])
set(gca,'FontSize',18);

% draw upslope area - slope relationship
subplot(2,1,2)
scatter(upstreamAreaProf,slopePerDist(:,chosenSlopeForAreaSlope) ...
    ,'Marker','*');
set(gca,'XScale','log','YScale','log');
grid off
title('Upstream Contributing Area - Slope')
xlabel('Area [m^2]')
ylabel('Slope')
xlim([upstreamAreaProf(1),upstreamAreaProf(end)])
set(gca,'FontSize',18);


% Distance - area curve

figure(25); clf;
set(gcf, 'Color',[1,1,1]);

plot(upstreamAreaProf./10^6,distFromInit)

grid on
set(gca,'FontSize',13,'fontWeight','bold' ...
    ,'XMinorTick','on','YMinorTick','on')

title('Upstream Contributing Area - Distance')
xlabel('Area [km^2]')
ylabel('Distance')
xlim([min(upstreamAreaProf(:)./10^6),max(upstreamAreaProf(:)./10^6)])
ylim([min(distFromInit(:)),max(distFromInit(:))])

%% Binned area-slope relationship: analyze slope-area relationship
% based on the vertically fixed interval method
% Note that errors in DEM produce many depressions and negative slope values
% along a stream profile. Therefore, this algorithm concentrates on
% making a stream profile without depression. Using this stream profile,
% it calculates new slope values and log binnined average slope values.
% By the way, using a smoothed stream proile as an input data makes analysis
% better. Because the smoothed stream profile has less depressions than
% raw stream profile data.

% Make a vertically fixed interval
ipElevMin = 10 * floor(0.1 * profElev(end));
ipElevMax = 10 * ceil(0.1 * profElev(1));
ipElev = ipElevMax:-contInterval:ipElevMin;

% Make a distance set interpolated from the correlation of elevation and
% distance. Note that, for avoiding errors from depression, depression are
% to be removed

tStrProfElev = profElev; % fixed elevation profile
% Remove depressions in an input stream profile.
nCells = numel(distFromInit);
for a = nCells:-1:2
    if tStrProfElev(a-1) <= tStrProfElev(a)
        tStrProfElev(a-1) = tStrProfElev(a) + 0.000001;
    end
end
% Interpolated distance from channel head
ipDistFromInit = interp1(tStrProfElev,distFromInit,ipElev);

% Draw stream profiles against vertically fixed interval and interpolated
% distance

figure(26); clf;

subplot(5,1,1)
% for debug
% plot(distFromInit,tStrProfElev,'o',ipDistFromInit,ipElev);
plot(ipDistFromInit,ipElev);
grid on
title('Longitudinal Stream Profile')
xlabel('Interpolated Distance From Channel Head [m]')
ylabel('Elevation [m]')
xlim([0 max(ipDistFromInit)])
ylim([min(ipElev),max(ipElev)])

% calculate new slope values using interpolated elevation and distance data
% note: basically new slope means local slope rather than 'trend slope'.
nCells2 = numel(ipElev);
newSlope = zeros(nCells2,1);
for a = 1
    newSlope(a) = (ipElev(a) - ipElev(a+1)) ...
        / (ipDistFromInit(a+1) - ipDistFromInit(a));
end
for a = 2:nCells2-1
    newSlope(a) = (ipElev(a-1) - ipElev(a+1)) ...
        / (ipDistFromInit(a+1) - ipDistFromInit(a-1));
end
for a = nCells2
    newSlope(a) = (ipElev(a-1) - ipElev(a)) ...
        / (ipDistFromInit(a) - ipDistFromInit(a-1));
end

% draw newly calculated slope vs distance relationship
subplot(5,1,2)
plot(ipDistFromInit,newSlope,'*')

% grid on
xlabel('Interpolated Distance From Channel Head [m]')
ylabel('Slope from Interpolated Elevation and Distance')
xlim([0 max(ipDistFromInit)])
ylim([min(newSlope),max(newSlope)])

% draw newly calculated upstream area vs elevation relationship
ipUpstreamArea = interp1(tStrProfElev,upstreamAreaProf,ipElev);

subplot(5,1,3)
plot(ipDistFromInit,ipUpstreamArea);
set(gca,'YScale','log');
grid on
title('Interpolated Upstream Area Profile')
xlabel('Interpolated Distance From Channel Head [m]')
ylabel('Interpolated Upstream Area [m^2]')
xlim([0 max(ipDistFromInit)])
ylim([min(ipUpstreamArea),max(ipUpstreamArea)])

% draw newly calculated upstream area vs slope relationship
subplot(5,1,4)
scatter(ipUpstreamArea,newSlope,'m+')
set(gca,'XScale','log','YScale','log');
grid on
title('Interpolated Upstream Area and Slope')
xlabel('Interpolated Upstream Area [m^2]')
ylabel('Slope from Interpolated Elevation and Distance')
xlim([min(ipUpstreamArea),max(ipUpstreamArea)])

% determine logarithmic-binned average slopes
minPower = floor(log10(min(ipUpstreamArea))/logBinSize) * logBinSize;
maxPower = ceil(log10(max(ipUpstreamArea))/logBinSize) * logBinSize;
diffPower = maxPower - minPower;

cellsNo3 = floor(diffPower/logBinSize);

tempArea = zeros(1,cellsNo3);
tempSlope = zeros(1,cellsNo3);

for a = 1:cellsNo3
    
    logAreaMin = minPower + (a-1) * logBinSize;
    logAreaMax = logAreaMin + logBinSize;
    
    tempArea(a) = 10^(logAreaMin + logBinSize * 0.5);
    
    j = find(ipUpstreamArea >= 10^logAreaMin ...
        & ipUpstreamArea < 10^logAreaMax);
    k = size(nonzeros(j));
    if k(1) > 0
        tempSlope(a) = 10^(mean(log10(newSlope(j))));
    end
end

% extract non-zero elements
area_bin = tempArea(find(tempSlope));
slope_bin = tempSlope(find(tempSlope));

subplot(5,1,5)
loglog(area_bin,slope_bin,'rs','MarkerFaceColor','r','MarkerSize',3)
% grid on
title('Logarithmically binned Upstream Area and Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([min(ipUpstreamArea),max(ipUpstreamArea)])