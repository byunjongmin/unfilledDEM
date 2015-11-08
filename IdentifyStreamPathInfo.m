function [streamPathInfo,streamNetDepth,streamGradMap,streamPathMap] ...
    = IdentifyStreamPathInfo(noNodeBtwUpDown,streamNet,streamNetElement ...
    ,streamNetNodeInfo,streamNetWithInterval,DEM,upstreamCellsNo)
%
% function
% Identify major stream paths and calculate gradients of them and make a
% stream gradient map

% Constants
[mRows,nCols] = size(streamNetElement);
ROOT_NODE = 1;

% Output Variables
% Every stream path information
% Note: terminal leaf node index, the number of all path nodes, distance
% from divide at a terminal leaf node, indicies on a path, elevation on a
% path, distance from divide on a path, gradient on a path
streamPathInfo = [];
% Depth pf a stream network tree data structure
streamNetDepth = streamNet.depthtree;
streamGradMap = nan(mRows,nCols);

% Find all indices with leaf node
HEAD_NODE = 4;
[leafNodeY,leafNodeX] = find(streamNetElement == HEAD_NODE);
[nLeafNode,tmp] = size(leafNodeY);

% Iterate
for ithPath = 1:nLeafNode
    
    sprintf('%d/%d in IdentifyStreamPathInfo function\n',ithPath,nLeafNode) % for debug
    
    % Identify the ith leaf node in the stream network
    ithLeafNodeIdx ...
        = sub2ind([mRows,nCols],leafNodeY(ithPath),leafNodeX(ithPath));
    ithLeafNode = find(streamNetNodeInfo(:,2) == ithLeafNodeIdx);
    
    % Calculate the ith stream path node number    
    ithStreamPathNodeNo ...
        = ceil(streamNetDepth.get(ithLeafNode) / noNodeBtwUpDown);
    
    % Initialize ith stream path information variables
    ithStreamPathIdx = zeros(1,ithStreamPathNodeNo);    
    ithStreamPathElev = zeros(1,ithStreamPathNodeNo);    
    ithStreamPathGradient = zeros(1,ithStreamPathNodeNo);    
    ithStreamPathDistBtwUpDown = zeros(1,ithStreamPathNodeNo);    
    ithStreamPathDistFromDivide = zeros(1,ithStreamPathNodeNo);
    ithStreamPathUpStreamCellsNo =  zeros(1,ithStreamPathNodeNo);
    
    % Identify the information of the ith stream profile beginning from the
    % channel head (or terminal leaf node in a stream network)
    upStreamNode = ithLeafNode;
    for ithNode = 1:ithStreamPathNodeNo
        
        % Upstream index
        upStreamIdx = streamNet.get(upStreamNode);
        ithStreamPathIdx(ithNode) = upStreamIdx;
    
        % Upstream elevation
        upStreamElev = DEM(upStreamIdx);
        ithStreamPathElev(ithNode) = upStreamElev;
        
        % Distance from the upstream cell to the pre-determined downstream
        % cell
        passedNodeNo = 0;
        accDistanceToDownStream = 0;
        
        tempDownStreamNode = zeros(1,noNodeBtwUpDown);
        
        while passedNodeNo ~= noNodeBtwUpDown
        
            % If upstream node is not a main outlet, follows the belows
            if upStreamNode ~= ROOT_NODE
            
                % Each distance
                distToDownStream = streamNetWithInterval.get(upStreamNode);
                % Accumlated distance
                accDistanceToDownStream ...
                    = accDistanceToDownStream + distToDownStream;
                % Check the passed node number
                passedNodeNo = passedNodeNo + 1;
                % Identify the downstream node
                downStreamNode = streamNet.getparent(upStreamNode);
                
                tempDownStreamNode(passedNodeNo) = downStreamNode;
                
                % Replace the upstream node
                upStreamNode = downStreamNode;
            
            % If upstream node approach a main outlet, end the loop
            else
                
                passedNodeNo = noNodeBtwUpDown;
                
            end            
        end
        
        % Distance from the upstream node to the pre-determined downstream
        % node
        ithStreamPathDistBtwUpDown(ithNode) = accDistanceToDownStream;
        
        % Downstream node index and elevation
        downStreamIdx = streamNet.get(downStreamNode);
        downStreamElev = DEM(downStreamIdx);
            
        % Stream gradient btw upstream node and downstream node
        ithStreamPathGradient(ithNode) = (upStreamElev - downStreamElev) ...
            / accDistanceToDownStream;
        
        % Make a stream gradient map
        passedNodeNo = 0;
        while passedNodeNo ~= noNodeBtwUpDown
            
            passedNodeNo = passedNodeNo + 1;
            
            if tempDownStreamNode(passedNodeNo) ~= 0                
                
                tempDownStreamIdx ...
                    = streamNet.get(tempDownStreamNode(passedNodeNo));
                if isnan(streamGradMap(tempDownStreamIdx))
                    
                    streamGradMap(tempDownStreamIdx) ...
                        = ithStreamPathGradient(ithNode);
                    
                else
                    
                    if ithStreamPathGradient(ithNode) ...
                            > streamGradMap(tempDownStreamIdx)
                        
                        streamGradMap(tempDownStreamIdx) ...
                            = ithStreamPathGradient(ithNode);
                    
                    end
                end
            else
                
                passedNodeNo = noNodeBtwUpDown;
                
            end           
        end
        
        % Distance from divide
        if ithNode ~= 1
        
            ithStreamPathDistFromDivide(ithNode) ...
                = ithStreamPathDistFromDivide(ithNode-1) ...
                + accDistanceToDownStream;
            
        else
            
            ithStreamPathDistFromDivide(ithNode) = accDistanceToDownStream;
            
        end
        
        % Upstream Area
        ithStreamPathUpStreamCellsNo(ithNode) = upstreamCellsNo(upStreamIdx);
    
    end
    
    % Record the information of the ith stream path
    streamPathInfo = [streamPathInfo;{ithLeafNodeIdx} ...
        ,{ithStreamPathNodeNo} ...
        ,{ithStreamPathDistFromDivide(ithStreamPathNodeNo)} ...
        ,{ithStreamPathIdx},{ithStreamPathElev} ...
        ,{ithStreamPathDistFromDivide},{ithStreamPathGradient} ...
        ,{ithStreamPathUpStreamCellsNo}];
        
end

% Make stream path map
streamPathMap = zeros(mRows,nCols) - 10;

for ithPath = 1:nLeafNode
    
    % Identify the ith leaf node in the stream network
    ithLeafNodeIdx ...
        = sub2ind([mRows,nCols],leafNodeY(ithPath),leafNodeX(ithPath));
    ithLeafNode = find(streamNetNodeInfo(:,2) == ithLeafNodeIdx);
    
    upStreamNode = ithLeafNode;
    while upStreamNode ~= ROOT_NODE
        
        % Upstream Nbr index
        upStreamIdx = streamNet.get(upStreamNode);
        
        if streamPathMap(upStreamIdx) == -10
            
            streamPathMap(upStreamIdx) = ithPath;
            
        end
        
        % Identify the downstream node
        downStreamNode = streamNet.getparent(upStreamNode);
        
        % Replace downstream with upstream
        upStreamNode = downStreamNode;
        
    end
    
end