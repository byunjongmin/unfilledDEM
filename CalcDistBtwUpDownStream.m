function streamNetWithInterval ...
    = CalcDistBtwUpDownStream(streamNet,streamNetElement ...
    ,streamNetNodeInfo,dY,dX)
%
% function
% Calculate distances btw every upstream and downstream

% Output Variale
streamNetWithInterval = tree(streamNet,'clear');

% Constants
[mRows,nCols] = size(streamNetElement);
NOT_CALCULATED = 0;
ROOT_NODE = 1;

% Find all indices with leaf node
HEAD_NODE = 4;
[leafNodeY,leafNodeX] = find(streamNetElement == HEAD_NODE);
[nLeafNode,tmp] = size(leafNodeY);

for ithPath = 1:nLeafNode

    % Identify the ith leaf node in the stream network
    ithLeafNodeIdx ...
        = sub2ind([mRows,nCols],leafNodeY(ithPath),leafNodeX(ithPath));
    ithLeafNode = find(streamNetNodeInfo(:,2) == ithLeafNodeIdx);
    
    upStreamNode = ithLeafNode;
    
    pathNotEnd = true; 
    while pathNotEnd
    
        % 상류 셀과의 거리를 구하지 않았다면, 아래의 과정을 거쳐 거리를 구함
        if numel(streamNetWithInterval.get(upStreamNode)) == NOT_CALCULATED

            downStreamNode = streamNetWithInterval.getparent(upStreamNode);
            
            % upStream Node가 outlet일 경우는 다음 leaf node로 이동함
            if upStreamNode ~= ROOT_NODE

                % Coordinate of upstream cell
                upStreamIdx = streamNet.get(upStreamNode);
                [upStreamY,upStreamX] = ind2sub([mRows,nCols],upStreamIdx);
                
                % Coordinate of downstream cell
                downStreamIdx = streamNet.get(downStreamNode);
                [downStreamY,downStreamX] ...
                    = ind2sub([mRows,nCols],downStreamIdx);

                % Calculate the distance btw them
                distanceBtwThem ... 
                    = sqrt(double(((downStreamY-upStreamY)*dY)^2 ...
                    +((downStreamX-upStreamX)*dX)^2));

                % Insert the distance value
                streamNetWithInterval ...
                = streamNetWithInterval.set(upStreamNode,distanceBtwThem);
                
                % Change upstream node with downstream node
                upStreamNode = downStreamNode;

            % 상류 셀과의 거리를 구했다면, 다음 leaf node로 이동함
            else
                
                pathNotEnd = false;
                
            end
            
        else
            
            pathNotEnd = false;
                
        end
    end
end
