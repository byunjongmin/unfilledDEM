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
    
        % ��� ������ �Ÿ��� ������ �ʾҴٸ�, �Ʒ��� ������ ���� �Ÿ��� ����
        if numel(streamNetWithInterval.get(upStreamNode)) == NOT_CALCULATED

            downStreamNode = streamNetWithInterval.getparent(upStreamNode);
            
            % upStream Node�� outlet�� ���� ���� leaf node�� �̵���
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

            % ��� ������ �Ÿ��� ���ߴٸ�, ���� leaf node�� �̵���
            else
                
                pathNotEnd = false;
                
            end
            
        else
            
            pathNotEnd = false;
                
        end
    end
end
