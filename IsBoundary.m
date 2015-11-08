function isTrue = IsBoundary(y,x,Y_TOP_BND,Y_BOTTOM_BND,X_LEFT_BND,X_RIGHT_BND)
%
% function IsBoundary
%
% =========================================================================
%> @section INTRO IsBoundary
%>
%> - �ԷµǴ� ���� ��ǥ(y,x)�� ���� �ܰ� ��迡 ��ġ�ϴ����� Ȯ���ϴ� �Լ�
%>
%> @version 0.8
%> @callgraph
%> @callergraph
%>
%> @retval isTrue       : �ԷµǴ� ���� ��ǥ�� �ܰ� ��迡 ��ġ�ϴ°��� ǥ���ϴ� ����
%>
%> @param y             : �ԷµǴ� ���� Y ��ǥ��
%> @param x             : �ԷµǴ� ���� X ��ǥ��
%> @param Y_TOP_BND     : ���� �ܰ� �� ��� Y ��ǥ��
%> @param Y_BOTTOM_BND  : ���� �ܰ� �Ʒ� ��� Y ��ǥ��
%> @param X_LEFT_BND    : ���� �ܰ� �� ��� X ��ǥ��
%> @param X_RIGHT_BND   : ���� �ܰ� �� ��� X ��ǥ��
% =========================================================================

if (x == X_RIGHT_BND) || (x==X_LEFT_BND)
    
    isTrue = true;

elseif (y == Y_BOTTOM_BND) || (y==Y_TOP_BND)
    
    isTrue = true;

else
    isTrue = false;

end