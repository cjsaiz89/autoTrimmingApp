%======================================================
% This function insures that all of the indices in retIndex
% exist in the input array(s)
%======================================================
function [retIndex] = valIndex(index, array)
    retIndex = index(index <= min((length(array))));
end
