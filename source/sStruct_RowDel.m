% This function removes rows, from fields in a structure
% It's assumed that all fields have the same  number oof rows.
% "IndextoRemove" is a binary vector of the same size as the filed, or the
% indices of the rows to be removed.
% 
% Written by Pejman Mohammadi
% pejman.mohammadi@bsse.ethz.ch
%

function Struct = sStruct_RowDel(Struct, IndextoRemove)

Fs = fieldnames(Struct);
for i = 1:length(Fs)
    Struct.(Fs{i})(IndextoRemove,:) = [];
end
end