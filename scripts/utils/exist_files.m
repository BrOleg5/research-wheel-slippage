function res = exist_files(file)
% EXIST_FILES  Check existence of files.
%   res = EXIST_FILES(dirs);
% See also EXIST
arguments
    file {mustBeText}
end
res = false(size(file));
for i = 1:length(file)
   res(i) = isfile(file(i));
end
end