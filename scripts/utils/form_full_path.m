function full_path = form_full_path(paths)
% FORM_FULL_PATH Form full path of files and directories from dir() output struct.
%   full_path = FORM_FULL_PATH(paths) Form full path of files and
%   directories from path struct. Return full_path string array.

arguments
    paths struct {mustBeNonempty}
end

n = length(paths);
full_path = strings(1, n);
for i = 1:n
    full_path(i) = fullfile(paths(i).folder, paths(i).name);
end
end

