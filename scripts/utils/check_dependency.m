function tf = check_dependency(dependency, func_handle)
% CHECK_DEPENDENCIES Check for existance of files and if its are not, then
% run function by handle.

arguments
    dependency {mustBeText, mustBeNonempty, mustBeNonzeroLengthText},
    func_handle {validateattributes(func_handle, {'cell', 'function_handle'}, {'vector', 'nonempty'})}
end
existance = exist_files(dependency);
if(any(~existance))
    fprintf("'%s' doesn't exist.\n", dependency(~existance));
    for i = 1:length(func_handle)
        if(iscell(func_handle))
            func = func_handle{i};
        elseif(length(func_handle) == 1)
            func = func_handle;
        end
        assert(isa(func, "function_handle"), "func isn't function handle.");
        s = functions(func);
        fprintf("Run %s() function.\n", s.function)
        func();
    end
    tf = false;
    return;
else
    disp("All dependencies exist.");
    tf = true;
end
end

