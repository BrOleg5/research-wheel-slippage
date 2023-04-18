function ru_surface = translate_surface(en_surface)
% TRANSLATE_SURFACE

arguments
    en_surface {mustBeVector, mustBeText, mustBeNonzeroLengthText, mustBeNonempty}
end
surface_dict = containers.Map(["table", "gray", "green"], ["гладкая", "серая", "зелёная"]);
assert(all(isKey(surface_dict, cellstr(en_surface))), ...
       "All elements of en_surface must be key of surface_dict Map container.");
ru_surface = string(values(surface_dict, cellstr(en_surface)));
end