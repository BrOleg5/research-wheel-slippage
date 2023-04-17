function lgd = legend_translate(txt_dict, language, varargin)
% LEGEND_TRANSLATE Wrapper of ylabel function with legend translation

assert(isa(txt_dict, "containers.Map"), "txt_dict must be containers.Map");
mustBeTextScalar(language);
assert(isKey(txt_dict, language), "language must be key of txt_dict Map container.");

if(isempty(varargin))
    lgd = legend(txt_dict(language));
else
    lgd = legend(txt_dict(language), varargin{:});
end
end