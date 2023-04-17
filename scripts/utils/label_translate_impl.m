function hh = label_translate_impl(txt_dict, language, label_func, varargin)
% LABEL_TRANSLATE_IMPL Wrapper of label functions (xlabel, ylabel, zlabel) with label translation

assert(isa(txt_dict, "containers.Map"), "txt_dict must be containers.Map");
mustBeTextScalar(language);
assert(isKey(txt_dict, language), "language must be key of txt_dict Map container.");
assert(isa(label_func, "function_handle"), "label_func must be functio handler");
s = functions(label_func);
assert(any(strcmp(s.function, ["xlabel", "ylabel", "zlabel"])), "label_func must be function handler of xlabel, ylabel or zlabel");

if(isempty(varargin))
    hh = label_func(txt_dict(language));
else
    hh = label_func(txt_dict(language), varargin{:});
end
end