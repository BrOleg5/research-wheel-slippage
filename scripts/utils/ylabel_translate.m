function hh = ylabel_translate(txt_dict, language, varargin)
% YLABEL_TRANSLATE Wrapper of ylabel function with label translation

hh = label_translate_impl(txt_dict, language, @ylabel, varargin{:});
end