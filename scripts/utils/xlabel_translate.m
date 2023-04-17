function hh = xlabel_translate(txt_dict, language, varargin)
% XLABEL_TRANSLATE Wrapper of xlabel function with label translation

hh = label_translate_impl(txt_dict, language, @xlabel, varargin{:});
end