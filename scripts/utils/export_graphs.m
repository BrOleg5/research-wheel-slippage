function export_graphs(obj, filename, formats, options)
% EXPORT_GRAPHS Save graphics in several formats
% 
%   See also EXPORTGRAPHICS

arguments
    obj {iptcheckhandle(obj, {'axes', 'figure'}, 'export_graphs', 'obj', 1)},
    filename {mustBeText, mustBeNonzeroLengthText, mustBeNonempty},
    formats {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, mustBeMember(formats, ...
        ["jpg", "jpeg", "png", "tif", "tiff", "gif", "eps", "emf", "pdf"])} = "pdf"
    options.ContentType {mustBeTextScalar, mustBeMember(options.ContentType, ["auto", "vector", "image"])} ...
        = "vector",
    options.Resolution (1, 1) {mustBeNumeric, mustBeInteger, mustBePositive} = 150,
    options.BackgroundColor = "white"
    options.Append (1, 1) {mustBeNumericOrLogical} = false,
    options.Colorspace {mustBeTextScalar, mustBeMember(options.Colorspace, ["rgb", "gray", "cmyk"])} ...
        = "rgb"
end

background_color_rgb = validatecolor(options.BackgroundColor, 'multiple');
for ext = formats
    filename_with_ext = strcat(filename, ".", ext);
    exportgraphics(obj, filename_with_ext, 'ContentType', options.ContentType, ...
                   'Resolution', options.Resolution, 'BackgroundColor', background_color_rgb, ...
                   'Append', options.Append, 'Colorspace', options.Colorspace);
end
end

