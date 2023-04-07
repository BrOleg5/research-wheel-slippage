function plot_movement_delta_angle_graphs()
% PLOT_MOVEMENT_DELTA_ANGLE_GRAPHS  Entry point

export_graph = true;
input_dataset_file = "output_files/datasets/for_research_movement_direction/one_surface_type/averaged_data.csv";
output_dir = "output_files/graphs/for_research_movement_direction";

check_dependency(input_dataset_file, @average_one_surface_type_dataset);

opts = detectImportOptions(input_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
T = readtable(input_dataset_file, opts);

idx = doublecmp(T.speedamp , 0.1);
T = T(idx, :);

plot_func = localfunctions;
for i = 1:length(plot_func)
    plot_func{i}(T, 'ExportGraphs', export_graph, 'ExportGraphFolder', output_dir);
end
end
%% Local functions
function plot_delta_angle_between_set_point_and_camera_direction(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "vx", "vy"]);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    resT = varfun(@mean, sample_T, 'GroupingVariables', 'movedir', 'InputVariables', ["vx", "vy"]);
    resT.camdir = wrapTo360(rad2deg(atan2(resT.mean_vy, resT.mean_vx)));
    resT.deltaAngle = wrapTo180(resT.movedir - resT.camdir);
    fig = figure("Name", surf_str);
    plot(resT.movedir, resT.deltaAngle, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
         "MarkerEdgeColor", "black");
    grid on;
    ax = gca;
    ax.XTick = 0:10:350;
    ax.XLimitMethod = 'tight';
    ax.YLimitMethod = 'padded';
    xlabel("Movement direction, ^\circ", "Interpreter", "tex");
    ylabel("Delta angle, ^\circ", "Interpreter", "tex");
    ax.FontSize = 12;
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "delta_movement_direction", ...
                              "between_set_point_and_camera");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, surf_str);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end

fig = figure("Name", "Plot delta angle between set point and camera movement direction");
grid on;
hold on;
surface_color_dict = get_surface_color_dict();
for i = 1:length(surftype)
    idx = T.surftype == surftype(i);
    sample_T = T(idx, :);
    resT = varfun(@mean, sample_T, 'GroupingVariables', 'movedir', 'InputVariables', ["vx", "vy"]);
    resT.camdir = wrapTo360(rad2deg(atan2(resT.mean_vy, resT.mean_vx)));
    resT.deltaAngle = wrapTo180(resT.movedir - resT.camdir);
    surf_str = string(surftype(i));
    plot(resT.movedir, resT.deltaAngle, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
         "MarkerEdgeColor", surface_color_dict(surf_str));
end
ax = gca;
ax.XTick = 0:10:350;
ax.XLimitMethod = 'tight';
ax.YLimitMethod = 'padded';
xlabel("Movement direction, ^\circ", "Interpreter", "tex");
ylabel("Delta angle, ^\circ", "Interpreter", "tex");
ax.FontSize = 12;
legend(string(surftype), 'Location', 'best', 'FontSize', 12);

if(options.ExportGraphs)
    output_dir = fullfile(options.ExportGraphFolder, "delta_movement_direction", ...
                          "between_set_point_and_camera");
    if(~isfolder(output_dir))
        mkdir(output_dir);
    end
    file_name = fullfile(output_dir, "common_graph");
    export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                  'BackgroundColor', 'white', 'Colorspace', 'rgb');
    close(fig);
end
end

function plot_delta_angle_between_set_point_and_current_direction(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "xcur", "ycur"]);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    resT = varfun(@mean, sample_T, 'GroupingVariables', 'movedir', 'InputVariables', ["xcur", "ycur"]);
    resT.curdir = wrapTo360(rad2deg(atan2(resT.mean_ycur, resT.mean_xcur)));
    resT.deltaAngle = wrapTo180(resT.movedir - resT.curdir);
    fig = figure("Name", surf_str);
    plot(resT.movedir, resT.deltaAngle, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
         "MarkerEdgeColor", "black");
    grid on;
    ax = gca;
    ax.XTick = 0:10:350;
    ax.XLimitMethod = 'tight';
    ax.YLimitMethod = 'padded';
    xlabel("Movement direction, ^\circ", "Interpreter", "tex");
    ylabel("Delta angle, ^\circ", "Interpreter", "tex");
    ax.FontSize = 12;
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "delta_movement_direction", ...
                              "between_set_point_and_current");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, surf_str);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end

function plot_delta_angle_between_camera_and_current_direction(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "vx", "vy", "xcur", "ycur"]);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);

    resT = varfun(@mean, sample_T, 'GroupingVariables', 'movedir', ...
                  'InputVariables', ["xcur", "ycur", "vx", "vy"]);
    resT.curdir = wrapTo360(rad2deg(atan2(resT.mean_ycur, resT.mean_xcur)));
    resT.camdir = wrapTo360(rad2deg(atan2(resT.mean_vy, resT.mean_vx)));
    resT.deltaAngle = wrapTo180(resT.camdir - resT.curdir);
    fig = figure("Name", surf_str);
    plot(resT.movedir, resT.deltaAngle, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
         "MarkerEdgeColor", "black");
    grid on;
    ax = gca;
    ax.XTick = 0:10:350;
    ax.XLimitMethod = 'tight';
    ax.YLimitMethod = 'padded';
    xlabel("Movement direction, ^\circ", "Interpreter", "tex");
    ylabel("Delta angle, ^\circ", "Interpreter", "tex");
    ax.FontSize = 12;
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "delta_movement_direction", ...
                              "between_camera_and_current");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, surf_str);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end

fig = figure("Name", "Plot delta angle between camera movement direction and current direction");
grid on;
hold on;
surface_color_dict = get_surface_color_dict();
for i = 1:length(surftype)
    idx = T.surftype == surftype(i);
    sample_T = T(idx, :);

    resT = varfun(@mean, sample_T, 'GroupingVariables', 'movedir', ...
                  'InputVariables', ["xcur", "ycur", "vx", "vy"]);
    resT.curdir = wrapTo360(rad2deg(atan2(resT.mean_ycur, resT.mean_xcur)));
    resT.camdir = wrapTo360(rad2deg(atan2(resT.mean_vy, resT.mean_vx)));
    resT.deltaAngle = wrapTo180(resT.camdir - resT.curdir);
    surf_str = string(surftype(i));
    plot(resT.movedir, resT.deltaAngle, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
         "MarkerEdgeColor", surface_color_dict(surf_str));
end
ax = gca;
ax.XTick = 0:10:350;
ax.XLimitMethod = 'tight';
ax.YLimitMethod = 'padded';
xlabel("Movement direction, ^\circ", "Interpreter", "tex");
ylabel("Delta angle, ^\circ", "Interpreter", "tex");
ax.FontSize = 12;
legend(string(surftype), 'Location', 'best', 'FontSize', 12);
if(options.ExportGraphs)
    output_dir = fullfile(options.ExportGraphFolder, "delta_movement_direction", ...
                          "between_camera_and_current");
    if(~isfolder(output_dir))
        mkdir(output_dir);
    end
    file_name = fullfile(output_dir, "common_graph");
    export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                  'BackgroundColor', 'white', 'Colorspace', 'rgb');
    close(fig);
end
end