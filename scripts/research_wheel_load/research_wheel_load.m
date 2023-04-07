function research_wheel_load()
% RESEARCH_WHEEL_LOAD  Entry point

export_graph = true;

input_nonuniform_dataset_file = "output_files/datasets/for_research_movement_direction/one_surface_type/averaged_data.csv";

prefix_path = "output_files/datasets/wheel_load_distribution/gray_surface";
no_average_uniform_dataset_file = fullfile(prefix_path, "no_averaged_data.csv");
input_uniform_dataset_file = fullfile(prefix_path, "averaged_data.csv");

output_dir = "output_files/graphs/wheel_load_distribution";

check_dependency(input_nonuniform_dataset_file, @average_one_surface_type_dataset);

% Process uniform dataset
union_uniform_data = @() union_one_surface_type_dataset("datasets/wheel_load_distribution", ...
                                                        no_average_uniform_dataset_file, ...
                                                        "gray");
process_uniform_data = @() process_one_surface_type_dataset(no_average_uniform_dataset_file);
average_uniform_data = @() average_one_surface_type_dataset(no_average_uniform_dataset_file, ...
                                                            input_uniform_dataset_file);

check_dependency(no_average_uniform_dataset_file, {union_uniform_data, process_uniform_data});
check_dependency(input_uniform_dataset_file, {average_uniform_data});

% Read uniform averaging dataset
opts = detectImportOptions(input_uniform_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
uniform_T = readtable(input_uniform_dataset_file, opts);

% Read nonuniform averaging dataset
opts = detectImportOptions(input_nonuniform_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
nonuniform_T = readtable(input_nonuniform_dataset_file, opts);
% Get the same surface type and speed amplitude as uniform dataset
idx = (nonuniform_T.surftype == uniform_T.surftype(1)) & (nonuniform_T.speedamp == uniform_T.speedamp(1));
nonuniform_T = nonuniform_T(idx, :);

plot_func = localfunctions;
for i = 1:length(plot_func)
    plot_func{i}(nonuniform_T, uniform_T, 'ExportGraphs', export_graph, 'ExportGraphFolder', output_dir);
end
end
%% Local functions
function plot_polar_current_half_second_average(nonuniform_T, uniform_T, options)
arguments
    nonuniform_T table,
    uniform_T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(nonuniform_T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                         "m3cur"]);
check_table_vars(uniform_T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                      "m3cur"]);

nonuniform_surftype = unique(nonuniform_T.surftype);
uniform_surftype = unique(uniform_T.surftype);
assert(length(nonuniform_surftype) == 1, "Nonuniform dataset must has one surface type.");
assert(length(uniform_surftype) == 1, "Uniform dataset must has one surface type.");
assert(nonuniform_surftype == uniform_surftype, ...
       "Surface type of uniform and nonuniform datasets must be the same.");

surf_str = string(uniform_surftype);
for m = ["1", "2", "3"]
    fig_name = surf_str + " motor " + m;
    fig = figure('Name', fig_name);
    mot_str = "m" + m + "cur";
    polarplot(deg2rad(nonuniform_T.movedir), nonuniform_T.(mot_str), ...
              'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', 'red');
    hold on;
    polarplot(deg2rad(uniform_T.movedir), uniform_T.(mot_str), ...
              'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', 'blue');
    pax = gca;
    % Add degree symbol to axis ticks
    pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
    lgd = legend(["nonuniform", "uniform"], 'Location', 'best');
    title(lgd, "Wheel load");
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "polar_current", ...
                                 "half_second_averaging");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, fig_name);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end

function plot_polar_current_full_average(nonuniform_T, uniform_T, options)
arguments
    nonuniform_T table,
    uniform_T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeNonzeroLengthText, mustBeText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(nonuniform_T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                         "m3cur"]);
check_table_vars(uniform_T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                      "m3cur"]);

nonuniform_surftype = unique(nonuniform_T.surftype);
uniform_surftype = unique(uniform_T.surftype);
assert(length(nonuniform_surftype) == 1, "Nonuniform dataset must has one surface type.");
assert(length(uniform_surftype) == 1, "Uniform dataset must has one surface type.");
assert(nonuniform_surftype == uniform_surftype, ...
       "Surface type of uniform and nonuniform datasets must be the same.");

nonuniform_movement_direction = unique(nonuniform_T.movedir);
uniform_movement_direction = unique(uniform_T.movedir);
assert(length(nonuniform_movement_direction) == length(uniform_movement_direction), ...
       "Movement direction vector of uniform and nonuniform datasets must have the same length.");
assert(all(nonuniform_movement_direction == uniform_movement_direction), ...
       "Movement direction vector of uniform and nonuniform datasets must have the same.");

movement_direction = uniform_movement_direction;

surf_str = string(uniform_surftype);
for m = ["1", "2", "3"]
    fig_name = surf_str + " motor " + m;
    fig = figure('Name', fig_name);
    mot_str = strcat("m", m, "cur");
    
    mean_current = group_mean(nonuniform_T.(mot_str), nonuniform_T.movedir, movement_direction);
    move_dir = [movement_direction; movement_direction(1)];
    mean_current = [mean_current; mean_current(1)];
    polarplot(deg2rad(move_dir), mean_current, 'LineStyle', '-', 'Color', "red", ...
              'Marker', '.', 'MarkerSize', 8, 'MarkerEdgeColor', "red");

    hold on;

    mean_current = group_mean(uniform_T.(mot_str), uniform_T.movedir, movement_direction);
    move_dir = [movement_direction; movement_direction(1)];
    mean_current = [mean_current; mean_current(1)];
    polarplot(deg2rad(move_dir), mean_current, 'LineStyle', '-', 'Color', "blue", ...
              'Marker', '.', 'MarkerSize', 8, 'MarkerEdgeColor', "blue");

    % Add degree symbol
    pax = gca;
    pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
    lgd = legend(["nonuniform", "uniform"], 'Location', 'best');
    title(lgd, "Wheel load");
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "polar_current", "full_averaging");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, fig_name);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end

function plot_robot_unsigned_delta_speed_vs_axis_current(nonuniform_T, uniform_T, options)
arguments
    nonuniform_T table,
    uniform_T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(nonuniform_T.Properties.VariableNames, ["surftype", "movedir", "xcur", "ycur", ...
                                                         "rotcur", "ddvx", "ddvy", "ddomega"]);
check_table_vars(uniform_T.Properties.VariableNames, ["surftype", "movedir", "xcur", "ycur", ...
                                                         "rotcur", "ddvx", "ddvy", "ddomega"]);

nonuniform_surftype = unique(nonuniform_T.surftype);
uniform_surftype = unique(uniform_T.surftype);
assert(length(nonuniform_surftype) == 1, "Nonuniform dataset must has one surface type.");
assert(length(uniform_surftype) == 1, "Uniform dataset must has one surface type.");
assert(nonuniform_surftype == uniform_surftype, ...
       "Surface type of uniform and nonuniform datasets must be the same.");

fig_names = ["x_axis", "y_axis", "rotation"];
varcur = ["xcur", "ycur", "rotcur"];
varspeed = ["ddvx", "ddvy", "ddomega"];

surf_str = string(uniform_surftype);
for j = 1:3
    fig = figure("Name", fig_names(j));
    grid on;
    hold on;
    switch j
        case {1, 2}
            nonuniform_ddspeed = nonuniform_T.(varspeed(j)) * 1e3;
            uniform_ddspeed = uniform_T.(varspeed(j)) * 1e3;
        case 3
            nonuniform_ddspeed = nonuniform_T.(varspeed(j));
            uniform_ddspeed = uniform_T.(varspeed(j));
    end
    nonuniform_cur = nonuniform_T.(varcur(j));
    uniform_cur = uniform_T.(varcur(j));
    plot(abs(nonuniform_cur), abs(nonuniform_ddspeed), "LineStyle", "none", "Marker", ".", ...
         "MarkerSize", 10, "MarkerEdgeColor", "red");
    plot(abs(uniform_cur), abs(uniform_ddspeed), "LineStyle", "none", "Marker", ".", ...
         "MarkerSize", 10, "MarkerEdgeColor", "blue");
    switch j
        case 1
            xlabel("Current along X axis, A");
            ylabel("v_{x enc} - v_{x cam}, mm/s", 'Interpreter', 'tex');
        case 2
            xlabel("Current along Y axis, A");
            ylabel("v_{y enc} - v_{y cam}, mm/s", 'Interpreter', 'tex');
        case 3
            xlabel("Rotational current, A");
            ylabel("\omega_{enc} - \omega_{cam}, rad/s", 'Interpreter', 'tex');
    end
    lgd = legend(["nonuniform", "uniform"], 'Location', 'best');
    title(lgd, "Wheel load");
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "robot_delta_speed_vs_axis_current", surf_str);
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, fig_names(j));
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end

function plot_wheel_slippage_vs_motor_current(nonuniform_T, uniform_T, options)
arguments
    nonuniform_T table,
    uniform_T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(nonuniform_T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                         "m3cur", "w1linslip", "w2linslip", ...
                                                         "w3linslip"]);

check_table_vars(uniform_T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                      "m3cur", "w1linslip", "w2linslip", "w3linslip"]);

nonuniform_surftype = unique(nonuniform_T.surftype);
uniform_surftype = unique(uniform_T.surftype);
assert(length(nonuniform_surftype) == 1, "Nonuniform dataset must has one surface type.");
assert(length(uniform_surftype) == 1, "Uniform dataset must has one surface type.");
assert(nonuniform_surftype == uniform_surftype, ...
       "Surface type of uniform and nonuniform datasets must be the same.");

surf_str = string(uniform_surftype);
for m = 1:3
    fig_name = surf_str + " motor " + num2str(m);
    fig = figure("Name", fig_name);
    grid on;
    hold on;
   
    mot_str = strcat("m", num2str(m));
    wh_str = strcat("w", num2str(m));

    nonuniform_slip = nonuniform_T.(strcat(wh_str, "linslip"));
    uniform_slip = uniform_T.(strcat(wh_str, "linslip"));

    nonuniform_cur = nonuniform_T.(strcat(mot_str, "cur"));
    uniform_cur = uniform_T.(strcat(mot_str, "cur"));

    plot(nonuniform_cur, nonuniform_slip, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
         "MarkerEdgeColor", "red");
    plot(uniform_cur, uniform_slip, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
         "MarkerEdgeColor", "blue");

    xlabel("Current, A");
    ylabel("Slippage");
    lgd = legend(["nonuniform", "uniform"], 'Location', 'best');
    title(lgd, "Wheel load");
    ax = gca;
    ax.XLim = [0, 1.6];
    ax.YLim = [0, 1];
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "wheel_slippage_vs_motor_current");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, fig_name);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end