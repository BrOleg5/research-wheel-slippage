function plot_robot_slippage_vs_axis_current_graphs()
% PLOT_ROBOT_SLIPPAGE_VS_AXIS_CURRENT_GRAPHS  Entry point

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
function plot_robot_signed_delta_speed_vs_axis_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "ddvx", "ddvy", "ddomega"]);

fig_names = ["x_axis", "y_axis", "rotation"];
varcur = ["xcur", "ycur", "rotcur"];
varspeed = ["ddvx", "ddvy", "ddomega"];
movedir = unique(T.movedir);
surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for j = 1:3
        fig = figure("Name", strcat(surf_str, " ", fig_names(j)), 'WindowState', 'maximized');
        grid on;
        hold on;
        switch j
            case {1, 2}
                ddspeed = sample_T.(varspeed(j)) * 1e3;
            case 3
                ddspeed = sample_T.(varspeed(j));
        end
        cur = sample_T.(varcur(j));
        n = length(movedir);
        legends = strings(n, 1);
        colors = generate_distrinct_colors(n);
        for k = 1:n
            idx = sample_T.movedir == movedir(k);
            plot(cur(idx), ddspeed(idx), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(k, :));
            legends(k) = string(movedir(k)) + char(176);
        end
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
        lgd = legend(legends, 'Location', 'bestoutside', 'NumColumns', 2, 'FontName', 'Arial', 'FontSize', 12);
        lgd.Title.String = "Directions";
        ax = gca;
        ax.FontName = "Arial";
        ax.FontSize = 14;
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "robot_delta_speed_vs_axis_current", ...
                                  "signed_delta_speed");
            output_dir = fullfile(output_dir, surf_str);
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
end

function plot_robot_delta_speed_vector_vs_current_vector(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "xcur", "ycur", "ddspeedvec"]);

movedir = unique(T.movedir);
surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    fig = figure("Name", surf_str, 'WindowState', 'maximized');
    grid on;
    hold on;
    cur = sqrt(sample_T.xcur.^2 + sample_T.ycur.^2);
    n = length(movedir);
    legends = strings(n, 1);
    colors = generate_distrinct_colors(n);
    for k = 1:n
        idx = sample_T.movedir == movedir(k);
        plot(cur(idx), sample_T.ddspeedvec(idx) * 1e3, "LineStyle", "none", "Marker", ".", ...
             "MarkerSize", 10, "MarkerEdgeColor", colors(k, :));
        legends(k) = string(num2str(movedir(k))) + char(176);
    end
    xlabel("Current, A");
    ylabel("v_{x enc} - v_{x cam}, mm/s", 'Interpreter', 'tex');
    lgd = legend(legends, 'Location', 'bestoutside', 'NumColumns', 2, 'FontName', 'Arial', 'FontSize', 12);
    lgd.Title.String = "Directions";
    ax = gca;
    ax.FontName = "Arial";
    ax.FontSize = 14;
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "robot_delta_speed_vector_vs_current_vector", ...
                              "separate_graphs");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, surf_str);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end

surface_color_dict = get_surface_color_dict();
fig = figure("Name", "Delta robot speed vector vs current vector");
grid on;
hold on;
for i = 1:length(surftype)
    idx = T.surftype == surftype(i);
    sample_T = T(idx, :);
    cur = sqrt(sample_T.xcur.^2 + sample_T.ycur.^2);
    surf_str = string(surftype(i));
    plot(cur, sample_T.ddspeedvec * 1e3, "LineStyle", "none", "Marker", ".", ...
         "MarkerSize", 10, "MarkerEdgeColor", surface_color_dict(surf_str));
end
xlabel("Current, A");
ylabel("v_{enc} - v_{cam}, mm/s", 'Interpreter', 'tex');
legend(string(surftype), 'Location', 'best');
ax = gca;
ax.FontName = "Arial";
ax.FontSize = 14;
if(options.ExportGraphs)
    output_dir = fullfile(options.ExportGraphFolder, "robot_delta_speed_vector_vs_current_vector", ...
                          "common_graphs");
    if(~isfolder(output_dir))
        mkdir(output_dir);
    end
    file_name = fullfile(output_dir, "all_surface");
    export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                  'BackgroundColor', 'white', 'Colorspace', 'rgb');
    close(fig);
end
end

function plot_robot_unsigned_delta_speed_vs_axis_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "xcur", "ycur", "rotcur", ...
                                              "ddvx", "ddvy", "ddomega"]);

fig_names = ["x_axis", "y_axis", "rotation"];
varcur = ["xcur", "ycur", "rotcur"];
varspeed = ["ddvx", "ddvy", "ddomega"];

movedir = unique(T.movedir);
surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for j = 1:3
        fig = figure("Name", strcat(surf_str, " ", fig_names(j)), 'WindowState', 'maximized');
        grid on;
        hold on;
        switch j
            case {1, 2}
                ddspeed = sample_T.(varspeed(j)) * 1e3;
            case 3
                ddspeed = sample_T.(varspeed(j));
        end
        cur = sample_T.(varcur(j));
        n = length(movedir);
        legends = strings(n, 1);
        colors = generate_distrinct_colors(n);
        for k = 1:n
            idx = sample_T.movedir == movedir(k);
            plot(abs(cur(idx)), abs(ddspeed(idx)), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(k, :));
            legends(k) = string(num2str(movedir(k))) + char(176);
        end
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
        lgd = legend(legends, 'Location', 'bestoutside', 'NumColumns', 2, 'FontName', 'Arial', 'FontSize', 12);
        lgd.Title.String = "Directions";
        ax = gca;
        ax.FontName = "Arial";
        ax.FontSize = 14;
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "robot_delta_speed_vs_axis_current", ...
                                  "unsigned_delta_speed", "separate_graphs");
            output_dir = fullfile(output_dir, surf_str);
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

surface_color_dict = get_surface_color_dict();
for j = 1:3
    fig = figure("Name", fig_names(j));
    grid on;
    hold on;
    for i = 1:length(surftype)
        idx = T.surftype == surftype(i);
        sample_T = T(idx, :);
        switch j
            case {1, 2}
                ddspeed = sample_T.(varspeed(j)) * 1e3;
            case 3
                ddspeed = sample_T.(varspeed(j));
        end
        cur = sample_T.(varcur(j));
        surf_str = string(surftype(i));
        plot(abs(cur), abs(ddspeed), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
             "MarkerEdgeColor", surface_color_dict(surf_str));
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
    end
    legend(string(surftype), 'Location', 'best');
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "robot_delta_speed_vs_axis_current", ...
                              "unsigned_delta_speed", "common_graphs");
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