function plot_wheel_slippage_vs_axis_current_graphs()
% PLOT_WHEEL_SLIPPAGE_VS_AXIS_CURRENT_GRAPHS  Entry point

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
plot_wheel_slippage_vs_motor_current(T, 'ExportGraphs', export_graph, ...
                                     'ExportGraphFolder', output_dir, 'ApplyLinearization', true);
end
%% Local functions
function plot_wheel_slippage_vs_motor_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.ApplyLinearization (1, 1) {mustBeNumericOrLogical} = false
end

if(options.ApplyLinearization)
    check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                  "m3cur", "w1linslip", "w2linslip", "w3linslip"]);
else
    check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                  "m3cur", "w1slip", "w2slip", "w3slip"]);
end

movedirs = unique(T.movedir);
surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for m = 1:3
        mot_str = strcat("m", num2str(m));
        wh_str = strcat("w", num2str(m));
        cur = sample_T.(strcat(mot_str, "cur"));
        fig = figure("Name", strcat(surf_str, " motor ", num2str(m)), 'WindowState', 'maximized');
        grid on;
        hold on;
        n = length(movedirs);
        legends = strings(n, 1);
        colors = generate_distrinct_colors(n);
        if(options.ApplyLinearization)
            slip = sample_T.(strcat(wh_str, "linslip"));
        else
            slip = sample_T.(strcat(wh_str, "slip"));
        end
        for j = 1:n
            idx = sample_T.movedir == movedirs(j);
            plot(cur(idx), slip(idx), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(j, :));
            legends(j) = string(movedirs(j)) + char(176);
        end
        xlabel("Current, A");
        ylabel("Slippage");
        lgd = legend(legends, 'Location', 'bestoutside', 'NumColumns', 2, 'FontName', 'Arial', 'FontSize', 12);
        lgd.Title.String = "Directions";
        ax = gca;
        ax.XLim(1) = 0;
        ax.YLim(1) = 0;
        ax.FontName = "Arial";
        ax.FontSize = 14;
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "wheel_slippage_vs_motor_current", ...
                                  "separate_graphs");
            if(options.ApplyLinearization)
                output_dir = fullfile(output_dir, "linearization");
            else
                output_dir = fullfile(output_dir, "not_linearization");
            end
            output_dir = fullfile(output_dir, surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, strcat("motor", num2str(m)));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'rgb');
            close(fig);
        end
    end
end

surface_color_dict = get_surface_color_dict();
for m = 1:3
    fig = figure("Name", strcat("Motor ", num2str(m)));
    grid on;
    hold on;
    for i = 1:length(surftype)
        idx = T.surftype == surftype(i);
        sample_T = T(idx, :);

        mot_str = strcat("m", num2str(m));
        wh_str = strcat("w", num2str(m));
        cur = sample_T.(strcat(mot_str, "cur"));
        if(options.ApplyLinearization)
            slip = sample_T.(strcat(wh_str, "linslip"));
        else
            slip = sample_T.(strcat(wh_str, "slip"));
        end
        surf_str = string(surftype(i));
        plot(cur, slip, "LineStyle", "none", "Marker", ".", "MarkerSize", 5, ...
             "MarkerEdgeColor", surface_color_dict(surf_str));
    end
    xlabel("Current, A");
    ylabel("Slippage");
    legend(string(surftype), 'Location', 'best');
    ax = gca;
    ax.XLim = [0, 1.6];
    ax.YLim = [0, 1];
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "wheel_slippage_vs_motor_current", ...
                              "common_graphs");
        if(options.ApplyLinearization)
            output_dir = fullfile(output_dir, "linearization");
        else
            output_dir = fullfile(output_dir, "not_linearization");
        end
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, strcat("motor", num2str(m)));
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end
%% Delta velocity vs current graphs
function plot_wheel_delta_velocity_vs_motor_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "w1ddvel", "w2ddvel", "w3ddvel"]);

movedirs = unique(T.movedir);
surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for m = 1:3
        mot_str = strcat("m", num2str(m));
        wh_str = strcat("w", num2str(m));
        domega = sample_T.(strcat(wh_str, "ddvel"));
        cur = sample_T.(strcat(mot_str, "cur"));
        fig = figure("Name", strcat(surf_str, " motor ", num2str(m)), 'WindowState', 'maximized');
        grid on;
        hold on;
        n = length(movedirs);
        legends = strings(n, 1);
        colors = generate_distrinct_colors(n);
        for j = 1:length(movedirs)
            idx = sample_T.movedir == movedirs(j);
            plot(cur(idx), abs(domega(idx)), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(j, :));
            legends(j) = string(num2str(movedirs(j))) + char(176);
        end
        xlabel("Current, A");
        ylabel("\omega_{enc} - \omega_{cam}, rad/s", 'Interpreter', 'tex');
        lgd = legend(legends, 'Location', 'bestoutside', 'NumColumns', 2, 'FontName', 'Arial', 'FontSize', 12);
        lgd.Title.String = "Directions";
        ax = gca;
        ax.XLim(1) = 0;
        ax.YLim(1) = 0;
        ax.FontName = "Arial";
        ax.FontSize = 14;
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "wheel_delta_velocity_vs_motor_current", ...
                                  "separate_graphs");
            output_dir = fullfile(output_dir, surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, strcat("motor", num2str(m)));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'rgb');
            close(fig);
        end
    end
end
end