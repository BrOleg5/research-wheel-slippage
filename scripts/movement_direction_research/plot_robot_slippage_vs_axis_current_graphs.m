function plot_robot_slippage_vs_axis_current_graphs()
% PLOT_ROBOT_SLIPPAGE_VS_AXIS_CURRENT_GRAPHS  Entry point

lang = "ru";
export_graph = true;
set(groot, "defaultAxesFontName", "Arial");
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
    plot_func{i}(T, 'ExportGraphs', export_graph, 'ExportGraphFolder', output_dir, 'Language', lang);
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
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
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
        colors = generate_distrinct_colors(n);
        for k = 1:n
            idx = sample_T.movedir == movedir(k);
            plot(cur(idx), ddspeed(idx), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(k, :));
        end
        switch j
            case 1
                xlabel_dict = containers.Map(["en", "ru"], ["Current along X axis, A", "Ток по оси X, А"]);
                ylabel_dict = containers.Map(["en", "ru"], ["v_{x enc} - v_{x cam}, mm/s", ...
                                                            "v_{x энк} - v_{x кам}, мм/с"]);
                xlabel_translate(xlabel_dict, options.Language);
                ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
            case 2
                xlabel_dict = containers.Map(["en", "ru"], ["Current along Y axis, A", "Ток по оси Y, А"]);
                ylabel_dict = containers.Map(["en", "ru"], ["v_{y enc} - v_{y cam}, mm/s", ...
                                                            "v_{y энк} - v_{y кам}, мм/с"]);
                xlabel_translate(xlabel_dict, options.Language);
                ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
            case 3
                xlabel_dict = containers.Map(["en", "ru"], ["Rotational current, A", "Ток вращения, А"]);
                ylabel_dict = containers.Map(["en", "ru"], ["\omega_{enc} - \omega_{cam}, rad/s", ...
                                                            "\omega_{энк} - \omega_{кам}, рад/с"]);
                xlabel_translate(xlabel_dict, options.Language);
                ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
        end
        lgd = legend(string(movedir), 'Location', 'bestoutside', 'NumColumns', 2, 'FontSize', 12);
        legend_title_dict = containers.Map(["en", "ru"], ["Movement direction, \circ", ...
                                                          "Направление движения, \circ"]);
        title(lgd, legend_title_dict(options.Language), 'Interpreter', 'tex');
        ax = gca;
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
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
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
    colors = generate_distrinct_colors(n);
    for k = 1:n
        idx = sample_T.movedir == movedir(k);
        plot(cur(idx), sample_T.ddspeedvec(idx) * 1e3, "LineStyle", "none", "Marker", ".", ...
             "MarkerSize", 10, "MarkerEdgeColor", colors(k, :));
    end
    xlabel_dict = containers.Map(["en", "ru"], ["Current vector, A", "Вектор тока, А"]);
    ylabel_dict = containers.Map(["en", "ru"], ["v_{enc} - v_{cam}, mm/s", ...
                                                "v_{энк} - v_{кам}, мм/с"]);
    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
    lgd = legend(string(movedir), 'Location', 'bestoutside', 'NumColumns', 2, 'FontSize', 12);
    legend_title_dict = containers.Map(["en", "ru"], ["Movement direction, \circ", ...
                                                      "Направление движения, \circ"]);
    title(lgd, legend_title_dict(options.Language), 'Interpreter', 'tex');
    ax = gca;
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
xlabel_dict = containers.Map(["en", "ru"], ["Current vector, A", "Вектор тока, А"]);
ylabel_dict = containers.Map(["en", "ru"], ["v_{enc} - v_{cam}, mm/s", ...
                                            "v_{энк} - v_{кам}, мм/с"]);
xlabel_translate(xlabel_dict, options.Language);
ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');

surf_str = string(surftype);
legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
title(lgd, legend_title_dict(options.Language));
ax = gca;
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
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
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
        colors = generate_distrinct_colors(n);
        for k = 1:n
            idx = sample_T.movedir == movedir(k);
            plot(abs(cur(idx)), abs(ddspeed(idx)), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(k, :));
        end
        switch j
            case 1
                xlabel_dict = containers.Map(["en", "ru"], ["Current along X axis, A", "Ток по оси X, А"]);
                ylabel_dict = containers.Map(["en", "ru"], ["v_{x enc} - v_{x cam}, mm/s", ...
                                                            "v_{x энк} - v_{x кам}, мм/с"]);
                xlabel_translate(xlabel_dict, options.Language);
                ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
            case 2
                xlabel_dict = containers.Map(["en", "ru"], ["Current along Y axis, A", "Ток по оси Y, А"]);
                ylabel_dict = containers.Map(["en", "ru"], ["v_{y enc} - v_{y cam}, mm/s", ...
                                                            "v_{y энк} - v_{y кам}, мм/с"]);
                xlabel_translate(xlabel_dict, options.Language);
                ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
            case 3
                xlabel_dict = containers.Map(["en", "ru"], ["Rotational current, A", "Ток вращения, А"]);
                ylabel_dict = containers.Map(["en", "ru"], ["\omega_{enc} - \omega_{cam}, rad/s", ...
                                                            "\omega_{энк} - \omega_{кам}, рад/с"]);
                xlabel_translate(xlabel_dict, options.Language);
                ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
        end
        lgd = legend(string(movedir), 'Location', 'bestoutside', 'NumColumns', 2, 'FontSize', 12);
        legend_title_dict = containers.Map(["en", "ru"], ["Movement direction, \circ", ...
                                                          "Направление движения, \circ"]);
        title(lgd, legend_title_dict(options.Language), 'Interpreter', 'tex');
        ax = gca;
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
    end
    switch j
        case 1
            xlabel_dict = containers.Map(["en", "ru"], ["Current along X axis, A", "Ток по оси X, А"]);
            ylabel_dict = containers.Map(["en", "ru"], ["v_{x enc} - v_{x cam}, mm/s", ...
                                                        "v_{x энк} - v_{x кам}, мм/с"]);
            xlabel_translate(xlabel_dict, options.Language);
            ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
        case 2
            xlabel_dict = containers.Map(["en", "ru"], ["Current along Y axis, A", "Ток по оси Y, А"]);
            ylabel_dict = containers.Map(["en", "ru"], ["v_{y enc} - v_{y cam}, mm/s", ...
                                                        "v_{y энк} - v_{y кам}, мм/с"]);
            xlabel_translate(xlabel_dict, options.Language);
            ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
        case 3
            xlabel_dict = containers.Map(["en", "ru"], ["Rotational current, A", "Ток вращения, А"]);
            ylabel_dict = containers.Map(["en", "ru"], ["\omega_{enc} - \omega_{cam}, rad/s", ...
                                                        "\omega_{энк} - \omega_{кам}, рад/с"]);
            xlabel_translate(xlabel_dict, options.Language);
            ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
    end
    surf_str = string(surftype);
    legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
    lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
    title(lgd, legend_title_dict(options.Language));
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