function plot_polar_graphs()
% PLOT_POLAR_GRAPHS  Entry point

lang = "ru";
export_graph = true;
set(groot, "defaultAxesFontName", "Arial");
set(groot, "defaultPolarAxesFontName", "Arial");
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
function plot_polar_current_half_second_average(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", "m3cur"]);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for m = ["1", "2", "3"]
        fig = figure('Name', strcat(surf_str, " motor ", m));
        mot_str = strcat("m", m, "cur");
        polarplot(deg2rad(sample_T.movedir), sample_T.(mot_str), 'LineStyle', 'none', ...
                  'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', 'black');
        if(isMATLABReleaseOlderThan("R2022a"))
            % Add degree symbol to axis ticks
            pax = gca;
            pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
        end
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "polar_current", ...
                                     "half_second_averaging", surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, strcat("motor", m));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'gray');
            close(fig);
        end
    end
end

surface_color_dict = get_surface_color_dict();
surftype = unique(T.surftype);
for m = ["1", "2", "3"]
    fig = figure('Name', strcat("motor ", m));
    for st = surftype.'
        idx = T.surftype == st;
        sample_T = T(idx, :);
        surf_str = string(st);
        mot_str = strcat("m", m, "cur");
        polarplot(deg2rad(sample_T.movedir), sample_T.(mot_str), 'LineStyle', 'none', ...
                  'Marker', '.', 'MarkerSize', 10, 'MarkerEdgeColor', surface_color_dict(surf_str));
        hold on;
        grid on;
    end
    if(isMATLABReleaseOlderThan("R2022a"))
        % Add degree symbol to axis ticks
        pax = gca;
        pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
    end
    surf_str = string(surftype);
    legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
    lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
    title(lgd, legend_title_dict(options.Language));
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "polar_current", ...
                                 "half_second_averaging");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, strcat("motor", m));
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end

function plot_polar_velocity_half_second_average(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "w1vel", "w2vel", "w3vel"]);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for w = ["1", "2", "3"]
        fig = figure('Name', strcat(surf_str, " wheel ", w));
        mot_str = strcat("w", w, "vel");
        polarplot(deg2rad(sample_T.movedir), abs(sample_T.(mot_str)), 'LineStyle', 'none', ...
                  'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', 'black');
        if(isMATLABReleaseOlderThan("R2022a"))
            % Add degree symbol to axis ticks
            pax = gca;
            pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
        end
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "polar_velocity", ...
                                     "half_second_averaging", surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, strcat("wheel", w));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'gray');
            close(fig);
        end
    end
end
end

function plot_polar_slippage_half_second_average(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeNonzeroLengthText, mustBeText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "w1slip", "w2slip", "w3slip"]);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for m = ["1", "2", "3"]
        fig = figure('Name', strcat(surf_str, " motor ", m));
        mot_str = strcat("w", m, "slip");
        polarplot(deg2rad(sample_T.movedir), sample_T.(mot_str), 'LineStyle', 'none', ...
                  'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', 'black');
        if(isMATLABReleaseOlderThan("R2022a"))
            % Add degree symbol to axis ticks
            pax = gca;
            pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
        end
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "polar_slippage", ...
                                     "half_second_averaging", surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, strcat("motor", m));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'gray');
            close(fig);
        end
    end
end
end

function plot_polar_robot_delta_speed_half_second_average(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeNonzeroLengthText, mustBeText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "ddvx", "ddvy", "ddomega"]);

surftype = unique(T.surftype);
varnames = ["ddvx", "ddvy", "ddomega"];
fig_names = ["x_axis", "y_axis", "rotation"];
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for j = 1:3
        fig = figure('Name', strcat(surf_str, " ", fig_names(j)));
        polarplot(deg2rad(sample_T.movedir), abs(sample_T.(varnames(j))), 'LineStyle', 'none', 'Marker', '.', ...
                  'MarkerSize', 5, 'MarkerEdgeColor', 'black');
        if(isMATLABReleaseOlderThan("R2022a"))
            % Add degree symbol to axis ticks
            pax = gca;
            pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
        end
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "polar_delta_speed", ...
                                     "half_second_averaging", surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, fig_names(j));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'gray');
            close(fig);
        end
    end
end
end

function plot_polar_current_full_average(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeNonzeroLengthText, mustBeText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", "m3cur"]);

surftype = unique(T.surftype);
movement_direction = unique(T.movedir);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for m = ["1", "2", "3"]
        fig = figure('Name', strcat(surf_str, " motor ", m));
        mot_str = strcat("m", m, "cur");
        mean_current = group_mean(sample_T.(mot_str), sample_T.movedir, movement_direction);
        move_dir = [movement_direction; movement_direction(1)];
        mean_current = [mean_current; mean_current(1)];
        polarplot(deg2rad(move_dir), mean_current, 'LineStyle', '-', 'Color', 'black', ...
                  'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', 'black');
        if(isMATLABReleaseOlderThan("R2022a"))
            % Add degree symbol to axis ticks
            pax = gca;
            pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
        end
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "polar_current", "full_averaging", surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, strcat("motor", m));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'gray');
            close(fig);
        end
    end
end

motor_color = ["black", "blue", "red"];
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    fig = figure('Name', surf_str);
    for m = 1:3
        mot_str = strcat("m", num2str(m), "cur");
        mean_current = group_mean(sample_T.(mot_str), sample_T.movedir, movement_direction);
        move_dir = [movement_direction; movement_direction(1)];
        mean_current = [mean_current; mean_current(1)];
        polarplot(deg2rad(move_dir), mean_current, 'LineStyle', '-', 'Color', motor_color(m), ...
                  'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', 'black');
        hold on;
        if(isMATLABReleaseOlderThan("R2022a"))
            % Add degree symbol to axis ticks
            pax = gca;
            pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
        end
    end
    lgd = legend(string(1:3), 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Motor", "Двигатель"]);
    title(lgd, legend_title_dict(options.Language));
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "polar_current", "full_averaging");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, surf_str);
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        saveas(fig, file_name + ".fig");
        close(fig);
    end
end

surface_color_dict = get_surface_color_dict();
for m = ["1", "2", "3"]
    fig = figure('Name', strcat("motor ", m));
    mot_str = strcat("m", m, "cur");
    for i = 1:length(surftype)
        idx = T.surftype == surftype(i);
        sample_T = T(idx, :);
        mean_current = group_mean(sample_T.(mot_str), sample_T.movedir, movement_direction);
        move_dir = [movement_direction; movement_direction(1)];
        mean_current = [mean_current; mean_current(1)];
        surf_str = string(surftype(i));
        polarplot(deg2rad(move_dir), mean_current, 'LineStyle', '-', 'Color', surface_color_dict(surf_str), ...
                  'Marker', '.', 'MarkerSize', 10, 'MarkerEdgeColor', surface_color_dict(surf_str));
        hold on;
    end
    if(isMATLABReleaseOlderThan("R2022a"))
        % Add degree symbol to axis ticks
        pax = gca;
        pax.ThetaTickLabel = string(pax.ThetaTickLabel) + char(176);
    end
    surf_str = string(surftype);
    legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
    lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
    title(lgd, legend_title_dict(options.Language));
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "polar_current", "full_averaging");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, strcat("motor", m));
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end