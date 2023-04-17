function plot_wheel_velocity_vs_motor_current()
% PLOT_WHEEL_VELOCITY_VS_MOTOR_CURRENT Entry point

lang = "ru";
export_graph = true;
set(groot, "defaultAxesFontName", "Arial");
prefix_path = "output_files/datasets/for_research_movement_direction/one_surface_type";
dataset_file = fullfile(prefix_path, ["no_averaged_data.csv", "averaged_data.csv"]);

prefix_path = "output_files/graphs/for_research_movement_direction/wheel_velocity_vs_motor_current";
output_dir = fullfile(prefix_path, ["no_averaged_data", "averaged_data"]);

check_dependency(dataset_file, @average_one_surface_type_dataset);

plot_func = localfunctions;
for i = 1:length(dataset_file)
    opts = detectImportOptions(dataset_file(i), "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
    opts = setvartype(opts, "surftype", "categorical");
    T = readtable(dataset_file(i), opts);
    
    idx = doublecmp(T.speedamp , 0.1);
    T = T(idx, :);

    for j = 1:length(plot_func)
        plot_func{j}(T, 'ExportGraphs', export_graph, 'ExportGraphFolder', output_dir(i), ...
                     'Language', lang);
    end
end
end
%% Local functions
function plot_wheel_velocity_vs_motor_current_separate_with_set_velocity(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "m1setvel", "m2setvel", "m3setvel", ...
                                              "m1cur", "m2cur", "m3cur", "w1vel", "w2vel", "w3vel"]);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for i = 1:3
        curstr = "m" + num2str(i) + "cur";
        velstr = "w" + num2str(i) + "vel";
        setvelstr = "m" + num2str(i) + "setvel";

        sample_T.(velstr) = abs(sample_T.(velstr));
        sample_T.(setvelstr) = abs(sample_T.(setvelstr));

        % Optimization size vector image: get only unique observation
        % Impact only not averaging dataset
        sample_T1 = unique(sample_T(:, [curstr, velstr, setvelstr]));

        fig = figure("Name", "motor " + num2str(i) + " " + surf_str, "WindowState", "maximized");
        grid on;
        hold on;
        setvel = unique(sample_T1.(setvelstr));
        n = length(setvel);
        legends = strings(n, 1);
        colors = generate_distrinct_colors(n);
        for k = 1:n
            idx = sample_T1.(setvelstr) == setvel(k);
            plot(sample_T1(idx, curstr).Variables, sample_T1(idx, velstr).Variables, ...
                 "LineStyle", "none", "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", colors(k, :));
            legends(k) = sprintf("%.5f", setvel(k));
        end
        xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
        ylabel_dict = containers.Map(["en", "ru"], ["Velocity, rad/s", "Угловая скорость, рад/с"]);
        xlabel_translate(xlabel_dict, options.Language);
        ylabel_translate(ylabel_dict, options.Language);
        lgd = legend(legends, 'Location', 'best', 'FontSize', 12);
        legend_title_dict = containers.Map(["en", "ru"], ["Set velocity, rad/s", ...
                                                          "Уставка угловой скорости, рад/с"]);
        title(lgd, legend_title_dict(options.Language));
        ax = gca;
        ax.FontSize = 14;
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "separate_with_set_velocity", surf_str);
            if(~isfolder(output_dir))
                mkdir(output_dir);
            end
            file_name = fullfile(output_dir, "motor" + num2str(i));
            export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                          'BackgroundColor', 'white', 'Colorspace', 'rgb');
            close(fig);
        end
    end
end
end

function plot_wheel_velocity_vs_motor_current_separate_with_set_vel_2(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "m1setvel", "m2setvel", "m3setvel", ...
                                              "m1cur", "m2cur", "m3cur", "w1vel", "w2vel", "w3vel"]);

surftype = unique(T.surftype);
surface_color_dict = get_surface_color_dict();
T1 = T;
for i = 1:3
    curstr = "m" + num2str(i) + "cur";
    velstr = "w" + num2str(i) + "vel";
    setvelstr = "m" + num2str(i) + "setvel";

    T1.(velstr) = abs(T1.(velstr));
    T1.(setvelstr) = abs(T1.(setvelstr));

    % Optimization size vector image: get only unique observation
    % Impact only not averaging dataset
    T2 = unique(T1(:, ["surftype", curstr, velstr, setvelstr]));

    fig = figure("Name", "motor " + num2str(i));
    grid on;
    hold on;
    for st = surftype.'
        idx = T2.surftype == st;
        sample_T = T2(idx, :);
        surf_str = string(st);
        plot(sample_T.(curstr), sample_T.(velstr), "LineStyle", "none", "Marker", ".", ...
             "MarkerSize", 10, "MarkerEdgeColor", surface_color_dict(surf_str));
    end
    xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
    ylabel_dict = containers.Map(["en", "ru"], ["Velocity, rad/s", "Угловая скорость, рад/с"]);
    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    surf_str = string(surftype);
    legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
    lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
    title(lgd, legend_title_dict(options.Language));
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "separate_with_surface_type");
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, "motor" + num2str(i));
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end