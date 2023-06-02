function plot_data_with_speed_amplitude()
% PLOT_SLIPPAGE_VS_CURRENT_WITH_SEVERAL_SPEED_GRAPHS  Entry point

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

plot_func = localfunctions;
for i = 1:length(plot_func)
    plot_func{i}(T, 'ExportGraphs', export_graph, 'ExportGraphFolder', output_dir, 'Language', lang);
end
plot_wheel_slippage_vs_motor_current(T, 'ExportGraphs', export_graph, 'ExportGraphFolder', ...
                                     output_dir, 'ApplyLinearization', true, 'Language', lang);
end
%% Local functions
function plot_robot_unsigned_delta_speed_vs_axis_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "xcur", "ycur", "rotcur", "ddvx", ...
                                              "ddvy", "ddomega"]);

fig_names = ["x_axis", "y_axis", "rotation"];
varcur = ["xcur", "ycur", "rotcur"];
varspeed = ["ddvx", "ddvy", "ddomega"];
colors = ["#FF7C7C", "#B9005B", "#820000", "#7DCE13", "#5BB318", "#2B7A0B"];

speedamp = unique(T.speedamp);
n = length(speedamp);
surftype = unique(T.surftype);
surftype(surftype == "table") = [];
surf_str = string(surftype);
h = length(surftype);
for j = 1:3
    fig = figure("Name", fig_names(j));
    grid on;
    hold on;
    en_legends = strings(3*h, 1);
    ru_legends = strings(3*h, 1);
    for i = 1:h
        for k = 1:n
            leg_idx = (i-1)*n + k;
            en_legends(leg_idx) = sprintf("%s surface %.1f m/s", surf_str(i), speedamp(k));
            ru_legends(leg_idx) = sprintf("%s поверхность %.1f м/с", translate_surface(surf_str(i)), ...
                                          speedamp(k));
            idx = (T.surftype == surftype(i)) & (T.speedamp == speedamp(k));
            sample_T = T(idx, :);
            switch j
                case {1, 2}
                    ddspeed = sample_T.(varspeed(j)) * 1e3;
                case 3
                    ddspeed = sample_T.(varspeed(j));
            end
            plot(abs(sample_T.(varcur(j))), abs(ddspeed), "LineStyle", "none", "Marker", ".", ...
                 "MarkerSize", 10, "MarkerEdgeColor", colors(leg_idx));
        end
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
    legend_dict = containers.Map(["en", "ru"], {en_legends, ru_legends});
    legend_translate(legend_dict, options.Language, 'Location', 'best');
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "robot_delta_speed_vs_axis_current", ...
                              "speed_amplitude_influence");
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

function plot_wheel_delta_velocity_vs_motor_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "m1cur", "m2cur", "m3cur", "w1ddvel", ...
                                              "w2ddvel", "w3ddvel"]);

colors = ["#FF7C7C", "#B9005B", "#820000", "#7DCE13", "#5BB318", "#2B7A0B"];

speedamp = unique(T.speedamp);
n = length(speedamp);
surftype = unique(T.surftype);
surftype(surftype == "table") = [];
h = length(surftype);
for m = 1:3
    fig_name = "motor " + num2str(m);
    fig = figure("Name", fig_name);
    grid on;
    hold on;
    en_legends = strings(3*h, 1);
    ru_legends = strings(3*h, 1);
    for i = 1:h
        surf_str = string(surftype(i));
        for k = 1:n
            leg_idx = (i-1)*n + k;
            en_legends(leg_idx) = sprintf("%s surface %.1f m/s", surf_str, speedamp(k));
            ru_legends(leg_idx) = sprintf("%s поверхность %.1f м/с", translate_surface(surf_str), ...
                                          speedamp(k));
            idx = (T.surftype == surftype(i)) & (T.speedamp == speedamp(k));
            sample_T = T(idx, :);
            cur = sample_T.("m" + num2str(m) + "cur");
            ddvel = sample_T.("w" + num2str(m) + "ddvel");
            plot(cur, ddvel, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(leg_idx));
        end
    end
    xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
    ylabel_dict = containers.Map(["en", "ru"], ["\omega_{enc} - \omega_{cam}, rad/s", ...
                                                "\omega_{энк} - \omega_{кам}, рад/с"]);
    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');

    legend_dict = containers.Map(["en", "ru"], {en_legends, ru_legends});
    legend_translate(legend_dict, options.Language, 'Location', 'best');
    ax = gca;
    ax.XLim = [0, 1.6];
    ax.YLim = [0, 6];
    ax.FontSize = 12;
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "wheel_delta_velocity_vs_motor_current", ...
                              "speed_amplitude_influence");
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

function plot_wheel_slippage_vs_motor_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.ApplyLinearization (1, 1) {mustBeNumericOrLogical} = false,
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

if(options.ApplyLinearization)
    check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                  "m3cur", "w1linslip", "w2linslip", "w3linslip"]);
else
    check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                                  "m3cur", "w1slip", "w2slip", "w3slip"]);
end

colors = ["#FF7C7C", "#B9005B", "#820000", "#7DCE13", "#5BB318", "#2B7A0B"];

speedamp = unique(T.speedamp);
n = length(speedamp);
surftype = unique(T.surftype);
surftype(surftype == "table") = [];
h = length(surftype);
for m = 1:3
    fig_name = "motor " + num2str(m);
    fig = figure("Name", fig_name);
    grid on;
    hold on;
    en_legends = strings(3*h, 1);
    ru_legends = strings(3*h, 1);
    for i = 1:h
        surf_str = string(surftype(i));
        for k = 1:n
            leg_idx = (i-1)*n + k;
            en_legends(leg_idx) = sprintf("%s surface %.1f m/s", surf_str, speedamp(k));
            ru_legends(leg_idx) = sprintf("%s поверхность %.1f м/с", translate_surface(surf_str), ...
                                          speedamp(k));
            idx = (T.surftype == surftype(i)) & (T.speedamp == speedamp(k));
            sample_T = T(idx, :);
            cur = sample_T.("m" + num2str(m) + "cur");
            if(options.ApplyLinearization)
                slip = sample_T.("w" + num2str(m) + "linslip");
            else
                slip = sample_T.("w" + num2str(m) + "slip");
            end
            plot(cur, slip, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(leg_idx));
        end
    end
    xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
    ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    legend_dict = containers.Map(["en", "ru"], {en_legends, ru_legends});
    legend_translate(legend_dict, options.Language, 'Location', 'best');
    ax = gca;
    ax.XLim(1) = 0;
    ax.YLim(1) = 0;
    ax.FontSize = 12;
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "wheel_slippage_vs_motor_current", ...
                              "speed_amplitude_influence");
        if(options.ApplyLinearization)
            output_dir = fullfile(output_dir, "linearization");
        else
            output_dir = fullfile(output_dir, "not_linearization");
        end
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