function plot_wheel_slippage_vs_current_graphs()
% PLOT_WHEEL_SLIPPAGE_VS_AXIS_CURRENT_GRAPHS  Entry point

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
plot_wheel_slippage_vs_motor_current(T, 'ExportGraphs', export_graph, ...
                                     'ExportGraphFolder', output_dir, 'Language', lang, ...
                                     'ApplyLinearization', true);
end
%% Local functions
function plot_shit_wheel_slippage_vs_motor_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.ApplyLinearization (1, 1) {mustBeNumericOrLogical} = false,
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

if(options.ApplyLinearization)
    check_table_vars(T.Properties.VariableNames, ["surftype", "m1cur", "m2cur", ...
                                                  "m3cur", "w1linslip", "w2linslip", "w3linslip"]);
else
    check_table_vars(T.Properties.VariableNames, ["surftype", "m1cur", "m2cur", ...
                                                  "m3cur", "w1slip", "w2slip", "w3slip"]);
end

fig = figure("Name", "First try");
grid on;
hold on;
for m = 1:3
    curstr = "m" + num2str(m) + "cur";
    slipstr = "w" + num2str(m) + "slip";
    plot(T.(curstr), T.(slipstr), "LineStyle", "none", "Marker", ".", "MarkerSize", 5, ...
         "MarkerEdgeColor", "black");
end
xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
xlabel_translate(xlabel_dict, options.Language);
ylabel_translate(ylabel_dict, options.Language);
ax = gca;
ax.YLim = [0, 1];
if(options.ExportGraphs)
    output_dir = fullfile(options.ExportGraphFolder, "shit_wheel_slippage_vs_motor_current");
    if(~isfolder(output_dir))
        mkdir(output_dir);
    end
    file_name = fullfile(output_dir, "first_try");
    export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                  'BackgroundColor', 'white', 'Colorspace', 'gray');
    close(fig);
end
end

function plot_wheel_slippage_vs_motor_current(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
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

movedir = unique(T.movedir);
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
        n = length(movedir);
        colors = generate_distrinct_colors(n);
        if(options.ApplyLinearization)
            slip = sample_T.(strcat(wh_str, "linslip"));
        else
            slip = sample_T.(strcat(wh_str, "slip"));
        end
        for j = 1:n
            idx = sample_T.movedir == movedir(j);
            plot(cur(idx), slip(idx), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(j, :));
        end
        xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
        ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
        xlabel_translate(xlabel_dict, options.Language);
        ylabel_translate(ylabel_dict, options.Language);
        lgd = legend(string(movedir), 'Location', 'bestoutside', 'NumColumns', 2, 'FontSize', 12);
        legend_title_dict = containers.Map(["en", "ru"], ["Movement direction, \circ", ...
                                                          "Направление движения, \circ"]);
        title(lgd, legend_title_dict(options.Language), 'Interpreter', 'tex');
        ax = gca;
        ax.XLim(1) = 0;
        ax.YLim(1) = 0;
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
    xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
    ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    surf_str = string(surftype);
    legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
    lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
    title(lgd, legend_title_dict(options.Language));
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
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "w1ddvel", "w2ddvel", "w3ddvel"]);

movedir = unique(T.movedir);
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
        n = length(movedir);
        colors = generate_distrinct_colors(n);
        for j = 1:length(movedir)
            idx = sample_T.movedir == movedir(j);
            plot(cur(idx), abs(domega(idx)), "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
                 "MarkerEdgeColor", colors(j, :));
        end
        xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
        ylabel_dict = containers.Map(["en", "ru"], ["\omega_{enc} - \omega_{cam}, rad/s", ...
                                                    "\omega_{энк} - \omega_{кам}, рад/с"]);
        xlabel_translate(xlabel_dict, options.Language);
        ylabel_translate(ylabel_dict, options.Language, 'Interpreter', 'tex');
        lgd = legend(string(movedir), 'Location', 'bestoutside', 'NumColumns', 2, 'FontSize', 12);
        legend_title_dict = containers.Map(["en", "ru"], ["Movement direction, \circ", ...
                                                          "Направление движения, \circ"]);
        title(lgd, legend_title_dict(options.Language), 'Interpreter', 'tex');
        ax = gca;
        ax.XLim(1) = 0;
        ax.YLim(1) = 0;
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

function plot_lin_wheel_slippage_vs_prod_cur_and_vel(T, options)
arguments
    T table,
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_table_vars(T.Properties.VariableNames, ["surftype", "movedir", "m1cur", "m2cur", ...
                                              "m3cur", "m1vel", "m2vel", "m3vel", ...
                                              "w1linslip", "w2linslip", "w3linslip"]);

movedir = unique(T.movedir);
surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for m = 1:3
        mot_str = strcat("m", num2str(m));
        wh_str = strcat("w", num2str(m));
        productCurrentAndVelocity = sample_T.(strcat(mot_str, "cur")) .* sample_T.(strcat(mot_str, "vel"));
        fig = figure("Name", strcat(surf_str, " motor ", num2str(m)), 'WindowState', 'maximized');
        grid on;
        hold on;
        n = length(movedir);
        colors = generate_distrinct_colors(n);
        slip = sample_T.(strcat(wh_str, "linslip"));
        for j = 1:n
            idx = sample_T.movedir == movedir(j);
            plot(productCurrentAndVelocity(idx), slip(idx), "LineStyle", "none", "Marker", ".", ...
                 "MarkerSize", 10, "MarkerEdgeColor", colors(j, :));
        end
        xlabel_dict = containers.Map(["en", "ru"], ["Current \times velocity, A \times rad/s", ...
                                                    "Ток \times скорость, А \times рад/с"]);
        ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
        xlabel_translate(xlabel_dict, options.Language, 'Interpreter', 'tex');
        ylabel_translate(ylabel_dict, options.Language);
        lgd = legend(string(movedir), 'Location', 'bestoutside', 'NumColumns', 2, 'FontSize', 12);
        legend_title_dict = containers.Map(["en", "ru"], ["Movement direction, \circ", ...
                                                          "Направление движения, \circ"]);
        title(lgd, legend_title_dict(options.Language), 'Interpreter', 'tex');
        ax = gca;
        ax.XLim(1) = 0;
        ax.YLim(1) = 0;
        ax.FontSize = 14;
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, "lin_wheel_slippage_vs_prod_cur_and_vel", ...
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
        productCurrentAndVelocity = sample_T.(strcat(mot_str, "cur")) .* sample_T.(strcat(mot_str, "vel"));
        slip = sample_T.(strcat(wh_str, "linslip"));
        surf_str = string(surftype(i));
        plot(productCurrentAndVelocity, slip, "LineStyle", "none", "Marker", ".", "MarkerSize", 5, ...
             "MarkerEdgeColor", surface_color_dict(surf_str));
    end
    xlabel_dict = containers.Map(["en", "ru"], ["Current \times velocity, A \times rad/s", ...
                                                "Ток \times скорость, А \times рад/с"]);
    ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
    xlabel_translate(xlabel_dict, options.Language, 'Interpreter', 'tex');
    ylabel_translate(ylabel_dict, options.Language);
    surf_str = string(surftype);
    legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
    lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
    title(lgd, legend_title_dict(options.Language));
    ax = gca;
    ax.XLim(1) = 0;
    ax.YLim = [0, 1];
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "lin_wheel_slippage_vs_prod_cur_and_vel", ...
                              "common_graphs");
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