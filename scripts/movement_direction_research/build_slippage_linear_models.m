function build_slippage_linear_models(options)
% BUILD_SLIPPAGE_LINEAR_MODELS Entry point

arguments
    options.ExportGraphs (1, 1) logical = true,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = ...
        "output_files/graphs/for_research_movement_direction/slippage_vs_current_linear_models",
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "ru",
    options.ExportModel (1, 1) logical = true,
    options.ExportModelFolder {mustBeText, mustBeNonempty} = ...
        "output_files/models",
end

set(groot, "defaultAxesFontName", "Arial");

input_dataset_file = "output_files/datasets/for_research_movement_direction/one_surface_type/averaged_data.csv";

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

lin_reg_model = impl_build_linear_models(T, 'ExportModel', options.ExportModel, ...
    'ExportModelFolder', options.ExportModelFolder);

plot_compare_data_with_models(T, lin_reg_model, 'ExportGraphs', options.ExportGraphs, ...
    'ExportGraphExtensions', options.ExportGraphExtensions, ...
    'ExportGraphFolder', options.ExportGraphFolder, 'Language', options.Language);

plot_compare_models(T, lin_reg_model, 'ExportGraphs', options.ExportGraphs, ...
    'ExportGraphExtensions', options.ExportGraphExtensions, ...
    'ExportGraphFolder', options.ExportGraphFolder, 'Language', options.Language);
end
%% Local functions
function lin_model = impl_build_linear_models(T, options)
arguments
    T table,
    options.ExportModel (1, 1) logical = false,
    options.ExportModelFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "m1cur", "m2cur", "m3cur", ...
    "w1linslip", "w2linslip", "w3linslip"]);

lin_model = struct;

log_folder = "output_files/logs";
if(~isfolder(log_folder))
    mkdir(log_folder);
end
log_file_name = fullfile(log_folder, "slippage_linear_models.txt");
logID = fopen(log_file_name, "w");

surftype = unique(T.surftype);

fit_func = 'poly1';
fitopts = fitoptions(fit_func, 'Normalize', 'off', 'Robust', 'LAR');

for m = 1:3
    for st = surftype.'
        idx = T.surftype == st;
        sample_T = T(idx, :);
        surf_str = string(st);
        cur_str = "m" + num2str(m) + "cur";
        slip_str = "w" + num2str(m) + "linslip";
        outliers_idx = isoutlier(sample_T.(slip_str));
        sample_T = sample_T(~outliers_idx, :);
        if(st == "table")
            % inverse polynominal
            [fitobj, gof] = fit(sample_T.(slip_str), sample_T.(cur_str), fit_func, fitopts);
            p = [1, -fitobj.p2] / fitobj.p1;
        else
            [fitobj, gof] = fit(sample_T.(cur_str), sample_T.(slip_str), fit_func, fitopts);
            p = [fitobj.p1, fitobj.p2];
        end
       
        % Debug graph
%         plot(reg_mdl);
%         waitforbuttonpress;

        % Write log
        fprintf(logID, "Motor %d, %s surface\n\n", m, surf_str);
        gof_txt = formattedDisplayText(gof, "NumericFormat", "short", "LineSpacing", "compact", ...
            "SuppressMarkup", true, "UseTrueFalseForLogical", true);
        fprintf(logID, gof_txt);
        fprintf(logID, "\nLinear model: s = %.4g * I + %.4g\n", p(1), p(2));
        fprintf(logID, "-----------------------------------------------------------------------\n");
        lin_model(m).(surf_str) = p;
    end
end
fclose(logID);

if(options.ExportModel)
    if(~isfolder(options.ExportModelFolder))
        mkdir(options.ExportModelFolder);
    end
    save(fullfile(options.ExportModelFolder, "slippage_vs_motor_linear_regression_models.mat"), ...
        "lin_model");
end
end

function plot_compare_data_with_models(T, lin_reg_model, options)
arguments
    T table,
    lin_reg_model struct,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = ...
        "output_files/graphs/for_research_movement_direction/slippage_vs_current_linear_models",
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

surftype = unique(T.surftype);
surftype_leg = string(surftype);
colors = get_surface_color_dict();
for m = 1:3
    fig = figure("Name", "Wheel slippage vs current");
    grid on;
    hold on;
    hl = [];
    for st = surftype.'
        idx = T.surftype == st;
        sample_T = T(idx, :);
        surf_str = string(st);
        cur_str = "m" + num2str(m) + "cur";
        slip_str = "w" + num2str(m) + "linslip";
        plot(sample_T.(cur_str), sample_T.(slip_str), "LineStyle", "none", ...
            "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", colors(surf_str));

        cur = [0, max(sample_T.(cur_str))];
        h = plot(cur, polyval(lin_reg_model(m).(surf_str), cur), "Color", "green", "LineWidth", 2);
        hl = [hl; h];
    end
    uistack(hl, "top");
    xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
    ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    ax = gca;
    ax.YLim = [0, 1];
    ax.XLim = [0, 1.6];
    ax.FontSize = 12;
    legend_dict = containers.Map(["en", "ru"], {[surftype_leg; "model"], ...
        [translate_surface(surftype_leg); "модель"]});
    legend_translate(legend_dict, options.Language, 'Location', 'best');
    if(options.ExportGraphs)
        if(~isfolder(options.ExportGraphFolder))
            mkdir(options.ExportGraphFolder);
        end
        file_name = fullfile(options.ExportGraphFolder, "motor" + num2str(m));
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end

function plot_compare_models(T, lin_reg_model, options)
arguments
    T table,
    lin_reg_model struct,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = ...
        "output_files/graphs/for_research_movement_direction/slippage_vs_current_linear_models",
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

surftype = unique(T.surftype);
surftype_leg = string(surftype);
colors = get_surface_color_dict();

line_style = ["-", "--", "-."];

fig = figure("Name", "Wheel slippage vs motor current");
grid on;
hold on;
for m = 1:3
    i = 1;
    for st = surftype.'
        idx = T.surftype == st;
        sample_T = T(idx, :);
        surf_str = string(st);
        cur_str = "m" + num2str(m) + "cur";

        cur = [0, max(sample_T.(cur_str))];
        hpl(m, i) = plot(cur, polyval(lin_reg_model(m).(surf_str), cur), ...
            "Color", colors(surf_str), "LineWidth", 1, "LineStyle", line_style(m));
        i = i + 1;
    end
end
xlabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
ylabel_dict = containers.Map(["en", "ru"], ["Slippage", "Проскальзывание"]);
xlabel_translate(xlabel_dict, options.Language);
ylabel_translate(ylabel_dict, options.Language);
ax = gca;
ax.YLim = [0, 1];
ax.XLim = [0, 1.6];
ax.FontSize = 12;
legend_dict = containers.Map(["en", "ru"], {surftype_leg, translate_surface(surftype_leg)});
surftype_title = containers.Map(["en", "ru"], ["Surface Type", "Тип поверхности"]);
wheel_title = containers.Map(["en", "ru"], ["Wheel number", "Номер колеса"]);
[hl(1).leg, hl(1).obj, hl(1).hout, hl(1).mout] = legendflex(hpl(1, 1:3), ...
    cellstr(legend_dict(options.Language)), 'anchor', {'nw','nw'}, 'buffer', [10 -10], ...
    'fontsize', 11, 'title', surftype_title(options.Language));
[hl(2).leg, hl(2).obj, hl(2).hout, hl(2).mout] = legendflex(hpl(1:3, 3), {'1', '2', '3'}, ...
    'ref', hl(1).leg, 'anchor', {'ne','nw'}, 'buffer', [1 0],'fontsize', 11, ...
    'title', wheel_title(options.Language));
if(options.ExportGraphs)
    if(~isfolder(options.ExportGraphFolder))
        mkdir(options.ExportGraphFolder);
    end
    file_name = fullfile(options.ExportGraphFolder, "compare_models");
    export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                  'BackgroundColor', 'white', 'Colorspace', 'rgb');
    close(fig);
end
end