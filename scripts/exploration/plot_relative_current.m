function plot_relative_current(input_dataset_file, output_graph_dir, options)
arguments
    input_dataset_file {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/datasets/for_research_movement_direction/one_surface_type/averaged_data.csv",
    output_graph_dir {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/graphs/explorations/relative_current"
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
end

check_dependency(input_dataset_file, @average_one_surface_type_dataset);

opts = detectImportOptions(input_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
T = readtable(input_dataset_file, opts);

if(~isfolder(output_graph_dir))
    mkdir(output_graph_dir);
end

surftype = unique(T.surftype);
T.sumcur = sum(T(:, ["m1cur", "m2cur", "m3cur"]).Variables, 2);
surface_color_dict = get_surface_color_dict();
for m = 1:3
    fig = figure("Name", strcat("Motor ", num2str(m)));
    grid on;
    hold on;
    for i = 1:length(surftype)
        idx = T.surftype == surftype(i);
        sample_T = T(idx, :);

        mot_str = strcat("m", num2str(m));
        rel_cur = sample_T.(mot_str + "cur") ./ sample_T.sumcur;
        movedir = sample_T.movedir;
        surf_str = string(surftype(i));
        plot(movedir, rel_cur, "LineStyle", "none", "Marker", ".", "MarkerSize", 10, ...
             "MarkerEdgeColor", surface_color_dict(surf_str));
    end
    xlabel_dict = containers.Map(["en", "ru"], ...
        ["Movement direction, \circ", "Направление движения, \circ"]);
    ylabel_dict = containers.Map(["en", "ru"], ["Relative current", "Отностельный ток"]);
    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    surf_str = string(surftype);
    legend_dict = containers.Map(["en", "ru"], {surf_str, translate_surface(surf_str)});
    lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
    legend_title_dict = containers.Map(["en", "ru"], ["Surface type", "Тип поверхности"]);
    title(lgd, legend_title_dict(options.Language));
    ax = gca;
    ax.XTick = 0:10:350;
    ax.XLimitMethod = 'tight';
    ax.YLimitMethod = 'padded';
    file_name = fullfile(output_graph_dir, strcat("motor", num2str(m)));
    savefig(fig, file_name);
    close(fig);
end
end
