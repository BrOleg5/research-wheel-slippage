function calc_correlation_between_motor_current_and_wheel_velocity(options)
% CALC_CORRELATION_BETWEEN_MOTOR_CURRENT_AND_WHEEL_VELOCITY  Entry point

arguments
    options.ExportGraphs (1, 1) logical = true,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = ...
        "output_files/graphs/for_research_movement_direction/correlation_between_motor_current_and_wheel_velocity"
end

dataset_file = "output_files/datasets/for_research_movement_direction/one_surface_type/averaged_data.csv";

check_dependency(dataset_file, @average_one_surface_type_dataset);

opts = detectImportOptions(dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
T = readtable(dataset_file, opts);

idx = doublecmp(T.speedamp , 0.1);
T = T(idx, :);

surftype = unique(T.surftype);
for st = surftype.'
    idx = T.surftype == st;
    sample_T = T(idx, :);
    surf_str = string(st);
    for i = 1:3
        curstr = "m" + num2str(i) + "cur";
        velstr = "w" + num2str(i) + "vel";
        sample_T.(velstr) = abs(sample_T.(velstr));
        fig = figure("Name", "Correlation motor " + num2str(i) + " " + surf_str);
        corrplot(sample_T, "DataVariables", [curstr, velstr], "Type", "Pearson", "TestR", "on");
        if(options.ExportGraphs)
            output_dir = fullfile(options.ExportGraphFolder, surf_str);
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