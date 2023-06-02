function build_linear_regression_models(options)
% BUILD_LINEAR_REGRESSION_MODELS Entry point

arguments
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = ...
        "output_files/graphs/for_research_movement_direction/slippage_vs_current_linear_models",
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en",
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

lin_reg_model = impl_build_linear_regression_models(T, 'ExportModel', options.ExportModel, ...
    'ExportModelFolder', options.ExportModelFolder);

end
%% Local functions
function lin_reg_model = impl_build_linear_regression_models(T, options)
arguments
    T table,
    options.ExportModel (1, 1) logical = false,
    options.ExportModelFolder {mustBeText, mustBeNonempty}
end

check_table_vars(T.Properties.VariableNames, ["surftype", "m1cur", "m2cur", "m3cur", ...
    "w1linslip", "w2linslip", "w3linslip"]);

lin_reg_model = struct;

surftype = unique(T.surftype);
for m = 1:3
    for st = surftype.'
        idx = T.surftype == st;
        sample_T = T(idx, :);
        surf_str = string(st);
        cur_str = "m" + num2str(m) + "cur";
        slip_str = "w" + num2str(m) + "linslip";
        if(st == "table")
            % inverse polynominal
            reg_mdl = fitlm(sample_T, "linear", "PredictorVars", slip_str, "ResponseVar", cur_str);
            p = flip(reg_mdl.Coefficients.Estimate);
            p = [1, -p(2)] / p(1);
        else
            reg_mdl = fitlm(sample_T, "linear", "PredictorVars", cur_str, "ResponseVar", slip_str);
            p = flip(reg_mdl.Coefficients.Estimate);
        end

        % Debug graph
%         plot(reg_mdl);
%         waitforbuttonpress;

        lin_reg_model(m).(surf_str) = p;
    end
end

if(options.ExportModel)
    if(~isfolder(options.ExportModelFolder))
        mkdir(options.ExportModelFolder);
    end
    save(fullfile(options.ExportModelFolder, "slippage_vs_motor_linear_regression_models.mat"), ...
        "lin_reg_model");
end
end