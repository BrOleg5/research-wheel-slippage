function process_wheel_load_distribution()
% PROCESS_WHEEL_LOAD_DISTRIBUTION  Entry point

input_dataset_folder = "datasets/wheel_load_distribution";
output_dataset_folder = "output_files/datasets/wheel_load_distribution";
file_name = ["nonuniform_load.csv", "uniform_load.csv"];
input_dataset_file = fullfile(input_dataset_folder, file_name);
output_dataset_file = fullfile(output_dataset_folder, file_name);

assert(length(input_dataset_file) == length(output_dataset_file), ...
       "input_dataset_file and output_dataset_file must have the same length.");

for i = 1:length(input_dataset_file)
    opts = detectImportOptions(input_dataset_file(i), "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
    T = readtable(input_dataset_file(i), opts);

    processed_T = process_wheel_load_data(T);
    
    output_file_dir = fileparts(output_dataset_file(i));
    if(~isfolder(output_file_dir))
        mkdir(output_file_dir);
    end
    writetable(processed_T, output_dataset_file(i), "WriteRowNames", true);
end
end
%% Local functions
function T_out = process_wheel_load_data(T_in)
arguments
    T_in table
end

check_table_vars(T_in.Properties.VariableNames, ["num", "w1load", "w2load", "w3load"]);

T_out = table('RowNames', ["wheel1", "wheel2", "wheel3"]);

T1 = varfun(@mean, T_in, 'InputVariables', ["w1load", "w2load", "w3load"]);
T_out.mean_load = T1.Variables';

T2 = varfun(@std, T_in, 'InputVariables', ["w1load", "w2load", "w3load"]);
T_out.std_load = T2.Variables';

std_mean = @(x)(std(x)/sqrt(length(x)));
T3 = varfun(std_mean, T_in, 'InputVariables', ["w1load", "w2load", "w3load"]);
T_out.std_mean_load = T3.Variables';

t = tinv(0.95, height(T_in));
T4 = varfun(@(x)(t*x), T3);
T_out.conf_bond_load = T4.Variables';

g = 9.81;
T_out.mean_norm_force = g .* T_out.mean_load ./  1e3;
T_out.std_norm_force = g .* T_out.std_load ./ 1e3;
T_out.std_mean_norm_force = g .* T_out.std_mean_load ./ 1e3;
T_out.conf_bond_norm_force = g .* T_out.conf_bond_load ./ 1e3;
end