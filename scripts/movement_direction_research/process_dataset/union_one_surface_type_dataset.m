function union_one_surface_type_dataset(input_dataset_folder, output_dataset_file, surf_type)
% UNION_SURFACE_TYPE_DATASET Entry point

arguments
    input_dataset_folder {mustBeTextScalar, mustBeNonzeroLengthText} = "datasetsfor_research_movement_direction/one_surface_type",
    output_dataset_file {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/datasets/for_research_movement_direction/one_surface_type/no_averaged_data.csv",
    surf_type {mustBeText, mustBeNonempty, mustBeNonzeroLengthText} = ["table", "gray", "green"]
end

suffix = "_surface";
input_dirs = fullfile(input_dataset_folder, strcat(surf_type, suffix));

movement_direction = 0:5:355; % degree
speed_amplitude = 0.1:0.1:0.3; % m/s
one_speed_nums = length(movement_direction);

output_table = table;
exp_num = 0;

for i = 1:length(input_dirs)
    assert(isfolder(input_dirs(i)), "Directory %s not exist!", input_dirs(i));
    list_of_files = dir(strcat(input_dirs(i), "/exp*.csv"));
    if(isempty(list_of_files))
        fprintf("Directory %s is empty.", input_dirs(i));
        continue;
    end
    idx = sort_filenames_by_number({list_of_files.name});
    list_of_files = list_of_files(idx);
    n = length(list_of_files);
    file_nums = 0:n-1;
    idx = mod(file_nums, one_speed_nums) + 1;
    md = movement_direction(idx);
    idx = idivide(int32(file_nums), int32(one_speed_nums)) + 1;
    sp_amp = speed_amplitude(idx);
    for j = 1:n
        file_path = fullfile(list_of_files(j).folder, list_of_files(j).name);
        T = readtable(file_path, "VariableDescriptionsLine", 1, "VariableUnitsLine", 2, ...
                      "VariableNamesLine", 3, "LeadingDelimitersRule", "error", ...
                      "TrailingDelimitersRule", "ignore", "ConsecutiveDelimitersRule", "error", ...
                      "MissingRule", "error", "ImportErrorRule", "error", "EmptyLineRule", "skip");
        
        assert(validate_dataset(T), "Validation of dataset '%s' is fault.", file_path);

        % Remove start zeros
        idx = (T.m1vel == 0) & (T.m2vel == 0) & (T.m3vel == 0);
        T(idx, :) = [];
        
        % Reset to zero start experiment time
        T.t = T.t - T.t(1);

        assert(all(T.t >= 0), "Time must be positive");

        h = size(T, 1);
        movedir = repelem(md(j), h, 1);
        speedamp = repelem(sp_amp(j), h, 1);
        surftype = repelem(surf_type(i), h, 1);
        expnum = repelem(exp_num, h, 1);
        T = addvars(T, expnum, 'After', 't');
        T = addvars(T, surftype, 'After', 'expnum');
        T = addvars(T, speedamp, 'After', 'surftype');
        T = addvars(T, movedir, 'After', 'speedamp');

        output_table = [output_table; T];
        exp_num = exp_num + 1;
    end
end

[output_dir, ~, ~] = fileparts(output_dataset_file);
if(~isfolder(output_dir))
    mkdir(output_dir);
end

% Read headlines of file
fid = fopen(file_path);
column_description = fgetl(fid);
column_units = fgetl(fid);
fclose(fid);

% Remove last delimeter
column_description = column_description(1:end-1);
column_units = column_units(1:end-1);
% Remove first description and units - time variable
column_description = remove_item(column_description, 1);
column_units = remove_item(column_units, 1);

% Write headerlines in file
fid = fopen(output_dataset_file, "w");
fprintf(fid, "Time;Experiment number;Surface type;Movement direction;Speed amplitude;%s\nms;none;none;degree;m/s;%s", ...
        column_description, column_units);
fclose(fid);

writetable(output_table, output_dataset_file, 'Delimiter', ';', 'WriteMode', 'Append', ...
           'WriteVariableNames', true);
end
%% Local function
function idx = sort_filenames_by_number(filenames)
% SORT_FILENAMES_BY_NUMBER  Sort array of filenames.
% 
% Filename form is *n.* where n is number. Eg. exp05.csv, video_52.csv etc. 
% Correctly sort filenames without leading zeros. Eg. array of exp15.csv, exp1.csv and exp150.csv 
% will be sorted like this: exp1.csv, exp15.csv and exp150.csv.
% 
% Implementation inspered by <a href="matlab:web('https://uk.mathworks.com/matlabcentral/answers/360531-how-do-i-sort-filenames-containing-text-and-numbers-in-numerical-order-in-matlab#answer_285085')">
% this answer<\a>.

arguments
    filenames string {mustBeVector, mustBeNonzeroLengthText, mustBeNonempty}
end
str_number = regexp(filenames, '\d+', 'match', 'once');
[~, idx] = sort(str2double(str_number));
end