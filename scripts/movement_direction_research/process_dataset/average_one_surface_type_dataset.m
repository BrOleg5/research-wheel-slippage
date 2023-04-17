function average_one_surface_type_dataset(input_dataset_file, output_dataset_file, Ts_ms)
% AVERAGE_ONE_SURFACE_TYPE_DATASET Entry point

arguments
    input_dataset_file {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/datasets/for_research_movement_direction/one_surface_type/no_averaged_data.csv",
    output_dataset_file {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/datasets/for_research_movement_direction/one_surface_type/averaged_data.csv",
    Ts_ms (1, 1) {mustBeNumeric, mustBePositive} = 500
end

if check_dependency(input_dataset_file, @process_one_surface_type_dataset)
    % Check data is processed
    opts = detectImportOptions(input_dataset_file, "VariableDescriptionsLine", 1, ...
                               "VariableUnitsLine", 2, "VariableNamesLine", 3);
%     if(length(opts.VariableNames) ~= 48)
    if(length(opts.VariableNames) ~= 51)
        fprintf("'%s' isn't processed.\nRun %s() function.\n", input_dataset_file, "process_one_surface_type_dataset");
        process_one_surface_type_dataset();
    end
end

% Read headlines of file
fid = fopen(input_dataset_file);
column_description = fgetl(fid);
column_units = fgetl(fid);
fclose(fid);

opts = detectImportOptions(input_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
T = readtable(input_dataset_file, opts);

% Convert time from double to duration with rounding
% because milliseconds must be integer in Matlab
T.t = milliseconds(floor(T.t));
% summary(T);

expnum = unique(T.expnum, 'stable');
output_T = timetable;
ts = milliseconds(Ts_ms);
for n = expnum.'
    idx = T.expnum == n;
    processed_T = table2timetable(T(idx, :), "RowTimes", "t");

    desired_vars_order = processed_T.Properties.VariableNames;

    % Save surftype from current experiment and remove it, because it isn't
    % averaged
    surftype = processed_T.surftype(1);
    processed_T = removevars(processed_T, "surftype");

    %%%%%%%%%%%%%%%%%%%%%%%
    % Averaging in 500 ms %
    %%%%%%%%%%%%%%%%%%%%%%%s

    % Durty hack for remove start averaging bin
    % Shift time on 10 ms
    dummyRow = processed_T(1, :);
    timeShift = milliseconds(10);
    processed_T.t = processed_T.t + timeShift;
    processed_T = [dummyRow; processed_T];

    newTimes = [processed_T.t(1), processed_T.t(2):ts:processed_T.t(end)];
    % Special averaging for motor positions (take last value in interval)
    pos_vals = ["m1pos", "m2pos", "m3pos"];
    pos_T = processed_T(:, pos_vals);
    processed_pos_T = retime(pos_T, newTimes, "lastvalue", "IncludedEdge", "right");

    time = processed_T.t;

    processed_T = removevars(processed_T, pos_vals);
    processed_T = retime(processed_T, newTimes, "mean", "IncludedEdge", "right");

    processed_T = join(processed_T, processed_pos_T);

    % Remove time shift
    processed_T(1:2, :) = [];
    processed_T.t = processed_T.t - timeShift;
    
    % Add surftype in timetable
    h = size(processed_T, 1);
    surftype = repelem(surftype, h, 1);
    processed_T = addvars(processed_T, surftype, 'Before', 'speedamp');

    processed_T = resortvars(processed_T, desired_vars_order);

    output_T = [output_T; processed_T];
end
% Change time unit from ms to s
% Remove first description and units - time variable
column_description = remove_item(column_description, 1);
column_units = remove_item(column_units, 1);

% Write headlines in file
fid = fopen(output_dataset_file, "w");
fprintf(fid, "Time;%s\ns;%s", column_description, column_units);
fclose(fid);

% Convert from timetable to table, because writetimetable write duration
% with unit in file
output_T = timetable2table(output_T);
output_T.t = seconds(output_T.t);

% Write table in file
writetable(output_T, output_dataset_file, 'Delimiter', ';', 'WriteMode', 'Append', ...
           'WriteVariableNames', true);
end
%% Local functions
function T2 = resortvars(T1, desired_vars_order)
arguments
    T1 timetable,
    desired_vars_order {mustBeText, mustBeVector}
end
current_vars_order = T1.Properties.VariableNames;
assert(length(current_vars_order) == length(desired_vars_order), ...
       "Length current_vars_order and desired_vars_order must have the same length.");

[~, varOrder] = ismember(current_vars_order, desired_vars_order); 
[~, resortOrder] = sort(varOrder); 
T2 = T1(:, resortOrder);
end