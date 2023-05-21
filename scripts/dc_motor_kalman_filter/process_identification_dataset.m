function out_T = process_identification_dataset(in_T)
%PROCESS_IDENTIFICATION_DATASET

arguments
    in_T table
end

check_table_vars(in_T.Properties.VariableNames, ["t", "m1cur", "m2cur", "m3cur", ...
                                                 "m1vel", "m2vel", "m3vel"]);

out_T = in_T;
% Conver time from ms to s
out_T.t = out_T.t / 1000;

% Shift time by 0.03 s so that time start from 0
out_T.t = out_T.t + 0.03;

% Add start row to table. Fill this row with zeros
start_row = zeros(1, width(out_T));
out_T = [array2table(start_row, "VariableNames", out_T.Properties.VariableNames); out_T];

% Add sign to current values. Sign taken from motor/wheel velocity
for m = 1:3
    mstr = "m" + num2str(m);
    out_T.(mstr + "cur") = sign(out_T.(mstr + "vel")) .* out_T.(mstr + "cur");
end
end