function out_T = process_test_dataset(in_T)
% ADD_CURRENT_SIGN

arguments
    in_T table
end

check_table_vars(in_T.Properties.VariableNames, ["t", "m1cur", "m2cur", "m3cur", ...
                                                 "m1vel", "m2vel", "m3vel"]);

out_T = in_T;
% Conver time from ms to s
dt = diff(out_T.t);
if(mean(dt) > 1)
    out_T.t = out_T.t / 1e3;
end

% Add sign to current values. Sign taken from motor/wheel velocity
for m = 1:3
    mstr = "m" + num2str(m);
    out_T.(mstr + "cur") = sign(out_T.(mstr + "vel")) .* out_T.(mstr + "cur");
end
end