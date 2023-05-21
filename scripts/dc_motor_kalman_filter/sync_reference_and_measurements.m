function [out_T, tau] = sync_reference_and_measurements(in_T, options)
%SYNC_REFERENCE_AND_MEASUREMENTS

arguments
    in_T table,
    options.DeadZone = 3
end

check_table_vars(in_T.Properties.VariableNames, ["t", "m1setvel", "m2setvel", "m3setvel", ...
                                                 "m1vel", "m2vel", "m3vel"]);

out_T = in_T;
% Find indexes of reference step changing
ref_step_idx = find(diff(abs(out_T.m1setvel) < options.DeadZone));
% Find indexes of velocity step changing
vel_step_idx = find(diff(abs(out_T.m1vel) < options.DeadZone));

assert(length(ref_step_idx) == length(vel_step_idx), ...
       "ref_step_idx and vel_step_idx have the same length.");
assert(rem(length(ref_step_idx), 2) == 0, "ref_step_idx and vel_step_idx must have even length.");

for i = 1:length(ref_step_idx)/2
    rise_gap_idx = ref_step_idx(2*i-1):vel_step_idx(2*i-1);
    fall_gap_idx = ref_step_idx(2*i):vel_step_idx(2*i);
    for m = 1:3
        setvelstr = "m" + num2str(m) + "setvel";
        out_T.(setvelstr)(rise_gap_idx) = 0;
        out_T.(setvelstr)(fall_gap_idx) = out_T.(setvelstr)(fall_gap_idx(1)-1);
    end
end

tau = mean(out_T.t(vel_step_idx) - out_T.t(ref_step_idx));
end