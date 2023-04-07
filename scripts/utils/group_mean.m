function group_mean = group_mean(arr, prop, group)
    arguments
        arr {mustBeNumeric, mustBeVector},
        prop {mustBeVector},
        group {mustBeVector}
    end
    assert(length(arr) == length(prop), "Arr and prop must be the same length");
    group_mean = zeros(size(group));
    for i = 1:length(group)
        idx = prop == group(i);
        group_mean(i) = mean(arr(idx));
    end
end