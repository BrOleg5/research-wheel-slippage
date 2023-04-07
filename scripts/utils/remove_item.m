function out_str = remove_item(in_str, num, delimeter)
    arguments
        in_str string {mustBeNonempty},
        num (1, 1) {mustBeNumeric, mustBeInteger, mustBePositive},
        delimeter (1, 1) char = ';'
    end
    split_str = split(in_str, delimeter);
    n = length(split_str);
    assert(num <= n, "Out of range");
    if(num == 1)
        out_str = join(split_str(2:end), delimeter);
    elseif(num == n)
        out_str = join(split_str(1:end-1), delimeter);
    else
        out_str = join([split_str(1:num-1); split_str(num+1:end)], delimeter);
    end
end