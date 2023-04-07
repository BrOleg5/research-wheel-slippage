function out_str = replace_item(in_str, num, new_str, delimeter)
    arguments
        in_str {mustBeTextScalar},
        num (1, 1) {mustBeNumeric, mustBeInteger, mustBePositive},
        new_str {mustBeTextScalar}
        delimeter (1, 1) char = ';'
    end
    split_str = split(in_str, delimeter);
    n = length(split_str);
    assert(num <= n, "Out of range");
    split_str(num) = new_str;
    out_str = join(split_str, delimeter);
end