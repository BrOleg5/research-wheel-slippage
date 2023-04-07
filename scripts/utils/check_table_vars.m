function check_table_vars(table_vars, required_vars)
% CHECK_TABLE_VARS 
%   CHECK_TABLE_VARS(table_vars, desired_vars)

arguments
    table_vars {mustBeText},
    required_vars {mustBeText}
end

tf = ismember(required_vars, table_vars);
assert(all(tf), "Table hasn't%s variables", sprintf(" %s", required_vars(~tf)));
end