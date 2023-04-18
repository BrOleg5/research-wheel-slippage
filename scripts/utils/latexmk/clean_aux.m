function clean_aux(target)
% CLEAN_AUX  Clean auxiliary files after compiling target report

arguments
    target {mustBeTextScalar}
end

target_path = fullfile("reports", target);
assert(isfolder(target_path), "Path '%s' doesn't exist.", target_path);

old_current_path = cd(target_path);
out_dir = fullfile("../../output_files/reports", target);
command = sprintf("latexmk -c -output-directory=%s", out_dir);
system(command);
cd(old_current_path);
end