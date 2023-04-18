function build_report(target, rc_file)
% BUILD_REPORT  Build target report

arguments
    target {mustBeTextScalar},
    rc_file {mustBeTextScalar} = "../common/.latexmkrc"
end

target_path = fullfile("reports", target);
assert(isfolder(target_path), "Path '%s' doesn't exist. Target '%s' may not exist or may be located in different path.", ...
       target_path, target);
assert(isfile(fullfile(target_path, target + ".tex")), "File '%s.tex' doesn't exist or locates in different path.", ...
       target);

old_current_path = cd(target_path);
out_dir = fullfile("../../output_files/reports", target);
command = sprintf("latexmk -r %s -output-directory=%s %s.tex", rc_file, out_dir, target);
system(command);
cd(old_current_path);
end