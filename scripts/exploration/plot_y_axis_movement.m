function plot_y_axis_movement(input_dataset_file, output_graph_dir)
arguments
    input_dataset_file {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/datasets/for_research_movement_direction/one_surface_type/no_averaged_data.csv",
    output_graph_dir {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/graphs/explorations/y_axis_movement"
end

check_dependency(input_dataset_file, @process_one_surface_type_dataset);

opts = detectImportOptions(input_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
T = readtable(input_dataset_file, opts);

if(~isfolder(output_graph_dir))
    mkdir(output_graph_dir);
end

movedir = 90;
surftype = "green";
motnum = 2;
mstr = "m" + num2str(motnum);
for speedamp = 0.1:0.1:0.3
    idx = (T.movedir == movedir) & (T.surftype == surftype) & doublecmp(T.speedamp, speedamp);
    T1 = T(idx, ["t", mstr + "cur", mstr + "vel", "xpos", "ypos"]);
    T1 = sortrows(T1, "t");
    dx = [0; diff(T1.xpos)] * 1e3;
    dy = [0; diff(T1.ypos)] * 1e3;
    T1.t = T1.t / 1e3;
    figname = sprintf("%s_motor_%d_%d_deg_speed_%.1f_m_s.fig", surftype, motnum, movedir, speedamp);
    fig = figure("Name", figname, "WindowState", "maximized");
    plot_data(fig, T1.t, T1.(mstr + "cur"), T1.(mstr + "vel"), dx, dy);
    figpath = fullfile(output_graph_dir, figname);
    savefig(fig, figpath);
    close(fig);
end

surftype = "gray";
speedamp = 0.1;
idx = (T.movedir == movedir) & (T.surftype == surftype) & doublecmp(T.speedamp, speedamp);
T1 = T(idx, ["t", mstr + "cur", mstr + "vel", "xpos", "ypos"]);
T1 = sortrows(T1, "t");
dx = [0; diff(T1.xpos)] * 1e3;
dy = [0; diff(T1.ypos)] * 1e3;
T1.t = T1.t / 1e3;
figname = sprintf("%s_motor_%d_%d_deg_speed_%.1f_m_s.fig", surftype, motnum, movedir, speedamp);
fig = figure("Name", figname, "WindowState", "maximized");
plot_data(fig, T1.t, T1.(mstr + "cur"), T1.(mstr + "vel"), dx, dy);
figpath = fullfile(output_graph_dir, figname);
savefig(fig, figpath);
close(fig);

surftype = "table";
speedamp = 0.1;
idx = (T.movedir == movedir) & (T.surftype == surftype) & doublecmp(T.speedamp, speedamp);
T1 = T(idx, ["t", mstr + "cur", mstr + "vel", "xpos", "ypos"]);
T1 = sortrows(T1, "t");
dx = [0; diff(T1.xpos)] * 1e3;
dy = [0; diff(T1.ypos)] * 1e3;
T1.t = T1.t / 1e3;
figname = sprintf("%s_motor_%d_%d_deg_speed_%.1f_m_s.fig", surftype, motnum, movedir, speedamp);
fig = figure("Name", figname, "WindowState", "maximized");
plot_data(fig, T1.t, T1.(mstr + "cur"), T1.(mstr + "vel"), dx, dy);
figpath = fullfile(output_graph_dir, figname);
savefig(fig, figpath);
close(fig);
end
%% Local functions
function plot_data(obj, time_s, current_A, velocity_rad_s, dx_mm, dy_mm)
arguments
    obj {iptcheckhandle(obj, {'axes', 'figure'}, 'export_graphs', 'obj', 1)},
    time_s {mustBeReal, mustBeVector},
    current_A {mustBeReal, mustBeVector},
    velocity_rad_s {mustBeReal, mustBeVector},
    dx_mm {mustBeReal, mustBeVector},
    dy_mm {mustBeReal, mustBeVector}
end
assert(length(time_s) == length(current_A), "time_s and current_A must be the same length.");
assert(length(time_s) == length(velocity_rad_s), "time_s and velocity_rad_s must be the same length.");
assert(length(time_s) == length(dx_mm), "time_s and dx_mm must be the same length.");
assert(length(time_s) == length(dy_mm), "time_s and dy_mm must be the same length.");

data = [current_A, velocity_rad_s, dx_mm, dy_mm];

y_label_txt = ["Current, A", "Velocity, rad/s", "dx, mm", "dy, mm"];
figure(obj);
tiledlayout(2, 2, "TileSpacing", "tight", "Padding", "compact");
for i = 1:4
    nexttile;
    grid on;
    hold on;
    plot(time_s, data(:, i), "Color", "black", "LineWidth", 0.8);
    xlabel("Time, s");
    ylabel(y_label_txt(i));
end
end