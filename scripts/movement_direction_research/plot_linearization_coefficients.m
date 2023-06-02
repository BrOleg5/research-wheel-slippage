function plot_linearization_coefficients(options)
%PLOT_LINEARIZATION_COEFFICIENTS Entry point

arguments
    options.ExportGraphs (1, 1) logical = true,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = "output_files/graphs/for_research_movement_direction",
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "ru"
end

set(groot, "defaultPolarAxesFontName", "Arial");

[~, robotKinematics] = get_robot_parameters();
dirs = 0:5:360;
K = calc_linearization_coefficient(deg2rad(dirs), robotKinematics);

colors = ["black", "blue", "red"];
fig = figure("Name", "Wheel coeffs vs movement direction");
grid on;
for m = 1:3
    polarplot(deg2rad(dirs), squeeze(K(m, :, :)), 'Color', colors(m), 'Marker', '.', 'MarkerSize', 5, ...
              'MarkerEdgeColor', colors(m));
    hold on;
end
lgn = legend(string(1:3), "Location", "best");
legend_title_dict = containers.Map(["en", "ru"], ["Wheel", "Колесо"]);
title(lgn, legend_title_dict(options.Language));
if(options.ExportGraphs)
    if(~isfolder(options.ExportGraphFolder))
        mkdir(options.ExportGraphFolder);
    end
    file_name = fullfile(options.ExportGraphFolder, ...
        "linearization_coefficients_vs_movement_directions");
    export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                  'BackgroundColor', 'white', 'Colorspace', 'rgb');
    close(fig);
end
end