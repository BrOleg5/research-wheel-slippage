function test_kalman_filter(options)

arguments
    options.ExportGraphs (1, 1) logical = true,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = "output_files/graphs/dc_motor_kalman_filter",
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "ru",
    options.ModelPath {mustBeText, mustBeNonempty} = "output_files/model_parameters/motor_models.mat",
end

load(options.ModelPath, "motor_model");

% Datasets
dataset_path = "datasets\dc_motor_model_identification";
dataset_file = ["data_for_identification_without_integral.csv", "validation_data.csv", ...
    "green_gray_validation_data.csv", "table_green_validation_data.csv", "gray_0.1_175_validation.csv"];
dataset_name = ["train_data", "validation_data", "green_gray_validation_data", ...
    "table_green_validation_data", "gray_0.1_175_validation"];
dataset_path = fullfile(dataset_path, dataset_file);

for i = 1:length(dataset_path)
    % Load dataset
    dataset = readtable(dataset_path(i), "NumHeaderLines", 2);
    dataset = process_identification_dataset(dataset);
    impl_test_kalman_filter(dataset, motor_model, 'ExportGraphs', options.ExportGraphs, ...
        'ExportGraphExtensions', options.ExportGraphExtensions, ...
        'ExportGraphFolder', options.ExportGraphFolder, 'Language', options.Language, ...
        'ExportFileName', dataset_name(i));
end
end
%% Local functions
function impl_test_kalman_filter(dataset, motor_model, options)

arguments
    dataset table,
    motor_model struct,
    options.ExportGraphs (1, 1) logical = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"],
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = ...
    "output_files/graphs/dc_motor_kalman_filter",
    options.ExportFileName {mustBeTextScalar},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "ru",
end

P = [3^2,      0,       0;
       0,  300^2,       0;
       0,      0,   0.2^2];
x_hat = zeros(size(P, 1), 1);

for motor_num = 1:3
    mstr = "m" + num2str(motor_num);
    motor_data = dataset(:, ["t", mstr + "setvel", mstr + "cur", mstr + "vel"]);
    motor_data = renamevars(motor_data,[mstr + "setvel", mstr + "cur", mstr + "vel"], ...
        ["setvel", "cur", "vel"]);

    [x_hat_history, P_history] = kalman_filter_process(x_hat, P, motor_data, motor_model(motor_num));

    ci = calc_confidence_interval(x_hat_history, P_history);
    ci = [squeeze(ci(:, 1, :)), squeeze(ci(:, 2, end:-1:1))];

    [state_fig, uncertancy_fig] = plot_compare(x_hat_history, P_history, motor_data, ci);

    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, "motor" + num2str(motor_num));
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, options.ExportFileName);
        export_graphs(state_fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(state_fig);

        if(motor_num == 1)
            file_name = fullfile(output_dir, "uncertancy");
            export_graphs(uncertancy_fig, file_name, options.ExportGraphExtensions, ...
                'ContentType', 'vector', 'BackgroundColor', 'white', 'Colorspace', 'rgb');
        end
        close(uncertancy_fig);
    end
end
end

function [state_fig, uncertancy_fig] = plot_compare(x_hat_history, P_history, dataset, ci)
arguments
    x_hat_history (3, 1, :) {mustBeNumeric, mustBeReal},
    P_history (3, 3, :) {mustBeNumeric, mustBeReal},
    dataset table,
    ci {mustBeNumeric, mustBeReal}
end

n_x = size(x_hat_history, 1);

check_table_vars(dataset.Properties.VariableNames, ["t", "setvel", "cur", "vel"]);

state_fig = figure("Name", "Kalman filter", "WindowState", "maximized");
y_name = ["Ток, А", "Скорость, рад/с", "Момент, Нм"];
t_name = ["cur", "vel"];
tiledlayout(n_x, 1, "TileSpacing", "tight", "Padding", "compact");
for i = 1:n_x
    nexttile;
    grid on;
    hold on;
    fill([dataset.t; dataset.t(end:-1:1)], ci(i, :), "green", "FaceAlpha", 0.4, "EdgeColor", "none");
    if(i ~= 3)
        plot(dataset.t, dataset.(t_name(i)), 'Color', 'blue', 'LineWidth', 1);
    end
    plot(dataset.t, squeeze(x_hat_history(i, :, :)), 'Color', 'red', 'LineWidth', 1);
    xlabel("Время, с");
    ylabel(y_name(i));
    switch i
        case {1, 2}
            legend(["Доверительный интервал 95%", "Измерения", "Отфильтрованное"], ...
                "Location", "best");
        case 3
            legend(["Доверительный интервал 95%", "Оценка"], ...
                "Location", "best");
    end
    ax = gca;
    ax.FontSize = 13;
    ax.XLimitMethod = 'tight';
    ax.YLimitMethod = 'tickaligned';
end

uncertancy_fig = figure("Name", "Kalman filter uncertancy", "WindowState", "maximized");
uns_name = ["Неопределённость тока", "Неопределённость скорости", "Неопределённость момента"];
tiledlayout(n_x, 1, "TileSpacing", "tight", "Padding", "compact");
for i = 1:n_x
    nexttile;
    grid on;
    hold on;
    plot(dataset.t, squeeze(P_history(i, i, :)), 'Color', 'black', 'LineWidth', 2);
    xlabel("Время, с");
    ylabel(uns_name(i));
    ax = gca;
    ax.FontSize = 13;
    ax.XLimitMethod = 'tight';
    ax.YLimitMethod = 'padded';
end
end

% Confidence interval
function ci = calc_confidence_interval(x_hat, P)
arguments
    x_hat (:, :, :) {mustBeNumeric, mustBeNonempty}
    P (:, :, :) {mustBeNumeric, mustBeNonempty}
end
    for i = 1:size(x_hat, 3)
        for j = 1:size(x_hat, 1)
            if(P(j, j, i) ~= 0)
                ci(j, :, i) = norminv([0.025, 0.975], x_hat(j, 1, i), sqrt(P(j, j, i)));
            else
                ci(j, :, i) = x_hat(j, 1, i) * [1, 1];
            end
        end
    end
end