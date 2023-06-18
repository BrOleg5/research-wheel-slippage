function [x_hat_history, P_history] = kalman_filter_process(x_hat, P, dataset, model, options)
%KALMAN_FILTER_PROCESS

arguments
    x_hat (3, 1) {mustBeNumeric, mustBeReal},
    P (3, 3) {mustBeNumeric, mustBeReal},
    dataset table,
    model struct,
    options.ReferenceDelay = 140e-3
end

check_table_vars(dataset.Properties.VariableNames, ["t", "setvel", "cur", "vel"]);

i_prev = 1;
x_hat_history = x_hat;
P_history = P;
for i = 1:height(dataset)
    if((dataset.t(i) - dataset.t(i_prev)) > options.ReferenceDelay)
        i_prev = i_prev + 1;
    end
    [x_hat, P] = predict(x_hat, dataset.setvel(i_prev), P, model);
    z = [dataset.cur(i); dataset.vel(i)];
    [x_hat, P] = correct(x_hat, z, P, model);
    x_hat_history(:, :, i) = x_hat;
    P_history(:, :, i) = P;
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

% Kalman filter function
function [x_hat, P] = predict(x_hat, u, P, model)
arguments
    x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
    u (:, 1) {mustBeNumeric, mustBeNonempty},
    P (:, :) {mustBeNumeric, mustBeNonempty},
    model struct {mustBeNonempty}
end
    assert(size(P, 1) == size(P, 2), "P matrix must be square.");
    assert(size(model.F, 1) == size(model.F, 2), "F matrix must be square.");
    assert(size(model.Q, 1) == size(model.Q, 2), "Q matrix must be square.");
    assert(length(x_hat) == size(P, 1), "Number of state x_hat must be equal size of P matrix.");
    assert(length(x_hat) == size(model.Q, 1), "Number of state x_hat must be equal size of Q matrix.");
    assert(length(x_hat) == size(model.F, 1), "Number of state x_hat must be equal size of F matrix.");
    assert(length(x_hat) == size(model.G, 1), "Number of state x_hat  must be equal row number of G matrix.");
    assert(length(u) == size(model.G, 2), "Number of control variable u must be equal column number of G matrix.");
    
    x_hat = model.F * x_hat + model.G * u;
    P = model.F * P * model.F.' + model.Q;
end

function [x_hat, P, K] = correct(x_hat, z, P, model)
arguments
    x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
    z (:, 1) {mustBeNumeric, mustBeNonempty},
    P (:, :) {mustBeNumeric, mustBeNonempty},
    model struct {mustBeNonempty}
end
    assert(size(P, 1) == size(P, 2), "P matrix must be square.");
    assert(size(model.R, 1) == size(model.R, 2), "R matrix must be square.");
    assert(length(x_hat) == size(P, 1), "Number of state x_hat must be equal size of P matrix.");
    assert(length(x_hat) == size(model.H, 2), "Number of state x_hat must be equal column number of H matrix.");
    assert(length(z) == size(model.R, 1), "Number of measurements z must be equal size of R matrix.");
    assert(length(z) == size(model.H, 1), "Number of measurements z must be equal row number of H matrix.");

    K = P * model.H.' / (model.H * P * model.H.' + model.R);
    x_hat = x_hat + K * (z - model.H * x_hat);
    I = eye(length(x_hat));
    Z = I - K * model.H;
%     P = Z * P * Z.' + K * model.R * K.';
    P = Z * P;
end