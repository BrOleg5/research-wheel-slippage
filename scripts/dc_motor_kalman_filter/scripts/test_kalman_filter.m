%% Configuration
n_x = 3;
n_u = 1;
n_z = 2;
model = struct('F', F, 'G', G, 'H', H, 'Q', Q_d, 'R', R_d);
% model = struct('F', F, 'G', G, 'H', H, 'Q', Qest1, 'R', Rest1);
x_hat = zeros(n_x, 1);
P = [3^2,      0,       0;
       0,  300^2,       0;
       0,      0,   0.2^2];
% P = zeros(n_x);
z = zeros(n_z);
%% Load experimental data
% T = readtable("data\data_for_identification_without_integral.csv", "NumHeaderLines", 2);
% T = readtable("data\validation_data.csv", "NumHeaderLines", 2);
T = readtable("data\green_gray_validation_data.csv", "NumHeaderLines", 2);
% T = readtable("data\table_green_validation_data.csv", "NumHeaderLines", 2);
T.t = T.t / 1000;
T.t = T.t + 0.03;
start_row = zeros(1, width(T));
T = [array2table(start_row, "VariableNames", T.Properties.VariableNames); T];
T.m1cur = sign(T.m1vel) .* T.m1cur;
T.m2cur = sign(T.m2vel) .* T.m2cur;
T.m3cur = sign(T.m3vel) .* T.m3cur;
%% Simulation
i_prev = 1;
x_hat_history = x_hat;
P_history = P;
for i = 1:height(T)
    if((T.t(i) - T.t(i_prev)) > motor1.tau)
        i_prev = i_prev + 1;
    end
    [x_hat, P] = predict(x_hat, T.m1setvel(i_prev), P, model);
    z = [T.m1cur(i); T.m1vel(i)];
%     if((abs(z(1)) < 1e-2) && (abs(z(2)) < 3))
%         model.R = R_d / 10;
%     else
%         model.R = R_d;
%     end
    [x_hat, P] = correct(x_hat, z, P, model);
    x_hat_history(:, :, i) = x_hat;
    P_history(:, :, i) = P;
end
ci = calc_confidence_interval(x_hat_history, P_history);
ci = [squeeze(ci(:, 1, :)), squeeze(ci(:, 2, end:-1:1))];
%% Plot result
fig = figure("Name", "Kalman filter", "WindowState", "maximized");
y_name = ["Current, A", "Velocity, rad/s", "Load torque, N*m"];
t_name = ["m1cur", "m1vel"];
uns_name = ["Current uncertancy", "Velocity uncertancy", "Torue uncertancy"];
tiledlayout(n_x, 2, "TileSpacing", "tight", "Padding", "compact");
for i = 1:n_x
%     subplot(n_x, 2, 2*i-1);
    nexttile;
    fill([T.t; T.t(end:-1:1)], ci(i, :), "green", "FaceAlpha", 0.2, "EdgeColor", "none");
    grid on;
    hold on;
    if(i ~= 3)
        plot(T.t, T.(t_name(i)), 'Color', 'blue', 'LineWidth', 0.75);
%     else
%         plot(T.t, -k_t.*squeeze(x_hat_history(1, :, :)) + b.*squeeze(x_hat_history(2, :, :)));
    end
    plot(T.t, squeeze(x_hat_history(i, :, :)), 'Color', 'red', 'LineWidth', 0.75);
    xlabel("Time, s");
    ylabel(y_name(i));
    if(i == 2)
        legend(["Confidence interval 95%", "Measure", "Filtered"], "Location", "best", "FontSize", 11);
%     else
%         legend(["Confidence interval 95%", "Filetred"], "Location", "best", "FontSize", 11);
%         legend(["Confidence interval 95%", "Naive", "Filetred"], "Location", "best");
    end
    ax = gca;
    ax.FontSize = 11;
%     subplot(n_x, 2, 2*i);
    nexttile;
    plot(T.t, squeeze(P_history(i, i, :)), 'Color', 'black', 'LineWidth', 0.75);
    grid on;
    hold on;
    xlabel("Time, s");
    ylabel(uns_name(i));
    ax = gca;
    ax.FontSize = 11;
end
% subplot(n_x, 2, 1);
% ax = gca;
% ax.YLim = [-1.1, 1.1];
% subplot(n_x, 2, 5);
% ax = gca;
% ax.YLim = [-0.022, 0.022];
%% Kalman filter function
% function [x_hat, P] = predict(x_hat, u, P, model, Ts)
% arguments
%     x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
%     u (:, 1) {mustBeNumeric, mustBeNonempty},
%     P (:, :) {mustBeNumeric, mustBeNonempty},
%     model struct {mustBeNonempty},
%     Ts (1, 1) {mustBeNumeric}
% end
%     assert(size(P, 1) == size(P, 2), "P matrix must be square.");
%     assert(size(model.F(1), 1) == size(model.F(1), 2), "F matrix must be square.");
%     assert(size(model.Q(1), 1) == size(model.Q(1), 2), "Q matrix must be square.");
%     assert(length(x_hat) == size(P, 1), "Number of state x_hat must be equal size of P matrix.");
%     assert(length(x_hat) == size(model.Q(1), 1), "Number of state x_hat must be equal size of Q matrix.");
%     assert(length(x_hat) == size(model.F(1), 1), "Number of state x_hat must be equal size of F matrix.");
%     assert(length(x_hat) == size(model.G(1), 1), "Number of state x_hat  must be equal row number of G matrix.");
%     assert(length(u) == size(model.G(1), 2), "Number of control variable u must be equal column number of G matrix.");
%     
%     F = model.F(Ts);
%     G = model.G(Ts);
%     Q = model.Q(Ts);
%     x_hat = F * x_hat + G * u;
%     P = F * P * F.' + Q;
% end
% 
% function [x_hat, P, K] = correct(x_hat, z, P, model, Ts)
% arguments
%     x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
%     z (:, 1) {mustBeNumeric, mustBeNonempty},
%     P (:, :) {mustBeNumeric, mustBeNonempty},
%     model struct {mustBeNonempty},
%     Ts (1, 1) {mustBeNumeric}
% end
%     assert(size(P, 1) == size(P, 2), "P matrix must be square.");
%     assert(size(model.R(1), 1) == size(model.R(1), 2), "R matrix must be square.");
%     assert(length(x_hat) == size(P, 1), "Number of state x_hat must be equal size of P matrix.");
%     assert(length(x_hat) == size(model.H, 2), "Number of state x_hat must be equal column number of H matrix.");
%     assert(length(z) == size(model.R(1), 1), "Number of measurements z must be equal size of R matrix.");
%     assert(length(z) == size(model.H, 1), "Number of measurements z must be equal row number of H matrix.");
% 
%     R = model.R(Ts);
%     K = P * model.H.' / (model.H * P * model.H.' + R);
%     x_hat = x_hat + K * (z - model.H * x_hat);
%     I = eye(size(x_hat));
%     Z = I - K * model.H;
%     P = Z * P * Z.' + K * R * K.';
% end

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
%% Confidence interval
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