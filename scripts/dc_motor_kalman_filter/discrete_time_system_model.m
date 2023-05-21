%% State space model
% x = [i; omega; T_l];
A = [-motor1.R_a/motor1.L_a, -(motor1.k_p + motor1.k_e)/motor1.L_a,           0; 
        motor1.k_t/motor1.J,                    -motor1.b/motor1.J,  1/motor1.J;
                          0,                                     0,           0];
B = [motor1.k_p/motor1.L_a; 0; 0];
C = [1, 0, 0;
     0, 1, 0];
%% Discretization
% syms Ts positive;
% Big = expm([A, B; zeros(1, 4)]*Ts);
% matlabFunction(vpa(Big(1:3, 1:3)), Ts, 'File', 'scripts/dc_motor_kalman_filter/utils/F_matrix');
% matlabFunction(vpa(Big(1:3, 4)), Ts, 'File', 'scripts/dc_motor_kalman_filter/utils/G_matrix');
% % Observation matrix
% H = C;

Ts = 30e-3; % s
Big = expm([A, B; zeros(1, 4)]*Ts);
F = Big(1:3, 1:3);
G = Big(1:3, 4);
% Observation matrix
H = C;
%% Check observability of discrete system
O = rank([H; H*F; H*(F)^2]);
if(O == size(A, 1))
    disp("Discrete system is fully observable.");
else
    warning("Discrete system isn't fully observable.");
end
%% Constructing the process noise matrix Q
% Q = [1e-7,      0,       0;
%         0, 2.2e-1,       0;
%         0,      0,    1e-9];
% 
% syms tau;
% Q_s = vpa(int(expm(A*tau) * Q * expm(A.'*tau), tau, [0, Ts]));
% % matlabFunction(Q_s, Ts, 'File', 'scripts/dc_motor_kalman_filter/utils/Q_matrix');
% Q_d = double(Q_s);

% Q_d = [4.83e-4,  1.02e-1,  -1.28e-5;
%        1.02e-1,    26.47, -3.31e-3;
%        -1.28e-5, -3.31e-3,  4.14e-7];

Q_d = [4.83e-4,  1.02e-2,  -1.28e-5;
       1.02e-2,    26.47, -3.31e-4;
       -1.28e-5, -3.31e-4,  4.14e-7];

if(~issymmetric(Q_d))
    error("Process noise covariance matrix Q doesn't symmetric.");
end
%% Constructing the measurement noise matrix R
% Encoder measurement noise
% link: https://www.researchgate.net/figure/Diagram-of-speed-measurement-noise-caused-by-encoder-pulse-counting-rotor-speed-50-rpm_fig1_347615973
% encoder_resolution = 500 * 2 * 2; % ticks per revolition (including two channel and rising and falling)
% f_mot_s = 900; % Hz
% T_mot_s = 1 / f_mot_s;
% omega_noise = 2 * pi / (encoder_resolution * T_mot_s); % rad/s
% omega_noise_sigma = omega_noise / 6;
% omega_variance = omega_noise_sigma.^2;
% 
% current_variance = (1.5/156/6).^2;
% 
% R = [current_variance,                0;
%                     0,   omega_variance];
% 
% % matlabFunction(R/Ts, Ts, 'File', 'scripts/dc_motor_kalman_filter/utils/R_matrix');
% R_d = R/Ts;

% R_d = [4.80e-3,     0;
%           0,   29.57];

R_d = [4.80e-3,     0;
            0,   29.57];

if(~issymmetric(R_d))
    error("Measurement noise covariance matrix R doesn't symmetric.");
end
%% Load experimental data
T = readtable("data\data_for_identification_without_integral.csv", "NumHeaderLines", 2);
T.t = T.t / 1000;
T.t = T.t + 0.03;
start_row = zeros(1, width(T));
T = [array2table(start_row, "VariableNames", T.Properties.VariableNames); T];
T.m1cur = sign(T.m1vel) .* T.m1cur;
T.m2cur = sign(T.m2vel) .* T.m2cur;
T.m3cur = sign(T.m3vel) .* T.m3cur;
%% Correlation matrix
idx = (T.m1vel ~= 0) | (T.m1cur ~= 0);
TT = T(idx, ["m1vel", "m1cur"]);
idx = find(diff(sign(TT.m1vel)) ~= 0);
idx = [0; idx; height(TT)];
intervals = zeros(height(TT), 1);
for i = 1:length(idx)-1
    intervals(idx(i)+1:idx(i+1)) = i;
end

idx = intervals == 2;
X = table2array(TT(idx, :));
plotmatrix(X);
[Rr, PValue, RL, RU] = corrcoef(X);
%% Estimate covariance
idx = diff(T.m1setvel) ~= 0;
interval_t = T.t(idx) + motor1.tau;
interval = zeros(size(interval_t));
for i = 1:length(interval_t)
    interval(i) = find(T.t > interval_t(i), 1);
end
C_history = zeros(2, 2);
for i = 1:round(length(interval)/2)
    start_idx = interval(2*i-1)+1;
    end_idx = interval(2*i)-1; 
    cur = T.m1cur(start_idx:end_idx);
    vel = T.m1vel(start_idx:end_idx);
    C_history(:, :, i) = cov(cur, vel);
end
C_mean = mean(C_history, 3);
C_mode = mode(C_history, 3);