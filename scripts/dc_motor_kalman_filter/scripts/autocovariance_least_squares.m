%% Load experimental data
T = readtable("data\data_for_identification_without_integral.csv", "NumHeaderLines", 2);

T.t = T.t / 1000;
T.t = T.t + 0.03;
start_row = zeros(1, width(T));
T = [array2table(start_row, "VariableNames", T.Properties.VariableNames); T];
T.m1cur = sign(T.m1vel) .* T.m1cur;
T.m2cur = sign(T.m2vel) .* T.m2cur;
T.m3cur = sign(T.m3vel) .* T.m3cur;

% idx = (T.m1vel ~= 0) | (T.m1cur ~= 0);
% T = T(idx, :);
% idx = find(diff(sign(T.m1vel)) ~= 0);
% idx = [0; idx; height(T)];
% intervals = zeros(height(T), 1);
% for i = 1:length(idx)-1
%     intervals(idx(i)+1:idx(i+1)) = i;
% end
% 
% idx = intervals == 8;
% T = T(idx, :);

u = [];
z = [];
% Remove time delay
i_prev = 1;
for i = 1:height(T)
    if((T.t(i) - T.t(i_prev)) > motor1.tau)
        i_prev = i_prev + 1;
    end
    u(i) = T.m1setvel(i_prev);
    z(:, i) = [T.m1cur(i); T.m1vel(i)];
end
%% Setup Autocovariance Least-Squares problem

% initial guesses
% Encoder measurement noise
% link: https://www.researchgate.net/figure/Diagram-of-speed-measurement-noise-caused-by-encoder-pulse-counting-rotor-speed-50-rpm_fig1_347615973
encoder_resolution = 500 * 2 * 2; % ticks per revolition (including two channel and rising and falling)
f_mot_s = 900; % Hz
T_mot_s = 1 / f_mot_s;
omega_noise = 2 * pi / (encoder_resolution * T_mot_s); % rad/s
omega_noise_sigma = omega_noise / 6;
omega_variance = omega_noise_sigma.^2;

current_variance = (1.5/156/6).^2;

G_hat = eye(3);
Qw_hat = [1e-7,      0,       0;
             0,   2e-1,       0;
             0,      0,    1e-9];
Rv_hat= [current_variance,                0;
                        0,   omega_variance];

L = dlqe(F, G_hat, H, Qw_hat, Rv_hat);

model.A = F;
model.B = G;
model.C = H;
model.G = G_hat;
model.xhat0 = [z(:, 1); 0];

data.datapts = length(u);
data.yk = z;
data.uk = u;

N = 15;

estimator.L = L;
estimator.Q = Qw_hat;
estimator.R = Rv_hat;

[Qest_cell, Rest_cell] = als_sdp_mrQ(data, N, model, estimator);
Qest1 = Qest_cell{1}
Rest1 = Rest_cell{1}