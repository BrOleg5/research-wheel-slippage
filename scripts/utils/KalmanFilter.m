classdef KalmanFilter
    %KALMANFILTER Implementation Kalman Filter
    
    properties
        F {mustBeNumeric},
        G {mustBeNumeric},
        H {mustBeNumeric},
        Q {mustBeNumeric},
        R {mustBeNumeric},
        state_number (1, 1) {mustBeNumeric, mustBeInteger},
        control_number (1, 1) {mustBeNumeric, mustBeInteger},
        measurement_number (1, 1) {mustBeNumeric, mustBeInteger}
    end
    
    methods
        function obj = KalmanFilter(F, G, H, Q, R)
            %KALMANFILTER Construct an instance of this class
        
        arguments
            F {mustBeNumeric},
            G {mustBeNumeric},
            H {mustBeNumeric},
            Q {mustBeNumeric},
            R {mustBeNumeric}
        end
            assert(size(F, 1) == size(F, 2), "F matrix must be square.");
            assert(size(Q, 1) == size(Q, 2), "Q matrix must be square.");
            obj.state_number = size(F, 1);
            assert(obj.state_number == size(G, 1), "State number and G matrix row number must be the same.");
            assert(obj.state_number == size(H, 2), "State number and H matrix column number must be the same.");
            assert(obj.state_number == size(Q, 1), "State number and Q matrix size must be the same.");
            assert(size(R, 1) == size(R, 2), "R matrix must be square.");

            obj.control_number = size(G, 2);
            obj.measurement_number = size(H, 1);

            assert(size(R, 1) == obj.measurement_number, "Measurement number and R matrix size must be the same.");
        
            obj.F = F;
            obj.G = G;
            obj.H = H;
            obj.Q = Q;
            obj.R = R;
        end
        
        function [x_hat, P] = predict(obj, x_hat, u, P)
        arguments
            obj KalmanFilter,
            x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
            u (:, 1) {mustBeNumeric, mustBeNonempty},
            P (:, :) {mustBeNumeric, mustBeNonempty},
        end
            assert(size(P, 1) == size(P, 2), "P matrix must be square.");
            assert(length(x_hat) == obj.state_number, "Length of x_hat must be equal to state number");
            assert(size(P, 1) == obj.state_number, "P matrix size and state number must be the same.");
            assert(length(u) == obj.control_number, "Number of control variable u and control number must be the same.");
            
            x_hat = obj.F * x_hat + obj.G * u;
            P = obj.F * P * obj.F.' + obj.Q;
        end

        function [x_hat, P, K] = correct(obj, x_hat, z, P)
        arguments
            obj KalmanFilter,
            x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
            z (:, 1) {mustBeNumeric, mustBeNonempty},
            P (:, :) {mustBeNumeric, mustBeNonempty}
        end
            assert(size(P, 1) == size(P, 2), "P matrix must be square.");
            assert(length(x_hat) == obj.state_number, "Length of x_hat must be equal to state number");
            assert(size(P, 1) == obj.state_number, "P matrix size and state number must be the same.");
            assert(length(z) == obj.measurement_number, "Length of z and measurement number must be the same.");
        
            K = P * obj.H.' / (obj.H * P * obj.H.' + obj.R);
            x_hat = x_hat + K * (z - obj.H * x_hat);
            I = eye(length(x_hat));
            Z = I - K * obj.H;
        %     P = Z * P * Z.' + K * obj.R * K.';
            P = Z * P;
        end
    end
end

