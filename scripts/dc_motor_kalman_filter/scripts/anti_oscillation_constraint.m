function W_cl = anti_oscillation_constraint()
% ANTI_OSCILLATION_CONSTRAINT

syms s k_p L_a R_a k_t J b;
W_a = 1/(L_a * s) / (1 + R_a / (L_a * s));
W_m = 1/(J * s) / (1 + b / (J * s));
W_1 = (W_a * k_t * W_m) / (1 + W_a * k_t * W_m * k_t);
W_cl = (k_p * W_1) / (1 + k_p * W_1);
W_cl = simplify(W_cl);
end