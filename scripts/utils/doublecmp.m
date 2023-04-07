function tf = doublecmp(d1, d2, esp)
% DOUBLECMP Compare doubles
%   tf = DOUBLECMP(d1, d2, esp);

arguments
    d1 {mustBeNumeric, mustBeReal},
    d2 {mustBeNumeric, mustBeReal},
    esp (1, 1) {mustBeNumeric, mustBeReal} = 1e-6
end

tf = abs(d1 - d2) < esp;
end

