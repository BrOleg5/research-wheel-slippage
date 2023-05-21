%% Configuration
show_plot = true;
anderson_darling_test = true; 
kolmogorov_smirnov_test = true;
%% Normal distribution test
T = readtable("data\data_for_identification_without_integral.csv", "NumHeaderLines", 2);

m = 1;
mstr = "m" + num2str(m);

vel_step_idx = find(abs(diff(T.(mstr + "vel"))) > 18);

i = 2;
assert(i <= (length(vel_step_idx)-1)/2, "Index out of range");
idx = vel_step_idx(2*i-1)+1:vel_step_idx(2*i);
time = T.t(idx);
value = T.(mstr + "vel")(idx);

for m = 1
    if(show_plot)
        fig = figure("Name", sprintf("Motor %d", m));
        histogram(value);
%         plot(time, value);
%         xlabel("Current, A");
        xlabel("Velocity, rad/s");
        ylabel("Sample count");
    end
    if(anderson_darling_test)
        fprintf("Anderson–Darling test\n");
        [h,p,adstat,cv] = adtest(value);
        fprintf("\tMotor %d\n\t\tHypothesis: %d\n\t\tp-value: %f\n\t\tcritical value: %f\n", m, h, p, cv);
    end
    if(kolmogorov_smirnov_test)
        fprintf(" Kolmogorov-Smirnov test\n");
        mu = mean(value);
        S = std(value);
        x = (value - mu)/S;
        [h,p,ksstat,cv] = kstest(x);
        fprintf(strcat("\tMotor %d\n\t\tMean: %f\n\t\tStandard deviation: %f\n\t\t", ...
                       "Hypothesis: %d\n\t\tp-value: %f\n\t\tcritical value: %f\n"), m, mu, S, h, p, cv);
    end
end