fig = figure("Name", "Kalman filter", "WindowState", "maximized");
y_name = ["Current, A", "Velocity, rad/s", "Load torque, N*m"];
t_name = ["m1cur", "m1vel"];
font_size = 12;
h = gobjects(n_x, 2);
tiledlayout(n_x, 1, "TileSpacing", "tight", "Padding", "compact");

for i = 1:n_x
    nexttile;
    if(i ~= 3)
        h(i, 1) = plot(T.t, T.(t_name(i)), 'LineWidth', 2, 'Color', 'blue');
        hold on;
    end
    h(i, 2) = plot(T.t, squeeze(x_hat_history(i, :, :)), 'LineWidth', 2, 'Color', 'red');
    grid on;
    xlabel("Time, s");
    ylabel(y_name(i));
    if(i ~= 3)
        legend(["Measure", "Filtered"], "Location", "best", "FontSize", font_size);
    else
        legend(["Estimated"], "Location", "southeast", "FontSize", font_size);
    end
end

for i = 1:n_x
    nexttile(i);
    xlim([0, 20]);
    ax = gca;
    ax.FontSize = font_size;
    lim = ax.YLim;
    ylim("manual");
    ylim(lim);
end

filename = "output_files/graphs/green_gray_kalman_filter.gif";
for j = 1:height(T)
    for i = 1:n_x
        if(i ~= 3)
            x = T.t(1:j);
            y = T.(t_name(i))(1:j);
            h(i, 1).XData = x;
            h(i, 1).YData = y;
            refreshdata(h(i, 1));
        end
        x = T.t(1:j);
        y = squeeze(x_hat_history(i, :, 1:j));
        h(i, 2).XData = x;
        h(i, 2).YData = y;
        refreshdata(h(i, 2));
    end
    drawnow;
    F = getframe(fig);
    im = frame2im(F);
    [imind,cm] = rgb2ind(im,256);
    if j == 1
      imwrite(imind, cm, filename, 'gif', 'Loopcount', 1, 'DelayTime', 30e-3);
    else
      imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 30e-3);
    end
end
close(fig);