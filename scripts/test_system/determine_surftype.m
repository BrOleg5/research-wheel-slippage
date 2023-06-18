function determine_surftype(top_view_image_path, input_dataset_folder)

arguments
    top_view_image_path {mustBeText, mustBeFile},
    input_dataset_folder {mustBeTextScalar, mustBeFolder}
end

surfType = categorical(["table", "green", "gray"]);

frame = imread(top_view_image_path);
catMap = findSurfaceBorder(frame, surfType);

list_of_files = dir(input_dataset_folder + "/exp*.csv");
for j = 1:length(list_of_files)
    input_dataset_file = fullfile(list_of_files(j).folder, list_of_files(j).name);
    write_surftype(catMap, input_dataset_file);
end
end
%% Local functions
function write_surftype(catMap, input_dataset_file)
arguments
    catMap categorical,
    input_dataset_file {mustBeTextScalar, mustBeFile}
end

surfType = categorical(["table", "green", "gray"]);

[robotParameters, ~] = get_robot_parameters();

fid = fopen("config\camera_params.json");
raw = fread(fid, inf);
fclose(fid);
camera_param = jsondecode(char(raw'));
cameraParameters = struct('x', camera_param.pixel_resolution_x/1e3, ...
                          'y', camera_param.pixel_resolution_y/1e3);

T = readtable(input_dataset_file, "NumHeaderLines", 2);

% Read headlines of file
fid = fopen(input_dataset_file);
column_description = fgetl(fid);
column_units = fgetl(fid);
fclose(fid);

T.w1surf = repelem(surfType(1), height(T), 1);
T.w2surf = T.w1surf;
T.w3surf = T.w1surf;
pixelSize = mean([cameraParameters.x, cameraParameters.y]);
for j = 1:height(T)
    wheel_coord = caclWheelCoordinate([T.xpos(j); T.ypos(j)], T.ang(j), ...
        robotParameters.robot_radius / pixelSize);
    wheel_coord = round(wheel_coord);
    for k = 1:3
        T(j, "w" + num2str(k) + "surf").Variables = catMap(wheel_coord(k, 2), wheel_coord(k, 1));
    end
end

column_description = sprintf("%s;1st wheel surface type; 2nd wheel surface type;3rd wheel surface type", column_description);
column_units = sprintf("%s;none;none;none", column_units);

% Write headlines in file
fid = fopen(input_dataset_file, "w");
fprintf(fid, "%s\n%s", column_description, column_units);
fclose(fid);
% Write table in file
writetable(T, input_dataset_file, 'Delimiter', ';', 'WriteMode', 'Append', ...
           'WriteVariableNames', true);
end

function catMap = findSurfaceBorder(frame, surface_type)
arguments
    frame (:, :, 3) {mustBeNonempty},
    surface_type (1, :) categorical {mustBeNonempty}
end
    app.fig = uifigure("Name", "Find surface border", "Position", [0, 0, 1920*0.7, 1080*0.7]);
    widthBorder = 0;
    heightBorder = 30;
    axesPos = [widthBorder, heightBorder, app.fig.InnerPosition(3) - 2*widthBorder, ...
               app.fig.InnerPosition(4) - 2*heightBorder];
    app.ax = uiaxes(app.fig, "Position", axesPos);
    hold(app.ax, "on");
    imshow(frame, 'parent', app.ax, 'Border', 'tight');

    app.positionLabel = uilabel(app.fig, "Text", "Current position:" , "Position", [120, 90, 150, 30], ...
                                "FontSize", 14);
    app.positionEdit = uieditfield(app.fig, "Editable", "off", "Value", "", ...
                                   "Position", [225, 90, 110, 30], "FontSize", 14, ...
                                   'HorizontalAlignment', 'center');

    app.PointLabel(1) = uilabel(app.fig, "Text", "First border point:" , ...
                                  "Position", [400, 90, 150, 30], ...
                                  "FontSize", 14);
    app.PointEdit(1) = uieditfield(app.fig, "Editable", "off", "Value", "", ...
                                   "Position", [515, 90, 110, 30], "FontSize", 14, ...
                                   'HorizontalAlignment', 'center');

    app.PointLabel(2) = uilabel(app.fig, "Text", "Second border point:" , ...
                                "Position", [640, 90, 150, 30], ...
                                "FontSize", 14);
    app.PointEdit(2) = uieditfield(app.fig, "Editable", "off", "Value", "", ...
                                   "Position", [775, 90, 110, 30], "FontSize", 14, ...
                                   'HorizontalAlignment', 'center');

    app.confirmButton = uibutton(app.fig, "FontSize", 14, "Text", "Confirm", ...
                                 "Position", [900, 90, 110, 30]);
    app.resetButton = uibutton(app.fig, "FontSize", 14, "Text", "Reset", ...
                                 "Position", [900, 55, 110, 30]);
    
    % link: https://uk.mathworks.com/matlabcentral/answers/583454-how-draw-lines-between-two-coordinates-saved-with-mouse-s-click-in-uiaxes-appdesigner
    app.pm = iptPointerManager(app.fig);
    app.pm.enterFcn = [];
    app.pm.exitFcn = @(~, ~) set(app.positionEdit, "Value",  "");
    app.pm.traverseFcn = @(~, ~) set(app.positionEdit, "Value", ...
                                     sprintf("%d, %d", round(app.ax.CurrentPoint(1, 1:2))));
    iptSetPointerBehavior(app.ax, app.pm);
    app.fig.WindowButtonMotionFcn = @(~, ~) NaN;

    app.fig.UserData.borderCoordinate = [NaN, NaN;
                                         NaN, NaN];
    app.fig.UserData.pointCount = 1;

    app.fig.WindowButtonDownFcn = @(~, ~) mouseClickCallback(app);

    app.confirmButton.ButtonPushedFcn = @(~, ~) confirmButtonCallback(app);
    app.resetButton.ButtonPushedFcn = @(~, ~) resetButtonCallback(app);

    uiwait(app.fig);
    X = app.fig.UserData.borderCoordinate(:, 1);
    Y = app.fig.UserData.borderCoordinate(:, 2);
    [slope, intersect] = estimateSlopeIntersect(X, Y);
    frameWidth = size(frame, 2);
    frameHeight = size(frame, 1);
    [imageColumn, imageRow] = meshgrid(1:frameWidth, 1:frameHeight);
    splitImageByLine = @(X, Y, slope, intersect) X > ((Y - intersect) ./ slope);
    outMap = splitImageByLine(imageColumn, imageRow, slope, intersect);
    string_surface = string(surface_type);
    ddPosition = [50, 400, 110, 30];
    app.surfaseDropDown(1) = uidropdown(app.fig, "Items", string_surface, "Value", string_surface(1), ...
                                        "Position", ddPosition);
    ddPosition = [1000, 400, 110, 30];
    app.surfaseDropDown(2) = uidropdown(app.fig, "Items", string_surface, "Value", string_surface(1), ...
                                        "Position", ddPosition);

    uiwait(app.fig);
    surface_categoty = categorical({app.surfaseDropDown(1).Value, app.surfaseDropDown(2).Value}, ...
                                   categories(surface_type));
    catMap = repelem(surface_categoty(1), frameHeight, frameWidth);
    catMap(outMap) = surface_categoty(2);
    close(app.fig);
end

function mouseClickCallback(app)
arguments
    app struct
end
    border = app.fig.UserData;
    if(border.pointCount < 3)
        border.borderCoordinate(border.pointCount, :) = round(app.ax.CurrentPoint(1, 1:2));
        app.PointEdit(border.pointCount).Value = sprintf("%d, %d", ...
                                                         border.borderCoordinate(border.pointCount, :));
        if(border.pointCount == 1)
            plot(app.ax, border.borderCoordinate(border.pointCount, 1), ...
                 border.borderCoordinate(border.pointCount, 2), ...
                 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8, 'MarkerEdgeColor', 'red');
        elseif(border.pointCount == 2)
            plot(app.ax, border.borderCoordinate(border.pointCount, 1), ...
                 border.borderCoordinate(border.pointCount, 2), ...
                 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8, 'MarkerEdgeColor', 'red');
            plot(app.ax, border.borderCoordinate(:, 1), ...
                 border.borderCoordinate(:, 2), 'Color', 'red', 'LineWidth', 1);
        end
        border.pointCount = border.pointCount + 1;
        app.fig.UserData = border;
    end
end

function resetButtonCallback(app)
    app.UserData.borderCoordinate = [NaN, NaN;
                                     NaN, NaN];
    app.fig.UserData.pointCount = 1;
    for i=1:2
        app.PointEdit(i).Value = "";
    end
    ClearLinesFromAxes(app.ax);
end

function ClearLinesFromAxes(ax)
    axesHandlesToChildObjects = findobj(ax, 'Type', 'line');
    if ~isempty(axesHandlesToChildObjects)
      delete(axesHandlesToChildObjects);
    end  
end

function confirmButtonCallback(app)
        uiresume(app.fig);
end

function [funcSlope, funcIntersect] = obtainSlopeIntersect()
    syms a b real;
    syms X Y [2, 1] real;
    eq = Y == a*X + b;
    [a, b] = solve(eq, [a, b]);
    a = simplify(a);
    b = simplify(b);
    funcSlope = matlabFunction(a, 'Vars', [X, Y]);
    funcIntersect = matlabFunction(b, 'Vars', [X, Y]);
end

function [slope, intersect] = estimateSlopeIntersect(X, Y)
arguments
    X (2, 1) {mustBeNumeric, mustBeNonempty},
    Y (2, 1) {mustBeNumeric, mustBeNonempty}
end
    [funcSlope, funcIntersect] = obtainSlopeIntersect();
    slope = funcSlope(X(1), X(2), Y(1), Y(2));
    intersect = funcIntersect(X(1), X(2), Y(1), Y(2));
end

function outputType = majorCategory(inType)
arguments
    inType categorical {mustBeNonempty}
end
    n = size(inType, 2);
    outputType = repelem(inType(1, 1), n, 1);
    for i = 1:n
        [N, categories] = histcounts(inType(:, i));
        [~, idx] = max(N);
        outputType(i) = categories(idx);
    end
end

function findWheelDebugApp(vidObj, T, robotParameters, pixelSize, heightCamera)
arguments
    vidObj audiovideo.internal.IVideoReader,
    T (:, 25) table,
    robotParameters struct,
    pixelSize (1, 1) double,
    heightCamera (1, 1) double
end
    app.fig = uifigure("Name", "Find wheel", "Position", [0, 0, 1920*0.7, 1080*0.7]);
    widthBorder = 0;
    heightBorder = 30;
    axesPos = [widthBorder, heightBorder, app.fig.InnerPosition(3) - 2*widthBorder, ...
               app.fig.InnerPosition(4) - 2*heightBorder];
    app.ax = uiaxes(app.fig, "Position", axesPos);
    hold(app.ax, "on");

    app.frameLabel = uilabel(app.fig, "Text", "Current frame:" , "Position", [120, 50, 150, 30], ...
                                "FontSize", 14);
    app.frameEdit = uieditfield(app.fig, "numeric", "Editable", "on", "Value", 1, ...
                                   "Position", [225, 50, 110, 30], "FontSize", 14, ...
                                   'HorizontalAlignment', 'center');
    app.frameSlider = uislider(app.fig, "Value", 1, "Limits", [1, vidObj.NumFrames], ...
                               "Position", [125, 115, 1050, 3]);

    app.tableLabel = uilabel(app.fig, "Text", "Current row:" , "Position", [420, 50, 150, 30], ...
                                "FontSize", 14);
    app.tableEdit = uieditfield(app.fig, "numeric", "Editable", "on", "Value", 1, ...
                                   "Position", [500, 50, 110, 30], "FontSize", 14, ...
                                   'HorizontalAlignment', 'center');
    app.tableSlider = uislider(app.fig, "Value", 1, "Limits", [1, height(T)], ...
                               "Position", [125, 35, 1050, 3]);

    frame = read(vidObj, app.frameSlider.Value);
    imshow(frame, 'parent', app.ax, 'Border', 'tight');

    app.confirmButton = uibutton(app.fig, "FontSize", 14, "Text", "Exit", ...
                                 "Position", [1200, 45, 110, 30]);

    frameWidth = size(frame, 2);
    frameHeight = size(frame, 1);
    [imageColumn, imageRow] = meshgrid(1:frameWidth, 1:frameHeight);

    inEllipse = @(X, Y, R, alpha, beta) (X.^2./cos(alpha) + Y.^2/cos(beta))./R.^2 <= 1;
    alpha = atan2(pixelSize * T.xpos, heightCamera);
    beta = atan2(pixelSize * T.ypos, heightCamera);
        
    
    cicleBoundary = getRobotBoundary(app.tableEdit.Value);
    hp = plot(app.ax, cicleBoundary(:, 2), cicleBoundary(:, 1), "Color", "white", "LineWidth", 1);
    cp = plot(app.ax, T.xpos(app.tableEdit.Value), T.ypos(app.tableEdit.Value), "LineStyle", "none", ...
              "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", "green");
    wheel_coord = caclWheelCoordinate([T.xpos(app.tableEdit.Value); T.ypos(app.tableEdit.Value)], ...
                                      T.ang(app.tableEdit.Value), robotParameters.R / pixelSize);
    wh = plot(app.ax, wheel_coord(:, 1), wheel_coord(:, 2), "LineStyle", "none", ...
              "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", "red");

    app.frameSlider.ValueChangedFcn = @(~, ~) updateFrameSlider();
    app.frameSlider.ValueChangingFcn = @(~, event) sliderValueChanging(event);
    app.frameEdit.ValueChangedFcn = @(~, ~) updateFrameEdit();

    app.tableSlider.ValueChangedFcn = @(~, ~) updateTableSlider();
    app.tableSlider.ValueChangingFcn = @(~, event) tableSliderValueChanging(event);
    app.tableEdit.ValueChangedFcn = @(~, ~) updateTableEdit();

    app.confirmButton.ButtonPushedFcn = @(~, ~) confirmButtonCallback(app);

    uiwait(app.fig);
    close(app.fig);

    function updateTableSlider()
        app.tableEdit.Value = round(app.tableSlider.Value);
        boundary = getRobotBoundary(app.tableEdit.Value);
        set(hp, 'XData', boundary(:, 2), 'YData', boundary(:, 1));
        set(cp, 'XData', T.xpos(app.tableEdit.Value), 'YData', T.ypos(app.tableEdit.Value));
        wheel_coord = caclWheelCoordinate([T.xpos(app.tableEdit.Value); T.ypos(app.tableEdit.Value)], ...
                                      T.ang(app.tableEdit.Value), robotParameters.R / pixelSize);
        set(wh, 'XData', wheel_coord(:, 1), 'YData', wheel_coord(:, 2));
    end
    
    function updateTableEdit()
        app.tableSlider.Value = app.tableEdit.Value;
        boundary = getRobotBoundary(app.tableEdit.Value);
        set(hp, 'XData', boundary(:, 2), 'YData', boundary(:, 1));
        set(cp, 'XData', T.xpos(app.tableEdit.Value), 'YData', T.ypos(app.tableEdit.Value));
        wheel_coord = caclWheelCoordinate([T.xpos(app.tableEdit.Value); T.ypos(app.tableEdit.Value)], ...
                                      T.ang(app.tableEdit.Value), robotParameters.R / pixelSize);
        set(wh, 'XData', wheel_coord(:, 1), 'YData', wheel_coord(:, 2));
    end
    
    function tableSliderValueChanging(event)
        app.tableEdit.Value = round(event.Value);
    end
    
    function boundary = getRobotBoundary(index)
        circlePixel = inEllipse(imageColumn - T.xpos(index), imageRow - T.ypos(index), ...
                                robotParameters.R / pixelSize, alpha(index), beta(index));
        boundary = bwboundaries(circlePixel, "noholes");
        boundary = boundary{1};
    end

    function updateFrameSlider()
        app.frameEdit.Value = round(app.frameSlider.Value);
        frame = read(vidObj, app.frameEdit.Value);
        imshow(frame, 'parent', app.ax, 'Border', 'tight');
        cicleBoundary = getRobotBoundary(app.tableEdit.Value);
        hp = plot(app.ax, cicleBoundary(:, 2), cicleBoundary(:, 1), "Color", "white", "LineWidth", 1);
        cp = plot(app.ax, T.xpos(app.tableEdit.Value), T.ypos(app.tableEdit.Value), "LineStyle", "none", ...
                  "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", "green");
        wheel_coord = caclWheelCoordinate([T.xpos(app.tableEdit.Value); T.ypos(app.tableEdit.Value)], ...
                                      T.ang(app.tableEdit.Value), robotParameters.R / pixelSize);
        wh = plot(app.ax, wheel_coord(:, 1), wheel_coord(:, 2), "LineStyle", "none", ...
              "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", "red");
    end
    
    function updateFrameEdit()
        frame = read(vidObj, app.frameEdit.Value);
        app.frameSlider.Value = app.frameEdit.Value;
        imshow(frame, 'parent', app.ax, 'Border', 'tight');
        cicleBoundary = getRobotBoundary(app.tableEdit.Value);
        hp = plot(app.ax, cicleBoundary(:, 2), cicleBoundary(:, 1), "Color", "white", "LineWidth", 1);
        cp = plot(app.ax, T.xpos(app.tableEdit.Value), T.ypos(app.tableEdit.Value), "LineStyle", "none", ...
                  "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", "green");
        wheel_coord = caclWheelCoordinate([T.xpos(app.tableEdit.Value); T.ypos(app.tableEdit.Value)], ...
                                      T.ang(app.tableEdit.Value), robotParameters.R / pixelSize);
        wh = plot(app.ax, wheel_coord(:, 1), wheel_coord(:, 2), "LineStyle", "none", ...
              "Marker", ".", "MarkerSize", 10, "MarkerEdgeColor", "red");
    end
    
    function sliderValueChanging(event)
        app.frameEdit.Value = round(event.Value);
        cicleBoundary = getRobotBoundary(app.tableEdit.Value);
        hp = plot(app.ax, cicleBoundary(:, 2), cicleBoundary(:, 1), "Color", "white", "LineWidth", 1);
    end
end

function wheel_coordinate = caclWheelCoordinate(position, orientation, robotRadius)
arguments
    position (2, 1) {mustBeNumeric, mustBeNonempty},
    orientation (1, 1) {mustBeNumeric},
    robotRadius (1, 1) {mustBeNumeric}
end
    wheel_angle = [pi/3, pi, 5*pi/3];
    wheel_coordinate = zeros(3, 2);
    wheel_coordinate(:, 1) = position(1) - robotRadius * sin(deg2rad(orientation) + wheel_angle);
    wheel_coordinate(:, 2) = position(2) - robotRadius * cos(deg2rad(orientation) + wheel_angle);
end