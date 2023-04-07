% link: https://www.quora.com/How-do-I-generate-n-visually-distinct-RGB-colours-in-Python
function rgb_colors = generate_distrinct_colors(n)
    arguments
        n (1, 1) {mustBeNumeric}
    end
    hue_partition = 1 / (n + 1);
    hue_array = 0:hue_partition:1;
    hue_array = hue_array';
    sat_array = ones(size(hue_array));
    val_array = ones(size(hue_array));
    hsv = [hue_array, sat_array, val_array];
    rgb_colors = hsv2rgb(hsv);
end