function surface_color_dict = get_surface_color_dict()
surface_type = ["table", "gray", "green"];
color = ["black", "red", "blue"];
surface_color_dict = containers.Map(surface_type, color);
end