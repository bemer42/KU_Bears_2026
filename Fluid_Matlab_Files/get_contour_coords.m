function [x_coords, y_coords] = get_contour_coords(x_vec, y_vec, Z_matrix, target_level)
try
    C = contourc(double(x_vec), double(y_vec), double(Z_matrix), [target_level, target_level]);
catch
    x_coords = NaN;
    y_coords = NaN;
    return;
end

x_coords = [];
y_coords = [];

i = 1;
total_cols = size(C, 2);
while i < total_cols
    num_points = C(2, i);
    x_seg = C(1, i + 1 : i + num_points);
    y_seg = C(2, i + 1 : i + num_points);
    
    x_coords = [x_coords, x_seg, NaN];
    y_coords = [y_coords, y_seg, NaN]; 
    i = i + num_points + 1;
end

% Return a single NaN instead of empty array if no contour was generated
if isempty(x_coords)
    x_coords = NaN;
    y_coords = NaN;
else
    x_coords(end) = [];
    y_coords(end) = [];
end
end