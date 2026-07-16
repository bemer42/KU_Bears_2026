% Contour Creator
function [x_coords, y_coords] = get_contour_coords(x_vec, y_vec, Z_matrix, target_level)


% Step 1: Compute contour data behind the scenes
try
    C = contourc(double(x_vec), double(y_vec), double(Z_matrix), [target_level, target_level]);
catch
    % Return empty arrays if contourc fails or target is out of bounds
    x_coords = [];
    y_coords = [];
    return;
end

% Pre-allocate coordinate arrays
x_coords = [];
y_coords = [];

% Step 2: Loop through the contour matrix and unpack every segment
i = 1;
total_cols = size(C, 2);

while i < total_cols
    % C(1, i) is the contour level value
    % C(2, i) is the number of points in this specific segment
    num_points = C(2, i);

    % Extract coordinates for this segment
    x_seg = C(1, i + 1 : i + num_points);
    y_seg = C(2, i + 1 : i + num_points);

    % Append the segment coordinates followed by a NaN separator
    x_coords = [x_coords, x_seg, NaN];
    y_coords = [y_coords, y_seg, NaN]; 

    % Move the pointer to the next segment header
    i = i + num_points + 1;
end

% Strip the trailing NaN if we successfully found segments
if ~isempty(x_coords)
    x_coords(end) = [];
    y_coords(end) = [];
end
end