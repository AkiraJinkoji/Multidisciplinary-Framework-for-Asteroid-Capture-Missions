
% this function would be used to calculate the initial attitude of the
% asteroid so that the spacecraft points its capture system toward the
% asteroid spin axis

function quaternion = compute_quaternion_attitude_from_2vec(v1, v2)
% Normalize the vectors
v1 = v1 / norm(v1);
v2 = v2 / norm(v2);

% Calculate the axis of rotation (cross product)
axis = cross(v1, v2);
axis = axis / norm(axis); % Normalize the axis

% Calculate the angle of rotation (dot product)
cos_theta = dot(v1, v2);
angle = acos(cos_theta); % Angle in radians

%Convert angle to degrees (optional)
angle_deg = rad2deg(angle);

% Calculate quaternion components
qw = cos(angle / 2);
qx = sin(angle / 2) * axis(1);
qy = sin(angle / 2) * axis(2);
qz = sin(angle / 2) * axis(3);

% Combine into quaternion
quaternion = [qx; qy; qz; qw];

% Display the results
% disp('Axis of rotation:');
% disp(axis);
% disp('Angle of rotation (radians):');
% disp(angle_deg);
% disp('Quaternion:');
% disp(quaternion);

end