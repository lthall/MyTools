function [ Earth_X, Earth_Y, Earth_Z, DCM_BE ] = Vector_Body2Earth( Roll_A, Pitch_A, Yaw_A, Body_X, Body_Y, Body_Z )
%VECTOR_BODY2EARTH Summary of this function goes here
%   Detailed explanation goes here

% Earth_X = Body_X * cos(Pitch_A) * cos(Yaw_A) ...
%     + Body_Y * (sin(Roll_A) * sin(Pitch_A) * cos(Yaw_A) - cos(Roll_A) * sin(Yaw_A)) ...
%     + Body_Z * (cos(Roll_A) * sin(Pitch_A) * cos(Yaw_A) + sin(Roll_A) * sin(Yaw_A));
% 
% Earth_Y = Body_X * cos(Pitch_A) * sin(Yaw_A) ...
%     + Body_Y * (sin(Roll_A) * sin(Pitch_A) * sin(Yaw_A) + cos(Roll_A) * cos(Yaw_A)) ...
%     + Body_Z * (cos(Roll_A) * sin(Pitch_A) * sin(Yaw_A) - sin(Roll_A) * cos(Yaw_A));
% 
% Earth_Z = - Body_X * sin(Pitch_A) ...
%     + Body_Y * sin(Roll_A) * cos(Pitch_A) ...
%     + Body_Z * cos(Roll_A) * cos(Pitch_A);

DCM_BE = [cos(Pitch_A) * cos(Yaw_A) ...
    (sin(Roll_A) * sin(Pitch_A) * cos(Yaw_A) - cos(Roll_A) * sin(Yaw_A)) ...
    (cos(Roll_A) * sin(Pitch_A) * cos(Yaw_A) + sin(Roll_A) * sin(Yaw_A)); ...
    ...
    cos(Pitch_A) * sin(Yaw_A) ...
    (sin(Roll_A) * sin(Pitch_A) * sin(Yaw_A) + cos(Roll_A) * cos(Yaw_A)) ...
    (cos(Roll_A) * sin(Pitch_A) * sin(Yaw_A) - sin(Roll_A) * cos(Yaw_A)); ...
    ...
    - sin(Pitch_A) ...
    sin(Roll_A) * cos(Pitch_A) ...
    cos(Roll_A) * cos(Pitch_A)];
    
    Earth_XYZ = DCM_BE * [Body_X; Body_Y; Body_Z];
    
    Earth_X = Earth_XYZ(1);
    Earth_Y = Earth_XYZ(2);
    Earth_Z = Earth_XYZ(3);
end

