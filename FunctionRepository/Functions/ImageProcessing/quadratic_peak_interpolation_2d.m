function [move_sub_x,move_sub_y] = quadratic_peak_interpolation_2d(parabolic_fit_values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the shift on sub-pixel level required to center the maximum of a 3x3
% matrix using parabolic curve fitting and mean square approximation
%
% Assume that the input values are centered around (x,y)=(0,0) with the
% maximum value in the center and the input matrix is the center pixel with
% the surrounding pixels. 
% Find the least square solution to the overdetermined system assuming that
% b=c1+c2x+c3y+c4x^2+c5y^2+c6xy
%
% Output: 
%     [move_sub_x, move_sub_y] - Shift necessary to move the maximum of the
%     fitted parabola to the center of the matrix. (Will be <=1 or >=-1)
% Inputs:
%     parabolic_fit_values - Values to be used when calculating a fitted
%     quadratic. (must be a 3x3 matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[X,Y]=meshgrid(-1:1,-1:1);

b=parabolic_fit_values(:);

A=[ones(length(X(:)),1) X(:) Y(:) X(:).^2 Y(:).^2 X(:).*Y(:)];

c=pinv(A)*b;

% By setting derivatives to zero I got x_max=(2*c2*c5-c3*c6)/(c6^2-4*c4*c5)
% and y_max=(2*c3*c4-c2*c6)/(c6^2-4*c4*c5) Use this to find the sub pixel
% movement

move_sub_x=(2*c(2)*c(5)-c(3)*c(6))/(c(6)^2-4*c(4)*c(5));
move_sub_y=(2*c(3)*c(4)-c(2)*c(6))/(c(6)^2-4*c(4)*c(5));


% show_images_for_poster(X,Y,b,c,move_sub_x,move_sub_y)
% Make sure the algorithm doesn't fail too much----

% End Conditions --- The image is no allowed to move more than 1 pixel even
% if the peak is centered outside our area correlation
if move_sub_x>1
    move_sub_x=1;
end
if move_sub_x<-1
    move_sub_x=-1;
end
if move_sub_y>1
    move_sub_y=1;
end
if move_sub_y<-1
    move_sub_y=-1;
end