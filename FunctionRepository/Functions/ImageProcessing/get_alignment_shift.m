function [move_total_x,move_total_y] = get_alignment_shift(in_pic1,in_pic2, max_shift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs the total x,y shifts required to align two images (Only works if
% initial alignment disparity is less than max_shift)
%
% ***Usage Instructions***: To align two images use:
% alligned_pic = shift_image(in_pic2,move_total_x,move_total_y);
%
% Output: 
%     [move_total_x, move_total_y] - Shift necessary to move the second
%     input image to align with the first input image
% Inputs:
%     in_pic1 - Image with which to align second input image
%     in_pic2 - Image to align with the first input image
%     max_shift - maximum allowable shift that can take place between the two images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of steps to check in each direction during correlation for various images
ds4_steps = floor(max_shift/4);
ds2_steps = 1;
orig_steps = 1;
sub_steps = 1;

% Create a pyramid of lower resolution images------------------

% lpf = low pass filtered, ds'x' = downsampled by factor x
lpf=lowpass(11,0.5);
lpf_pic1 = filter2(lpf, in_pic1);
lpf_pic2 = filter2(lpf, in_pic2);

ds2_pic1 = lpf_pic1(1:2:end,1:2:end);
ds2_pic2 = lpf_pic2(1:2:end,1:2:end);

lpf_ds2_pic1 = filter2(lpf, ds2_pic1);
lpf_ds2_pic2 = filter2(lpf, ds2_pic2);

ds4_pic1 = lpf_ds2_pic1(1:2:end,1:2:end);
ds4_pic2 = lpf_ds2_pic2(1:2:end,1:2:end);


% Correlate with downsample factor x4 ------------------------------------

total_ds4_steps=2*ds4_steps+1;
% From 1+num_steps --> end-num_steps
compare_pic1 = ds4_pic1(1+ds4_steps:end-ds4_steps,1+ds4_steps:end-ds4_steps); % Remove outer borders (depending on number of steps to take)

[row_corr,col_corr]=size(compare_pic1);
% Get Average (To remove background) (Assume: Background is ~ avg. intensity of image)
compare_avg1 = sum(compare_pic1(:))/(row_corr*col_corr);

for x=1:total_ds4_steps
    for y=1:total_ds4_steps
        % Get Average (To remove background)        
        ds4_pic2_avg = sum(sum(ds4_pic2(y:row_corr+y-1,x:col_corr+x-1)))/(row_corr*col_corr);
        % Correlation matrix
        gamma_4(y,x)=sum(sum((ds4_pic2(y:row_corr+y-1,x:col_corr+x-1)-ds4_pic2_avg).*(compare_pic1-compare_avg1)));
    end
end

% Find maximum correlation value
[max_junk,max_ind]=max(gamma_4(:));
% Get coordinates of max correlation value (b/w 1 and total_ds4_steps)
[I,J]=ind2sub([total_ds4_steps,total_ds4_steps],max_ind);
% Calculate shift ex. for ds4_steps = 2, gives 5x5 matrix and [3,3] is
% center corresponding to no movement
move_4_x=(J-ds4_steps-1);
move_4_y=(I-ds4_steps-1);

% Correlate with downsample factor x2 ------------------------------------

%Remove Border to allow for shift (include maximum possible movement from
%ds4 correlation (2*ds4_steps+1) (the 2 corresponding to the sample factor difference) 
%and maximum possible movement from ds2 correlation (ds2_steps)
compare_pic1=ds2_pic1(1+(2*ds4_steps+ds2_steps):end-(2*ds4_steps+ds2_steps),1+(2*ds4_steps+ds2_steps):end-(2*ds4_steps+ds2_steps));

[row_corr,col_corr]=size(compare_pic1);
x_pos=1;

% Perform downsample x2 correlation around best correlation from downsample x4 pictures
compare_avg2=sum(compare_pic1(:))/(row_corr*col_corr);

for x=(1+(2*ds4_steps+1))+2*move_4_x-ds2_steps:(1+(2*ds4_steps+1))+2*move_4_x+ds2_steps
    y_pos=1;
    for y=(1+(2*ds4_steps+1))+2*move_4_y-ds2_steps:(1+(2*ds4_steps+1))+2*move_4_y+ds2_steps
        
        ds2_pic2_avg=sum(sum(ds2_pic2(y:row_corr+y-1,x:col_corr+x-1)))/(row_corr*col_corr);
        gamma_2(y_pos,x_pos)=sum(sum((ds2_pic2(y:row_corr+y-1,x:col_corr+x-1)-ds2_pic2_avg).*(compare_pic1-compare_avg2)));
        y_pos=y_pos+1;
        
    end
    x_pos=x_pos+1;
end

[max_junk,max_ind]=max(gamma_2(:));
[I,J]=ind2sub([2*ds2_steps+1,2*ds2_steps+1],max_ind);
move_2_x=(J-ds2_steps-1);
move_2_y=(I-ds2_steps-1);

% Correlate with original image-------------------------------------------------------

compare_pic1=in_pic1(1+(4*ds4_steps+2*ds2_steps+orig_steps):end-(4*ds4_steps+2*ds2_steps+orig_steps),1+(4*ds4_steps+2*ds2_steps+orig_steps):end-(4*ds4_steps+2*ds2_steps+orig_steps));
[row_corr,col_corr]=size(compare_pic1);
x_pos=1;
compare_avg1=sum(compare_pic1(:))/(row_corr*col_corr);
        
for x=1+(4*ds4_steps+2*ds2_steps+orig_steps)+4*move_4_x+2*move_2_x-1:(1+(4*ds4_steps+2*ds2_steps+orig_steps))+4*move_4_x+2*move_2_x+1
    y_pos=1;
    for y=(1+(4*ds4_steps+2*ds2_steps+orig_steps))+4*move_4_y+2*move_2_y-1:(1+(4*ds4_steps+2*ds2_steps+orig_steps))+4*move_4_y+2*move_2_y+1
        
        in_pic2_avg=sum(sum(in_pic1(y:row_corr+y-1,x:col_corr+x-1)))/(row_corr*col_corr);
        gamma_1(y_pos,x_pos)=sum(sum((in_pic2(y:row_corr+y-1,x:col_corr+x-1)-in_pic2_avg).*(compare_pic1-compare_avg1)));
        y_pos=y_pos+1;
    end
    x_pos=x_pos+1;
end

[max_junk,max_ind]=max(gamma_1(:));
[I,J]=ind2sub([2*orig_steps+1,2*orig_steps+1],max_ind);
move_1_x=(J-orig_steps-1);
move_1_y=(I-orig_steps-1);

%  Correlate for subpixel accuracy-----------------------------------------------------------

compare_pic1=in_pic1(1+(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps):end-(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps),...
    1+(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps):end-(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps));
[row_corr,col_corr]=size(compare_pic1);
x_pos=1;
        compare_avg1=sum(compare_pic1(:))/(row_corr*col_corr);
        
for x=(1+(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps))+4*move_4_x+2*move_2_x+move_1_x-1:(1+(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps))+4*move_4_x+2*move_2_x+move_1_x+1
    y_pos=1;
    for y=(1+(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps))+4*move_4_y+2*move_2_y+move_1_y-1:(1+(4*ds4_steps+2*ds2_steps+orig_steps+sub_steps))+4*move_4_y+2*move_2_y+move_1_y+1
        
        in_pic2_avg=sum(sum(in_pic2(y:row_corr+y-1,x:col_corr+x-1)))/(row_corr*col_corr);
        gamma_sub(y_pos,x_pos)=sum(sum((in_pic2(y:row_corr+y-1,x:col_corr+x-1)-in_pic2_avg).*(compare_pic1-compare_avg1)));
        y_pos=y_pos+1;
    end
    x_pos=x_pos+1;
end

[move_sub_x,move_sub_y]=quadratic_peak_interpolation_2d(gamma_sub);

move_total_x=4*move_4_x+2*move_2_x+move_1_x+move_sub_x;
move_total_y=4*move_4_y+2*move_2_y+move_1_y+move_sub_y;