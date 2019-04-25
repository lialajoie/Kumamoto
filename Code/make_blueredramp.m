% Lia Lajoie
% CSM
% 23 Aug 2017

% blue-red color ramp 
% inputs are data value range of interest (1x2 vector) and the color ramp
% mode, either "zero" or "center". "zero" will distribute the colors such
% that a value of zero is white, with positive values in red and negative
% values in blue. "center" will distribute the colors such that white is at
% the center of the input data range
function [blueredramp] = make_blueredramp(c_range,varargin)
if nargin > 1
    mode = varargin;
else
    mode = 'center'; % if user doesn't specify a mode, default is "center"
end

num_bar = 256; % 256 color ramp
if strcmp(mode,'zero') == 1 % if user wants zero to separate red and blue ramps
    numcolors_end1 = ceil((range(c_range) - abs(c_range(2)))/range(c_range)*num_bar);
    numcolors_end2 = floor((range(c_range) - abs(c_range(1)))/range(c_range)*num_bar);
elseif strcmp(mode,'center') == 1 % if user wants red and blue to be equally represented
    numcolors_end1 = ceil(num_bar/2);
    numcolors_end2 = floor(num_bar/2);
end
color1_ones = ones(numcolors_end1,1);
color1_grad = linspace(0,1,numcolors_end1)';
color2_ones = ones(numcolors_end2,1);
color2_grad = linspace(0,1,numcolors_end2)';
blue_color = [color1_grad color1_grad color1_ones]; % blue color ramp
red_color = [color2_ones color2_grad color2_grad];  % create your own red color map
red_flip = flipud(red_color);
blueredramp = vertcat(blue_color,red_flip);
end