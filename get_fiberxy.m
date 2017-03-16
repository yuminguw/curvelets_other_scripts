
function  fiber_xy = get_fiberxy(imagename,pathname,fz,midpointEST,saveon,markeron)
% get_fiberxy get the fiber coordinates
% input: 
%    imagename: name of the image
%    pathname: path of the image directory
%    fz: font size of the fiber label
%    midpointEST: method to estimate the coordinates of the middle point associated with a fiber
%        1: based on the coordinates of the end points
%        2: based on the fiber length
%    saveon: 1: save the results into a csv file
%              0: do not save the results
%    markeron: 1: show marker of the center point
%              0: do not show the marker of the center point

% output:
%   fiber_xy: include 7 columes
%       c1: fiber index, c2-c3: xy of fiber start point, c4-c5:xy of fiber end point, c6-c7:xy of fiber center 
% usage:
%    example1: no input
%        xy = get_fiberxy;
%    example2: use an image in current directory,do not save result but indicate the middle point position, 
%        xy = get_fiberxy('testimage1.tif','.',8,2,0,1);

% started on Jan, 2017, Yuming Liu, LOCI, UW-Madison
close all; clc;home;
if nargin == 0
    imagename = 'testimage1.tif';
    pathname = pwd;
    fz = 8;  % fontsize of the fiber label
    midpointEST = 2;  % use interpolated coordinates to estimate the center of a fiber
    saveon = 0;% do not save the coordinates
    markeron = 0;
end
[~,imagenameNE,~] = fileparts(imagename);
matfilename = ['ctFIREout_' imagenameNE '.mat'];
% path and name of the csv file of the xy coordinates
saveon = 0 ;  
xy_savename = fullfile(pwd,'ctFIREout',[imagenameNE '_fiberXY_start_end_center.csv']);

if midpointEST == 1
    disp('Fiber middle point estimation is based on the coordinates of the two end points.')
elseif midpointEST == 2
    disp('Fiber middle point estimation is based on the fiber length.')
end
    

img_data = imread(fullfile(pathname,imagename));
pixh = size(img_data,1);
pixw = size(img_data,2);
mat_data = load(fullfile(pathname,'ctFIREout',matfilename));
length_threshold = mat_data.cP.LL1;  % get length threshold from the .mat file

data = mat_data.data;
FN = find(data.M.L>length_threshold);  % get fibers who are longer than the length threshold
LFa = length(FN);
rng(1001) ;
clrr2 = rand(LFa,3); % set random color

figure('Pos', [300 200 768 768]);
imshow(fullfile(pathname, imagename));
title(sprintf('%s',imagename));
hold on
start_xy = nan(LFa,2);
end_xy = nan(LFa,2);
center_xy = nan(LFa,2);
for LL = 1:LFa
    VFa_LL = data.Fa(1,FN(LL)).v;
    
    if midpointEST == 1   %use fiber end points to estimate the center position of a fiber
        fsp = data.Fa(FN(LL)).v(1);
        fep = data.Fa(FN(LL)).v(end);
        sp = data.Xa(fep,:);
        ep = data.Xa(fsp,:);
        cen = round(mean([sp; ep]));
        center_xy(LL,1:2) = cen(1,1:2); 

    elseif midpointEST == 2 %use interpolated coordinates to estimate the center position of a fiber
        VFai_LL = data.Fai(1,FN(LL)).v;
        center_indx = round(length(VFai_LL)/2);
        center_xy(LL,1:2) = round(data.Xai(VFai_LL(center_indx),1:2));
        if center_xy(1) > pixw ||center_xy(2)>pixh || center_xy(1) < 1 || center_xy(2) < 1
            vertex_indices = data.Fa(LL).v;
            s2 = size(vertex_indices,2);
            cen(1) = round(data.Xa(vertex_indices(round(s2/2)),1));
            cen(2) = round(data.Xa(vertex_indices(round(s2/2)),2));
            center_xy(LL,1:2) = [cen(1) cen(2)];
            fprintf('Interpolated coordinates of fiber %d is out of boundary,orignial coordinates is used for length-based fiber middle point estimation. \n',LL)
        end
    end
    XFa.LL = data.Xa(VFa_LL,:);
    plot(XFa.LL(:,1),XFa.LL(:,2), '-','color',clrr2(LL,1:3),'linewidth',1);  % plot a fiber
    if XFa.LL(1,1)> XFa.LL(end,1)
        start_xy(LL,1:2) = round([XFa.LL(end,1), XFa.LL(end,2)]);
        end_xy(LL,1:2) = round([XFa.LL(1,1), XFa.LL(1,2)]);
    else
        start_xy(LL,1:2) = round([XFa.LL(1,1), XFa.LL(1,2)]);
        end_xy(LL,1:2) = round([XFa.LL(end,1), XFa.LL(end,2)]);
    end
    if markeron == 1
        plot(center_xy(LL,1),center_xy(LL,2),'*', 'MarkerSize',7,'color',clrr2(LL,1:3))
    end
    text(center_xy(LL,1),center_xy(LL,2),sprintf('%d',LL), 'fontsize',fz,'color',clrr2(LL,1:3))
    if mod(LL-1,20) == 0
        disp(imagename)
        fprintf('%5s   %6s   %6s   %6s   %6s   %6s   %6s \n', 'Fiber#','X_end1',...
            'Y_end1', 'X_end2','Y_end2', 'X_mid','Y_mid')
    end
    fprintf('%5d   %6d   %6d  %6d   %6d   %6d   %6d \n',...
       LL,start_xy(LL,1),start_xy(LL,2),end_xy(LL,1),end_xy(LL,2),...
        center_xy(LL,1),center_xy(LL,2));
end
hold off
fib_indx = [1:LFa]';
fiber_xy = horzcat(fib_indx,start_xy,end_xy,center_xy);  % c1: fiber index, c2-c3: xy of start, c4-c5:xy of end, c6-c7:xy of center
if saveon == 1
    csvwrite(xy_savename,fiber_xy)
    fprintf('fiber xy coordinates are saved into %s \n',xy_savename)
end

end
%end of get_fiberxy

