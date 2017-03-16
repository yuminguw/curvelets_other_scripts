% get_fiberxy get the fiber coordinates
% Jan, 2017, Yuming Liu, LOCI, UW-Madison

function  fiber_xy = get_fiberxy(imagename,pathname,fz)
close all
fz = 8;  % fontsize of the fiber label
imagename = 'testimage1.tif';
[~,imagenameNE,~] = fileparts(imagename);
pathname = pwd;
matfilename = ['ctFIREout_' imagenameNE '.mat'];
% path and name of the csv file of the xy coordinates
saveon = 0 ;  % 1: to save the coordinates
xy_savename = fullfile(pwd,'ctFIREout',[imagenameNE '_fiberXY_start_end_center.csv']);


img_data = imread(fullfile(pathname,imagename));
pixh = size(img_data,1);
pixw = size(img_data,2);
mat_data = load(fullfile(pathname,'ctFIREout',matfilename));
length_threshold = mat_data.cP.LL1;  % get length threshold from the .mat file

data = mat_data.data;
FN = find(data.M.L>length_threshold);  % get fibers who are longer than the length threshold
FLout = data.M.L(FN);
LFa = length(FN);
rng(1001) ;
clrr2 = rand(LFa,3); % set random color

OL_fig = figure('Pos', [300 200 768 768]);
imshow(fullfile(pathname, imagename));
title(sprintf('%s',imagename));
hold on

for LL = 1:LFa
    VFa_LL = data.Fa(1,FN(LL)).v;
    %Use interpolated coordinates to estimate the center of fiber
    VFai_LL = data.Fai(1,FN(LL)).v;
    center_indx = round(length(VFai_LL)/2);
    center_xy(LL,1:2) = round(data.Xai(VFai_LL(center_indx),1:2));
    
    XFa.LL = data.Xa(VFa_LL,:);
    plot(XFa.LL(:,1),XFa.LL(:,2), '-','color',clrr2(LL,1:3),'linewidth',1);  % plot a fiber
    if XFa.LL(1,1)> XFa.LL(end,1)
        start_xy(LL,1:2) = round([XFa.LL(end,1), XFa.LL(end,2)]);
        end_xy(LL,1:2) = round([XFa.LL(1,1), XFa.LL(1,2)]);
    else
        start_xy(LL,1:2) = round([XFa.LL(1,1), XFa.LL(1,2)]);
        end_xy(LL,1:2) = round([XFa.LL(end,1), XFa.LL(end,2)]);
    end
    text(center_xy(LL,1),center_xy(LL,2),sprintf('%d',LL), 'fontsize',fz,'color',clrr2(LL,1:3))
    %     disp(imagename)
    %     disp(sprintf('%d, FN = %d, start = [%d  %d] end = [%d  %d],center = [%d  %d], ',...
    %         LL,FN(LL),start_xy(LL,1),start_xy(LL,2),end_xy(LL,1),end_xy(LL,2), center_xy(LL,1),center_xy(LL,2)));
    %     disp('Press any key to continue ...')
    %     pause
   
end
hold off
fib_indx = [1:LFa]';
fiber_xy = horzcat(fib_indx,start_xy,end_xy,center_xy);  % c1: fiber index, c2-c3: xy of start, c4-c5:xy of end, c6-c7:xy of center
if saveon == 1
    csvwrite(xy_savename,fiber_xy)
    disp(sprintf('fiber xy coordinates are saved into %s',xy_savename))
end

return

%end of get_fiberxy

