function [M Movie] = ellipse_plot_henderson(xyzs_id,xyzs_id_columns,filename,image_mat_in,nuc,dead_inds)
% xyzs_id is the output from tracking code
% [ xyzs_idN, xyzs_id_columns, filename, framerate] = run_tracking_contour2();
% filename is the original image you want to compare it to
frameindx = xyzs_id_columns-1;
cellindx = xyzs_id_columns;
nframes = max(xyzs_id(:,frameindx));
ncells = max(xyzs_id(:,cellindx));
Movie = [];
filename = [filename '.tif'];
if nargin < 5
    nuc = 0;
end
if nargin < 6
   dead_inds = [(1:ncells).' nframes*ones(ncells,1)]; 
end
if nargin > 2
    tiffInfo = imfinfo(filename);
    no_frame = numel(tiffInfo);
    Movie = cell(no_frame,1);
    for iFrame = 1:no_frame
        Movie{iFrame} = uint16(imread(filename,'Index',iFrame,'Info',tiffInfo));
    end
    M = struct('cdata',[],'colormap',[]);
end
%disp([num2str(nargin) '  ' num2str(nuc)])
clr = jet(ncells);
figure
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
for i = 1:nframes
    hold on
    %if nargin > 2
    axis([0,size(Movie{1},2),0,size(Movie{1},1)]);
    hold on
    im = imagesc(imadjust(Movie{i},image_mat_in(i,:)));
    colormap(gray);
    set(im, 'AlphaData', 0.5);
    %else
    %         hold on
    %         axis([0,1000,0,1000]);
    %end
    X = xyzs_id(xyzs_id(:,frameindx) == i,:);
    if nargin > 3 && nuc > 0
    X = X(X(:,cellindx) == nuc,:);
    end
    
    %if size(X,1)/nframes > 0.8
    for j = 1:size(X,1)
        hold on
        ellipse(X(j,3),X(j,4),X(j,5),X(j,1),X(j,2),clr(X(j,cellindx),:));
        hold on
        text(X(j,1)+5,X(j,2),num2str(X(j,cellindx)),'color','r')
    end
    title(['frame  ' num2str(i)])
    %end
    M(i) = getframe;
    if i ~= nframes
        clf;
    end
end
% if nuc == 0
%     savename = filename(1:end-4);
%     movie2avi(M,[savename '_nuc.avi'],'FPS',3);
% else
%     savename = [filename(1:end-4) '_' num2str(nuc)];
%     movie2avi(M,savename,'FPS',3);
% end
movie2avi(M,'ellipse_video.avi','FPS',3)