col = size(xyzs_idG,2);
M(length(Movie)) = struct('cdata',[],'colormap',[]);
clr = jet(numel(unique(xyzs_idG(:,col))));
figure
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
for i = 1:size(clusterCell,1)
    xyzs_temp = xyzs_idG(xyzs_idG(:,col-1) == i, :);
    xyzs_temp = sortrows(xyzs_temp, col-2);
    hold on
    axis([0,size(Movie{1},2),0,size(Movie{1},1)])
    im = imagesc(im2bw(MovieO{i},imadj(i,1)));
    colormap(gray);
    set(im, 'AlphaData', 0.5);
    title(['Frame ' num2str(i) ' number of clusters ' num2str(length(clusterCell{i}))]);
    for  jj = 1:length(clusterCell{i})
        disp(['plotting frame ' num2str(i) ' cluster ' num2str(jj)])
        cmap = zeros(size(clusterCell{i}{jj},1),1);
        cmap(:,1) = clr(xyzs_temp(xyzs_temp(:,col-2) == jj,col),1);
        cmap(:,2) = clr(xyzs_temp(xyzs_temp(:,col-2) == jj,col),2);
        cmap(:,3) = clr(xyzs_temp(xyzs_temp(:,col-2) == jj,col),3);
        scatter(clusterCell{i}{jj}(:,1), clusterCell{i}{jj}(:,2), 10, cmap);
        text(xyzs_temp(jj,1) + 20, xyzs_temp(jj, 2), num2str(xyzs_temp(xyzs_temp(:,col-2) == jj,col)), 'color', 'r');
    end
    M(i) = getframe;
    clf;
end

%{
xx = [];
x = xyzs_idG(:,end);
for i = 1:max(x)
    xx(i) = sum(x == i);
end
hist(xx)
%}