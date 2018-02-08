function [golgi_id, nuclei_id] = golgiCombo(xyzs_idG, xyzs_idN,clusterCell,plot_option,filenameN,filenameG, maxD,image_mat_in,imadj)
%{
how to combine lists of cell ids...
we have two xyzs_id matricies
just match cell ids with centers of mass
re-check every so often

look at first instance of new cell id from nuclei data and compare it to
golgi locations at that timestep
%}
%%
try
    xyzs_idN = xyzs_id;
end
tempG = [xyzs_idG(:,1) xyzs_idG(:,2) xyzs_idG(:,11) xyzs_idG(:,12)];
tempN = [xyzs_idN(:,1) xyzs_idN(:,2) xyzs_idN(:,12) xyzs_idN(:,13)];
%xG = xyzs_idG;
nCells = max(tempN(:,4));
nGol = numel(unique(tempG(:,4)));
nFrames = max(tempN(:,3));
%nFrames = 20;
for i = 1:max(tempN(:,4))
    f = tempN(:,4) == i;
    if sum(f) < 1%5
        tempN(tempN(:,4)==i,:) = [];
    end
end

for i = 1:nGol
    tg = tempG(tempG(:,4)==i,3);
    if length(tg) < 1%nFrames/10
        tempG(tempG(:,4)==i,:) = [];
    end
    %     if rand(1) < 0.3
    %         tempG(tempG(:,4)==i,:) =[];
    %         xG(xG(:,12)==i,:) = [];
    %     end
end

%this list will have the id from the nuclei code in the first column and
%the golgi it corresponds to in the second, can then change the outputs to
%match based off of this
%go through all the cells initial positions and see which Golgi they most
%likely correlate with
% tempG = [xyzs_idG(:,1) xyzs_idG(:,2) xyzs_idG(:,11) xyzs_idG(:,12)];
% tempN = [xyzs_idN(:,1) xyzs_idN(:,2) xyzs_idN(:,11) xyzs_idN(:,12)];
distcell = cell(nCells,1);
distcell1 = cell(nCells,1);
for i = 1:nCells
    %disp(['cell ' num2str(i)])
    distcell{i} = zeros(nGol,nFrames);
    dd = zeros(nGol,2);
    posN = tempN(tempN(:,4) == i,:);
    for k = 1:nGol
        posG = tempG(tempG(:,4)== k,:);
        I = intersect(posN(:,3),posG(:,3));
        indn = ismember(posN(:,3),I);
        indg = ismember(posG(:,3),I);
        Dist = sqrt((posN(indn,1)-posG(indg,1)).^2+(posN(indn,2)-posG(indg,2)).^2);
        %         for j = 1:nFrames
        %
        %             try
        %                 distcell{i}(k,j) = dist(posG(posG(:,4)==k,1),posG(posG(:,4)==k,2),posN(posN(:,4)==i,1),posN(posN(:,4)==i,2));
        %             end
        %         end
        d1 = numel(Dist)*sum(Dist)/(numel(Dist(Dist~=0)).^2);
        %d1 = numel(distcell{i}(k,:))*sum(distcell{i}(k,:))/(numel(distcell{i}(k,distcell{i}(k,:) ~= 0)).^2);
        d2 = sum(Dist)/numel(Dist(Dist~=0));
        %d2 = sum(distcell{i}(k,:))/(numel(distcell{i}(k,distcell{i}(k,:) ~= 0)));
        %d1 = d2;
        %distcell{i}(k,1) = d1;
        %distcell{i}(k,2) = d2;
        dd(k,:) = [d1 d2];
    end
    distcell1{i} = distcell{i};
    %     cpy1 = distcell{i}(:,1);
    %     cpy2 = distcell{i}(:,2);
    distcell{i} = zeros(nGol,2);
    distcell{i}(:,1) = 1:nGol;
    %     distcell{i}(:,2) = cpy1;
    %     distcell{i}(:,3) = cpy2;
    distcell{i}(:,2:3) = dd;
    distcell{i} = sortrows(distcell{i},3);
end
dmat = zeros(nCells,1);
flist = [];
for i = 1:nCells
    check = 0;
    m = 1;
    while check == 0
        if distcell{i}(m,3) < maxD
            if abs(distcell{i}(m,3) - distcell{i}(m+1,3)) < maxD
                check1 = 0;
                mmm = 1;
                tempdistcell = sortrows(distcell{i}(distcell{i}(:,3)<maxD,:),2);
                while check1 == 0
                    if ismember(tempdistcell(mmm,1),flist) ~= 1
                        dmat(i) = tempdistcell(mmm,1);
                        check = 1;
                        check1 = 1;
                        flist = [flist dmat(i)];
                    end
                    if mmm == length(tempdistcell(:,1))
                        check1 = 1;
                    end
                    mmm = mmm + 1;
                end
            elseif abs(distcell{i}(m,3) - distcell{i}(m+1,3)) > maxD
                if ismember(distcell{i}(m,1),flist) ~= 1
                    dmat(i) = distcell{i}(m,1);
                    check = 1;
                    flist = [flist dmat(i)];
                end
            end
        end
        if m == length(distcell{i}(:,1))
            check = 1;
        end
        m = m + 1;
    end
end
indC = cell(nGol,1);
dG = cell(nGol,1);
for i = 1:nGol
    indC{i} = tempG(tempG(:,4) == i,3);
    for j = 1:nCells
        dG{i}(j,:) = [j distcell{j}(distcell{j}(:,1) == i,2) distcell{j}(distcell{j}(:,1) == i,3)];
    end
end

idList = zeros(nGol,3);
idList(:,3) = inf;
idList(:,1) = 1:nGol;
for i = 1:nCells
    j = dmat(i);
    if j ~= 0
        idList(j,2) = i;
    end
end
nn = nCells;
for i = 1:nGol
    dGtemp = sortrows(dG{i},3);
    if idList(i,2) == 0
        check = 0;
        for j = 1:size(dGtemp,1)
            if dGtemp(j,3) < maxD
                v1 = indC{i};
                v2 = indC{idList(idList(:,2) == dGtemp(j,1),1)};
                v1c = ismember(v1,v2);
                v2c = ismember(v2,v1);
                vf = [v1c.' v2c.'];
                if sum(vf) == 0 %&& dGtemp(j,3) < idList(idList(:,2) == dGtemp(j,1),3)
                    idList(i,2) = dGtemp(j,1);
                    idList(i,3) = dGtemp(j,3);
                    check = 1;
                    break
                end
                if check == 1
                    break
                end
            end
            if check == 1
                break
            end
        end
        if check == 0
            idList(i,2) = nn + 1;
            nn = nn + 1;
        end
    end
end

golgi_id = xyzs_idG;
for i = 1:size(golgi_id,1)
    golgi_id(i,12) = idList(golgi_id(i,12),2);
end
nuclei_id = xyzs_idN;

%%
if plot_option ~= 0

    tiffInfo = imfinfo(filenameN);
    no_frame = numel(tiffInfo);
    MovieN = cell(no_frame,1);
    for iFrame = 1:no_frame
        MovieN{iFrame} = uint16(imread(filenameN,'Index',iFrame,'Info',tiffInfo));
    end
    M(no_frame) = struct('cdata',[],'colormap',[]);
    clr = jet(nn);
    
    tiffInfo = imfinfo(filenameG);
    no_frame = numel(tiffInfo);
    MovieG = cell(no_frame,1);
    for iFrame = 1:no_frame
        MovieG{iFrame} = uint16(imread(filenameG,'Index',iFrame,'Info',tiffInfo));
    end
    %M(nFrames) = struct('cdata',[],'colormap',[]);
    %clr = jet(nn);
    
end
%%
if plot_option == 1
    disp('Golgi shown as blue circles, nuclei as green. Identification number for that cell in red. Only cells with both golgi and nuclei shown.')
    f = figure;
    pause(3);
    gmap = [0,0,0;
        0,255,0;
        0,0,255;];
    frame_h = get(f,'JavaFrame');
    set(frame_h,'Maximized',1);
    pause(3);
    for i = 1:length(MovieN)
        Ntem = nuclei_id(nuclei_id(:,12)==i,:);
        Gtem = golgi_id(xyzs_idG(:,11)==i,:);
        Ntem = sortrows(Ntem,13);
        Gtem = sortrows(Gtem,12);
        hold on
        axis([0,size(MovieN{1},2),0,size(MovieN{1},1)])
        hold on
        colormap(gmap./255);
        imn = uint8(im2bw(imadjust(MovieN{i},image_mat_in(i,:)))*200);
        img = uint8(im2bw(imadjust(MovieG{i},imadj(i,:)))*100);
        im = imagesc(imn+img);
        %im = imagesc(imadjust(Movie{i},image_mat_in(i,:)));
        %set(im, 'AlphaData', 1);
        for j = 1:size(Ntem,1)
            if ismember(Ntem(j,13),Gtem(:,12)) == 1
                %if length(nuclei_id(nuclei_id(:,13)==Ntem(j,13),1)) > nFrames/10
                %scatter(Ntem(j,1),Ntem(j,2),100,'r');%clr(Ntem(j,13),:));
                text(Ntem(j,1) + 10, Ntem(j, 2), num2str(Ntem(j,13)), 'color', 'r');
                try
                    q = quiver(Ntem(j,1),Ntem(j,2),-(Ntem(j,1)-mean(Gtem(Gtem(:,12)==Ntem(j,13),1))), -(Ntem(j,2)-mean(Gtem(Gtem(:,12)==Ntem(j,13),2))),0.5,'r');%'r','LineWidth',2);
                    %set(q,'MaxHeadSize',1)
                end
            end
        end
        
        for j = 1:size(Gtem,1)
            if ismember(Gtem(j,12),Ntem(:,13)) == 1
                %if length(golgi_id(golgi_id(:,12)==Gtem(j,12),1)) > nFrames/10
                %scatter(Gtem(j,1),Gtem(j,2),100,'r');%clr(Gtem(j,12),:));
                %text(Gtem(j,1) + 10, Gtem(j, 2), num2str(Gtem(j,12)), 'color', 'b');
            end
        end
        %set(gca,'ydir','reverse')
        M(i) = getframe;
        if i ~= no_frame
            clf;
        end
    end
    try
        movie2avi(M, 'combo.avi', 'FPS', 3);
    end
end
col = 12;
if plot_option == 2
    for i = 1:nFrames
        xyzs_temp = xyzs_idG(xyzs_idG(:,col-1) == i, :);
        xyzs_temp = sortrows(xyzs_temp, col);
        hold on
        axis([0,size(Movie{1},1),0,size(Movie{1},2)])
        im = image(Movie{i});
        colormap(gray);
        set(im, 'AlphaData', 0.5);
        for k = 1:length(clusterCell{i})
            disp(['plotting frame ' num2str(i) ' cluster ' num2str(k)])
            ind = xyzs_temp(k,12);
            scatter(clusterCell{i}{k}(:,1), clusterCell{i}{k}(:,2), 10, clr(ind,:));
            text(xyzs_temp(k,1) + 50, xyzs_temp(k, 2), num2str(xyzs_temp(xyzs_temp(:,col-2) == k,col)), 'color', 'r');
        end
        Ntem = xyzs_idN(xyzs_idN(:,11)==i,:);
        for k = 1:size(Ntem,1)
            ellipse(Ntem(j,3),Ntem(j,4),Ntem(j,5),Ntem(j,1),Ntem(j,2),clr(Ntem(j,12),:));
        end
        M(i) = getframe;
        clf;
    end
end
%disp('done with combo')
%%
    function d = dist(x1,y1,x2,y2)
        d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
    end
end