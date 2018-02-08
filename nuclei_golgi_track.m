%{ 
This Matlab script contains the necessary steps to idenfity and track
nuclei and Golgi bodies and then combine their identification lists to link
nuclei and Golgi bodies that belong to the same cell. 

Steps: 
1. Image import from tiff stacks. A user interface will appear prompting
you to select two tiff video files of tracked nuclei and Golgi bodies (in
seperate channels). If your videos are not in this format you can convert
it to a tiff stack in ImageJ or modify this Matlab code directly to suit
your needs. 
2. Image pre-processing. Images can often contain backround noise or have
insufficient brightness and contrast for objects to be picked up by the
tracking code. Please modify your inputs for this section until you see
images where you can identify objects by eye. 
3. Nuclei tracking. We use the ACTIVE nuclei tracking software package
developed by the Henderson lab. 
4. Golgi tracking. This code identifies and tracks Golgi bodies and
constructs an identification list similar to the ACTIVE package. 
5. Combination of ID lists. We look at cell positions of both nuclei and
Golgi and pair them appropriately to create one master ID list. 

Following this analysis: 
The primary tool we used to look at mean and standard deviation of
orientations was the AngleSpread2 function with the following syntax: 

[FinalStd, MaxAngle]  = AngleSpread2(theta);

where theta is a vector list of orientations. This function can be found in
the Nuclear_Alignment_Video_Analysis folder. As a note to the user this
code is more to provide the initial imaging and tracking software necessary
to identify and track irregularly shaped objects such as the Golgi body.
Any futher analysis is largely based on the user's needs. 

Developed by Giuseppe Passucci from 2013 to 2017 as part of the Manning
Group at Syracuse University in collaboration with the Henderson (SU) and
Turner (SUNY Upstate) groups. 
%}

%{
Needed external files: 
- ACTIVE tracking package
- ellipse_plot_henderson.m
- golgiFind.m
- quick_golgi_movie.m
- golgiCombo.m
%}
[filename, pathname] = uigetfile({'*.tif'}, 'Select a nuclei file');
full_nameN = [pathname, filename];
[filename, pathname] = uigetfile({'*.tif'}, 'Select a golgi file');
full_nameG = [pathname, filename];
choice = questdlg('Would you like to see plots for visual verification of tracking?', ...
	'Plot option', ...
	'Yes','No','No');
switch choice
    case 'Yes'
        disp('Plots and movies will be shown.')
        plot_toggle = 1;
    case 'No'
        disp('Plots and movies will not be shown.')
        plot_toggle = 0;
end

filename = full_nameN;
tiffInfo = imfinfo(filename);
no_frame = numel(tiffInfo);
MovieN = cell(no_frame,1);
tmax = length(MovieN);
for iFrame = 1:no_frame
    %MovieN{iFrame} = imcomplement(uint8(imread(filename,'Index',iFrame,'Info',tiffInfo)));
    MovieN{iFrame} = uint16(imread(filename,'Index',iFrame,'Info',tiffInfo));
end
% m2 = MovieN(2:2:end,1);
% m1 = MovieN(1:2:end,1);
% MovieN = cellfun(@(x,y) x+y,m1,m2,'UniformOutput',0);
disp('Determing image adjustment values for nuclei.')
%%
prompt = {'First cutoff, sigma ','Second cutoff, sigma'};
dlg_title = 'Nuclei intensity histogram cutoff';
num_lines = 1;
def = {'3','4'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
sig1 = str2num(answer{1});
sig2 = str2num(answer{2});

if plot_toggle == 1
M = struct('cdata',[],'colormap',[]);
fig = figure;
end
image_mat_in = zeros(no_frame,2);
for i = 1:no_frame
    x = MovieN{i};
    x = double(x(:));
    [N,X] = hist(x,10000);
    X = X(N~=0);
    N = N(N~=0);
    N = N/trapz(X,N);
    f = find(N<max(N)*0.05);
    fm = find(N == max(N));
    f = f(f>fm(end));
    f = f(1);
    [mu, sigma] = normfit(x(x<X(f)));
    y = normpdf(X,mu,sigma);
%     scatter(X,N)
%     hold on
%     plot(X,y,'r')
%     xlim([0 2*10^4])
%     set(gca,'FontSize',24)
%     ylabel('PDF','FontSize',30)
%     xlabel('Pixel intensity','FontSize',30)
    c = [0 mu sigma]/2^16;
    c1 = c(2)+c(3)*sig1;
    c2 = c(2)+c(3)*sig2;
    if sig1==0 && sig2==0
        c1 = 0;
        c2 = 1;
    end
    image_mat_in(i,:) = [c1 c2];
    if plot_toggle == 1
        subplot(1,2,1);
        s1 = imagesc(MovieN{i});
        subplot(1,2,2);
        s2 = imagesc(imadjust(MovieN{i},[c1 c2]));
        M(i) = getframe(fig);
        if i ~= no_frame
            clf;
        end
    end
end
if plot_toggle == 1
%movie2avi(M,'nuclei_adjust.avi')
end
%%
disp('Determing image adjustment values for golgi.')
prompt = {'First cutoff, sigma ','Second cutoff, sigma'};
dlg_title = 'Golgi intensity histogram cutoff';
num_lines = 1;
def = {'4','5'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
sig1 = str2num(answer{1});
sig2 = str2num(answer{2});
filename = full_nameG;
tiffInfo = imfinfo(filename);
no_frame = numel(tiffInfo);
MovieG = cell(no_frame,1);
tmax = length(MovieG);
for iFrame = 1:no_frame
    MovieG{iFrame} = uint16(imread(filename,'Index',iFrame,'Info',tiffInfo));
    %MovieG{iFrame} = imcomplement(uint8(imread(filename,'Index',iFrame,'Info',tiffInfo)));
end
if plot_toggle == 1
    close all
    M = struct('cdata',[],'colormap',[]);
    fig = figure;
end

imadj = zeros(no_frame,2);
for i = 1:no_frame
    try
        c1 = 0;
        c2 = 0;
        x = MovieG{i};
        x = double(x(:));
        [mu, sigma] = normfit(x);
        c = [0 mu sigma]/2^16;
        c1 = c(2)+c(3)*sig1;
        c2 = c(2)+c(3)*sig2;
        if sig1==0 && sig2==0
            c1 = 0.1;
            c2 = 1;
        end
        imadj(i,:) = [c1 c2];
        if plot_toggle == 1
        subplot(1,2,1);
        s1 = imagesc(MovieG{i});
        subplot(1,2,2);
        s2 = imagesc(imadjust(MovieG{i},[c1 c2]));
        M(i) = getframe(fig);
        if i ~= no_frame
            clf;
        end
        end
    catch
        continue
    end
end
%%
if plot_toggle == 1
%movie2avi(M,'golgi_adjust.avi')
end
%close all
for i = 1:length(imadj)
   if imadj(i,1) == 0 
       f = find(imadj(1:i,1));
       imadj(i,:) = imadj(f(end),:);
   end
end

disp('Nuclei tracking:')
figure
imagesc(imadjust(MovieN{1},image_mat_in(1,:)))
[ xyzs_idN, xyzs_id_columns, filename, framerate, new_dir] = run_tracking_contour2(image_mat_in, full_nameN);

if plot_toggle == 1
    ellipse_plot_henderson(xyzs_idN,xyzs_id_columns,full_nameN(1:end-4),image_mat_in);
    close all
end

disp('Golgi tracking:')
figure
imagesc(imadjust(MovieG{1},imadj(1,:)))
[xyzs_idG, Movie, MovieO, clusterCell, filename, avg_sep, min_area, maxD] = golgiFind(full_nameG,imadj);


disp('Combining identification lists')
[golgi_id, nuclei_id] = golgiCombo(xyzs_idG, xyzs_idN,clusterCell,0,full_nameN,full_nameG, maxD,image_mat_in,imadj);
disp(full_nameG)
save('data.mat','image_mat_in','xyzs_idN','xyzs_id_columns','xyzs_idG','clusterCell','golgi_id','nuclei_id','full_nameG','full_nameN')
