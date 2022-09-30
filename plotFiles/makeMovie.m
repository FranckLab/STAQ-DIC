
try close(v); catch; end

close all; 
%%%%%% Read saved images %%%%%%%
files = dir('*_vel_coneplot.jpg'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

%%%%%% Write frames to videos %%%%%%
v = VideoWriter('video_vel_coneplot.mp4','MPEG-4');
v.FrameRate = 5;
open(v);  % set(gcf, 'Position', [100 100 500 500]);
for tempk = [ 1 : 1 : length(im) ]
    
    myfig = figure;
    % imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow( imread( im{tempk}) );  title(['Frame #',num2str(tempk+1)]); 
    % text(830,100,['Frame #',num2str(tempk+1)]);
    % set(gcf, 'Position', [100 100 500 500]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    clf(myfig); close all; 
    % waitbar(tempk/length(files));
    
end

 
close(v);



%% y

try close(v); catch; end

close all; 
%%%%%% Read saved images %%%%%%%
files = dir('*_Disp_y.jpg'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

%%%%%% Write frames to videos %%%%%%
v = VideoWriter('video_disp_y.mp4','MPEG-4');
v.FrameRate = 5;
open(v);  % set(gcf, 'Position', [100 100 500 500]);
for tempk = [ 1 : 1 : length(im) ]
    
    myfig = figure;
    % imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow( imread( im{tempk}) );  title(['Frame #',num2str(tempk+1)]); 
    % text(830,100,['Frame #',num2str(tempk+1)]);
    % set(gcf, 'Position', [100 100 500 500]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    clf(myfig); close all; 
    % waitbar(tempk/length(files));
    
end

 
close(v);


%%  try close(v); catch; end

close all; 
%%%%%% Read saved images %%%%%%%
files = dir('*_strain_exx.jpg'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

%%%%%% Write frames to videos %%%%%%
v = VideoWriter('video_strain_exx.mp4','MPEG-4');
v.FrameRate = 5;
open(v);  % set(gcf, 'Position', [100 100 500 500]);
for tempk = [ 1 : 1 : length(im) ]
    
    myfig = figure;
    % imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow( imread( im{tempk}) );  title(['Frame #',num2str(tempk+1)]); 
    % text(830,100,['Frame #',num2str(tempk+1)]);
    % set(gcf, 'Position', [100 100 500 500]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    clf(myfig); close all; 
    % waitbar(tempk/length(files));
    
end

 
close(v);

%%

try close(v); catch; end

close all; 
%%%%%% Read saved images %%%%%%%
files = dir('*_strain_exy.jpg'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

%%%%%% Write frames to videos %%%%%%
v = VideoWriter('video_strain_exy.mp4','MPEG-4');
v.FrameRate = 5;
open(v);  % set(gcf, 'Position', [100 100 500 500]);
for tempk = [ 1 : 1 : length(im) ]
    
    myfig = figure;
    % imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow( imread( im{tempk}) );  title(['Frame #',num2str(tempk+1)]); 
    % text(830,100,['Frame #',num2str(tempk+1)]);
    % set(gcf, 'Position', [100 100 500 500]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    clf(myfig); close all; 
    % waitbar(tempk/length(files));
    
end

 
close(v);

%%

try close(v); catch; end

close all; 
%%%%%% Read saved images %%%%%%%
files = dir('*_strain_eyy.jpg'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

%%%%%% Write frames to videos %%%%%%
v = VideoWriter('video_strain_eyyf.mp4','MPEG-4');
v.FrameRate = 5;
open(v);  % set(gcf, 'Position', [100 100 500 500]);
for tempk = [ 1 : 1 : length(im) ]
    
    myfig = figure;
    % imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow( imread( im{tempk}) );  title(['Frame #',num2str(tempk+1)]); 
    % text(830,100,['Frame #',num2str(tempk+1)]);
    % set(gcf, 'Position', [100 100 500 500]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    clf(myfig); close all; 
    % waitbar(tempk/length(files));
    
end

 
close(v);
%%
% figure, plot(alpha,F11,'.'); 
