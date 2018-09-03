function [] = animate3dPlot(filename)
% make animation of current 3d plot
% Matt Churgin, August 2018

if ~strcmp(filename(end-2:end),'gif')
    filename=[filename '.gif'];
end
elevationAngle=80;
startAz=270; % starting azimuth angle
revolutions=1; % number of revolutions to complete

for az=startAz:2:(startAz+revolutions*360)
    view([az,elevationAngle])
    pause(0.05)
    set(gca,'XTick',[50:50:250])
    set(gca,'YTick',[50:50:250])
    set(gca,'ZTick',[2:2:10])
    drawnow
    
    % Capture the plot as an image
    frame = getframe;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if az == startAz
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end