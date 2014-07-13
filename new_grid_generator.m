function [grid_x,grid_y]=grid_generator(FileNameBase,PathNameBase);

% Code to generate the DIC analysis grid
% Completely rewritten by Chris
% Programmed first by Dan and Rob 
% 
% Last revision: 12/27/06

% customised by Aditya Vipradas

% The grid_generator function will help you create grids of markers. The
% dialog has different options allowing you to create a marker grid which is rectangular,
% circular, a line or two rectangels of a shape or contains only of two
% markers. After choosing one of the shapes you will be asked for the base
% image which is typically your first image. After opening that image you
% will be asked to click at the sites of interest and the markers will be
% plotted on top of your image. You can choose if you want to keep these
% markers or if you want to try again.
% It has to be noted that you can
% always generate your own marker positions. Therefore the marker position
% in pixel has to be saved as a text based format where the x-position is
% saved as grid_x.dat and the y-position saved as grid_y.dat.
%


% Prompt user for base image
if exist('FileNameBase')==0
[FileNameBase,PathNameBase] = uigetfile( ...
    {'*.bmp;*.tif;*.jpg;*.TIF;*.BMP;*.JPG','Image files (*.bmp,*.tif,*.jpg)';'*.*',  'All Files (*.*)'}, ...
    'Open base image for grid creation');

end
cd(PathNameBase)
im_grid = imread(FileNameBase);

[grid_x,grid_y,FileNameBase,PathNameBase] = gridtypeselection(FileNameBase, PathNameBase, im_grid);

close all

%-------------------------------
%
% Decide which type of grid you want to create

function [grid_x,grid_y,FileNameBase,PathNameBase] = gridtypeselection(FileNameBase, PathNameBase, im_grid);

hold off
imshow(im_grid,'truesize');

%gridselection = menu(sprintf('Which type of grid do you want to use'),...
%    'Rectangular','Circular','Two Markers','Line','Two Rectangles of Markers','Cancel');


%menu is reduced down to what is required
gridselection = menu(sprintf('Which type of grid do you want to use'),...
    'Rectangular','Two Markers','Set extreme coordinates','Cancel');

if gridselection==1
    [grid_x,grid_y,FileNameBase,PathNameBase, numXelem, numYelem] = rect_grid(FileNameBase, PathNameBase, im_grid);
    return
end

if gridselection==2
    [grid_x,grid_y,FileNameBase,PathNameBase] = twop_grid(FileNameBase, PathNameBase, im_grid);
    return
end

if gridselection==3
    [grid_x, grid_y, FileNameBase,PathNameBase] = set_coord(FileNameBase, PathNameBase, im_grid);
    return
end

if gridselection==4
    return;
end

if gridselection==5
    [grid_x,grid_y,FileNameBase,PathNameBase] = circ_grid(FileNameBase, PathNameBase, im_grid);
    return
end

if gridselection==6
    [grid_x,grid_y,FileNameBase,PathNameBase] = line_grid(FileNameBase, PathNameBase, im_grid);
    return
end



%-------------------------------
%
%extreme coordinates requires selection of two points at the top right
%corner and bottom left corner so that the domain coordinates can be stored in xcoord.dat and 
%ycoord.dat for plotting purposes
function [grid_x, grid_y, FileNameBase,PathNameBase] = set_coord(FileNameBase, PathNameBase, im_grid);

title(sprintf('Pick two extreme points on the sample for domain limit setting') )

[x(1,1),y(1,1)]=ginput(1);
hold on
plot(x(1,1),y(1,1),'+g')

[x(2,1),y(2,1)]=ginput(1);
plot(x(2,1),y(2,1),'+g')
xcoord = [x(1,1), x(2,1)]
ycoord = [y(1,1), y(2,1)]

grid_x = 0;
grid_y = 0;
% Accept the chosen markers, try again or give up 

confirmcircselection = menu(sprintf('Do you want to use these markers?'),...
    'Yes','No, try again','Go back to grid-type selection');

if confirmcircselection==2
    close all
    hold off
    imshow(im_grid,'truesize');
    set_coord(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==3
    close all
    gridtypeselection(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==1
    close all
    save xcoord.dat xcoord -ascii -tabs
    save ycoord.dat ycoord -ascii -tabs
    gridtypeselection(FileNameBase, PathNameBase, im_grid);
end

%-------------------------------
%
% Define line and create markers

function [grid_x,grid_y,FileNameBase,PathNameBase] = line_grid(FileNameBase, PathNameBase, im_grid);

title(sprintf('Pick two points on the sample.') )

[x(1,1),y(1,1)]=ginput(1);
hold on
plot(x(1,1),y(1,1),'+g')

[x(2,1),y(2,1)]=ginput(1);
plot(x(2,1),y(2,1),'+g')


linelength=sqrt((x(2,1)-x(1,1))*(x(2,1)-x(1,1))+(y(2,1)-y(1,1))*(y(2,1)-y(1,1)));
lineslope=(y(2,1)-y(1,1))/(x(2,1)-x(1,1));
intersecty=y(1,1)-lineslope*x(1,1);
ycalc=zeros(2,1);
ycalc=lineslope*x+intersecty;
plot(x(:,1),ycalc(:,1),'-b')


prompt = {'Enter the number of intersections between markers on the line:'};
dlg_title = 'Input for grid creation';
num_lines= 1;
def     = {'30'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
linediv = str2num(cell2mat(answer(1,1)));
linestep=((max(x)-min(x))/linediv);
grid_x(1:linediv+1)=min(x)+linestep*(1:linediv+1)-linestep;
grid_y=lineslope*grid_x+intersecty;

plot(grid_x,grid_y,'ob')
title(['Selected grid has ',num2str(linediv), ' rasterpoints'])    % plot a title onto the image

% Accept the chosen markers, try again or give up 

confirmcircselection = menu(sprintf('Do you want to use these markers?'),...
    'Yes','No, try again','Go back to grid-type selection');

if confirmcircselection==2
    close all
    hold off
    imshow(im_grid,'truesize');
    line_grid(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==3
    close all
    gridtypeselection(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==1
    save grid_x.dat grid_x -ascii -tabs
    save grid_y.dat grid_y -ascii -tabs
end

%-------------------------------
%
%
%for two markers, select one point along the centerline of the beam and
%then the program asks the pixels distance to the left and right of the
%selected marker
%the two markers are then chosen accordingly
function [grid_x,grid_y,FileNameBase,PathNameBase] = twop_grid(FileNameBase, PathNameBase, im_grid)

title(sprintf('Pick two points on the sample.') );
[a,b] = ginput(1);

%select two grid points
prompt = 'Enter distance in pixels to the left and right of the selected point';
dlg_title = 'Enter distance in pixels to the left and right of the selected point';
num_lines = 1;
def = {'5'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
pix = str2num(cell2mat(answer(1,1)));
disp(pix);
x(1,1) = a - pix;
x(2,1) = a + pix;
y(1,1) = b;
y(2,1) = b;
%[x(1,1),y(1,1)]=[a-pix,b];
hold on
plot(x(1,1),y(1,1),'+g')
disp(x(1,1))
%[x(2,1),y(2,1)]=[a+pix,b];
plot(x(2,1),y(2,1),'+g')
disp(x(2,1))
disp('Distance between the two markers in x-direction is:')
disp(x(2,1) - x(1,1))
disp('Distance between the two markers in y-direction is:')
disp(y(2,1) - y(1,1))
% Accept the chosen markers, try again or give up 

confirmcircselection = menu(sprintf('Do you want to use these two markers?'),...
    'Yes','Try again','Go back to grid-type selection');

grid_x = x;
grid_y = y;

if confirmcircselection==2
    close all
    hold off
    imshow(im_grid,'truesize');
    twop_grid(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==3
    close all
    gridtypeselection(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==1
    save grid_x.dat grid_x -ascii -tabs
    save grid_y.dat grid_y -ascii -tabs
    close all
    hold off
end
%-------------------------------
%
% Select a circular area

function [grid_x,grid_y,FileNameBase,PathNameBase] = circ_grid(FileNameBase, PathNameBase, im_grid);

title(sprintf('Pick three points on the circle in clockwise order at the upper boundary of the sample.') )

[x(1,1),y(1,1)]=ginput(1);
hold on
plot(x(1,1),y(1,1),'+g')

[x(2,1),y(2,1)]=ginput(1);
plot(x(2,1),y(2,1),'+g')

[x(3,1),y(3,1)]=ginput(1);
plot(x(3,1),y(3,1),'+g')

xnew=x;
ynew=y;

% Calculate center between the 3 sorted points and the normal slope of the vectors
slope12=-1/((ynew(2,1)-ynew(1,1))/(xnew(2,1)-xnew(1,1)));
slope23=-1/((ynew(3,1)-ynew(2,1))/(xnew(3,1)-xnew(2,1)));
center12(1,1)=(xnew(2,1)-xnew(1,1))/2+xnew(1,1);
center12(1,2)=(ynew(2,1)-ynew(1,1))/2+ynew(1,1);
center23(1,1)=(xnew(3,1)-xnew(2,1))/2+xnew(2,1);
center23(1,2)=(ynew(3,1)-ynew(2,1))/2+ynew(2,1);
% plot(center12(1,1),center12(1,2),'+b')
% plot(center23(1,1),center23(1,2),'+b')

if slope12==slope23
    return
end

% Calculate the crossing point of the two vectors
achsenabschnitt1=center12(1,2)-center12(1,1)*slope12;
achsenabschnitt2=center23(1,2)-center23(1,1)*slope23;
xdata=min(x):max(x);
ydata1=achsenabschnitt1+slope12*xdata;
ydata2=achsenabschnitt2+slope23*xdata;
% plot(xdata,ydata1,'-b')
% plot(xdata,ydata2,'-b')
xcross=(achsenabschnitt2-achsenabschnitt1)/(slope12-slope23);
ycross=slope12*xcross+achsenabschnitt1;
plot(xcross,ycross,'or')

% Calculate radius and plot circle
R=sqrt((xcross-xnew(1,1))*(xcross-xnew(1,1))+(ycross-ynew(1,1))*(ycross-ynew(1,1)));
% ydata=ycross-sqrt(R*R-(xdata-xcross).*(xdata-xcross));
% plot(xdata,ydata,'-b')

% Calculate angle between vectors
xvector=[1;0];
x1vec(1,1)=xnew(1,1)-xcross;x1vec(2,1)=ynew(1,1)-ycross
x3vec(1,1)=xnew(3,1)-xcross;x3vec(2,1)=ynew(3,1)-ycross
alpha13=acos((dot(x1vec,x3vec))/(sqrt(x1vec'*x1vec)*sqrt(x3vec'*x3vec)))*180/pi;
alpha01=acos((dot(xvector,x1vec))/(sqrt(x1vec'*x1vec)*sqrt(xvector'*xvector)))*180/pi;
alpha03=acos((dot(xvector,x3vec))/(sqrt(xvector'*xvector)*sqrt(x3vec'*x3vec)))*180/pi;
totalangle=alpha13;
minangle=alpha01;
maxangle=alpha03;
angldiv=abs(round(totalangle))*10;
anglstep=(totalangle/angldiv);
anglall(1:angldiv+1)=maxangle+anglstep*(1:angldiv+1)-anglstep;
xcircle(1:angldiv+1)=xcross+R*cos(-anglall(1:angldiv+1)/180*pi);
ycircle(1:angldiv+1)=ycross+R*sin(-anglall(1:angldiv+1)/180*pi);
plot(xcircle,ycircle,'-b')
drawnow

title(['Segment of circle spreads over ',num2str(totalangle),'°'])


% Accept the chosen circle, try again or give up 

confirmcircselection = menu(sprintf('Do you want to use this circle as basis?'),...
    'Yes','No, try again','Go back to grid-type selection');

if confirmcircselection==2
    close all
    imshow(im_grid,'truesize');
    circ_grid(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==3
    close all
    gridtypeselection(FileNameBase, PathNameBase, im_grid);
end

if confirmcircselection==1
    
    prompt = {'Enter the number of intersections between markers on the circle:'};
    dlg_title = 'Input for grid creation';
    num_lines= 1;
    def     = {'30'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    angldiv = str2num(cell2mat(answer(1,1)));
    
    anglstep=(totalangle/angldiv);
    anglall(1:angldiv+1)=maxangle+anglstep*(1:angldiv+1)-anglstep;
    
    markerxpos(1:angldiv+1)=xcross+R*cos(-anglall(1:angldiv+1)/180*pi);
    markerypos(1:angldiv+1)=ycross+R*sin(-anglall(1:angldiv+1)/180*pi);
    
    plot(markerxpos,markerypos,'ob');
    
    % Pick the lower bound in the image
    title(sprintf('Pick three points lying on the circle in clockwise order. The first and last one define the width of the raster') )
    
    [x(4,1),y(4,1)]=ginput(1);
    hold on
    plot(x(1,1),y(1,1),'+r')
    
    lowboundx=x(4,1);
    lowboundy=y(4,1);
    
    R2=sqrt((xcross-lowboundx(1,1))*(xcross-lowboundx(1,1))+(ycross-lowboundy(1,1))*(ycross-lowboundy(1,1)));
    markerxposlb(1:angldiv+1)=xcross+R2*cos(-anglall(1:angldiv+1)/180*pi);
    markeryposlb(1:angldiv+1)=ycross+R2*sin(-anglall(1:angldiv+1)/180*pi);
    
    plot(markerxposlb,markeryposlb,'ob');
    
    prompt = {'Enter the number of intersections between the upper and lower bound:'};
    dlg_title = 'Input for grid creation';
    num_lines= 1;
    def     = {'5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    Rdiv = str2num(cell2mat(answer(1,1)));
    
    Rstep=((R-R2)/Rdiv);
    Rall(1:Rdiv+1)=R2+Rstep*(1:Rdiv+1)-Rstep;
    
    grid_x=ones(Rdiv+1,angldiv+1)*xcross;
    grid_y=ones(Rdiv+1,angldiv+1)*ycross;
    A=Rall;
    B=cos(-anglall(1:angldiv+1)/180*pi);
    C=A'*B;
    grid_x=grid_x+Rall'*cos(-anglall(1:angldiv+1)/180*pi);
    grid_y=grid_y+Rall'*sin(-anglall(1:angldiv+1)/180*pi);
    
    close all
    imshow(im_grid,'truesize');
    hold on
    plot(grid_x,grid_y,'.b')    
    
    title(['Selected grid has ',num2str(angldiv*Rdiv), ' rasterpoints'])    % plot a title onto the image

    
    % Do you want to keep the grid?
    confirmselection = menu(sprintf('Do you want to use this grid?'),...
        'Yes','No, try again','Go back to grid-type selection');
    
    if confirmselection==1
        % Save settings and grid files in the image directory for visualization/plotting later
        %         save settings.dat xspacing yspacing xmin_new xmax_new ymin_new ymax_new -ascii -tabs
        save grid_x.dat grid_x -ascii -tabs
        save grid_y.dat grid_y -ascii -tabs
    end
    
    if confirmselection==2
        close all
        hold off
        imshow(im_grid,'truesize');
        circ_grid(FileNameBase, PathNameBase, im_grid);
    end
    
    if confirmselection==3
        gridtypeselection(FileNameBase, PathNameBase, im_grid);
    end
    
end


return



%-------------------------------
%
%rectangular selection requires selection of two points along the
%centerline of the beam and then asks for number of pixels to the left and
%right of the chosen points
%the rectangle is formed accordingly
function [grid_x,grid_y,FileNameBase,PathNameBase, numXelem, numYelem] = rect_grid(FileNameBase, PathNameBase, im_grid);
clc;
title(sprintf('Define the region of interest.  Pick (single click) a point in the LOWER LEFT region of the gage section.\n  Do the same for a point in the UPPER RIGHT portion of the gage section.'))

[a, b] = ginput(1);
hold on;
plot(a,b,'+b')
hold on;
[c, d] = ginput(1);
plot(c,d,'+b')

%select two grid points
prompt = 'Enter distance in pixels to the left and right of the selected point';
dlg_title = 'Enter distance in pixels to the left and right of the selected point';
num_lines = 1;
def = {'5'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
pix = str2num(cell2mat(answer(1,1)));
disp(pix);

x(1,1) = a - pix
y(1,1) = b
x(2,1) = c + pix
y(2,1) = d
hold on;
plot(x(1,1),y(1,1),'+b')
hold on;
plot(x(2,1),y(2,1),'+b')
hold on;

drawnow

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

lowerline=[xmin ymin; xmax ymin];
upperline=[xmin ymax; xmax ymax];
leftline=[xmin ymin; xmin ymax];
rightline=[xmax ymin; xmax ymax];

plot(lowerline(:,1),lowerline(:,2),'-b')
plot(upperline(:,1),upperline(:,2),'-b')
plot(leftline(:,1),leftline(:,2),'-b')
plot(rightline(:,1),rightline(:,2),'-b')

% closereq

cd(PathNameBase)

% Prompt user for grid spacing/resolution
prompt = {'Enter horizontal (x) resolution for image analysis [pixels]:', ...
        'Enter vertical (y) resolution for image analysis [pixels]:'};
dlg_title = 'Input for grid creation';
num_lines= 1;
def     = {'50','50'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
xspacing = str2num(cell2mat(answer(1,1)));
yspacing = str2num(cell2mat(answer(2,1)));

% Round xmin,xmax and ymin,ymax "up" based on selected spacing
numXelem = ceil((xmax-xmin)/xspacing)-1;
numYelem = ceil((ymax-ymin)/yspacing)-1;
num = [numXelem, numYelem+1]
sprintf('Enter in first dialog box(number of elements in one row) = %d',numXelem)
sprintf('Enter in second dialog box(number of elements in one column) = %d',numYelem+1)
save num.dat num -ascii -tabs

xmin_new = (xmax+xmin)/2-((numXelem/2)*xspacing);
xmax_new = (xmax+xmin)/2+((numXelem/2)*xspacing);
ymin_new = (ymax+ymin)/2-((numYelem/2)*yspacing);
ymax_new = (ymax+ymin)/2+((numYelem/2)*yspacing);

% Create the analysis grid and show user
[x,y] = meshgrid(xmin_new:xspacing:xmax_new,ymin_new:yspacing:ymax_new);
[rows columns] = size(x);
zdummy = 200.*ones(rows,columns);
imshow(FileNameBase);
title(['Selected grid has ',num2str(rows*columns), ' rasterpoints'])    % plot a title onto the image
hold on;
plot(x,y,'+b')

grid_x=x;
grid_y=y;
% Do you want to keep the grid?
confirmselection = menu(sprintf('Do you want to use this grid?'),...
    'Yes','No, try again','Go back to grid-type selection');

if confirmselection==1
    % Save settings and grid files in the image directory for visualization/plotting later
    save settings.dat xspacing yspacing xmin_new xmax_new ymin_new ymax_new -ascii -tabs
    save grid_x.dat x -ascii -tabs
    save grid_y.dat y -ascii -tabs
    close all
    hold off
end

if confirmselection==2
    close all
    hold off
    imshow(im_grid,'truesize');
    rect_grid(FileNameBase, PathNameBase, im_grid);
end

if confirmselection==3
    close all
    hold off
    gridtypeselection(FileNameBase, PathNameBase, im_grid);
end
