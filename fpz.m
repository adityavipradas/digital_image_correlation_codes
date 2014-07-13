%DISPLACEMENTS IN X_DIRECTION VS DEPTH PLOTS

% Initialize data
% written by Chris and Dan

% customised by Aditya Vipradas

% Displacement.m allows you to analyze the data you aquiered with the
% correlation, fitting or mean routine. It only needs the validx and
% validy and can calculate strain from it. Before you start you should 
% consider cleaning up the data as described in the guide. After that step
% you can analyze parts of your data, or the full set. Try to use also the
% console command, e.g. if you want to analyze only image 100-110 since
% something really interesting happend there, load validx and validy into
% your workspace and call
% displacement(validx(:,100:110),validy(:,100:110));
% In this case displacement only loads the important images and you can
% clean this part of your data set.

% Changed 3. February 2008


function [validx,validy]=displacement(validx,validy);

%load data in case you did not load it into workspace yet
if exist('validx')==0
    [validxname,Pathvalidx] = uigetfile('*.dat','Open validx.dat');
    if validxname==0
        disp('You did not select a file!')
        return
    end
    cd(Pathvalidx);
    validx=importdata(validxname,'\t');
end
if exist('validy')==0
    [validyname,Pathvalidy] = uigetfile('*.dat','Open validy.dat');
    if validyname==0
        disp('You did not select a file!')
        return
    end
    cd(Pathvalidy);
    validy=importdata(validyname,'\t');
end

%define the size of the data set
sizevalidx=size(validx);
sizevalidy=size(validy);

%calculate the displacement relative to the first image in x and y
%direction

validxfirst=zeros(size(validx));
validxfirst=mean(validx(:,1),2)*ones(1,sizevalidx(1,2));
displx=validx-validxfirst;
validyfirst=zeros(size(validy));
validyfirst=mean(validy(:,1),2)*ones(1,sizevalidy(1,2));
disply=validy-validyfirst;

save displx.dat displx -ascii -tabs
save disply.dat disply -ascii -tabs

% update temporary matrices
validxtemp=validx;
validytemp=validy;

%pixels to mm conversion
prompt = 'Number of pixels corresponding to 1mm';
dlg_title = 'Number of pixels corresponding to 1mm';
num_lines = 1;
def = {'5.6'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
onemm = str2num(cell2mat(answer(1,1)));

displxtemp = displx./onemm;

same_x = size(load('grid_x.dat'));
same_x = same_x(1,1);
xcoord = [];
xdisp = [];
count = 0;

for i = same_x:same_x:sizevalidx(1,1)
    count = count + 1;
    xcoord(count,1) = validxtemp(i,1);
end

if rem(count, 2)~=0
    origin = xcoord((count+1)/2, 1);
else
    origin = (xcoord(count/2, 1) + xcoord((count/2)+1, 1))/2;
end

sizexcoord = size(xcoord);

xcoordplot = [];

for i = 1:1:sizexcoord(1,1)
    xcoordplot(i,1) = (xcoord(i,1) - origin)/onemm;
end

cnt1 = 0;
for i = 1:1:sizevalidx(1,2)
    for j = same_x:same_x:sizevalidx(1,1)
        cnt1 = cnt1 + 1;
        xdisp(cnt1,i) = displxtemp(j,i);
    end
    cnt1 = 0;
end

%least count on x-axis
prompt = 'Enter least count for x-axis';
dlg_title = 'Enter least count for x-axis';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
lsx = str2num(cell2mat(answer(1,1)));

%least count on y-axis
prompt = 'Enter least count for y-axis';
dlg_title = 'Enter least count for y-axis';
num_lines = 1;
def = {'0.05'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
lsy = str2num(cell2mat(answer(1,1)));

%specimen name
prompt = 'Enter specimen name';
dlg_title = 'Enter specimen name';
num_lines = 1;
def = {'CSRE-300-18.5-0.20d-D'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
specimen = cell2mat(answer(1,1));
sprintf('%s',specimen);


mkdir('fpz')
cd('fpz')
Vid = 'fpz'
Vid1 = 'fpzfig'

%individual
for i=1:1:sizevalidx(1,2)
    h = figure();
    plot(xcoordplot, xdisp(:,i), '-o','MarkerSize',2);
    set(gca, 'FontSize', 13)
    ylim([min(min(xdisp)) max(max(xdisp))+lsy])
    grid on;
    title(sprintf('Displacements along line MN for image %d (%s)',i,specimen));
    xlabel('line MN (mm)');
    ylabel('displacements in x-direction (mm)');
    u = i;
    ustr = num2str(u);
    videoname = [Vid ustr 'jpg'];
    videoname1 = [Vid1 ustr 'fig'];
    saveas(h, videoname, 'jpg');
    saveas(h, videoname1, 'fig');
    set(clf,'visible','off');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
