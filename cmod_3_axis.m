%3axes CMOD plots
%image load CMOD

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
clear displx;
validxfirst=zeros(size(validx));
validxfirst=mean(validx(:,1),2)*ones(1,sizevalidx(1,2));
displx=validx-validxfirst;
clear validxfirst
clear disply;
validyfirst=zeros(size(validy));
validyfirst=mean(validy(:,1),2)*ones(1,sizevalidy(1,2));
disply=validy-validyfirst;
clear validyfirst

save displx.dat displx -ascii -tabs
save disply.dat disply -ascii -tabs

[validx, validy,displx,disply]=mouth_crack(validx, validy,displx,disply);

%---------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by Aditya Vipradas
function [validx, validy,displx,disply] = mouth_crack(validx, validy,displx,disply) ;

clear xplot
clear sizevalidx
clear selection1
clear selection2
clear badpoints
sizevalidx=size(validx);

% update temporary matrices
displxtemp=displx;
validxtemp=validx;
validytemp=validy;

%epsilonMat holds cmod values
epsilonMat = zeros(sizevalidx(1,1)/2,sizevalidx(1,2));

%pixels to mm conversion
prompt = 'Number of pixels corresponding to 1mm';
dlg_title = 'Number of pixels corresponding to 1mm';
num_lines = 1;
def = {'5.6'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
onemm = str2num(cell2mat(answer(1,1)));

for badpointnum=sizevalidx(1,1)/2:-1:1
    badpointnum2 = badpointnum + sizevalidx(1,1)/2;
    epsilon1D = (displxtemp(badpointnum,:)-displxtemp(badpointnum2,:))/onemm
        
%positive for tensile and negative for compressive deformations i.e. give sign conventions
    for i=1:1:length(epsilon1D)
        %stretch for badpointnum and stretch for badpointnum2...tensile
        if displxtemp(badpointnum, i) < 0 & displxtemp(badpointnum2, i) > 0
            epsilon1D(1, i) = abs(epsilon1D(1,i))
        %press for badpointnum and press for badpointnum2...compressive
        elseif displxtemp(badpointnum, i) > 0 & displxtemp(badpointnum2, i) < 0 & epsilon1D(1, i) > 0
            epsilon1D(1,i) = -1 * epsilon1D(1,i)
        %press for badpointnum and stretch for badpointnum2 but disp of badpointnum is
        %greater than badpointnum2...compressive
        elseif displxtemp(badpointnum, i) > 0 & displxtemp(badpointnum2, i) > 0 & abs(displxtemp(badpointnum, i)) > abs(displxtemp(badpointnum2, i))
            if epsilon1D(1,i) > 0
                epsilon1D(1,i) = -1 * epsilon1D(1,i)
            end
        %press for badpointnum and stretch for badpointnum2 but disp of badpointnum is
        %less than badpointnum2...tensile
        elseif displxtemp(badpointnum, i) > 0 & displxtemp(badpointnum2, i) > 0 & abs(displxtemp(badpointnum, i)) < abs(displxtemp(badpointnum2, i))
            epsilon1D(1, i) = abs(epsilon1D(1,i))
        %stretch for badpointnum and press for badpointnum2 but disp of badpointnum is
        %greater than badpointnum2 ...tensile
        elseif displxtemp(badpointnum, i) < 0 & displxtemp(badpointnum2, i) < 0 & abs(displxtemp(badpointnum, i)) > abs(displxtemp(badpointnum2, i))
            epsilon1D(1, i) = abs(epsilon1D(1,i))
        %stretch for badpointnum and press for badpointnum2 but disp of badpointnum is
        %less than badpointnum2 ...compressive
        elseif displxtemp(badpointnum, i) < 0 & displxtemp(badpointnum2, i) < 0 & abs(displxtemp(badpointnum, i)) < abs(displxtemp(badpointnum2, i))
            if epsilon1D(1,i) > 0
                epsilon1D(1,i) = -1 * epsilon1D(1,i)
            end
        %badpointnum remains as it is but badpointnum2 is
        %stretched...tensile
        elseif displxtemp(badpointnum, i) == 0 & displxtemp(badpointnum2, i) > 0
            epsilon1D(1, i) = abs(epsilon1D(1,i))
        %badpointnum remains as it is but badpointnum2 is
        %pressed...compressive
        elseif displxtemp(badpointnum, i) == 0 & displxtemp(badpointnum2, i) < 0
            if epsilon1D(1,i) > 0
                epsilon1D(1,i) = -1 * epsilon1D(1,i)
            end
        %badpointnum is pressed and badpointnum2 remains as it
        %is...compressive
        elseif displxtemp(badpointnum, i) > 0 & displxtemp(badpointnum2, i) == 0
             if epsilon1D(1,i) > 0
                epsilon1D(1,i) = -1 * epsilon1D(1,i)
             end
        %badpointnum is stretched but badpointnum2 remains as it
        %is...tensile
        elseif displxtemp(badpointnum, i) < 0 & displxtemp(badpointnum2, i) == 0
            epsilon1D(1, i) = abs(epsilon1D(1,i))
        end
    end
    epsilonMat(badpointnum,:) = epsilon1D;
end

%read data for DIC loads
%excel datasheet
prompt = 'Enter name of the excel datasheet';
dlg_title = 'Enter name of the excel datasheet';
num_lines = 1
def = {'D3 Curves.xls'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
excel = cell2mat(answer(1,1));
disp(excel);

[num] = xlsread(excel)
remove = num(isfinite(num(:,7)), :)
load = remove(:,7)
disp(epsilonMat)
disp(load)
load = 10 * load


%perform separate interpolations for data above and below the max load
%because each load value has two displacement values on the curve(bell
%shaped)
maxrow = find(ismember(num(:,5), max(num(:,5))), 1);
maxdic = find(ismember(load, max(load)), 1);
disp(maxrow);
disp(maxdic);

%curve to the left of the experimental maximum load value
fleft = csapi(num(1:maxrow, 5), num(1:maxrow, 4))
for i=1:1:maxdic
    dispinter(i,1) = fnval(fleft, load(i,1));
    image(i,1) = i;
end

%curve to the right of the experimental maximum load value
fright = csapi(num(maxrow:length(num(:,5)), 5), num(maxrow:length(num(:,5)), 4))
for i=maxdic+1:1:length(load)
    dispinter(i,1) = fnval(fright, load(i,1));
    image(i,1) = i;
end

%least count on x-axis
prompt = 'Enter least count for image';
dlg_title = 'Enter least count for image';
num_lines = 1
def = {'1'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
lsx = str2num(cell2mat(answer(1,1)));
disp(lsx);

%least count on load y-axis
prompt = 'Enter least count for load y-axis';
dlg_title = 'Enter least count for load y-axis';
num_lines = 1
def = {'100'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
Plsy = str2num(cell2mat(answer(1,1)));
disp(Plsy);

%least count on disp y-axis
prompt = 'Enter least count for cmod y-axis';
dlg_title = 'Enter least count for cmod y-axis';
num_lines = 1
def = {'1'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
lsy = str2num(cell2mat(answer(1,1)));
disp(lsy);

%specimen name
prompt = 'Enter specimen name';
dlg_title = 'Enter specimen name';
num_lines = 1
def = {'CSRE-300-18.5-0.20d-D'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
specimen = cell2mat(answer(1,1));
disp(specimen);
sprintf('%s',specimen)

%cross check interpolation
figure();
plot(num(:,5), num(:,4))
hold on;
plot(load, dispinter, 'red')
grid on;
legend('experimental data', 'cubic interpolation');
title('relation between experimental and interpolated data');
xlabel('load(N)');
ylabel('CMOD(mm)');

%plot 2 y axes in one graph
figure();
[ax1, h1, h2] = plotyy(image, load, image, epsilonMat(1,:));
hold on;
[ax2, h3, h4] = plotyy(image, load, image, dispinter(:,1));
grid on;

%set axis limits and scale
if min(epsilonMat(1,:)) < min(dispinter(:,1))
    set(ax1(2), 'Ylim', [0, ceil(max(epsilonMat(1,:))) + 5]);
    set(ax2(2), 'Ylim', [0, ceil(max(epsilonMat(1,:))) + 5]);
    set(ax1(2), 'YTick', [0: lsy : ceil(max(epsilonMat(1,:))) + 5]);
    set(ax2(2), 'YTick', [0: lsy : ceil(max(epsilonMat(1,:))) + 5]);

else
    set(ax1(2), 'Ylim', [0, ceil(max(dispinter(:,1))) + 5]);
    set(ax2(2), 'Ylim', [0, ceil(max(dispinter(:,1))) + 5]);
    set(ax1(2), 'YTick', [0: lsy : ceil(max(dispinter(:,1))) + 5]);
    set(ax2(2), 'YTick', [0: lsy : ceil(max(dispinter(:,1))) + 5]);
end
set(ax1(1), 'Xlim', [1, max(image)]);
set(ax1(2), 'Xlim', [1, max(image)]);
set(ax2(1), 'Xlim', [1, max(image)]);
set(ax2(2), 'Xlim', [1, max(image)]);
set(ax1(1), 'XTick', [1: lsx : max(image)]);
set(ax1(2), 'XTick', [1: lsx : max(image)]);
set(ax2(1), 'XTick', [1: lsx : max(image)]);
set(ax2(2), 'XTick', [1: lsx : max(image)]);
set(ax1(1), 'YTick', [0: Plsy : max(load) + Plsy]);
set(ax2(1), 'YTick', [0: Plsy : max(load) + Plsy]);

% %set line type
% set(h2, 'Linestyle',':');
% set(h4, 'Linestyle','--');
set(get(ax1(1), 'Xlabel'),'String','image');
set(get(ax1(1), 'Ylabel'),'String','load(N)');
set(get(ax1(2), 'Ylabel'),'String','load-point displacement(mm)');
set(h1, 'Marker','*');
set(h2, 'Marker','o', 'MarkerSize',2);
set(h4, 'Marker','+');
set(h1, 'color', 'black');
set(h2, 'color', 'magenta');
set(h4, 'color', 'green');
legend([h1;h2;h4], 'load Vs image','DIC Vs image', 'experimental Vs image');
title(sprintf('CMOD (%s)',specimen));

%choose colors for plots
ask = menu(sprintf('Do you want to choose different colours for the 2 y-axes plot?'),'Yes','No');
if ask == 1
    reply = 2
    while reply==2
        c1 = uisetcolor('color for load Vs image');
        set(h1, 'color', c1);
        set(h3, 'color', c1);
        c2 = uisetcolor('color for DIC Vs image');
        set(h2, 'color', c2);
        c3 = uisetcolor('color for experimental Vs image');
        set(h4, 'color', c3);
        reply = menu(sprintf('Are the colors OK?'),...
        'Yes','No, try again');
        if reply==1
            break;
        end
    end
else
    set(h1, 'color', 'black');
    set(h2, 'color', 'magenta');
    set(h4, 'color', 'green');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%