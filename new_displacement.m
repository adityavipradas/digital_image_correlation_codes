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
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);

%calculate the displacement relative to the first image in x and y
%direction
clear displx;
validxfirst=zeros(size(validx));
validxfirst=mean(validx(:,1),2)*ones(1,sizevalidx(1,2));
displx=validx-validxfirst;
disp(displx)
clear validxfirst
clear disply;
validyfirst=zeros(size(validy));
validyfirst=mean(validy(:,1),2)*ones(1,sizevalidy(1,2));
disply=validy-validyfirst;
clear validyfirst

save displx.dat displx -ascii -tabs
save disply.dat disply -ascii -tabs

[validx, validy,displx,disply]=strain_1D_2Points_func(validx, validy,displx,disply);
   
%---------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strain between two markers
% written by Chris
function [validx, validy,displx,disply] = strain_1D_2Points_func(validx, validy,displx,disply) ; % 1D strain calculation

sizevalidx=size(validx);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);
defaultimage=loopimages;
numberbadpoints=1;


clear xplot
clear sizevalidx
clear selection1
clear selection2
clear badpoints
sizevalidx=size(validx);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);

% update temporary matrices
displxtemp=displx;
validxtemp=validx;
validytemp=validy;

epsilonMat = zeros(sizevalidx(1,1)/2,sizevalidx(1,2));
ycoord = zeros(sizevalidx(1,1)/2, 1);

%pixels to mm conversion
prompt = 'Number of pixels corresponding to 1mm';
dlg_title = 'Number of pixels corresponding to 1mm';
num_lines = 1;
def = {'5.6'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
onemm = str2num(cell2mat(answer(1,1)));

for badpointnum=sizevalidx(1,1)/2:-1:1
    y = (validy(sizevalidx(1,1),1) - validy(badpointnum, 1))/onemm;
    ycoord(badpointnum, 1) = y;
    badpointnum2 = badpointnum + sizevalidx(1,1)/2;
    epsilon1D=(displxtemp(badpointnum,:)-displxtemp(badpointnum2,:))/onemm;
    disp(epsilon1D)
    
%positive for tensile and negative for compressive strains i.e. give sign conventions
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
disp(epsilonMat)


%least count on x-axis
prompt = 'Enter least count for x-axis';
dlg_title = 'Enter least count for x-axis';
num_lines = 1
def = {'50'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
lsx = str2num(cell2mat(answer(1,1)));
disp(lsx);


%least count on y-axis
prompt = 'Enter least count for y-axis';
dlg_title = 'Enter least count for y-axis';
num_lines = 1
def = {'0.05'};
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


mkdir('disp_against_depth')
cd('disp_against_depth')
Vid = 'dispdepth'

%plot
cc = hsv(sizevalidx(1,2));
h = figure;
set(gca, 'YTick', min(min(epsilonMat)):lsy: max(max(epsilonMat))+lsy);
%set(gca, 'XTick', ycoord(sizevalidx(1,1)/2, 1):lsx:ycoord(1,1) + lsx);
hold on;
grid on;
title(sprintf('COD Vs depth plots for all images (%s)',specimen));
xlabel('depth(mm)');
ylabel('displacement(mm)');

%all in one
for i=1:1:sizevalidx(1,2)
    plot(ycoord, epsilonMat(:,i), '-o','MarkerSize',2,'color',cc(i,:))
end
u = 0
ustr = num2str(u)
videoname = [Vid ustr 'jpg']
saveas(h, videoname, 'jpg')
set(clf,'visible','off')

%individual plots
for i=1:1:sizevalidx(1,2)
    h = figure();
    p = polyfit(ycoord, epsilonMat(:,i), 1);
    disp(p)
    for j=1:1:length(ycoord)
        line_fit(j,1) = p(1) * ycoord(j,1) + p(2)
    end
    plot(ycoord, epsilonMat(:,i), '-o','MarkerSize',2)
    hold on;
    plot(ycoord, line_fit, 'red')
    grid on;
    legend('actual displacements', 'linear regression');
    title(sprintf('COD Vs depth for image %d (%s)',i,specimen));
    xlabel('depth(mm)');
    ylabel('COD(mm)');
    u = i
    ustr=num2str(u);
    videoname=[Vid ustr 'jpg']
    saveas(h,videoname,'jpg')
    set(clf,'visible','off')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
