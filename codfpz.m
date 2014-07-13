% Initialize data
% written by Chris and Dan

% customised by Aditya Vipradas
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
clear validxfirst;
clear disply;
validyfirst=zeros(size(validy));
validyfirst=mean(validy(:,1),2)*ones(1,sizevalidy(1,2));
disply=validy-validyfirst;
clear validyfirst;

save displx.dat displx -ascii -tabs
save disply.dat disply -ascii -tabs

% update temporary matrices
displxtemp=displx;
validxtemp=validx;
validytemp=validy;

%pixels to mm conversion
prompt = 'Number of pixels corresponding to 1mm';
dlg_title = 'Number of pixels corresponding to 1mm';
num_lines = 1;
def = {'5.6'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
onemm = str2num(cell2mat(answer(1,1)));

%notch dimensions
prompt = 'Notch width in mm';
dlg_title = 'Notch width in mm';
num_lines = 1;
def = {'2'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
nwidth = str2num(cell2mat(answer(1,1)));

prompt = 'Notch height in mm';
dlg_title = 'Notch height in mm';
num_lines = 1;
def = {'60'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
nheight = str2num(cell2mat(answer(1,1)));

same_x = size(load('grid_x.dat'));
same_x = same_x(1,1);

leftmode = [];
rightmode = [];
xcoord = [];
ycoord = [];
count = 0;
originy = validytemp(same_x, 1);
cod = [];

for i = same_x:same_x:sizevalidx(1,1)
    count = count + 1;
    xcoord(count,1) = validxtemp(i,1);
end

if rem(count, 2)~=0
    originx = xcoord((count+1)/2, 1);
else
    originx = (xcoord(count/2, 1) + xcoord((count/2)+1, 1))/2;
end

for i = 1:1:sizevalidx(1,2)
    for j = 1:1:same_x
        countl = 0;
        countr = 0;
        lm = [];
        rm = [];
        for k = j:same_x:sizevalidx(1,1)
            if validxtemp(k,1) <= originx
                countl = countl + 1;
                lm(1, countl) = displxtemp(k,i);
            else
                countr = countr + 1;
                rm(1, countr) = displxtemp(k,i);
            end
        end
        leftmode(j,i) = mean(lm)/onemm;
        rightmode(j,i) = mean(rm)/onemm;
        %COD
        tempdisp = leftmode(j,i) - rightmode(j,i);
        %positive for tensile and negative for compressive deformations i.e. give sign conventions
        %stretch for badpointnum and stretch for badpointnum2...tensile
        if leftmode(j,i) < 0 & rightmode(j,i) > 0
            tempdisp = abs(tempdisp);
        %press for badpointnum and press for badpointnum2...compressive
        elseif leftmode(j,i) > 0 & rightmode(j,i) < 0 & tempdisp > 0
            tempdisp = -1 * tempdisp;
        %press for badpointnum and stretch for badpointnum2 but disp of badpointnum is
        %greater than badpointnum2...compressive
        elseif leftmode(j,i) > 0 & rightmode(j,i) > 0 & abs(leftmode(j,i)) > abs(rightmode(j,i))
            if tempdisp > 0
                tempdisp = -1 * tempdisp
            end
        %press for badpointnum and stretch for badpointnum2 but disp of badpointnum is
        %less than badpointnum2...tensile
        elseif leftmode(j,i) > 0 & rightmode(j,i) > 0 & abs(leftmode(j,i)) < abs(rightmode(j,i))
            tempdisp = abs(tempdisp)
        %stretch for badpointnum and press for badpointnum2 but disp of badpointnum is
        %greater than badpointnum2 ...tensile
        elseif leftmode(j,i) < 0 & rightmode(j,i) < 0 & abs(leftmode(j,i)) > abs(rightmode(j,i))
            tempdisp = abs(tempdisp)
        %stretch for badpointnum and press for badpointnum2 but disp of badpointnum is
        %less than badpointnum2 ...compressive
        elseif leftmode(j,i) < 0 & rightmode(j,i) < 0 & abs(leftmode(j,i)) < abs(rightmode(j,i))
            if tempdisp > 0
                tempdisp = -1 * tempdisp
            end
        %badpointnum remains as it is but badpointnum2 is
        %stretched...tensile
        elseif leftmode(j,i) == 0 & rightmode(j,i) > 0
            tempdisp = abs(tempdisp)
        %badpointnum remains as it is but badpointnum2 is
        %pressed...compressive
        elseif leftmode(j,i) == 0 & rightmode(j,i) < 0
            if tempdisp > 0
                tempdisp = -1 * tempdisp
            end
        %badpointnum is pressed and badpointnum2 remains as it
        %is...compressive
        elseif leftmode(j,i) > 0 & rightmode(j,i) == 0
             if tempdisp > 0
                tempdisp = -1 * tempdisp
             end
        %badpointnum is stretched but badpointnum2 remains as it
        %is...tensile
        elseif leftmode(j,i) < 0 & rightmode(j,i) == 0
            tempdisp = abs(tempdisp)
        end
        cod(j,i) = tempdisp;
    end
end

disp(leftmode);
disp(rightmode);

for i = 1:1:same_x
    ycoord(i,1) = ((originy - validytemp(i,1))/onemm) + nheight;
end

mkdir('codfpz')
cd('codfpz')
Vid = 'codfpz'
Vid1 = 'codplot'

prompt = 'State equal positive and negative limits on x-axis';
dlg_title = 'State positive and negative limits on x-axis';
num_lines = 1;
def = {'50'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
xlimit = str2num(cell2mat(answer(1,1)));

prompt = 'State upper limit on the y-axis';
dlg_title = 'State upper limit on the y-axis';
num_lines = 1;
def = {'500'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
ylimit = str2num(cell2mat(answer(1,1)));

%specimen name
prompt = 'Enter specimen name';
dlg_title = 'Enter specimen name';
num_lines = 1;
def = {'CSRE-300-18.5-0.20d-D'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
specimen = cell2mat(answer(1,1));
sprintf('%s',specimen);

%plot
for i = 1:1:sizevalidx(1,2)
    h = figure();
    plot([-nwidth/2, -nwidth/2, nwidth/2, nwidth/2], [0, nheight, nheight, 0])
    hold on;
    for j = 1:1:same_x
        plot([leftmode(j,i), rightmode(j,i)],[ycoord(j,1), ycoord(j,1)])
        xlim([-xlimit xlimit]);
        ylim([0 ylimit]);
        hold on;
        grid on;
    end
    title(sprintf('Development of FPZ for image %d (%s)',i,specimen));
    xlabel('beam span(mm)');
    ylabel('beam depth(mm)');
    u = i
    ustr=num2str(u);
    videoname=[Vid ustr 'jpg']
    saveas(h,videoname,'jpg')
    set(clf,'visible','off')
    h1 = figure();
    for j = 1:1:same_x
        plot([cod(j,i), 0],[ycoord(j,1), ycoord(j,1)])
        hold on;
        grid on;
    end
    videoname1 = [Vid1 ustr 'jpg']
    saveas(h1, videoname1,'jpg')
    set(clf,'visible','off')
end