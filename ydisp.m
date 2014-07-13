%load point displacement and fracture energy for a beam
%LOAD VS LPD PLOTS

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

[validx, validy,displx,disply]=lpd(validx, validy,displx,disply);

%---------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load point displacement at beam notch
% written by Aditya Vipradas
function [validx, validy,displx,disply] = lpd(validx, validy,displx,disply) ;

clear xplot
clear sizevalidx
clear selection1
clear selection2
clear badpoints
sizevalidx=size(validx);

%pixels to mm conversion
prompt = 'Number of pixels corresponding to 1mm';
dlg_title = 'Number of pixels corresponding to 1mm';
num_lines = 1;
def = {'5.6'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
onemm = str2num(cell2mat(answer(1,1)));

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
disp(load)
load = 10 * load

%least count on x-axis
prompt = 'Enter least count for x-axis';
dlg_title = 'Enter least count for x-axis';
num_lines = 1
def = {'0.05'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
lsx = str2num(cell2mat(answer(1,1)));
disp(lsx);

%find cubic interpolation function of the load against lpd curve
% remove duplicate data for executing spline function and form in dvar for disply and lvar for load

r = 1;
count = 0;

for i=1:1:sizevalidx(1,2)
    for j=i+1:1:sizevalidx(1,2)
        if disply(1,i) ~= disply(1,j)
            count = 0;
        else
            count = 1;
            break;
        end
    end
    if count == 0
        dvar(1,r) = disply(1,i)/onemm;
        lvar(1,r) = load(i,1);
        r = r + 1;
    end
    count = 0;
end

%perform sorting in ascending order for dvar and then sort lvar accordingly
%in sortdvar and sortlvar respectively
sortdvar = sort(dvar);
for i=1:1:length(sortdvar)
    if dvar(1,i) ~= sortdvar(1,i)
        for j=1:1:length(sortdvar)
            if sortdvar(1,i) == dvar(1,j)
                sortlvar(1,i) = lvar(1,j);
                break;
            end
        end
    else
        sortlvar(1,i) = lvar(1,i);
    end
end
disp(dvar)
disp(sortdvar)
disp(lvar)
disp(sortlvar)

%cubic spline interpolation
pp = csapi(sortdvar, sortlvar);
func = fnint(pp);
fracture_energy = fnval(func, sortdvar(1,length(sortdvar)) - fnval(func, sortdvar(1,1)));

%area of the beam above the notch in mm2
prompt = 'Enter area of the beam above the notch in mm2';
dlg_title = 'Enter area of the beam above the notch in mm2';
num_lines = 1
def = {'18000'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
notch_area = str2num(cell2mat(answer(1,1)));
disp(notch_area);
%fracture energy in N/m
fracture_energy = (fracture_energy/notch_area) * 1000;
disp(sprintf('Fracture energy is %d N/m',fracture_energy))

%specimen name
prompt = 'Enter specimen name';
dlg_title = 'Enter specimen name';
num_lines = 1
def = {'CSRE-300-18.5-0.20d-D'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
specimen = cell2mat(answer(1,1));
disp(specimen);
sprintf('%s',specimen)

row = 1;
for i=sortdvar(1,1): 0.01 :sortdvar(1,length(sortdvar))
    x(row,1) = i;
    y(row,1) = fnval(pp, i);
    row = row + 1;
end

%consider left marker of the two selected markers for plotting
%plot
figure();
hold on;
plot(disply(1,:)/onemm, load, '-o')
grid on
plot(num(:,1), num(:,2), 'red')
title(sprintf('load Vs load-point displacement (%s)',specimen))
set(gca, 'XTick', min(disply(1,:)):lsx: max(disply(1,:))+lsx);
legend('DIC','experimental');
xlabel('load-point displacement(mm)')
ylabel('load(N)')

%cross check
figure();
hold on;
plot(disply(1,:)/onemm, load, '-o')
grid on
plot(x, y, 'red')
legend('DIC', 'cubic interpolation');
title('relation between actual and interpolated curve');
set(gca, 'XTick', min(disply(1,:)):lsx: max(disply(1,:))+lsx);
xlabel('load-point displacement(mm)')
ylabel('load(N)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%