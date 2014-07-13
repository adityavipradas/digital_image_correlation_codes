%split tensile load Vs displacement

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
def = {'4.5'};
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

%excel datasheet
prompt = 'Enter name of the excel datasheet';
dlg_title = 'Enter name of the excel datasheet';
num_lines = 1
def = {'D2 Curves.xls'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
excel = cell2mat(answer(1,1));
disp(excel);

[num] = xlsread(excel)
remove = num(isfinite(num(:,7)), :)
load = remove(:,7)
load = transpose(load);
disp(epsilonMat)
disp(load)
load = 10 * load

%least count on x-axis
prompt = 'Enter least count for x-axis';
dlg_title = 'Enter least count for x-axis';
num_lines = 1
def = {'0.005'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
lsx = str2num(cell2mat(answer(1,1)));
disp(lsx);

%specimen name
prompt = 'Enter specimen name';
dlg_title = 'Enter specimen name';
num_lines = 1
def = {'CSRE-300-18.5-0.20d-D'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
specimen = cell2mat(answer(1,1));
disp(specimen);
sprintf('%s',specimen)

%curve fit
prompt = 'Enter degree of fit';
dlg_title = 'Enter degree of fit';
num_lines = 1
def = {'3'};
answer = inputdlg(prompt, dlg_title,num_lines,def);
cfit = str2num(cell2mat(answer(1,1)));
disp(cfit);

pp = polyfit(load, epsilonMat, cfit);
disp(pp);
curvefit = zeros(1, sizevalidx(1,2));
coeff = 1
for image=1:1:sizevalidx(1,2)
    for i=cfit:-1:0
        curvefit(1, image) = curvefit(1, image) + (load(1, image)^i) * pp(coeff);
        coeff = coeff + 1;
    end
    coeff = 1;
end

%result of interest
sprintf('x-displacement corresponding to maximum load is %d mm',curvefit(1, sizevalidx(1,2)))

%plot
figure();
plot(curvefit, load, 'color', 'red', 'Marker', '*');
hold on;
plot(epsilonMat, load, '-o')
grid on
title(sprintf('load Vs displacement (%s)', specimen))
%set(gca, 'XTick', min(min(epsilonMat)):lsx: max(max(epsilonMat))+lsx);
legend(sprintf('Regression of order %d',cfit),'DIC');
xlabel('x-displacement(mm)')
ylabel('load(N)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%