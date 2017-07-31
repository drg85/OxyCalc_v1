




function Ox_Calc_Vectorized

%OX-CALC - Oxygen distribution tool. Version 1.0 
%Authors: Dr. David Robert Grimes & Dr. Daniel R. Warren

%This code allows the estimation of tissue oxygen from known vessel
%networks. Inputs are text files, formatted in the style outlined by Secomb
%et al. 

%Queries: E-mail davidrobert.grimes@oncology.ox.ac.uk or @drg1985

%resolution
choiceint = menu('Choose a grid resolution','5um - Rapidly computes, some loss of accuracy','2um - quick to compute, good accuracy','1um - Very accurate but can be slow to compute');

if choiceint == 1; 
    res = 5;
elseif choiceint == 2
    res = 2;
elseif choiceint == 3;
    res = 1;
else
    disp('Error! No resolution selected!')
end

clear choiceint

%Constants
D = 2*10^- 9;  %Diffusion of O2 in water m^2/s
om = 3.0318*10^7; %Omega, arising from Henry's law (See http://rsif.royalsocietypublishing.org/content/11/92/20131124.full for full derivation!) 
ros = 5*10^-6; % 'default' vessel radius 


prompt = {'Enter typical diffusion distance in microns from 5 micron radius source vessel:','Enter partial pressure at wall of typical 5 micron radius vessel in mmHg:'};
dlg_title = 'Standard inputs - 3D oxygen solver with 5um grid spacing';
num_lines = 1;
defaultans = {'150','40'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

RN = str2double(answer{1,1});
rn = RN*10^-6;

po = str2double(answer{2,1});

asl = ((3*D*po)/om)*((((sqrt(rn^2 - ros^2))/3)*(2*ros^2 - 8*rn^2) + 2*rn^3*log((rn + sqrt(rn^2 - ros^2))/ros))^-1);

%[filename,~]=uigetfile('*.txt','Select Vessel Network file');


[stmfile, stmpath] = uigetfile('*.txt', 'Select Vessel Network file');
pullin = importdata(fullfile(stmpath, stmfile));



%pullin = importdata(filename);

Grid_data = pullin.data; 

Segnum = size(Grid_data);
A = Segnum(1);
Segnum = A; 
clear A pullin filename pathname prompr dlg_title num_line defaultans answer 

%Now if correct vessel matrix selected, we jettison the first column
%(vessel numbers) and node data (2 + 3rd column when applicable) and round the values to the nearest micron! 

%VesselData format columns = X1 Y1 Z1 X2 Y2 Z2 radius Length! 

VesselData = zeros(Segnum,8); 

VesselData(:,1) = Grid_data(:,6); %get X1
VesselData(:,2) = Grid_data(:,7); %get Y1
VesselData(:,3) = Grid_data(:,8); %get Z1
VesselData(:,4) = Grid_data(:,9); %get X2
VesselData(:,5) = Grid_data(:,10); %get Y2
VesselData(:,6) = Grid_data(:,11); %get Z2
VesselData(:,7) = (Grid_data(:,4))./2; %get radius
VesselData(:,8) = sqrt((VesselData(:,1) - VesselData(:,4)).^2 + (VesselData(:,2) - VesselData(:,5)).^2 + (VesselData(:,3) - VesselData(:,6)).^2); %get length

clear Grid_data 

VesselData = round(VesselData); 

offset = max(VesselData(:,7)) + 1; %Offset to avoid out of bounds error
VesselData(:,1) = VesselData(:,1) + offset; 
VesselData(:,2) = VesselData(:,2) + offset;
VesselData(:,3) = VesselData(:,3) + offset;
VesselData(:,4) = VesselData(:,4) + offset; 
VesselData(:,5) = VesselData(:,5) + offset; 
VesselData(:,6) = VesselData(:,6) + offset;


xmax = [max(VesselData(:,1)) max(VesselData(:,4))];
xmax = max(xmax)+ 2*offset; 
ymax = [max(VesselData(:,2)) max(VesselData(:,5))];
ymax = max(ymax) + 2*offset; 
zmax = [max(VesselData(:,3)) max(VesselData(:,6))];
zmax = max(zmax) + 2*offset; 

prompt = {'In some data sets, vessel maps are confined to a single plane. Would you like to add a vertical off-set (in microns) to this? Enter value:'};
dlg_title = 'Add vertical offset?';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

voff = str2double(answer{1}); 



%Now we set out grid limits in selected micron spacing! 
xmax = round(xmax/res); 
ymax = round(ymax/res); 
zmax = round(zmax/res); 
zmax = round(zmax + voff/res); 
%%%%%%%%%%%%%%%%%%%%




%This part normalises the O2 partial pressure at the wall for small vessel
%sections to avoid over-estimation. Full details here : http://rsif.royalsocietypublishing.org/content/13/116/20160070

Cons = asl*om/(3*D);
rt = VesselData(:,7).*10^-6;
zm = (VesselData(:,8)/2).*10^-6;

Pnorms = Cons.*(zm.^3/3 + 2.*rn^3.*log((sqrt(zm.^2 + rt.^2) + zm)./(rt)) + (rt.^2 - 3*rn^2).*zm);

%This part finds the surface normal in 3D. Full details in supplementary material here : http://rsif.royalsocietypublishing.org/content/13/116/20160070
% VECTORIZED DW 08/04/16

RT = rt.*10^6;

A = VesselData(:,1) - VesselData(:,4);
B = VesselData(:,2) - VesselData(:,5);
C = VesselData(:,3) - VesselData(:,6);

theta = acos(C./(sqrt(A.^2+B.^2+C.^2)));
phi = atan2(B,A);

Cx = (VesselData(:,1) + VesselData(:,4))/2;
Cy = (VesselData(:,2) + VesselData(:,5))/2;
Cz = (VesselData(:,3) + VesselData(:,6))/2;

dCx = RT.*cos(theta).*cos(phi);
dCy = RT.*cos(theta).*sin(phi);
dCz = - RT.*sin(theta);

% xperp = round(Cx+dCx);
% yperp = round(Cy+dCy);
% zperp = round(Cz+dCz);

ceilabs = @(x) (x>0).*ceil(x) + (x<=0).*floor(x); %ceiling function on magnitude, keep sign
xperp = round(Cx + ceilabs(dCx));
yperp = round(Cy + ceilabs(dCy));
zperp = round(Cz + ceilabs(dCz));




%This part find the surface normal in res micron spacing!
xperp = round(xperp/res); 
yperp = round(yperp/res); 
zperp = round(zperp/res); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%edge limit correction 
% VECTORIZED DW 08/04/16

xperp(xperp<1) = 1;
yperp(yperp<1) = 1;
zperp(zperp<1) = 1;






%Introducing new constants for kernel estimation 

a0 = 1*10^-7;
omega = 3.0318*10^7;
D = 2*10^-9;
imdist = res*10^-6;   %gives pixel length; in this case, one pixel is res micron 

C = asl*omega/(3*D);

%This part scales the oxygen diffusion distances for different size
%vessels.Full details here : http://rsif.royalsocietypublishing.org/content/13/116/20160070
% VECTORIZED DW 08/04/16

rov = VesselData(:,7);
rop = rov*10^-6;
func = @(rop_i,Pnorms_i) fzero(@(rnf) ...
           real(C*((2*rop_i.^2 -8.*rnf.^2).*(sqrt(rnf.^2 - rop_i^2))./3 ... 
           + 2.*(rnf.^3).*log((rnf + sqrt(rnf.^2 - rop_i^2))./rop_i)) ...
           - Pnorms_i), [1 150]*rop_i);
rn = arrayfun(func,rop,Pnorms);
rov = rop;
value = ceil(rn./imdist);
h = 2*value + 1;

%assignin('base', 'valop', value)  %spit to workspace

%This part of the code finds the X Y Z position of kernel elements.Full details here : http://rsif.royalsocietypublishing.org/content/13/116/20160070

X = [VesselData(:,1) VesselData(:,4)];
Y = [VesselData(:,2) VesselData(:,5)];
Z = [VesselData(:,3) VesselData(:,6)];

[YMin,Minvy] = min(Y,[],2);
[YMax,Maxvy] = max(Y,[],2);
XMin = X(sub2ind(size(Y),1:size(Y,1),Minvy'))';
XMax = X(sub2ind(size(Y),1:size(Y,1),Maxvy'))';

angle = atan2(YMax - YMin,XMax - XMin);
angle(Maxvy == Minvy,:) = 0;

Length = VesselData(:,8);
Length = round(Length/2); % 2 micron correction! 

Ry = sqrt((X(:,1)-X(:,2)).^2 + (Y(:,1)-Y(:,2)).^2);
stepsydom = Ry./Length;

Az = arrayfun(@(x,b) b*(1:x),Length,stepsydom,'UniformOutput',false);

Xpair = cell(1,numel(Az));
Ypair = cell(1,numel(Az));
Zpair = cell(1,numel(Az));
for i = 1:numel(Az)
    if angle(i) == 0
        Xpair{i} = round(X(i,Minvy(i)) + Az{i}*cos(angle(i)));
        Ypair{i} = ones(size(Az{i}))*VesselData(i,2);
    else
        Xpair{i} = round(X(i,Minvy(i)) + Az{i}*cos(angle(i)));
        Ypair{i} = round(Y(i,Minvy(i)) + Az{i}*sin(angle(i)));
    end
    Zpair{i} = round(linspace(Z(i,Minvy(i)),Z(i,Maxvy(i)),numel(Az{i})));
end

%resize to res micron! 
vl = length(Xpair);

for hh = 1:vl
Xpair{hh} = round(Xpair{hh}/res); 
Ypair{hh} = round(Ypair{hh}/res);
Zpair{hh} = round(Zpair{hh}/res);
end



%This part adds each kernel, normalises to vessel when complete then moves
%to next vessel.Full details here : http://rsif.royalsocietypublishing.org/content/13/116/20160070
% VECTORIZED DW 21/04/16

finalgrid = zeros(xmax,ymax,zmax);

handle = waitbar(0);

%tic;

for i = 1:Segnum 

    handle = waitbar((i-1)/Segnum,handle,['Processing Segment ' num2str(i) ' of ' num2str(Segnum)]);
    
    [m, n, l] = meshgrid(1:h(i),1:h(i),1:h(i));   % should correct to 2 micron! 
    r = sqrt(((m - value(i)).*(imdist)).^2 + ((n - value(i)).*imdist).^2 + ((l - value(i)).*imdist).^2); 
    matrix = (a0*omega/(6*D)).*(r.^2 + (2*(rn(i).^3))./r - 3.*rn(i).*rn(i));
    matrix(r >= rn(i)) = 0;
    matrix(r <= rov(i)) = (a0*omega/(6*D)).*(rov(i).^2 + (2*(rn(i).^3))./rov(i) - 3.*rn(i).*rn(i)); % matrix is the point kernel
    roNaN = 1e-6*sqrt((Cx(i)-res.*xperp(i)).^2+(Cy(i)-res.*yperp(i)).^2+(Cz(i)-res.*zperp(i)).^2); %res micron correction! 
    matrix(r < (roNaN-2e-6)) = NaN; % centre of vessels to NaN

    
    vessgrid = zeros(xmax,ymax,zmax);

    for j = 1:Length(i)
        posnoffset = [Xpair{i}(j)-value(i)-1 Ypair{i}(j)-value(i)-1 Zpair{i}(j)-value(i)-1];
        %assignin('base', 'posoff', posnoffset)  %spit to workspace
        vessgrid = offsetadd(vessgrid,matrix,posnoffset);
    end
    
    val = vessgrid(xperp(i),yperp(i),zperp(i));
    %fix for NaNs - test now
    if isnan(val) == 1
        val = max(max(max(vessgrid))); 
    end
    
    scale = Pnorms(i)/val;
    finalgrid = finalgrid + scale*vessgrid;


end

%toc;

close(handle);

clear val vessgrid

%Statistics for section

for a = 1 : zmax
    
    meangrid(a) = nanmean(nanmean(finalgrid(:,:,a)));

    stdgrid(a) = nanstd(nanstd(finalgrid(:,:,a)));
    
    anoxia(a) = length(find(finalgrid(:,:,a) <= 0));
    severe(a) = length(find(finalgrid(:,:,a) <= 0.5));
    halfmax(a) = length(find(finalgrid(:,:,a) <= 2.5));
    OEReffect(a) = length(find(finalgrid(:,:,a) <= 20));
end

perc = 100/(xmax*ymax*zmax);

meangrid = nanmean(meangrid);
stdgrid = nanstd(stdgrid); 
anoxia_per = (sum(anoxia))*perc;
severe_per = (sum(severe))*perc;
halfmax_per = (sum(halfmax))*perc;
OEReffect_per = (sum(OEReffect))*perc;
relvol = sum(isnan(finalgrid(:)))/numel(finalgrid);

uiwait(msgbox(['For the section analysed, vessel volume is ' num2str(100*relvol) '%, the mean oxygen tension is ' num2str(meangrid) ' +/- ' num2str(stdgrid) ' mmHg. Of the entire volume, ' num2str(anoxia_per) ' % is completely anoxic, ' num2str(severe_per) ' % is below the oxygen threshold for mitosis (< 0.5 mmHg), ' num2str(halfmax_per) ' % of the volume is below the half-max OER threshold (2.5 mmHg) and in total ' num2str(OEReffect_per) ' of the volume is below OER saturation (taken as < 20 mmHg).'   ],'Result Statistics'));

%assignin('base', 'Grid3D', finalgrid);  %spit to workspace

%Visualisation stuff

z1 = mean(VesselData(:,3));
z2 =  mean(VesselData(:,6));
cent = mean([z1 z2]);
cent = round(cent);
clear z1 z2 



% prompt = {'Enter a vertical cut to visualise:'};
% dlg_title = 'Visualisation';
% num_lines = 1;
% defaultans = {num2str(cent)};
% answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
%vislevel = str2double(answer{1});
vislevel = round(cent/res); 

Slice = finalgrid(:,:,vislevel);
spitout = size(Slice);
xdata = 0:res:res.*(spitout(1)-1); %5 micron spaced
ydata = 0:res:res.*(spitout(2)-1); %5 micron spaced
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
surf(ydata,xdata,Slice); shading interp; colorbar
title(['Vertical slice through vertical (Z =  ' num2str(res*vislevel) ' microns.)'] );
xlabel('X-axis (microns)'); ylabel('Y-axis (microns)'); zlabel('O2 partial pressure (mmHg)')

tt = size(finalgrid); 
zsteps = tt(3); 
zdata = 0:res:((zsteps-1)*(res));
zstring = num2str(zdata);  

% Create pop-up menu
    popup = uicontrol('Style', 'popup',...
           'String', {zdata},...
           'Position', [20 340 100 50],...
           'Callback', @setmap);      
       
       
           function setmap(source,~)
       
        
        % For R2014a and earlier: 
         val = get(source,'Value');
        % maps = get(source,'String'); 
        poson = zdata(val);
        
        close all
        Slice = finalgrid(:,:,val);
spitout = size(Slice);
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
surf(ydata,xdata,Slice); shading interp; colorbar
title(['Vertical slice through vertical (Z =  ' num2str(res*(val-1)) ' microns.)'] );
xlabel('X-axis (microns)'); ylabel('Y-axis (microns)'); zlabel('O2 partial pressure (mmHg)')
          popup = uicontrol('Style', 'popup',...
           'String', {zdata},...
           'Position', [20 340 100 50],...
           'Callback', @setmap);      
        
        

        
           end


jj = size(finalgrid); 
jj1 = jj(1);
jj2 = jj(2);
jj3 = jj(3); 

preamble = ['OX-CALC - Oxygen distribution tool. Version 1.1. Authors: Dr. David Robert Grimes & Dr. Daniel R. Warren. This code allows the estimation of tissue oxygen from known vessel networks from the kernel method outlined here -http://rsif.royalsocietypublishing.org/content/13/116/20160070 - Code presented as is, and users are welcome to edit as desired.']; 
 
queries = ['Queries: E-mail davidrobert.grimes@oncology.ox.ac.uk / grimesd2@gmail.com or @drg1985. Please cite as: Estimating oxygen distribution from vasculature in three-dimensional tumour tissue David Robert Grimes, Pavitra Kannan, Daniel R. Warren, Bostjan Markelc, Russell Bates, Ruth Muschel, Mike Partridge. J. R. Soc. Interface 2016 13 20160070; DOI: 10.1098/rsif.2016.0070.'];

%resolution']


sizeinfo = ['The resolution of this simulation is ' num2str(res) ' micron, and the dimensions are ' num2str(res*(jj1 -1)) 'um x ' num2str(res*(jj2 -1)) 'um x ' num2str(res*(jj3 -1)) 'um'];       
sectioninfo = ['For the section analysed, vessel volume is ' num2str(100*relvol) '%, the mean oxygen tension is ' num2str(meangrid) ' +/- ' num2str(stdgrid) ' mmHg. Of the entire volume, ' num2str(anoxia_per) ' % is completely anoxic, ' num2str(severe_per) ' % is below the oxygen threshold for mitosis (< 0.5 mmHg), ' num2str(halfmax_per) ' % of the volume is below the half-max OER threshold (2.5 mmHg) and in total ' num2str(OEReffect_per) ' of the volume is below OER saturation (taken as < 20 mmHg).'   ];
vessinfo = ['There are ' num2str(Segnum) ' vessel segments in this sample.' ];

fid = fopen('OutputInfo.txt','wt');
fprintf(fid, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', preamble, '',vessinfo,'',sizeinfo, '', sectioninfo, '', queries);
fclose(fid);


save('GridOutput.mat','finalgrid');

%assignin('base', 'Results1', finalgrid)  %spit to workspace
end

function a = offsetadd(a,b,offset)

% a = offsetadd(a,b,offset)
%   Adds two n-dimensional matrices of different sizes with offset. The
%   vector 'offset' specifies the offset of matrix 'b' subscripts with
%   respect to matrix 'a' subscripts. Size of output is specified by matrix
%   a.
%
%   Daniel Warren
%   Department of Oncology, University of Oxford
%   April 2016

ndims = numel(offset);
for j = 1:ndims
	bi = 1:size(b,j);
	bio = bi+offset(j);
	bmask = (bio >= 1 & bio <= size(a,j));

	bi_min(j) = bi(find(bmask,1,'first'));
	bi_max(j) = bi(find(bmask,1,'last'));
	ai_min(j) = bio(find(bmask,1,'first'));
	ai_max(j) = bio(find(bmask,1,'last'));
end

a_inds = arrayfun(@colon,ai_min,ai_max,'UniformOutput',false);
b_inds = arrayfun(@colon,bi_min,bi_max,'UniformOutput',false);

a(a_inds{:}) = a(a_inds{:}) + b(b_inds{:});

end