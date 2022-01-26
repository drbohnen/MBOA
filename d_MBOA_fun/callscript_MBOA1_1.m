% call Script for MBOA Release 1.0  
% DRB (NCSU) and JKH (USC) 
% Original Public Release: 28 June 2011 
% 30 June 2011 - Added option to fill small holes in grid 
% 05 Nov 2011 - Fixed bug in H/R check & updated plotm lines for R2011b
%
disp('..............')
disp('..............')
disp('Modified Basal Outlining Algorithm by DRB (NCSU) & JKH (USC)')
disp('Release MBOA1.1, June 2011')
disp('All Matlab Windows Closed, All Variables Cleared')
disp('..............')
disp('..............')

tic


%% BLOCK 1, for Loading the Grid 
map=strcat(direc,mapgrid); 
[x,y,z]=grdread2(map);z=double(z); x=double(x); y=double(y); 
nx=length(x); ny=length(y); 
if max(x(:)) > 180; x=x-360; end  % note contourm always returns negative west values 

minx=min(x); maxx=max(x); dx=mode(diff(x)); %increment of x
miny=min(y); maxy=max(y); dy=mode(diff(y)); %increment of y
Olat=median(y); Olon=median(x); 
lz=min(z(:));hz=max(z(:));  [lat,lon] = meshgrat([miny maxy],[minx maxx],[ny nx]); 

intp_test=sum(isnan(z(:))); 

if intp_test > 0 && intpflag==1; 
    display('Filling small gaps in the topography...') 
    F=TriScatteredInterp(lon(~isnan(z)),lat(~isnan(z)),z(~isnan(z)));
    z=F(lon,lat);

        if contype == 1 
        direc='./';   mapgrid=strcat('fill_',mapgrid); 
        map=strcat(direc,mapgrid); 
        disp('Writing filled GMT grid to the local directory...') 
        disp('Use this as the input grid on future runs and reset intpflag = 0 ...') 
        grdwrite2mod(x,y,z,mapgrid); 
        end
end


% set tick intervals (uses round2) 

if (maxx-minx)/4 < 0.25 ||  (maxy-miny)/4 < 0.25
   lonlab=round2((maxx-minx)/4,.01); latlab=round2((maxy-miny)/4,.01);  
elseif (maxx-minx)/4 < 0.05 || (maxy-miny)/4 < 0.05 
    lonlab=round2((maxx-minx)/4,.001); latlab=round2((maxy-miny)/4,.001);  
elseif (maxx-minx)/4 < 0.025 || (maxy-miny)/4 < 0.025 
    lonlab=round2((maxx-minx)/4,.0001); latlab=round2((maxy-miny)/4,.0001);
else 
    lonlab=round2((maxx-minx)/4,.25); latlab=round2((maxy-miny)/4,.25); 
end
  
    
% Make initial map 
if plotlevel > 0 
    figure; axesm('MapProjection','mercator','ParallelLabel','on','MeridianLabel','on',...
    'LabelUnits','dm','fontsize',12,'MapLatLimit',[miny maxy],'MapLonLimit',[minx maxx]); 
    setm(gca,'MLabelLocation',lonlab,'PLabelLocation',latlab,'MlineLocation',lonlab,'PLineLocation',latlab); 
    pcolorm(lat,lon,z); view(0,90); caxis([lz,hz]); colorbar('vert'); tightmap;
    v=floor(lz):clstep*4:floor(hz); [C,h]=contourm(lat,lon,z,v, 'k','Color',[0.35 0.5 0.35],'LineStyle','--','LineWidth',.1);
end 

R = makerefmat(lon(1,1), lat(1,1), dx, dy); % refine the referencing matrix 


%% BLOCK to call closed contour routine 
if ptsmin < 8; ptsmin=8; disp('....'); disp('Resetting ptsmin = 8'); disp('....'); end

[ccontour_t, stats_t, cl]=cc_MBOA1_1(lz,hz, clstep,clstart,contype, lat, lon, z,...
dy, dx, Olat, Olon, ptsmin, elp_ratio,Efit,minHR,maxHR, maxdim,direc,mapgrid); 


if ~isempty(ccontour_t);   % only continue there are closed-contours returned. 
%%
%% BLOCK 2, for reformatting  contours 
disp('Reformatting basal contours...') 
stat=[];
for i=1:length(cl); eval([' stat=cat(2,stat, stats_t.S_' num2str(abs(cl(i))) ');']); end 

[~,c]=size(stat); 
Nrows=max(stat(1,:)); 
ccont=nan(Nrows,c,2); % preallocate space  

cc=1; %counter index 
for i=1:length(cl); 
    eval(['acont=ccontour_t.C_' num2str(abs(cl(i))) ';']); 
    [r,c,l]=size(acont); 
        if r > Nrows; acont=acont(1:Nrows,:,:); [r,c,l]=size(acont);  end;  %Might have some NaN rows extending past where Nrows 
    ccont(1:r,cc:cc+c-1,1:2)=acont; 
    cc=cc+c; % skip ahead the proper number of columns 
end 

% reorder largest to smallest 
[J,I]=sortrows(stat',-1); % I is the sort order 
stat=J'; ccont=ccont(:,I,:); 
disp('Done...') 


%% BLOCK 3, for call A/O modification routine  
if AdthresR > 5;  AdthresR=5; fprintf('Warning.... Setting AdthresR to its maximum = 5'); end 
if MaxHdiff > 5;  MaxHdiff=5; fprintf('Warning... Setting MaxHdiff to its maximum = 5'); end 

other_cc=ccont;  
[~,c,~]=size(ccont); 
for j=1:c
    clat=ccont(~isnan(ccont(:,j,2)),j,2); 
    clon=ccont(~isnan(ccont(:,j,2)),j,1);
    other_cc(:,j,:)=0; % null the current closed contour for inpolygon search in modccbase 
    CMOD=modcc_MBOA1_1(clat,clon,R,z,lat,lon,dx,dy,other_cc,MaxHdiff,AdthresR,SlopePortion,SlopeRatioThres,Naz,Rstep,FT,FL,AOcrit,plotlevel,verb); 
    other_cc(:,j,:)=NaN; %change back to NaN
    other_cc(1:length(CMOD),j,:)=CMOD; % replace with CMOD, which will always be smaller in length than the original. 
end

%% BLOCK 4, replace the trailing zeros  
[r,c,l]=size(other_cc); 
for k=1:c 
    a=find(other_cc(:,k,1)==0 & other_cc(:,k,2)==0); % assumes no contour has a point passing through exactly lat=lon=0; 
        if ~isempty(a); other_cc(a,k,:)=NaN; end 
end

numpts=find(~isnan(other_cc(:,1,1))==1, 1, 'last' );  % should have (Naz*2)+1 points 
MBC=other_cc(1:numpts,:,:);  % modified basal contours 
BC=ccont;      % basal closed contours 

%% Block 5, for calling the stats program (which also writes output to ascii) 
foldername=['d_' mapgrid(1:end-4) datestr(now,'yyyymmmdd_HHMM_run') '/'];
mkdir(foldername)
warning('off','all');
[statsMBC,MBCxyz]=stats_MBOA1_1(lat, lon, z, R, Olat,Olon,MBC,plotlevel,writestats,mapgrid,'mbc',foldername,figsavetype);  % do the stats for MBC (Modified closed contour)
[statsBC,BCxyz]=stats_MBOA1_1(lat, lon, z, R, Olat,Olon,BC,plotlevel,writestats,mapgrid,'bc',foldername,figsavetype);  % do the stats for BC (closed contour)
warning('on','all');

if plotlevel > 0 
    fhan=figure(1); 
    plotm(double(statsMBC(:,7)),double(statsMBC(:,8)),'k.','MarkerSize',6)
    setm(gca,'grid','on')
end
% 
if figsavetype==0
    stxt=['print -dpng ' foldername 'map.png']; eval(stxt);  % never save mat figure of map. Too big a file 
elseif  figsavetype==1
    stxt=['print -dpng ' foldername 'map.png']; eval(stxt);     
elseif  figsavetype==2
    stxt=['print -depsc2 ' foldername 'map.ps']; eval(stxt);    
end
%   
   
disp('......')
disp('......')
toc

if writestats==1
    disp('....')
    disp('saving parameter file...'); disp('')
    
txt=['para_'  datestr(now,'yyyymmmdd_HHMM_') mapgrid(1:end-3) 'txt'];
disp(txt)
fdP=fopen(strcat(foldername,txt),'w');
fprintf(fdP,'%% date & time run end %s\n%% mapgrid %s\n%% direc %s \n', datestr(now,'yyyymmmdd_HHMM'), num2str(mapgrid), num2str(direc));
fprintf(fdP,'%% clstep %s\n%% maxdim %s\n%% elp_ratio %s \n', num2str(clstep), num2str(maxdim), num2str(elp_ratio));
fprintf(fdP,'%% Efit %s\n%% minHR %s\n%% maxHRo %s \n', num2str(Efit), num2str(minHR), num2str(maxHR));
fprintf(fdP,'%% Naz %s\n%% AOcrit %s\n%% AdthresR %s Rstep %s\n', num2str(Naz), num2str(AOcrit), num2str(AdthresR), num2str(Rstep));
fprintf(fdP,'%% MaxHdiff %s\n%% Slope Portion %s\n%% Slope Ratio Thres %s \n', num2str(MaxHdiff), num2str(SlopePortion), num2str(SlopeRatioThres));
fprintf(fdP,'%% Filter Type %s\n%% Filter Length %s\n', num2str(FT), num2str(FL));
fclose(fdP);
disp(['your files are in directory...' foldername ''])
end 

disp('thanks'); disp('') 
finsound 

end % end for the only continue if there are closed contour statement 
