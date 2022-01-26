function [statsMBC,MBCxyz]=stats_MBOA1_1(lat, lon, z, R, Olat,Olon,MBC,plotlevel,writestats,mapgrid,fileext,foldername,figsavetype)
% lat lon z are the grid inputs 
% R is the reference matrix 
% Olat & Olon are the origin used in projecting to equal area 
% MBC is the contour file in 3-D matrix form with NaN padding 
% plotlevel == 2 will make a 3-D plot of each mound 
% writestats==1 write a stats file and contour base file with the results  
% mapgrid is the name of the grid file for naming the output stats 
% fileext is the file extension used to denote bc or mbc contours 
% foldername is the folder to write stats and images 
% filesavetype is the image file type (0 matlab fig, 1 png, 2 eps-color) 
% 
% DRB NCSU March 26, 2011 
% Version 1.1 June 2011 

%% BLOCK 1 Project the grid
[XX,YY]=grn2eqa(lat,lon,[Olat, Olon],almanac('earth','wgs84','meters')); %project to equal area

%% BLOCK 2  -- preallocate 
[r,c,~]=size(MBC); %number of rows of largest contour, number of contours, dimension
statsMBC=NaN(11,c);
MBCxyz=nan(r,c,3);  MBCxyz(:,:,1:2)=MBC; 


%% BLOCK 3 -- LOOK THRU ALL THE CONTOURS 
for i=1:c  % each column in the basal matrix 
    mbclat=MBC(:,i,2);mbclat=mbclat(~isnan(mbclat));  % basal lats 
    mbclon=MBC(:,i,1);mbclon=mbclon(~isnan(mbclon));  % basal lons 
    mbcz=ltln2val(z,R,mbclat,mbclon,'bilinear'); 
    if var(mbcz) < 1; mbcz=ones(length(mbcz),1).*round(mean(mbcz)); end  % assume traditional closed contour. Remove small variability
    MBCxyz(1:length(mbcz),i,3)=mbcz;  % output the contour elevation 
  
    %project base to equal area
    [mbcx,mbcy]=grn2eqa(mbclat,mbclon,[Olat, Olon],almanac('earth','wgs84','meters')); %project base lat lon to equal area
    mbc_per=sum(sqrt(diff(mbcx).^2 + diff(mbcy).^2)); 
    G=fit_ellipse(mbcx,mbcy); %fit the ellipse in x-y space to the active base
    mbc_area = polyarea(mbcx,mbcy);  % find area of active base         
            
    if isempty(G.a); statsMBC(7:11,i)=[NaN,NaN];
    else
        [latelp,lonelp]=eqa2grn(G.X0_in,G.Y0_in,[Olat, Olon],almanac('earth','wgs84','meters'));
        statsMBC(7:11,i)=[latelp,lonelp, G.long_axis/1000,G.short_axis/1000,rad2deg(G.phi)]; %write long and short axis to stats file
    end
     
    %find points inside ellipse
    IN=inpolygon(XX,YY,mbcx,mbcy); %find points inside base 
    Zinside=z(IN);                 %get z values inside base 
    [cent_lat,cent_lon]=eqa2grn(XX(IN),YY(IN),[Olat, Olon],almanac('earth','wgs84','meters'));  % project in meters 
    statsMBC(1:2,i)=meanm(cent_lat,cent_lon); 
    %
    Fbase = TriScatteredInterp(mbcx(1:end-1),mbcy(1:end-1),mbcz(1:end-1)); % fits a plane through base points 
    temp=Fbase(XX(IN),YY(IN)); % predicted z values on the bottom area of mound 
    ht=Zinside-temp;  hpos=ht(ht>0);
    mbc_vol=mbc_area*mean(hpos); 
    statsMBC(3:6,i)=[max(ht)/1000, mbc_per/1000, mbc_area/(1000^2),mbc_vol/(1000^3)];
    %
    if plotlevel>=2  && strcmp('mbc',fileext)
        fhan=figure; % figure handle for saving  
        Zingrid=nan(size(XX)); Zingrid(IN)=Zinside; surf(XX,YY,Zingrid,'LineStyle','none') % plot the surface 
        hold on; plot3(mbcx,mbcy,mbcz,'.r','MarkerSize',12); 
        dim=1.25*max([max(mbcx)-min(mbcx),max(mbcy)-min(mbcy)]);  % this is a pad for plot axis 
        axis([min(mbcx)-(0.125*dim),dim+(min(mbcx)-0.125*dim),min(mbcy)-(0.125*dim),dim+(min(mbcy)-0.125*dim)]); 
        [la]=degrees2dms(statsMBC(1,i)); [lo]=degrees2dms(statsMBC(2,i));
        
        ttxt=['title(''smt # ' num2str(i) ' lat ' num2str(la(1)) '  ' num2str(la(2)) '  ' num2str(round2(la(3),.1)) '      lon ' num2str(lo(1)) '  ' num2str(lo(2)) '  ' num2str(round2(lo(3),.1)) ''')'];
         eval(ttxt);                      % title text 
        stxt=['saveas(fhan,''' foldername 'smt' num2str(i) '.fig'')']; eval(stxt); % save text 
        close(fhan) 
    end
  
    % if you are plotting the closed contour bases, plot them ontop of the
    % modified bases 
    if plotlevel>=2  && strcmp('bc',fileext);
        ltxt=['fhan=open(''' foldername 'smt' num2str(i) '.fig'');']; eval(ltxt); 
        plot3(mbcx,mbcy,mbcz,'.w','MarkerSize',12); 
        h=legend('3-D surface','modified base','closed-contour case'); set(h,'Color',[.95 .95 .5]) 
    
        % figure save 
        if figsavetype==0
            stxt=['saveas(fhan,''' foldername 'smt' num2str(i) '.fig'')']; eval(stxt); 
            close(fhan)
        elseif  figsavetype==1  
            stxt=['print -dpng ' foldername 'smt' num2str(i) '.png']; eval(stxt);    
            dtxt=['delete(''' foldername 'smt' num2str(i) '.fig'')']; eval(dtxt);  
            close(fhan) 
        elseif  figsavetype==2
            stxt=['print -depsc2 ' foldername 'smt' num2str(i) '.ps']; eval(stxt);    
            dtxt=['delete(''' foldername 'smt' num2str(i) '.fig'')']; eval(dtxt);
            close(fhan)
        end
   end % end of figure save 
    
end

statsMBC=statsMBC'; 

%% BLOCK 4 - write stats files if wanted 
if writestats==1
    disp('....')
    disp('writing stats file...') 
    txt=['stats_'  datestr(now,'yyyymmmdd_HHMM_') mapgrid(1:end-3) fileext ];
    disp(txt)
    fdS=fopen(strcat(foldername,txt),'w');
    fprintf(fdS,'%%latcent loncent maxht(km) perim(km) area(km^2)  vol(km^3) latcentelps loncentelps lonaxis(km) shortaxis(km) elporient(deg_cc_from_east)\n');
    fprintf(fdS,'%02.6f %03.7f %1.3f %1.3f %1.6f %1.9f %02.6f %03.7f %1.3f %1.3f %1.2f\n', statsMBC');
    fclose(fdS);
end 

disp('writing ascii xyz bases file...') 
txt=['xyz_bases'  datestr(now,'yyyymmmdd_HHMM_') mapgrid(1:end-3) fileext ];
disp(txt)
fdC=fopen(strcat(foldername,txt),'w');
fprintf(fdC,'NaN NaN NaN\n');

for j=1:c
    tempz1=MBCxyz(:,j,3);tempz=tempz1(~isnan(tempz1));      % basal depths 
    templat=MBCxyz(:,j,2);templat=templat(~isnan(tempz1));  % basal lats 
    templon=MBCxyz(:,j,1);templon=templon(~isnan(tempz1));  % basal lons 
    tempout=[templon,templat,tempz]; 
    fprintf(fdS,'%1.7f %1.7f %1.2f\n', tempout');
    fprintf(fdC,'NaN NaN NaN\n');
end 
fclose(fdC);


disp('....')
disp('You may safely ignore duplicate data point warnings')




