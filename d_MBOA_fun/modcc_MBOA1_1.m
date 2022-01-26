function [COUT2]=modcc_MBOA1_1(clat,clon,R,z,lat,lon,dx,dy,other_cc,MaxHdiff,AdthresR,SlopePortion,SlopeRatioThres,Naz,Rstep,FT,FL,AOcrit,plotlevel,verb)
% INPUTS
% clat and clon are the contour lat and lon 
% R and z are the parameters and grid data  
% dx and dy are the grid dimensions 
% other_cc is all other closed contours 
% MaxHdiff - maximum percentage difference between height forward and reverse 
% AdthresR is the max percent increase in range along a profile 
% SlopePortion - percentage of the profile length to calculate slope at the base 
% SlopeRatioThres - ratio drop in slope at the edge 
% Naz - number of azimuthally distributed toographic profiles in a half circle 
% Rstep = range step in m to advance the topo profiles 
% FT,FL are the filter type and length 
% AOcrit check this criter is if = 1 
% plotlevel = plot level and verb = verbose setting 
%
% Modified Closed Contour v1.1 
% DRB (NCSU) and JKH (USC)
% Last updated 28 June 2011 
% drbohnen@ncsu.edu 

%  make a new figure for plotting the profiles of each seamount in a separate window. 
if plotlevel==3; figure; end 
 
% find peak within the closed contour base 
InPts=inpolygon(lon, lat,clon,clat);  % inpolygon works on x-y (lon-lat) inputs  
zsub=z(InPts); latsub=lat(InPts); lonsub=lon(InPts);  % elevations, lats & lons inside the contour 
[~,I]=max(zsub);  lonmaxz=lonsub(I); latmaxz=latsub(I); % I is position of the peak inside the contour 
dtemp=distance(clat,clon,latmaxz,lonmaxz); 
Pkslope=median((max(zsub)-min(zsub))./...
    (1000*deg2km(dtemp))); 
Pkslope=atand(Pkslope);  % median slope in degrees between peak and the base 


% find maximum range and number of points on the profile 

maxd=max(dtemp); mind=min(dtemp); % max & min distance from peak to base 
cellsize=mean([dx,dy]); % avg grid resolution 
rngdeg=min([(1+AdthresR)*(maxd+2*cellsize),5*(maxd+2*cellsize)]) ; % max range to extend in degrees, can be no more than 5 times

Rstep=km2deg(Rstep/1000); 
if Rstep < cellsize;  Rstep=cellsize;
    fprintf('Rstep too small, resetting Rstep to %s meters (avg grid size)\n', num2str(floor(deg2km(Rstep)*1000))); end 
if Rstep >10*cellsize;  Rstep=10*cellsize;
    fprintf('Rstep too large, resetting to %s meters (10 * avg grid size)\n', num2str(floor(deg2km(Rstep)*1000))); end 

if mind < 3*Rstep
    cellsize=mind*.3;  fprintf('reducing Rstep for this small edifice\n'); 
    else
    cellsize=Rstep; 
end

% parameters for extracting profiles 
N=floor(rngdeg/cellsize); % number of points in profile
Daz=180/Naz;  % angle between profile lines 
az_pos=0:Daz:180-Daz;  % azimuth of profiles (forward direction) [0, 90, 150] % 
COUT1=zeros(Naz,1,2); COUT2=zeros(Naz,1,2);  % preallocate 

for a=1:length(az_pos) ;  %5 look through each azimuth (Naz is set in input script)
    az=az_pos(a); % set current az to profile 
    rpeatflag=1;  % dummy variable 
    % get profiles forward and reverse along a given azimuth 
    [plat1,plon1]=track1(latmaxz,lonmaxz,az,rngdeg,[],'degrees',N);  % find profile forward sample points 
    plat1=flipud(plat1); plon1=flipud(plon1);  z1=ltln2val(z,R,plat1,plon1);  % flip from edge to peak orientation 
    plat1=plat1(~isnan(z1)); plon1=plon1(~isnan(z1)); z1=z1(~isnan(z1)); N1=length(z1); % eliminate NaN bathy points 
    [plat2,plon2]=track1(latmaxz,lonmaxz,az,-1*rngdeg,[],'degrees',N); % find profile reverse sample points 
    plat2(1)=[]; plon2(1)=[]; z2=ltln2val(z,R,plat2,plon2);
    plat2=plat2(~isnan(z2)); plon2=plon2(~isnan(z2)); z2=z2(~isnan(z2));N2=length(z2); % eliminate NaN bathy points 
    % find intersection with present basal closed contour in both directions 
    [~,~,I1]=polyxpoly(plon1,plat1,clon,clat); % find forward intersection, plon1 and plat1 are forward profiles  
    if isempty(I1); I1=1; I1_ind =1;   % happens if the contour is within one step of the NaN data edge. 
    else[~,I1_ind]=max(distance(plat1(I1(:,1)),plon1(I1(:,1)),latmaxz,lonmaxz)); % take most distant intersection 
        I1=I1(I1_ind,1)+1; %move to the inside of the contour base to start
    end  
    %
    [~,~,I2]=polyxpoly(plon2,plat2,clon,clat); %find reverse intersection, plon2 and plat2 are reverse profiles  
    if isempty(I2); I2=length(plon2);I2_ind=1; 
    else [~,I2_ind]=max(distance(plat2(I2(:,1)),plon2(I2(:,1)),latmaxz,lonmaxz));% take most distant intersection 
    end % happens if the contour is with one step of the NaN data edge. 
    I2=I2(I2_ind,1)+N1;  % add the offset from the forward profile 
    %
    plat=cat(1,plat1,plat2); plon=cat(1,plon1,plon2); zi=cat(1,z1,z2); % combine the forward and reverse profile passing through peak 
    ri=deg2km(distance(plat(1),plon(1),plat,plon))*1000; % calculate the distane along the profile in meters 
    diffr=mode(diff(ri)); %range step distance in meters 
    %
    if plotlevel==3  % plot ont the profile and position of basal contour  
        plot(ri, zi); hold on; plot(ri(I1),zi(I1),'ro'); plot(ri(I2),zi(I2),'rx'); 
        %figure(1); hold on; plotm(double(plat),double(plon),'k'); 
    end


    %% FORWARD END OF PROFILE (first pass) %%% 
    I1p=I1:-1:1;             %possible values, I1 is the index of the base to start with
    I1new=max([I1-1,1]);     % set as just outside the base initially 
    InitPkhght=zi(N1)-zi(I1);   % height of seamount in m, N1 is index @ profile center  
    Initrad1=ri(N1)-ri(I1);  % Initial Radius in m 

    for i=1:length(I1p) % each decreasing index 
        r_for=ri(I1p(i):I2); z_for=zi(I1p(i):I2); % get the ranges and elevations along the profile
        y_for=interp1([r_for(1),r_for(end)],[z_for(1),z_for(end)],r_for); % estimate the elevation of the bottom 
        h_for=z_for-y_for; A_for(i)=trapz(r_for,h_for);  % find height of seamount along profile and the crossectional area of the profile 
        %
        if AOcrit < 2; O_for(i)=sum(sqrt(diff(r_for).^2+diff(z_for).^2)); 
        else %+(r_for(end)-r_for(1)); % calculate the total outline "o" length of the profile
            O_for(i)=sum(sqrt(diff(r_for).^2+diff(z_for).^2))+(r_for(end)-r_for(1));  
        end 
        AO_for(i)=A_for(i)/O_for(i); % ratio of area to profile outline
        Curhght1=zi(N1)-zi(I1p(i));
        Currad1=ri(N1)-ri(I1p(i)); %Current height and radius
        SlpPts1=max([3,ceil(SlopePortion*(Currad1/diffr))])-1; % number of points to calculate slope over at base (min 2 points for slope calc)

        if N1-I1p(i) >= SlpPts1 
            BaseSlope1=slopefit(ri(I1p(i):I1p(i)+SlpPts1),zi(I1p(i):I1p(i)+SlpPts1));
        else 
            BaseSlope1=Pkslope; % force it to pass the slope test 
        end

        if i> 1  % check for the second point 
            if AO_for(i) - AO_for(i-1) < 0 && AOcrit > 0 ||...   % check for A/O increase 
                Curhght1 > (InitPkhght*MaxHdiff)+InitPkhght  ||... % check percentage deepening relative to initial peak  
                Currad1 > (Initrad1*AdthresR)+Initrad1 ||... % check if radius or depth increased by more than allowed 
                BaseSlope1/Pkslope < SlopeRatioThres         % check if the slope (degrees) fails to some fraction of the original 
                I1new=I1p(i); % keep current step forward if it fails 
                %
                if verb==1% 
                    if AO_for(i)- AO_for(i-1) < 0 && AOcrit > 0; fprintf('Forward profile %1.0f at step %1.0f :',a, i);disp('failed Ao'); 
                    else 
                        if BaseSlope1/Pkslope < SlopeRatioThres;fprintf('Forward profile %1.0f at step %1.0f :',a, i); disp('failed for base slope'); end 
                        if Currad1 > (Initrad1*AdthresR)+Initrad1;fprintf('Forward profile %1.0f at step %1.0f :',a, i); disp(['Radius increase too great, R = ' num2str(Currad1,'%1.1f') 'Initial R= ' num2str(Initrad1,'%1.1f') ]); end 
                        if Curhght1 > (InitPkhght*MaxHdiff)+InitPkhght; fprintf('Forward profile %1.0f at step %1.0f :',a, i);disp(['Height difference too large, H= ' num2str(Curhght1,'%1.1f') 'Initial H = ' num2str(InitPkhght,'%1.1f') ]) ; end 
                    end
                end % verbose 
            break
            else % test to make sure the limits don't enter another seamount 
                [ sumpoly] = ptinside(plon(I1p(i)),plat(I1p(i)),other_cc(:,:,1),other_cc(:,:,2)); 
                if sumpoly > 0 
                    I1new=I1p(i - 1); % back up outside if it fails 
                    if verb==1 
                    fprintf('Forward profile %1.0f at step %1.0f :',a, i);disp('HIT ANOTHER MOUND - YOU SHALL NOT PASS!'); end 
                break 
                end  
      
            end  % end check for all slope criteria 
        end    % end i > 1 
    end    % end of for 
    clear A_for O_for AO_for 


    %% REVERSE END OF PROFILE %%% 
    I2p=I2:1:N1+N2;           % Possible values 
    I2new=min([I2+1,N1+N2]);  % Set as the current value initially 
    Initrad2=ri(I2)-ri(N1);   % Initial Radius in m 

    for j=1:length(I2p);  % for each possible value 
        r_rev=ri(I1new:I2p(j)); z_rev=zi(I1new:I2p(j)); % get the ranges and elevations along the profile
        y_rev=interp1([r_rev(1),r_rev(end)],[z_rev(1),z_rev(end)],r_rev); % estimate the elevation of the bottom 
        h_rev=z_rev-y_rev; A_rev(j)=trapz(r_rev,h_rev);  % find height of seamount along profile and the area of the seamount
        if AOcrit < 2; O_rev(j)=sum(sqrt(diff(r_rev).^2+diff(z_rev).^2)); 
        else %+(r_rev(end)-r_rev(1)); % calculate the total outline "o" length of the profile
            O_rev(j)=sum(sqrt(diff(r_rev).^2+diff(z_rev).^2))+(r_rev(end)-r_rev(1)); 
        end 

        AO_rev(j)=A_rev(j)/O_rev(j); % ratio of area to profile outline
        Curhght2=zi(N1)-zi(I2p(j));
        Currad2=ri(I2p(j))-ri(N1); %Current height and radius
        SlpPts2=max([3,ceil(SlopePortion*(Currad2/diffr))]) -1 ;  % number of points to calculate slope over at base 

        if I2p(j)-N1 >= SlpPts2 % check if there are enough points in the 1/2 profile  
            BaseSlope2=-1*slopefit(ri(I2p(j)-SlpPts2:I2p(j)),zi(I2p(j)-SlpPts2:I2p(j)));
        else 
            BaseSlope2=Pkslope;  %  force it to pass the slope test 
        end

        Hdiff=abs(Curhght2-Curhght1)/max([Curhght2,Curhght1]);      % difference in height 

        if j> 1  % check for the second point 
            if AO_rev(j)- AO_rev(j-1) < 0 && AOcrit > 0||... 
                Currad2 > (Initrad2*AdthresR)+Initrad2 ||... % check if radius or depth of increased by more than allowed 
                BaseSlope2/Pkslope < SlopeRatioThres ||...         % check if the slope (degrees) fails to some fraction of the original         
                Hdiff > MaxHdiff  % check that one side isn't too much deeper than the other 
                I2new=I2p(j);  % keep current step forward if fail occurs 
                %
                if verb==1% 
                    if AO_rev(j)- AO_rev(j-1) < 0 && AOcrit > 0 ; fprintf('Reverse profile %1.0f at step %1.0f :',a, j); disp('failed Ao'); 
                    else 
                        if BaseSlope2/Pkslope < SlopeRatioThres; fprintf('Reverse profile %1.0f at step %1.0f :',a, j); disp('failed for base slope'); end 
                        if Currad2 > (Initrad2*AdthresR)+Initrad2; fprintf('Reverse profile %1.0f at step %1.0f :',a, j); disp(['Radius increase too great, R = ' num2str(Currad2,'%1.1f') 'Initial R= ' num2str(Initrad2,'%1.1f') ]); end 
                        if Hdiff > MaxHdiff; fprintf('Reverse profile %1.0f at step %1.0f :',a, j); disp('MaxHdiff failure, one side too deep relative to the other'); end 
                    end
                end % verbose 
  
            break
            else 
                [ sumpoly] = ptinside(plon(I2p(j)),plat(I2p(j)),other_cc(:,:,1),other_cc(:,:,2)); 
                if sumpoly > 0 
                    I2new=I2p(j-1);  % back up if fail occurs 
                    if verb==1 
                    fprintf('Reverse profile %1.0f at step %1.0f :',a, j);disp('HIT ANOTHER MOUND - YOU SHALL NOT PASS!'); end 
                break 
                end 
            % 
            end  % end check increasing AO_for 
        end  % end i > 1 
    end
    clear A_rev O_rev AO_rev i j 


    %% now try moving the first one again to extend deeper, if flag is set 
    if rpeatflag==1
        I1p=I1new:-1:1;  % possible values 

        for i=1:length(I1p) 
            r_for=ri(I1p(i):I2new); z_for=zi(I1p(i):I2new); % get the ranges and elevations along the profile
            y_for=interp1([r_for(1),r_for(end)],[z_for(1),z_for(end)],r_for); % estimate the elevation of the bottom 
            h_for=z_for-y_for; A_for(i)=trapz(r_for,h_for);  % find height of seamount along profile and the area of the seamount
            if AOcrit < 2; O_for(i)=sum(sqrt(diff(r_for).^2+diff(z_for).^2)); 
            else %+(r_for(end)-r_for(1)); % calculate the total outline "o" length of the profile
                O_for(i)=sum(sqrt(diff(r_for).^2+diff(z_for).^2))+(r_for(end)-r_for(1)); 
            end 
            AO_for(i)=A_for(i)/O_for(i); % ratio of area to profile outline
            Curhght1=zi(N1)-zi(I1p(i));
            Currad1=ri(N1)-ri(I1p(i)); %Current height and radius
            SlpPts1=max([3,ceil(SlopePortion*(Currad1/diffr))])-1;

            if I2new-I1p(i) >= SlpPts1 
                BaseSlope1=slopefit(ri(I1p(i):I1p(i)+SlpPts1),zi(I1p(i):I1p(i)+SlpPts1));
            else 
            BaseSlope1=Pkslope; % force to pass slope criteria 
            end

            Hdiff=abs(Curhght2-Curhght1)/max([Curhght2,Curhght1]);      % difference in height 

            if i> 1  % check for the second point 
                if AO_for(i)- AO_for(i-1) < 0 && AOcrit > 0  || ...  % check A/o increases and not going up hill  
                    Currad1 > (Initrad1*AdthresR)+Initrad1 ||... % check if radius or depth of increased by more than allowed 
                    BaseSlope1/Pkslope < SlopeRatioThres ||...         % check if the slope (degrees) fails to some fraction of the original 
                    Hdiff > MaxHdiff  % check that one side isn't too much deeper than the other 

                    I1new=I1p(i-1); % Don't automatically advance one if this fails. 
     
                    if verb==1% 
                        if AO_for(i)- AO_for(i-1) < 0 && AOcrit > 0 ; fprintf('Pass 2 forward profile %1.0f at step %1.0f :',a, i); disp('failed Ao'); 
                        else 
                            if BaseSlope1/Pkslope < SlopeRatioThres; fprintf('Pass 2 forward profile %1.0f at step %1.0f : ',a, i); disp('failed for base slope'); end 
                            if Currad1 > (Initrad1*AdthresR)+Initrad1; fprintf('Pass 2 forward profile %1.0f at step %1.0f : ',a, i); disp(['Radius increase too great, R = ' num2str(Currad1,'%1.1f') 'Initial R= ' num2str(Initrad1,'%1.1f') ]); end 
                            if Hdiff > MaxHdiff; fprintf('Pass 2 forward profile %1.0f at step %1.0f : ',a, i);disp('MaxHdiff failure, one side too deep relative to the other'); end 
                        end
                    end  % end verbose 
  
                break
                else       % test to make sure the limits don't inter another seamount 
                    [sumpoly] = ptinside(plon(I1p(i)),plat(I1p(i)),other_cc(:,:,1),other_cc(:,:,2)); 
                    if sumpoly > 0 
                        I1new=I1p(i - 1); % back up outside if it fails 
                        if verb==1 
                        fprintf('2nd Pass Forward profile %1.0f at step %1.0f :',a, i);disp('HIT ANOTHER MOUND - YOU SHALL NOT PASS!');end  
                    break 
                    end  
      
                end  % end check increasing AO_for 
            end  % end i > 1 
        end    % goes with the for statement 
    clear A_for O_for AO_for
    end % this end goes with the rpeatflag 


    %% Now try repeating the reverse profile again  
    I2p=I2new:1:N1+N2;           % Possible values 

    for j=1:length(I2p);  % for each possible value 
        r_rev=ri(I1new:I2p(j)); z_rev=zi(I1new:I2p(j)); % get the ranges and elevations along the profile
        y_rev=interp1([r_rev(1),r_rev(end)],[z_rev(1),z_rev(end)],r_rev); % estimate the elevation of the bottom 
        h_rev=z_rev-y_rev; A_rev(j)=trapz(r_rev,h_rev);  % find height of seamount along profile and the area of the seamount
        if AOcrit < 2; O_rev(j)=sum(sqrt(diff(r_rev).^2+diff(z_rev).^2)); 
        else %+(r_rev(end)-r_rev(1)); % calculate the total outline "o" length of the profile
            O_rev(j)=sum(sqrt(diff(r_rev).^2+diff(z_rev).^2))+(r_rev(end)-r_rev(1)); 
        end 
        AO_rev(j)=A_rev(j)/O_rev(j); % ratio of area to profile outline
        Curhght2=zi(N1)-zi(I2p(j));
        Currad2=ri(I2p(j))-ri(N1); %Current height and radius
        SlpPts2=max([3,ceil(SlopePortion*(Currad2/diffr))]) -1 ;  % number of points to calc slope over at base 
        %
        if I2p(j)-N1 >= SlpPts2 % check if there are enough points in the 1/2 profile  
            BaseSlope2=-1*slopefit(ri(I2p(j)-SlpPts2:I2p(j)),zi(I2p(j)-SlpPts2:I2p(j)));
        else 
            BaseSlope2=Pkslope;  %  force it to pass the slope test 
        end
        %
        Hdiff=abs(Curhght2-Curhght1)/max([Curhght2,Curhght1]);      % difference in height 
        if j> 1  % check for the second point 
            if AO_rev(j)- AO_rev(j-1) < 0 && AOcrit > 0||... 
                Currad2 > (Initrad2*AdthresR)+Initrad2 ||... % check if radius or depth of increased by more than allowed 
                BaseSlope2/Pkslope < SlopeRatioThres ||...         % check if the slope (degrees) fails to some fraction of the original         
                Hdiff > MaxHdiff  % check that one side isn't too much deeper than the other 
                I2new=I2p(j-1);  % keep current step forward if fail occurs 
                if verb==1% 
                    if AO_rev(j)- AO_rev(j-1) < 0 && AOcrit > 0 ; fprintf('Pass 2 Reverse profile %1.0f at step %1.0f :',a, j); disp('failed Ao'); 
                    else 
                        if BaseSlope2/Pkslope < SlopeRatioThres; fprintf('Pass 2 Reverse profile %1.0f at step %1.0f :',a, j); disp('failed for base slope'); end 
                        if Currad2 > (Initrad2*AdthresR)+Initrad2; fprintf('Pass 2 Reverse profile %1.0f at step %1.0f :',a, j); disp(['Radius increase too great, R = ' num2str(Currad2,'%1.1f') 'Initial R= ' num2str(Initrad2,'%1.1f') ]); end 
                        if Hdiff > MaxHdiff; fprintf('Pass 2 Reverse profile %1.0f at step %1.0f :',a, j); disp('MaxHdiff failure, one side too deep relative to the other'); end 
                    end
                end % verbose 
  
            break
            else 
                [sumpoly] = ptinside(plon(I2p(j)),plat(I2p(j)),other_cc(:,:,1),other_cc(:,:,2)); 
                if sumpoly > 0 
                    I2new=I2p(j-1);  % back up if fail occurs 
                    if verb==1 
                    fprintf('Pass 2 Reverse profile %1.0f at step %1.0f :',a, j);disp('HIT ANOTHER MOUND - YOU SHALL NOT PASS!'); end 
                break 
                end 
        % 
            end  % end check increasing AO_for 
        end  % end i > 1 
    end
    clear A_rev O_rev AO_rev i j 

    %% TRY MOVING BOTH ENDS OF THE PROFILE 
    I1p=I1new:-1:1;  % possible values 
    I2p=I2new:1:N1+N2;  % possible values 
    slen=min(length(I1p),length(I2p)); % possible number steps 
    InitDiam=Initrad2+Initrad1; %  Initial Diameter  

    for i=1:1:slen % for each step 
        if slen==1; break; end % if either point is already at the end, then nothing to do  
        r_for=ri(I1p(i):I2p(i)); z_for=zi(I1p(i):I2p(i)); % get the ranges and elevations along the profile 
        y_for=interp1([r_for(1),r_for(end)],[z_for(1),z_for(end)],r_for); % estimate the elevation of the bottom 
        h_for=z_for-y_for; A_for(i)=trapz(r_for,h_for);  % find height of seamount along profile and the area of the seamount
        if AOcrit < 2 O_for(i)=sum(sqrt(diff(r_for).^2+diff(z_for).^2)); 
        else %+(r_for(end)-r_for(1)); % calculate the total outline "o" length of the profile 
            O_for(i)=sum(sqrt(diff(r_for).^2+diff(z_for).^2))+(r_for(end)-r_for(1)); 
        end 
        AO_for(i)=A_for(i)/O_for(i); % ratio of area to profile outline 
        Curhght1=zi(N1)-zi(I2p(i)); %Current height
        Curhght2=zi(N1)-zi(I2p(i));
        Currad1=ri(N1)-ri(I1p(i)); %
        Currad2=ri(I2p(i))-ri(N1); %
        CurDiam=ri(I2p(i))-ri(I1p(i)); % current diameter now
        SlpPts1=max([3,ceil(SlopePortion*(Currad1/diffr))])-1;
        SlpPts2=max([3,ceil(SlopePortion*(Currad2/diffr))])-1;

        if I2p(i)-N1 >= SlpPts2
            BaseSlope2B=-1*slopefit(ri(I2p(i)-SlpPts2:I2p(i)),zi(I2p(i)-SlpPts2:I2p(i)));
        else 
            BaseSlope2B=Pkslope; 
        end
        if N1-I1p(i) >= SlpPts1 
            BaseSlope1B=slopefit(ri(I1p(i):I1p(i)+SlpPts1),zi(I1p(i):I1p(i)+SlpPts1));
        else 
            BaseSlope1B=Pkslope; 
        end  

        DiffH=abs(Curhght2-Curhght1)/max([Curhght2,Curhght1]);      % difference in height 

        if i> 1  % check for the second point 
            if AO_for(i)- AO_for(i-1) < 0 && AOcrit > 0  ||...
                DiffH > MaxHdiff || ... 
                CurDiam > (InitDiam*AdthresR)+InitDiam ||... 
                BaseSlope1B/Pkslope < SlopeRatioThres || ...
                BaseSlope2B/Pkslope < SlopeRatioThres    
                I1new=I1p(i-1);  %   set back on i if failed
                I2new=I2p(i-1);
                if verb==1% 
                    if AO_for(i)- AO_for(i-1) < 0 && AOcrit > 0 ; fprintf('Moving both ends %1.0f at step %1.0f :',a, i); disp('failed Ao'); 
                    else 
                        if BaseSlope1B/Pkslope < SlopeRatioThres || BaseSlope2B/Pkslope; fprintf('Moving both ends %1.0f at step %1.0f :',a, i); disp('failed for base slope'); end 
                        if CurDiam > (InitDiam*AdthresR)+InitDiam; fprintf('Moving both ends %1.0f at step %1.0f :',a, i); disp('Diameter too big'); end 
                        if Hdiff > MaxHdiff; fprintf('Moving both ends %1.0f at step %1.0f :',a, i); disp('MaxHdiff failure, one side too deep relative to the other'); end 
                    end
                end % verbose 
  
            break
            else         % test to make sure the limits don't inter another seamount 
                [sumpoly1] = ptinside(plon(I1p(i)),plat(I1p(i)),other_cc(:,:,1),other_cc(:,:,2));  % is the forward point now inside another larger on
                [sumpoly2] = ptinside(plon(I2p(i)),plat(I2p(i)),other_cc(:,:,1),other_cc(:,:,2)); % is the reverse end point now inside 
                if sumpoly1+sumpoly2 > 0  % then one of the points has moved within  
                    I1new=I1p(i-1);  %  set back on i if failed
                    I2new=I2p(i-1);
                    if verb==1 
                    fprintf('Moving both %1.0f at step %1.0f :',a, i);disp('HIT ANOTHER MOUND - YOU SHALL NOT PASS!');end  
                break
                end 

            end  % end check increasing AO loop  

        end  % end i > 1 
    end
         

    %%% WRITE OUT THE LAT LON 
    COUT1(a,1,1)=plon(I1new);  COUT1(a,1,2)=plat(I1new);
    COUT2(a,1,1)=plon(I2new);  COUT2(a,1,2)=plat(I2new);
    plot(ri(I1new),zi(I1new),'ko','LineWidth',2); plot(ri(I2new),zi(I2new),'go','LineWidth',2); 

    clear A_for O_for AO_for I1new I2new I1 I2 
end  %% end of the azimuth loop 
%% Filter Results
COUT=cat(1,COUT1,COUT2);  % combine new lat lon picks 
RAD=distance(COUT(:,2),COUT(:,1),latmaxz,lonmaxz); % distance to base in degrees 
RADpad=cat(1,RAD(end-FL+1:end),RAD,RAD(1:FL));  % cat three last points onto the beginning and three first onto end 
if FT==1 %median filter 
    fRAD=medfilt1(RADpad,FL);
    fRAD=fRAD(FL+1:end-FL); % gets rid of edges 
elseif FT==2  % max filte r
    disp('max filter')
    fRAD= minmaxfilt1(RADpad, FL,'max','full');
    kern=ones(1,FL)./FL;  
    fRAD=conv(fRAD,kern); % 3-pt mean filter 
    fRAD=fRAD(FL+ceil(FL/2):end-FL-floor(FL/2)); % gets rid of edges 
    %fRAD=fRAD(FL+ceil(FL/2):end-FL-floor(FL/2)); % gets rid of edges  
elseif FT==3  % max filter
    disp('min filter')
    fRAD= minmaxfilt1(RADpad, FL,'min','full');
    fRAD=fRAD(FL+ceil(FL/2):end-FL-floor(FL/2)); % gets rid of edges    
else % mean filter 
    kern=ones(1,FL)./FL;  
    fRAD=conv(RADpad,kern); % 3-pt mean filter 
    fRAD=fRAD(FL+ceil(FL/2):end-FL-floor(FL/2)); % gets rid of edges  
end

az_pos2=cat(2,az_pos,180+az_pos); az_pos2=az_pos2'; 

for k=1:length(az_pos2)
    [COUT2(k,1,2),COUT2(k,1,1)]=track1(latmaxz,lonmaxz,az_pos2(k),fRAD(k),[],'degrees',1); 
end
COUT2(k+1,1,:)=COUT2(1,1,:);  % make closed 

%[K,~]=convhull(COUT(:,1,1),COUT(:,1,2));  % calculate convex hull of new picks, which is a closed polygon! 
%COUT2=COUT(K,:,:);  % format to return to the seamount bases 

if plotlevel > 0 
    figure(1) % plots on the original map generated in CC_MBOA1_0.m
    plotm(double(COUT(:,1,2)),double(COUT(:,1,1)),'r','LineWidth',1.5)  % the A/O method points 
    plotm(double(COUT2(:,1,2)),double(COUT2(:,1,1)),'k','LineWidth',1.5) % new contour base cc
end
