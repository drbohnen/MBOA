function [ccontour_t, stats_t, cl]=cc_MBOA1_1(lz,hz,clstep,clstart, contype,lat, lon, z, dy, dx, Olat, Olon, ptsmin,elp_ratio,Efit,minHR,maxHR,maxdim,direc,mapgrid)
% INPUTS 
% lz, hz, clstep and clstart are the start, stop and interval for contouring 
% lat, lon, z, dy, dx are grid parameters 
% Olat, Olon are center of the grid 
% ptsmin = minimum points in contour line for consideration
% elp_ratio and Efit are elliptical fit shape and threshold misfit 
% minHR, maxHR, maxdim  limits for slopes (min & max, deg) and dimension(km) 
% 
% Modified Basal Outlining Algorithm MBOA 1.1 
% DRB (NCSU) and JKH (USC) 
% 
% modified on 23 December 2011 to accomodate case where lowest 
% two contours are at the same elevation (Added Block 5B). 
%
% modified on 05 Nov 2011 to fix bug in H/R ratio check for accepting
% a contour and to update plotm lines for Matlab R2011b. 
%
% modified on 01 July 2011 to allow multiple cfile directories to exist
%
% drbohnen@ncsu.edu 

disp('...')
disp('Applying streamlined version of Bohnenstiehl, Howell & Hey''s (G3 2008)') 
disp('closed contour algorithm') 
disp('...')

DEP=departure(min(lon(:)),min(lon(:))+1,mean([min(lat(:)), max(lat(:))]));  % used later for H/B ratio 

tic 
%% BLOCK 1
% GET CLOSED CONTOURS WITH > ptsmin;  STORE THEM IN ARRAYS NAMED BY LEVEL 
disp('Finding all closed contours...') 
cl_1=[];

if contype==1
disp('Using GMT contours...') 
cfiles=cat(2,'cfiles', num2str(floor(now*1000))); 
clvec=round2(lz,clstart):clstep:floor(hz); clvec=clvec';  
mkdir(cfiles); cd(cfiles)
fd=fopen('clvec.txt', 'w'); fprintf(fd,'%1.2f\n', clvec); fclose(fd);  
eval(['[status,result]=system(''grdcontour ../' direc mapgrid ' -Jm -Cclvec.txt -D >junk.ps'');']) 


for i=1:length(clvec) 
[CLmat,STmat,clout]=gmtcontours(clvec(i),ptsmin); 

if ~isempty(clout)
    cl_1=cat(1,cl_1,clout); % keep track of the CL you have now 
    eval(['C_' num2str(abs(clout)) '=CLmat;']) 
    eval(['S_' num2str(abs(clout)) '=STmat;']) 
end
end
cd ../
status=rmdir(cfiles,'s'); 

else 

% matlab contour routine 
disp('Using Matlab contours...')
for cl=round2(lz,clstart):clstep:floor(hz)
    v=[cl cl]; 
    [C]=contourm(lat,lon,z,v,'LineStyle','none');
    a=find(C(1,:)==cl & C(2,:)>=ptsmin); %at the current interval find the contours with more than ptsmin points
    atest=C(2,a+1)-C(2,a+C(2,a))==0 & C(1,a+1)-C(1,a+C(2,a))== 0 ; % equivalent to find, but more efficient 
    a=a(atest);   %get only those that are closed.
    
    if ~isempty(a)  %if there are close contours with > ptsmin points 
        cl_1=cat(1,cl_1,cl); % keep track of the CL you have now 
        
        %make a matrix of NaNs for contours: layer 1 lon,,layer 2 latitude
        eval(['C_' num2str(abs(cl)) '= NaN(max(C(2,a)),length(a), 2,''single'');']); 
        %make a matrix of NaNs for statistics 
        eval(['S_' num2str(abs(cl)) '= NaN(3,length(a),''single'');']);
        for b=1:length(a) %for each contour at cl 
                eval(['C_' num2str(abs(cl)) '(1:C(2,a(b)),b,1)= C(1,a(b)+1:a(b)+C(2,a(b)));']); %place longitude of each point in cc
                eval(['C_' num2str(abs(cl)) '(1:C(2,a(b)),b,2)= C(2,a(b)+1:a(b)+C(2,a(b)));']); %place latitude of each point in cc
                eval(['S_' num2str(abs(cl)) '(1,b)=C(2,a(b)); ']); %place number of points in first row
                %
                eval(['contlats=C_' num2str(abs(cl)) '(:,b,2);']); eval(['contlons=C_' num2str(abs(cl)) '(:,b,1);']); 
                contlats=contlats(~isnan(contlats)); contlons=contlons(~isnan(contlons)); 
               [meanlat,meanlong]=meanm(contlats,contlons); % input order is lat lon for mapping toolbox call 
                eval(['S_' num2str(abs(cl)) '(2:3,b)=[meanlong, meanlat];']); %longitude latitude of first point
        end
    end
end  % end of matlab contour routine 

end  % end of contype 

cl=cl_1; % reset to the remaining levels 

%% BLOCK 2 
% Remove those that don't have > 2 contours inside it. 
cl_2=cl; del_ind=[]; 
for i=1:1:length(cl)-2 
    eval(['Cact=C_' num2str(abs(cl(i))) ';']); %get active C level 
    eval(['Sact=S_' num2str(abs(cl(i))) ';']); %get active S level 
    eval(['Splus2=S_' num2str(abs(cl(i+2))) ';']); %get active S level + 2  
    t1=size(Cact); 

    del_ind=[]; 
    for j=1:t1(2); %for each closed contour at the active level 
        num_ind = sum(inpolygon(Splus2(2,:), Splus2(3,:),Cact(:,j,1), Cact(:,j,2))); 
        if num_ind < 1;
            del_ind=cat(1,del_ind,j); 
        end  
    end
    
    Cact(:,del_ind,:)=[]; Sact(:,del_ind,:)=[]; % deletes closed ones 
    eval(['C_' num2str(abs(cl(i))) '=Cact;']); %reassign C
    eval(['S_' num2str(abs(cl(i))) '=Sact;']); %reassign S  

    if isempty(Cact) 
        cl_2(i)=0; 
    end
    del_ind=[]; 
end
cl=cl_2(cl_2 ~= 0); % reset to the remaining levels  
disp(['Done, Elapsed time is ' num2str(ceil(toc)) ' seconds'])
disp('...')
%% 
%% BLOCK 3 
% now evaluate the shape of the ellipse 
disp('Testing their elliptical shape and maximum dimension...') 
cl_3=cl;  U=[]; del_ind=[];
for i=1:1:length(cl)
    eval(['Cact=C_' num2str(abs(cl(i))) ';']); %get active C level 
    eval(['Sact=S_' num2str(abs(cl(i))) ';']); %get active S level 
    t1=size(Cact); 

    % first layer of U is the x, second layer is the y
    % C stores lon lat, but input into grn2eqa is lat lon 
    [U(:,1:t1(2),1),U(:,1:t1(2),2)]=grn2eqa(Cact(:,:,2),Cact(:,:,1),[Olat, Olon],almanac('earth','wgs84','meters')); 

    for j=1:t1(2); %for each closed contour at the active level  
        
        G=fit_ellipse(U(~isnan(U(:,j,1)),j,1),U(~isnan(U(:,j,2)),j,2)); % fit ellipse 
        if ~isempty(G.a) % if an ellipse can be found 
            % note small efit is a better fit 0 is perfect 
            % save U.mat U Cact Olat Olon 
            e_gof=ellipse_gof(G,U(~isnan(U(:,j,1)),j,1),U(~isnan(U(:,j,2)),j,2)); % return the misfit value 
            MD=max(pdist([U(~isnan(U(:,j,1)),j,1),U(~isnan(U(:,j,2)),j,2)]))/1000; % max distance 
        else 
            e_gof=999999; %if there is no ellipse found 
        end

        if  e_gof > Efit || G.long_axis/G.short_axis > elp_ratio || MD > maxdim ; 
            del_ind=cat(1,del_ind,j); % delete column if fail 
            disp('Deleting contour due to failed ellipse shape criteria or max dimension'); 
        end  
    end
    Cact(:,del_ind,:)=[]; Sact(:,del_ind,:)=[]; % deletes closed ones 

    % now trim extra NaNs to keep the matrix small
    if ~isempty(Cact) 
        Cactmean=nanmean(Cact(:,:,1),2);
        [ro,~,~]=ind2sub(size(Cactmean),find(isnan(Cactmean)));
        if ~isempty(ro); Cact=Cact(1:min(ro)-1,:,:);  end
    end
 
    eval(['C_' num2str(abs(cl(i))) '=Cact;']); %reassign C
    eval(['S_' num2str(abs(cl(i))) '=Sact;']); %reassign S  

    if isempty(Cact) 
        cl_3(i)=0; 
    end
    U=[]; del_ind=[]; 
end 
cl=cl_3(cl_3 ~= 0); % reset to the remaining levels  
disp(['Done, Elapsed time is ' num2str(ceil(toc)) ' seconds'])
disp('...')

%% BLOCK 4 
% now test the Height-to-Base values
disp('Testing against slope criteria & maximum dimension...') 
cl_4=cl; del_ind=[];

for i=1:1:length(cl)
    eval(['Cact=C_' num2str(abs(cl(i))) ';']); %get active C level 
    eval(['Sact=S_' num2str(abs(cl(i))) ';']); %get active S level 
    t1=size(Cact); 

    del_ind=[]; 
    for j=1:t1(2); %for each closed contour at the active level 
        InPts=inpolygon(lon, lat,Cact(:,j,1), Cact(:,j,2)); 
        maxQ100=max(z(InPts==1)-cl(i)); minQ100=min(z(InPts==1)-cl(i)); 
        thr=maxQ100/sqrt((length(find(InPts(:)==1))*1000*deg2km(dx)*1000*deg2km(dy)*DEP)/pi);  %height:basal radius ratio (circular equivalent)
        % 
   
        if  maxQ100 < clstep*2 || minQ100 < -clstep || thr < minHR || thr > maxHR ;  % apply the test(s) 
            del_ind=cat(1,del_ind,j); % delete column if fail 
            disp('Deleting contour due to Min/Max Slope failure'); 
        end  
    end
    Cact(:,del_ind,:)=[]; Sact(:,del_ind,:)=[]; % deletes closed ones 

    % now trim extra NaNs to keep the matrix small
    if ~isempty(Cact) 
        Cactmean=nanmean(Cact(:,:,1),2);
        [ro,~,~]=ind2sub(size(Cactmean),find(isnan(Cactmean)));
        if ~isempty(ro); Cact=Cact(1:min(ro)-1,:,:);  end
    end
 
    eval(['C_' num2str(abs(cl(i))) '=Cact;']); %reassign C
    eval(['S_' num2str(abs(cl(i))) '=Sact;']); %reassign S  

    if isempty(Cact) 
        cl_4(i)=0; 
    end
    del_ind=[]; 
end
cl=cl_4(cl_4 ~= 0); % reset to the remaining levels  
disp(['Done, Elapsed time is ' num2str(ceil(toc)) ' seconds'])
disp('...')

%% BLOCK 5 
% anything remaining is a base, provided it doesn't lie within another
% contour 
disp('Finding all basal contours....') 
disp(''); 
cl=flipud(cl); % must start at the top and work down 
cl_5=cl; 
InPts=[]; del_ind=[]; 

for i=1:1:length(cl)
    eval(['Cact=C_' num2str(abs(cl(i))) ';']); %get active S level 
    eval(['Sact=S_' num2str(abs(cl(i))) ';']); %get active S level 
    t1=size(Sact); 
    cl(i); 
    for j=1:t1(2); % each contour represented in the S matrix  
        cnt=1;   % used to check against each contour at a higher level 
        for ii=i+1:length(cl) % loop through the higher levels 
            eval(['Cact2=C_' num2str(abs(cl(ii))) ';']); %get contours at next level(s) 
            t2=size(Cact2); 
            for jj=1:t2(2) % and all the contours at those levels 
                InPts(cnt)=inpolygon(Sact(2,j), Sact(3,j),Cact2(:,jj,1),Cact2(:,jj,2)); 
                cnt=cnt+1; 
            end
        end   
        if  sum(InPts) > 0 % if it lies in one higher level cont 
            del_ind=cat(1,del_ind,j); % delete column if fail 
        end  
    end 

    Cact(:,del_ind,:)=[]; Sact(:,del_ind,:)=[]; % deletes closed ones 

    % now trim extra NaNs to keep the matrix small
    if ~isempty(Cact) 
        Cactmean=nanmean(Cact(:,:,1),2);
        [ro,~,~]=ind2sub(size(Cactmean),find(isnan(Cactmean)));
        if ~isempty(ro); Cact=Cact(1:min(ro)-1,:,:);  end
    end
  
    eval(['C_' num2str(abs(cl(i))) '=Cact;']); %reassign C
    eval(['S_' num2str(abs(cl(i))) '=Sact;']); %reassign S  
    InPts=[]; % clear 
    
    if isempty(Cact); cl_5(i)=0; end
 
    del_ind=[];  % reset as you go to the next level Sact Cact 
end
cl=cl_5(cl_5 ~= 0); % reset to the remaining levels  


% BLOCK 5B % check the unusal case where there are two contours at the same
% level 
cl_5=cl; % use this cl_5 name again, 
InPts=[]; del_ind=[]; 

for i=1:1:length(cl)
    eval(['Cact=C_' num2str(abs(cl(i))) ';']); %get active S level 
    eval(['Sact=S_' num2str(abs(cl(i))) ';']); %get active S level 
    t1=size(Sact);  
    for j=1:t1(2); % each contour represented in the S matrix  
        cnt=1;     
        for jj=1:t1(2); % used to check against each contour at a same level
        InPts(cnt)=inpolygon(Sact(2,j), Sact(3,j),Cact(:,jj,1),Cact(:,jj,2)); 
                cnt=cnt+1; 
        if  sum(InPts) > 1 % if should be one, since it's inside itself
            del_ind=cat(1,del_ind,j); % this will contain the indices of the two contours at the same level 
        end  
        end 
    end
   
    % now lets keep the bigger one 
    if ~isempty(del_ind);  
        [~,I]=min(Sact(1,del_ind));  del_ind=del_ind(I);   % will keep the smaller one in the delete list  
    end  
    
    % repeat, if there are three at the same level (very unlikely) 
     if ~isempty(del_ind);  % still not empty   
        [~,I]=min(Sact(1,del_ind));  del_ind=del_ind(I); 
    end  
     
    Cact(:,del_ind,:)=[]; Sact(:,del_ind,:)=[]; % deletes closed ones 

    % now trim extra NaNs to keep the matrix small
    if ~isempty(Cact) 
        Cactmean=nanmean(Cact(:,:,1),2);
        [ro,~,~]=ind2sub(size(Cactmean),find(isnan(Cactmean)));
        if ~isempty(ro); Cact=Cact(1:min(ro)-1,:,:);  end
    end
  
    eval(['C_' num2str(abs(cl(i))) '=Cact;']); %reassign C
    eval(['S_' num2str(abs(cl(i))) '=Sact;']); %reassign S  
    InPts=[]; % clear 
    
    if isempty(Cact); cl_5(i)=0; end
 
    del_ind=[];  % reset as you go to the next level Sact Cact 
end


cl=cl_5(cl_5 ~= 0); % reset to the remaining levels  
cl=flipud(cl);      % return the order to low to high 
% 

%% BLOCK 6 - display final result 
C_row1=[]; C_row2=[]; ccount=0; % final number of bases 
for q=1:length(cl);
    eval(['Cact=C_' num2str(abs(cl(q))) ';']); %get active C level 
    eval(['Sact=S_' num2str(abs(cl(q))) ';']); %get active S level 
    tt=size(Cact);
    for qq=1:tt(2) 
        C_row1=cat(1,C_row1,Cact(1:Sact(1,qq),qq,1),NaN); C_row2=cat(1,C_row2,Cact(1:Sact(1,qq),qq,2),NaN); 
        ccount=ccount+1; 
    end 
end


if ~isempty(C_row1)
CNEW=cat(1,C_row1',C_row2'); plotm(double(CNEW(2,:)),double(CNEW(1,:)),'w','LineWidth',2);
%% 
%% BLOCK 7 - RETURN THE BASES 
% generate the structured array 
for k=1:length(cl);
    eval(['ccontour_t.C_' num2str(abs(cl(k))) '=C_' num2str(abs(cl(k))) ';']);
    eval(['stats_t.S_' num2str(abs(cl(k))) '=S_' num2str(abs(cl(k))) ';']);
end
disp(['Done, Located ' num2str(ccount) ' basal contours in ' num2str(ceil(toc)) ' seconds'])
disp('...')

else % if ~isempty(C_row1)
    disp('.......') 
    disp('.......') 
  disp('Sorry no closed contour bases, for the criteria specified') 
  ccontour_t=[]; 
  stats_t=[];
  cl=[];
end

    
%%


