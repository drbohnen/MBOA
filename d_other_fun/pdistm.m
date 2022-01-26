function p=pdistm(loc)
% [p]=distm(loc)
% input loc is a vector of [lat,lon]
% outputs 1xM vector that is equivalent 
% to the output of pdist, but calculated
% using greatcircle distance in deg 
% DRB June 2010 NCSU
% 

lat=loc(:,1); 
lon=loc(:,2); 
c=1; p=[]; 
while c < length(lon) % for each event
c2=(c+1):1:length(lon);
r=distance(lat(c),lon(c),lat(c2),lon(c2));  
p=cat(1,p,r);
c=c+1; 
end
p=p';