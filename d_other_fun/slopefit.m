function S=slopefit(r,h)
% where r is range in m 
% h is height in m
% S is the slope in Deg. 
% buily to work with modccbase function. 
% DRB (NC State University 2010)

[row,col]=size(r); 
if col > row; r=r'; end 
[row,col]=size(h); 
if col > row; h=h'; end 
A=[ones(length(r),1), r];
    x=A\h; 
    S= atand(x(2)); 
    