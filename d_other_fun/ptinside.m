function [ sumpoly] = ptinside(xo,yo,XV,YV)
% 
% function [sumpoly] = ptinside(xo,yo,XV,YV); 
% xo & yo are the coordinate of a point 
% XV and YV are closed polygon matrix 
% where each column is a new polygon 
% 
% This function test to see if the point xo yo 
% is within any of the contours 
% 
% OUTPUT 
% sumpoly gives the number of polygons the point lies within 
% if it does not lie within any polyon, then a zero is returned 
% DRB (NCSU) 2009

[r,c]=size(XV); 
CNT=zeros(1,c); 
    for i=1:c 
CNT(i)=inpolygon(xo,yo,XV(:,i),YV(:,i)); 
    end

sumpoly=sum(CNT); 