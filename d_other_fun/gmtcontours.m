function  [CLmat,STmat,clout]=gmtcontours(cl,ptsmin)
% 
% 
% 
% 
txt1=['[~,result_row] = system([''wc -l contour*' num2str(cl) '*i.xyz | grep xyz | minmax | sed ''''s/\// /g'''' | sed ''''s/>//g''''| awk ''''/ / {print $6}'''''']);']; 
eval(txt1);
txt2=['[~,result2] = system([''ls contour*' num2str(cl) '*i.xyz > level.list'']);']; 
eval(txt2);
[~,result_col] = system('wc -l level.list | awk ''/ / {print $1}''');
result_row=str2num(result_row); result_col=str2num(result_col); 

    CLmat=nan(result_row,result_col,2);
    STmat=nan(3,result_col);
    
if result_row > 0
     
    fd=fopen('level.list','r'); 
    S=textscan(fd,'%s'); f=S{1}; fclose(fd); 
    for i=1:length(f)
    dat=load(char(f(i,:))); 
    if length(dat) > ptsmin && length(dat) > 3; 
        CLmat(1:length(dat(:,1)),i,1)= dat(:,1); % lon 
        CLmat(1:length(dat(:,1)),i,2)= dat(:,2); % lat 
        STmat(1,i)=length(dat); 
        [meanlat,meanlong]=meanm(dat(:,2),dat(:,1)); % input order is lat lon for mapping toolbox call 
        STmat(2:3,i)=[meanlong,meanlat]'; % write out mean lat and lon 
    end
    
    end  
CLmat(:,isnan(CLmat(1,:,1)),:)=[]; % gets rid of any that didn't make the ptsmin cut 
STmat(:,isnan(STmat(1,:)))=[]; % gets rid of any that didn't make the ptsmin cut 
CLmat=single(CLmat); STmat=single(STmat); 
clout=cl; 

else  
    CLmat=[];
    clout=[];
    STmat=[]; 
end

txt3=['delete contour*' num2str(cl) '*i.xyz ']; 
eval(txt3);

