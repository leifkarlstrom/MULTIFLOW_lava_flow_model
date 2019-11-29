% sensitivity analysis
format longE
% a, b, c, FilteredWavelength, Area, Volume, Bi, Length 

%for w= 0.1:0.1:1.0
%    for q= 0.1:0.1:1.0
%        for r = 1:1:10
%            for FilteredWavelength=50:50:250

A = load('sensitivity_param.txt');
FilteredWavelength = 150;
b = 0.5;
c = 5;

for b=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 1]
for c = [1,2,3,4,5,6,7,8,9,10]
data= A( A(:,2)==b & A(:,3)==c & A(:,4)==FilteredWavelength, : );
a=0.1:0.1:1.0;
plotresults(a, data)
end 
end 

figure; 
for a=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 1]
for c= [1,2,3,4,5,6,7,8,9,10]
data= A( A(:,1)==a & A(:,3)==c & A(:,4)==FilteredWavelength, : );
b=0.1:0.1:1.0;
plotresults(b, data)
end 
end 

figure; 
for a=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 1]
for b=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 1]
data= A( A(:,1)==a & A(:,2)==b & A(:,4)==FilteredWavelength, : );
c=1:1:10;
plotresults(c, data)
end 
end 


function plotresults(x, data)
    %area
    y1=data(:,5);
    %length 
    y2=data(:,6);
    % bi
    y3=data(:,7);

    subplot(3,1,1)
    scatter(x,y1)
    ylabel('Area')
    hold on


    subplot(3,1,2)
    scatter(x,y2)
    ylabel("Length")
    hold on

    subplot(3,1,3)
    scatter(x,y3)
    ylabel("Branching Metric")

    hold on 
end 