% sensitivity analysis
format longE
% 1a, 2b, 3c, 4FilteredWavelength, 5Area, 6Volume, 7Bi, 8Length 

%for w= 0.1:0.1:1.0
%    for q= 0.1:0.1:1.0
%        for r = 1:1:10
%          for FilteredWavelength=50:50:250

data = load('param2.txt');
Area=data(:,5);
Bi=data(:,7);
L=data(:,8);


scatterplt(data,Area)
scatterplt(data,Bi)
scatterplt(data, L)



% area
function scatterplt(data, Area);

a=data(:,1);
b=data(:,2);
c= data(:,3);
fw= data(:,4);
figure;
subplot(4,1,1)

scatter(a,Area)

subplot(4,1,2)


scatter(b,Area)

subplot(4,1,3)
scatter(c,Area)

subplot(4,1,4)
scatter(fw, Area)
end 
% for b=[0.1, 0.5, 1]
% for c = [1,5,10]
% data= A( A(:,2)==b & A(:,3)==c & A(:,4)==FilteredWavelength, : );
% a=0.1:0.1:1.0;
% plotresults(a, data)
% end 
% end 
% 
% figure; 
% for a=[0.1, 0.5, 1]
% for c= [1,5,10]
% data= A( A(:,1)==a & A(:,3)==c & A(:,4)==FilteredWavelength, : );
% b=0.1:0.1:1.0;
% plotresults(b, data)
% end 
% end 
% 
% figure; 
% for a=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 1]
% for b=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 1]
% data= A( A(:,1)==a & A(:,2)==b & A(:,4)==FilteredWavelength, : );
% c=1:1:10;
% plotresults(c, data)
% end 
% end 


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