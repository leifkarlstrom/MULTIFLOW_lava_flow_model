% sensitivity analysis
format longE
close all
% 1a, 2b, 3c, 4FilteredWavelength, 5Area, 6Volume, 7Bi, 8Length 

%for w= 0.1:0.1:1.0
%    for q= 0.1:0.1:1.0
%        for r = 1:1:10
%          for FilteredWavelength=50:50:250

data1 = load('paramlong2.txt');
data2=load('paramlong.txt');
data1=data1( data1(:,3) ~= 0,:);


data= [ data1; data2];


Area=data(:,7);
Bi=data(:,9);
L=data(:,10);


% scatterplt(data,Area, 'Area')
% scatterplt(data,Bi, 'Branching Index')
% scatterplt(data, L, 'Length')

% scatterpltcloser(data,Area, 'Area')
% scatterpltcloser(data,Bi, 'Branching Index')
% scatterpltcloser(data, L, 'Length')


a=data(:,1);
b=data(:,2);
c= data(:,3);
fw= data(:,6);



%function oneatatime(data, metric, metricn)
    a= unique(data(:,1));
    b= unique(data(:,2));
    c= unique(data(:,3));
    fw = unique(data(:,6));

    figure;
    % a 
    for bb= 1 %[ b(1), b(ceil(end/2)), b(end)]
        for cc= 10 %[ c(1), c(ceil(end/2)), c(end)]
            for fww =  50 %[ fw(1), fw(ceil(end/2)), fw(end)]
                temp = data1( data1(:,2)== bb & data1(:,3)== cc & data1(:,6)==fww, : );
                plotresults(1,temp, 'a parameter')
                temp = data2( data2(:,2)== bb & data2(:,3)== cc & data2(:,6)==fww, : );
                plotresults(1,temp, 'a parameter')
            end 
        end 
    end

    figure;
    % b
    for aa= 1 %[ a(1), a(ceil(end/2)), a(end)]
        for cc= 10 %[ c(1), c(ceil(end/2)), c(end)]
            for fww = 50 %[ fw(1), fw(ceil(end/2)), fw(end)]
                temp = data( data(:,1)==aa & data(:,3)== cc & data(:,6)==fww, : );
                plotresults(2,temp, 'b parameter')
            end 
        end 
    end

    figure;
    %c
    for bb= 1 %[ b(1), b(ceil(end/2)), b(end)]
        for aa= 1 %[ a(1), a(ceil(end/2)), a(end)]
            for fww =  50 %[ fw(1), fw(ceil(end/2)), fw(end)]
                temp = data( data(:,1)== aa & data(:,2)== bb & data(:,6)==fww, : );
                plotresults(3,temp, 'c parameter')
            end 
        end 
    end

    figure;
    %fw
    for bb= 1 %[ b(1), b(ceil(end/2)), b(end)]
        for aa= 1 %[ a(1), a(ceil(end/2)), a(end)]
            for cc =  10 %[ c(1), c(ceil(end/2)), c(end)]
                temp = data( data(:,1)== aa & data(:,2)== bb & data(:,3)==cc, : );
                plotresults(6,temp, 'Filter Wavelength parameter')
            end 
        end 
    end

function plotresults(n, data, ntitle)
    nx=512; 
    ny=512;
    dx=10;

    totArea=ny*nx*dx*dx;
    xlimit=[ 15, 15, 100, 0, 0, 250];
    xlow=[ -1, -1, 10, 0, 0, 50];
    x= data(:,n);
    %area
    y1=data(:,7)/totArea;
    %length 
    y2=data(:,10)/max(data(:,10));
    % bi
    y3=data(:,9);

    % area
    subplot(3,1,1)
    sgtitle(ntitle)
    scatter(x,y1, 'k')
    hold on
    plot(x,y1, 'b')
    ylabel('Area')
    xlim([xlow(n), xlimit(n)])
    ylim([0,1])

    % length
    subplot(3,1,2)
    scatter(x,y2, 'k')
    hold on
    plot(x,y2, 'b')
    ylabel("Length")
    xlim([xlow(n), xlimit(n)])
    ylim([0,1])


    % branching
    subplot(3,1,3)
    scatter(x,y3, 'k')
    hold on
    plot(x,y3, 'b')
    ylabel("Branching Metric")
    xlim([xlow(n), xlimit(n)])
    ylim([0, 2])
end 