load flowbnd.mat
load('SierraNegraDEM.mat')
SNDEM=sierranegraX;

shp=shaperead('SNflowoutlines-polygon.shp');
[N,M]=size(SNDEM);

BOX=shp.BoundingBox;

X=linspace(BOX(1,1),BOX(2,1), M);
Y=linspace(BOX(1,2),BOX(2,2), N);
[Lon1,Lat1]=meshgrid(X,Y);
WholeFlow=zeros(N,M);
figure(1);
s=scatter(Lon1(:), Lat1(:));
hold off 
alpha(s, 0.1)

hold on 
for ii=1:length(flowbnd);
    x=flowbnd(ii).lon;
    y=flowbnd(ii).lat;
    mask=inpolygon(Lon1,Lat1, x,y);
    
    if flowbnd(ii).hole==0
        WholeFlow=WholeFlow+mask; 
        plot(x,y, 'r')
        hold on
    else 
        WholeFlow=WholeFlow-mask;
        plot(x,y, 'b')
        hold on
    end

end 


figure;
imshow(WholeFlow)

% SYL=zeros(1,length(flowbnd));
% 
% for ii=1:length(flowbnd)
%     if any(flowbnd(ii).lon>-91.116)
%     SYL(ii) = 1;
%     if flowbnd(ii).hole==0
%         %is internal boundary
%         if length(flowbnd(ii).lon)>100
%         figure(1)
%         %plot(flowbnd(ii).x,flowbnd(ii).y,'r')
%         %figure(2)
%         plot(flowbnd(ii).lon,flowbnd(ii).lat,'b')
%         %hold on
%         end
%     elseif flowbnd(ii).hole==1
%         if length(flowbnd(ii).lon)>100
%         %is internal boundary
%         figure(1)
%         %plot(flowbnd(ii).x,flowbnd(ii).y,'w')
%         %figure(2)
%         plot(flowbnd(ii).lon,flowbnd(ii).lat,'r')  
%         hold on
%         end
%     end
%     end
% end
% 
% syl=find(SYL==1);
% 
% for ii=1:length(syl)
%     SylviaFlow(ii).lon=flowbnd(syl(ii)).lon;
%     SylviaFlow(ii).lat=flowbnd(syl(ii)).lat;
%     SylviaFlow(ii).hole=flowbnd(syl(ii)).hole;
% end
% 
% 
% 
% 
% 
% SylviaBin=zeros(N,M);
% 
% for ii=1:length(SylviaFlow)
%     if SylviaFlow(ii).hole==0
%         x=SylviaFlow(ii).lon;
%         y=SylviaFlow(ii).lat; 
%         mask=inpolygon(Lon1,Lat1, x,y);
%         SylviaBin=SylviaBin+mask;
%     else
%         x=SylviaFlow(ii).lon;
%         y=SylviaFlow(ii).lat; 
%         mask=inpolygon(Lon1,Lat1, x,y);
%         SylviaBin=SylviaBin-mask;
%     end
%     
% end 
% 
% figure;
% imshow(SylviaBin)
% 
% 
% % 
% % cond1= isnan(x) | isnan(y); x(cond1)=[]; y(cond1)=[];
% % 
% % figure;
% % plot(x,y)
% % hold on
% % plot(Lon1, Lat1, 'ko')
% % 
% % mask=inpolygon(Lon1,Lat1,x,y);
% % figure; 
% % imshow(mask)
% % sum(mask(:))
% 
