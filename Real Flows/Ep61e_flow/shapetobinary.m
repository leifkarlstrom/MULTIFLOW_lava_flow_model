function mask=shapetobinary(s)
%s ='2006_ovest_outline.shp';
%s2='colata_1999_west.shp';
%s3='colata_2001_sud.shp';
shp=shaperead(s,'UseGeoCoords', true);
BOX=shp.BoundingBox;
x=(shp.Lon);
y=(shp.Lat);

[Lon1,Lat1]=meshgrid(BOX(1,1):BOX(2,1), BOX(1,2):BOX(2,2));

cond1= isnan(x) | isnan(y); x(cond1)=[]; y(cond1)=[];

figure;
plot(x,y)

mask=inpolygon(Lon1,Lat1,x,y);
figure; 
imshow(mask)
sum(mask(:))
end 

