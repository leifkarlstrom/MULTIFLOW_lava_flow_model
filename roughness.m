addpath("/Users/akubo/Desktop/DEMS/rainer");
addpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1")

filename='output_be.tif';
[A, R]=geotiffread('output_be.tif');
A(A<0)=0;
I = geotiffinfo(filename); 

[Ny, Nx] = size(A);

N= Ny*Nx;

figure;

n= 1000;
r=zeros(1,n);
area=zeros(1,n);
for j =1:1000
    x=randi([j+1, Nx-j-1])
    y=randi([j+1,Ny-j-1])

    y1=y-j
    y2=y+j 

    x1=x-j 
    x2=x+j
    window=A(y1:y2, x1:x2);
    area(j)= (x+2*j)*(y+2*j);
    r(j)=std(window(:));
    %[DEM_detrend, a,b] = Detrend2(temp);
    %krms=sqrt( sum(DEM_detrend(:).^2 )/i^2);
    %plot(log10(i^2), log10(krms), 'bo')
    
    %hold on 
    %temp=A(j:end, j:end);
    %[x,y]=size(temp);
    %N=x*y;
    %[DEM_detrend, a,b] = Detrend2(temp);
    %krms(j/10)=sqrt( sum(DEM_detrend(:).^2 )/N);
   
end
plot([1:1000]*2, r, 'b.')

