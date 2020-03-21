function DEMmirror = Multiflow_2Dspecfilter(DEM,Ny,Nx)
%mirror dem (to minimize edge effects) and perform spectral filtering
%Leif Karlstrom 10/2018

%assumes DEM is even in size
if rem(Ny,2)==1||rem(Nx,2)==1
    error('not even dimension')
end

%construct tiles that are 1/2 the DEM X and Y size, resulting DEM will be
%2*Ny by 2*Nx

LR = flip(DEM,2); %flip along vertical dimension
LS=LR(1:Ny,Nx/2+1:Nx); 
RS=LR(1:Ny,1:Nx/2); 

UD = flip(DEM,1); %flip along horizontal dimension
US=UD(1:Ny/2,1:Nx);
DS=UD(Ny/2+1:Ny,1:Nx);

%for corners, average each edge and taper
ULS = zeros(Ny/2,Nx/2); URS = zeros(Ny/2,Nx/2);
LLS = zeros(Ny/2,Nx/2); LRS = zeros(Ny/2,Nx/2);

taper=0.993; %amount to taper corners

            ULS(1,1:Nx/2)=LS(end,:)*taper;
            LLS(end,1:Nx/2)=LS(1,:)*taper;
            URS(1,1:Nx/2)=RS(end,:)*taper;
            LRS(end,1:Nx/2)=RS(1,:)*taper;
            
            ULS(:,end)=US(:,1)*taper;
            LLS(:,end)=DS(:,1)*taper;
            URS(:,1)=US(:,end)*taper;
            LRS(:,1)=DS(:,end)*taper;

for i=2:Ny/2
            ULS(i,1:Nx/2-i)=ULS(i-1,1:Nx/2-i)*taper;            
            LLS(end-i+1,1:Nx/2-i)=LLS(end-i+2,1:Nx/2-i)*taper;
            
            URS(i,i:Nx/2)=URS(i-1,i:Nx/2)*taper; 
            LRS(end-i+1,i:Nx/2)=LRS(end-i+2,i:Nx/2)*taper;            
end

for i=2:Nx/2
            ULS(i:Ny/2,end-i+1)=ULS(i:Ny/2,end-i+2)*taper;
            URS(i:Ny/2,i)=URS(i:Ny/2,i-1)*taper;
            
            LLS(1:Ny/2-i+1,end-i+1)=LLS(1:Ny/2-i+1,end-i+2)*taper;
            LRS(1:Ny/2-i+1,i)=LRS(1:Ny/2-i+1,i-1)*taper;
                      
end


DEMmirror1=cat(2,ULS,US,URS);
DEMmirror2=cat(2,LS,DEM,RS);
DEMmirror3=cat(2,LLS,DS,LRS);

DEMmirror = cat(1,DEMmirror3,DEMmirror2,DEMmirror1);
clear DEMmirror1 DEMmirror2 DEMmirror3

%multiple with cosine taper that is restricted to mirrored parts
xl = repmat(1:Nx/2,2*Ny,1);
DEMmirror(:,1:Nx/2) = DEMmirror(:,1:Nx/2).*(0.5*(1 + cos(pi*(Nx/2-xl)./(Nx/2-1)))).^2;

xr = repmat(3*Nx/2+1:2*Nx,2*Ny,1);
DEMmirror(:,3*Nx/2+1:2*Nx) = DEMmirror(:,3*Nx/2+1:2*Nx).*(0.5*(1 + cos(pi*(xr-(3*Nx/2+1))./(Nx/2)))).^2;

yd = repmat(transpose(1:Ny/2),1,2*Nx);
DEMmirror(1:Ny/2,:) = DEMmirror(1:Ny/2,:).*(0.5*(1 + cos(pi*(Ny/2-yd)./(Ny/2-1)))).^2;

yu = repmat(transpose(3*Ny/2+1:2*Ny),1,2*Nx);
DEMmirror(3*Ny/2+1:2*Ny,:) = DEMmirror(3*Ny/2+1:2*Ny,:).*(0.5*(1 + cos(pi*(yu-(3*Ny/2+1))./(Ny/2)))).^2;
