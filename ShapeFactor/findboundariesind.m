FlowMap=SylviaBin;
[N,M]=size(FlowMap);
FlowMap(FlowMap==0)=NaN;

[X,Y]=meshgrid(xpt,ypt);


flowX=X.*FlowMap;
flowY=Y.*FlowMap;

row1=find(xpt==min(flowX(:)));
row2=find(xpt==max(flowX(:)));

col1=find(ypt==min(flowY(:)));
col2=find(ypt==max(flowY(:)));