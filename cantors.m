function [N]= cantor(NN)

    [r,c] = size(NN);

    if r~= 2 && c~= 2
        error('cantor:not 2 columns', 'matrix must be r by 2')
    end 

    if c >r
        NN=transpose(NN);
    end 
    
    LENGTH = max(r,c);
    N= zeros(LENGTH, 1);
    for i = 1:LENGTH
        x= NN(i,1);
        y= NN(i,2);
        N(i)= 0.5*(x+y)*(x+y+1)+y;
    end
    

end  


