
function [NN]= cantinverse(N)
    [r,c]= size(N);
    NN= zeros(r, 2);
    for i= 1:r
        z= N(i);
        w= floor( (sqrt(8*z+1) - 1)/2 );
        t = w*(w+1)/2;
        y = z-t;
        x = w-y;

        NN(i,:)=[x,y];
    end 
end 