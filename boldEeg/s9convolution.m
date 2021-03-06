% Discrete time convolution takes two discrete time signals as input and gives a discrete time signal as output. 
% Syntax: 
% [y,n] = convolution(x1,n1,x2,n2); 
% where 
% x1 - values of the first input signal - should be a row vector 
% n1 - time index of the first input signal - should be a row vector 
% x2 - values of the second input signal - should be a row vector 
% n2 - time index of the second input signal - should be a row vector
function [y,yn]=s9convolution(x1,n1,x2,n2)
    
    if s10signalValidity(x1,n1)==0
        error('x1 & n1 must be of the same size');
    end 
    if s10signalValidity(x2,n2)==0
        error('x2 & n2 must be of the same size');
    end
    

    % Time reverse the second signal 
    [xr2,nr2]=s11timeReversal(x2,n2);
    
    N1=size(n1,2);
    N2=size(n2,2);
    
    % endpoints of signals
    n11=n1(1,1);
    n12=n1(1,N1);
    n21=n2(1,1);
    n22=n2(1,N2);
    
    % y-resulting signal ; yn-index set of the resulting signal
    y=zeros(1,N1+N2-1);
    yn=n11+n21:n12+n22;
    
    % convolution operation
    i=0;
    for n=n11+n21:n12+n22
        %         fprintf('\nconvolution step %d of %d', n, n12+n22);
        nrs2=s12timeShift(nr2,n);
        p=s13signalMult(x1,n1,xr2,nrs2);
        i=i+1;
        y(1,i)=sum(p);
    end
    
end