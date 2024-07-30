function x=fluweight(Em,Ex,opt,ram)
%x=fluweight(Em,Ex,opt,ram)
%
% INPUT
%  Em   Emission wavelengths (nm)
%  Ex   Excitation wavelengths (nm)
%  opt  First is for 1-order Rayleigh, Second is for Raman, 
%        third is for 2-order Rayleigh
%  ram  Distance from 1st order Rayleigh to the Raman scatter line
%        Default: 3600 cm^-1 (water)
% 
% OUTPUT
%  x    Matrix used for weighting the samples

%Finding where the Raman scatter line is centered.

%Converting from wavelength to wavenumber and back
if nargin<4 || isempty(ram)
    ram=3600;
end
exv=1./(Ex*1e-7)-ram;
emram=(1./exv)*1e7;

%Make a normal distribution curve
a=normpdf(-4:.01:4,0,1);

%Setting the weight landscape
x=zeros(length(Em),length(Ex));

%If all values are not given in the input, set to the default.
if nargin<3 || length(opt)<3
    opt=[opt zeros(1,3-length(opt))];
end

%Setting the weights
for i=1:length(Ex)
    
    %1st order Rayleigh scatter
    if opt(1)~=0
        [j,k]=sort(abs(Em-Ex(i)));
        m=sort(k(j<=opt(1)));
        b=-opt(1):2*opt(1)/800:opt(1);
        for n=1:length(m)
            [o,p]=min(abs(b-(Em(m(n))-Ex(i))));
            x(m(n),i)=a(p);
        end
        x(1:k(1),i)=a(b==0);
    end
    
    %2nd order Rayleigh scatter
    if opt(3)~=0
        [j,k]=sort(abs(Em-Ex(i)*2));
        m=sort(k(j<=opt(3)));
        b=-opt(3):2*opt(3)/800:opt(3);        
        for n=1:length(m)
            [o,p]=sort(abs(b-(Em(m(n))-Ex(i)*2)));
            x(m(n),i)=a(p(1))*3/4;
        end
        if k(1)==length(Em)
        else
            x(k(1):end,i)=a(b==0)*3/4;
        end
    end
    
    %Raman scatter
    if opt(2)~=0
        [j,k]=sort(abs(Em-emram(i)));
        m=sort(k(j<=opt(2)));
        b=-opt(2):2*opt(2)/800:opt(2);
        for n=1:length(m)
            [o,p]=sort(abs(b-(Em(m(n))-emram(i))));
%             x(m(n),i)=a(p(1))*.4;
            x(m(n),i)=a(p(1))*.2; %This is the "normal" one
        end
    end
end

%Setting the weights
x=1-x./max(a);