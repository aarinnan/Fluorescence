function [model,evolve]=flumod(X,Fac,Em,Ex,wid,op,cons,num)
% [model,evolve]=flumod(X,Fac,Em,Ex,wid,op,cons,num)
%
% Used to decompose an EEM (Excitation-Emission Matrix). The Rayleigh
% scatter line is modeled as a separate factor. 
% OBS! The PARAFAC used here is the PARAFAC from the NWAY-Toolbox
%      which must be first in the path if you have other parafac 
%      algorithms (download at www.models.kvl.dk)
%
% INPUT:
%  X      EEM (Should be Samples x Emission x Excitation)
%  Fac    Number of factors to extract (not in the Rayleigh modelling)
%  Em     Emission wavelengths
%  Ex     Excitation wavelengths
%  wid    Estimated width of Rayleigh peak. If set to empty fluwid will 
%          be used to estimate the width
%  op     Two binary #. 1 - Use (0) PARAFAC or (1) PCA for Rayleigh
%                       2 - Display (1) while calculating
%          Default = [0 1]
%  cons   Constraints for the fluorophores, as in parafac. Default = [0 0 0]
%  num    Number of factors to model Rayleigh. Default = 2
%
% OUTPUT:
%  model  A struct array containing different information concerning the
%          model
%  evolve Each iteration in the modeling step
% 
% Please refer to the paper
% "1st Order Rayleigh as a Separate Component in the 
% Decomposition of Fluorescence Landscapes", Rinnan, Booksh, Bro
% Anal.Chim.Acta, 2005, In Press.


% Copyright, 2005 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
% In case of doubt, contact the holder of the copyrights.
%
% Asmund Rinnan
% E-mail asmundrinnan@gmail.com

dim=size(X);
if length(Em)~=dim(2) || length(Ex)~=dim(3)
    id=[find(dim~=length(Em) & dim~=length(Ex)) find(dim==length(Em)) find(dim==length(Ex))];
    disp('The EEM should be in the format Samples x Emission x Excitation')
    disp(['Type ''X=reshape(nshape(X,[' num2str(id) ']),[' num2str(dim(id)) ']);'])
    disp('PS! Remember to change ''X'' to the variable name you are using!')
    disp(' ')
    error('Wrong dimension of X')
end

%If the width value is not given, it is found through the use of 'fluwid'
if nargin<5 || isempty(wid)
    wid=fluwid(X,Em,Ex);
end

% if length(find(Em>Ex(1)*2-wid))>0
%     disp(' ')
%     disp('This function does not handle 2. order Rayleigh. Please remove those parts of the EEM, and')
%     disp(['try again. Remove the ' num2str(length(find(Em>Ex(1)*2-wid)>0)) ' last points in the emission axis and the same in your x'])
%     error('Too large X data (see above)')
%     return
% end

%Search for the max-value in the matrix. If there are several numbers which
%are max, this is most likely because of saturation of the detector, and
%these numbers are set to missing.
if length(find(vec(X)==max(vec(X))))>1
    disp('It seems that the detector has been saturated. Max values set to missing.')
    X(X==max(vec(X)))=NaN;
end

%Sets the options for the function. The 'opt' is the vector used in
%PARAFAC.
opt=[0 0 0 0 -1 10000];
%Crit is the relative decrease in two consecutive iterations necessary to
%converge. maxit, is the number of iterations allowed in this algorithm.
%(10000 in the PARAFAC modeling - see above).
crit=1e-6;
maxit=1000;

%If some of the parameters are not given, they are set to the default
%values.
if nargin<6 || length(op)<3
    if nargin<6
        op=[0 1];
    else
        op=[op 1];
    end
end
if nargin<7 || length(cons)<3
    cons=[0 0 0];
end
if nargin<8
    num=2;
end

%Defining the rawdata as the X-matrix entered into the function
Xorg=X;
%Finds those Excitation wavelengths that are not influenced by the 1st
%ordered Rayleigh scatter line.
rm=length(find((Ex-Em(1)+wid)<0));
exrm=Ex(1:rm);
%Removes these wavelengths
Ex=Ex(rm+1:end);
X=X(:,:,rm+1:end);
if rm~=0
    model.info.Rayleigh=['Ex ' num2str(exrm) ' removed during Rayleigh model'];
    if op(2)==1
        disp(['The excitation wavelengths: ' num2str(exrm) ' have been removed during Rayleigh modelling'])
    end
end

%Prints the inserted values into the result structure
model.info.Options=num2str(op);
model.info.Constraints=num2str(cons);
model.info.EEMFac=num2str(Fac);
model.info.RayFac=num2str(num);

%Initalizing the algorithm
konv=nansum(vec(Xorg).^2);
olderr=konv;
it=0;
while konv>crit && it<maxit
    %Modeling the Rayleigh scatter
    [Xnew,ray,fac,iter]=raymod(X,Em,Ex,wid,Xorg,rm,op,num);
    
    %Error message revorded if the number of iterations to model the
    %Rayleigh scatter reaches the limit.
    if iter==opt(6)
        model.error{1,1}='Number of iterations in modelling Rayleigh reached max';
    end
    
    %Models the fluorophores
    [facpar,iter]=parafac(Xnew,Fac,opt,cons);
    
    %Error message recorded if the number of iterations in modeling the
    %fluorophores is reaches the limit.
    if iter==opt(6)
        model.error{2,1}='Number of iterations in modelling fluorophores reached max';
    end
    
    %Finds the scores and loadings for the fluorophores, and calculate the
    %fluorophore model.
    [ao,bo,co]=fac2let(facpar);
    mod=reshape(ao*seye(Fac,3)*kron(co,bo)',size(Xnew));
    
    %'err' is the sum of the squared residual, and find the relative change
    %in residual
    err=nansum(vec(Xnew-mod).^2);
    konv=(olderr-err)/olderr;
    olderr=err;
    it=it+1;
    if op(2)==1
        fprintf('%g   %12.8f \n',it,konv)
    end
    
    %Finding the residual X-matrix after the fluorophore model is
    %subtracted
    X=Xorg(:,:,rm+1:end)-mod(:,:,rm+1:end);
    
    %These variables are given if 'evolve' is chose as output
    a{it}=ao;
    b{it}=bo;
    c{it}=co;
    rayleigh{it}=ray;
    rayfac{it}=fac;
    res(it)=err;
end

%If the last iteration led to an increase in residual, the algorithm stops,
%and the previous iteration is the solution.
if konv<0
    ao=a{it-1};
    bo=b{it-1};
    co=c{it-1};
    mod=reshape(ao*seye(Fac,3)*kron(co,bo)',size(Xnew));
    ray=rayleigh{it-1};
    fac=rayfac{it-1};
    err=res(it-1);
    it=it-1;
    if op(2)==1
        disp('The second last iteration used for modelling')
    end
end

%Formating the output.
model.rayleigh=ray;
model.rayfac=fac;
model.mod=mod;
model.parafac={ao,bo,co};
model.err=err;
model.it=it;
if nargout==2
    evolve.parafac={a,b,c};
    evolve.rayfac=rayfac;
    evolve.raymod=rayleigh;
    evolve.res=res;
end

%-------------------------------------
function [Xnew,rayleigh,fac,it]=raymod(X,Em,Ex,wid,Xorg,rm,opt,num)
%The function that models the Rayleigh scatter.

op2=opt(2);
opt=opt(1);

%The options for the PARAFAC algorithm
op=[0 0 0 0 -1 10000];

%Building up the matrix later used to model the 1st order Rayleigh scatter
Xr=ones(size(X,1),wid*2+1,length(Ex))*NaN;
for i=1:length(Ex)
    k=find(Em<=Ex(i)+wid & Em>=Ex(i)-wid);
    j=find(abs(Em-Ex(i))==min(abs(Em-Ex(i))));
    m=find(j(1)==k);
    xray=squeeze(X(:,k,i));
    pos{i,1}=k;
    pos{i,2}=round(Em(k)-Ex(i))+wid+1;
    Xr(:,pos{i,2},i)=xray;
end
%Removing the 2nd mode slabs with only missing values.
nanval=zeros(size(Xr,2),1);
nanval(sum(isnan(squeeze(Xr(1,:,:))))==size(Xr,2))=1;    

%Modeling the Rayleigh peak - either by PCA or PARAFAC
if opt==1
    for i=1:size(Xr,1)
        x=squeeze(Xr(i,nanval==0,:));
        [t,p]=pcamiss(x,num);
        if sum(t)<0
            t=-t;
            p=-p;
        end
        ray(i,:,:)=t*p';
        fac{i}={t,p};
    end
    it=1;
else
    [fac,it]=parafac(Xr(:,nanval==0,:),num,op);
    [a,b,c]=fac2let(fac);
    ray=reshape(a*seye(num,3)*kron(c,b)',size(Xr(:,nanval==0,:)));
end

%Reshaping the model into the original data matrix.
xray=Xr;
Xr(:,nanval==0,:)=ray;
Xr(isnan(xray))=NaN;
rayleigh=zeros(size(Xorg));
for i=1:length(Ex)
    r=Xr(:,pos{i,2},i);
    rayleigh(:,pos{i,1},rm+i)=r;
end

%Subtracting the Rayleigh model from the raw data.
Xnew=Xorg-rayleigh;

%----------------------------------------------------------
function I=seye(n,dim)
%I=seye(n,dim)
%
% returns the superdiagonal identity array
% of dimensions n^dim

I=reshape(zeros(n^dim,1),ones(1,dim)*n);
for i=1:n
   e=num2str(i);
   for d=2:dim
      e=[e ',' num2str(i)];
   end
   eval(['I(',e,')=1;']);
end
I=reshape(I,n,n^(dim-1));

%---------------------------------------------------------
function xnew=vec(x)

xnew=x(:);

%---------------------------------------------------------
function [T,P,R] = pca(X,N);

% [T,P,R] = PCA(X,N)
%           Trekker ut N scorer og ladninger fra X
%           ved hjelp av SVD dekomponering.
%           R er en celle som inneholder N restmatriser
% [T,P] = PCA(X);        
%           gir 2 PCer
%
% Bj?rn Grung

if nargin < 2, N = 2; end

[r,k] = size(X);
if N>r | N>k
    N=min(r,k);
    disp(['The number of components has been reduced to ' num2str(N)])
end
if k < r
  covX = (X'*X)/(r-1);
  [U,S,V] = svd(covX);
else
  covX = (X*X')/(r-1);
  [U,S,V] = svd(covX);
  V = X'*V;
  for i = 1:r
    V(:,i) = V(:,i)/norm(V(:,i));
  end
end

P = V(:,1:N);
T = X*P;

for i=1:N
   R{i}=X-T(:,1:i)*P(:,1:i)';
end

%---------------------------------
function [T,P,a]=pcamiss(X,lv,opt)
% [T,P,a]=pcamiss(X,lv,opt)
%
% INPUT:
%  X   X-matrix
%  lv  Number of components
%  opt 1 - show progress
%
% OUTPUT:
%  T   Scores
%  P   Loadings
%  a   Number of iterations

% ?smund Rinnan 090105

if nargin<3
    opt=0;
end
C=isnan(X);
if sum(vec(C))==0
    [t,p]=pca(X,lv);
    a=1;
else
    X(find(C==1))=rand(sum(sum(C)'),1);
    for i=1:lv
        a=1;
        konv=1;
        while konv>1e-9
            [t,p]=pca(X,1);
            if a>1e5
                error('Convergence has not been reached within the iteration limits')
                break
            end
            Xnew=t*p';
            R=sum((vec(X(C==1))-vec(Xnew(C==1))).^2);
            Rn=sum(vec(Xnew(C==1).^2));
            konv=R/Rn;
            if opt==1
                fprintf('%12.8f \n',konv)
            end
            X(find(C==1))=Xnew(find(C==1));
            a=a+1;
        end
        X=X-t*p';
        T(:,i)=t;
        P(:,i)=p;
    end
end