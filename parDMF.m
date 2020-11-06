function[delay_margin]=parDMF(g,a,A,B,C)
% delay margin finder (DMF) for neutral systems
% A,B,C: square real matrices satisfying A+B is Hurwitz
%     g: grid step for alpha (coarse grid)
%     a: required accuracy of the solutions 

n = size(A,1);

if max(abs(eig(C)))>=1
error('Spectral radius of C is larger than unity')
end

if max(real(eig((eye(n)-C)\(A+B))))>=0
error('Delay-free system is unstable')
end

alpha_l=-1; alpha_u=1;
points = ceil(abs(alpha_l-alpha_u)/g)+1;
%delta = abs(alpha_l-alpha_u)/(points-1);
% delete(gcp('nocreate'))
% parpool('local',1)
p = gcp;                    % open new pool                    
poolsize = p.NumWorkers;

if 0==exist('Data','dir')       % check if Data folder already exists
    mkdir Data                  % create new directory to save the data
end

%==========================================================================
% Split alpha-grid (coarse) in chunks
%==========================================================================
pc=points-1;
alpha=zeros(1,points);
for kk=0:pc
    alpha(kk+1)=-cos(kk/pc*pi); % Chebyshev extrema points (or Chebyshev nodes)
end
[alpha_splitted]=split(alpha,poolsize);

alpha_chunks=alpha_splitted;
for i=1:poolsize-1
    alpha_chunks{i}=[alpha_splitted{i},alpha_splitted{i+1}(1)];   % extended intervals 
end

%==========================================================================
% First batch in parallel
%==========================================================================
alphaInter{poolsize,1}=[];  % allocate memory
lambdaInter{poolsize,1}=[];
parfor ii=1:poolsize    
    % Numbering
    %----------------------------------------------------------------------
    [FamLambda,alpha] = numbering(n,alpha_chunks{ii},A,B,C);
    ReFamLambda = real(FamLambda);
    
    parsave(sprintf('ReFamLambda%d.mat',ii),alpha,ReFamLambda,0)

    % alpha-solutions in the coarse grid
    %----------------------------------------------------------------------
    indices_zx = zeros(n,n);
    indices = find(sum(abs(ReFamLambda'))>g); % rule out lambda=0 eigenval 

    for kk = 1:length(indices)
        zx = zci(ReFamLambda(indices(kk),:));
        indices_zx(1:length(zx),indices(kk)) = zx;  
    end

    [q,m,index]=find(indices_zx);   % q-th alpha solution for the coarse grid
                                    % m-th factor          
    %alphaIntLocal{1,q}=[]
    if isempty(q)
       MinLocalTau(ii)=NaN; % jump to the "extraction of DM" part of the code
       MinLocalAlpha(ii)=NaN;
    else     
       MinLocalTau(ii)=NaN; %fill MinLocalTau vector with NaN
       MinLocalAlpha(ii)=NaN;
       alphaInter{ii} =  [alpha(index)',alpha(index+1)'];                 % alpha interval were the candidate solution may be found
       lambdaInter{ii} = [diag(FamLambda(m,index)),diag(FamLambda(m,index+1))];  % eigenvalue interval were the candidate solution may be found
    end
end

%==========================================================================
% Split alpha(interval)-solutions and the corresponding lambda values
%==========================================================================
alphaInter=alphaInter(~cellfun('isempty',alphaInter));   %rule out empty cells
lambdaInter=lambdaInter(~cellfun('isempty',lambdaInter));

alphaInter=cell2mat(alphaInter)';
lambdaInter=cell2mat(lambdaInter)';

nn=min(poolsize,size(alphaInter,2)); % with q empty, then nn=0
parsave('Indices.mat',poolsize,points,nn)

[alphaInter_splitted]=split(alphaInter,nn);
[lambdaInter_splitted]=split(lambdaInter,nn);


%==========================================================================
% Second batch in parallel
%==========================================================================
parfor ii=1:nn
    % increase accuracy of alpha-solutions 
    %----------------------------------------------------------------------
    alphaCRIT=zeros(size(alphaInter_splitted{ii},2),1);
    lambdaCRIT=zeros(size(alphaInter_splitted{ii},2),1);
    for jj=1:size(alphaInter_splitted{ii},2)
        [alphaCRIT(jj,:),lambdaCRIT(jj,:)] = intersections(a,alphaInter_splitted{ii}(1,jj), alphaInter_splitted{ii}(2,jj),lambdaInter_splitted{ii}(1,jj),lambdaInter_splitted{ii}(2,jj),A,B,C);
    end  

    % (omega,tau)-solutions in the refined grid
    %----------------------------------------------------------------------
    Omega = imag(lambdaCRIT);
    Alpha = alphaCRIT;

    indexOm = find(Omega); % consider only non-zero omegas 
    Omega = Omega(indexOm);
    Alpha = Alpha(indexOm);

    Numerator = atan2(sqrt(1-Alpha.^2),Alpha);
    indexNum = find(Numerator<0);
    Numerator(indexNum) = Numerator(indexNum)+2*pi;
             
    Tau = -Numerator./Omega;
    indexTau = find(Tau<0);
    Tau(indexTau) = Tau(indexTau)+2*pi./Omega(indexTau);
    Tau(Tau==0) = NaN; % rule out tau=0
    
    parsave(sprintf('AlphaSolns%d.mat',ii),Alpha,Omega,Tau)
      
    [MinTau,indexMinTau]=min(Tau);
    MinAlpha=Alpha(indexMinTau);

    MinLocalTau(ii)=MinTau;  
    MinLocalAlpha(ii)=MinAlpha; 
end

%==========================================================================
% extraction of DM
%==========================================================================
[DM,indexDM]=min(MinLocalTau,[],'omitnan'); % minimum excluding NaN values

AlphaDM=MinLocalAlpha(indexDM);
parsave('DM.mat',DM,AlphaDM,0)

indexDM(isnan(DM))=[];                     % in case MinLocalTau contains only NaN values then make indexDm empty
if isempty(indexDM)
       delay_margin = 'infinite';
       return
end
[delay_margin] = num2str(DM,15);

end


function parsave(fname, x,y,z)

savdir = '.\Data';
save(fullfile(savdir,fname),'x','y','z');

%save(fname, 'x', 'y', 'z')
end


function[zx]=zci(v)
% zci: zero-crossing indices

zx_temp = find(v(:).*circshift(v(:), [-1 0]) < 0);
    
    if isempty(zx_temp)
        zx=0;
        elseif  sign(v(1))*sign(v(end))<0
            zx=zx_temp(1:end-1); % rule out last element
        else
        zx=zx_temp;
    end
end


function[batch_splitted]=split(batch,numsectors)
    noPointsWorker=diff(fix(linspace(0, size(batch,2), numsectors+1)));  % size of the splited alpha vectors: to be send to the workers
    batch_splitted=mat2cell(batch,size(batch,1),noPointsWorker);           % split alpha vector into  (~) even chunks
end



function [alphaCRIT,lambdaCRIT] = intersections(a,alpha1,alpha2,lambda1,lambda2,A,B,C)
% Increasing accuracy of the zero-crossings
n = size(A,1);
iteration = 50;
tol = a;
cnt = 1;
flag1 = 0;
flag2 = 0;

% Zero-crossings
%--------------------------------------------------------------------------
    for i = 1:1
        if real(lambda1(i,1))*real(lambda2(i,1)) < 0
            cnt2 = i;
            break
        end
    end
    f2s = lambda1(cnt2,1);
    f2e = lambda2(cnt2,1);

% Iterations
%--------------------------------------------------------------------------
    while cnt<=iteration
          p = alpha1+(alpha2-alpha1)/2;
          tempF2 = eig( (eye(n)-C*(p + sqrt(-1)*sqrt(1-p^2)))\(A+B*(p + sqrt(-1)*sqrt(1-p^2))) ); % using the upper semi-circle: alpha + j*sqrt(1-alpha^2)
          [~,iii] = min(abs(real(tempF2)-(real(f2e)+real(f2s))/2));
          lambdaCRIT = tempF2(iii);

          if abs(alpha2-alpha1)/2 <= tol
             flag1 = 1;
          end
          if cnt==iteration
             flag2 = 1;
          end
          if (flag1==1) || (flag2==1)
              alphaCRIT = p;
              break
          end
          cnt = cnt+1;
          if f2s*real(lambdaCRIT) < 0
             alpha2 = p;
             f2e = lambdaCRIT;
          end
          if real(lambdaCRIT)*f2e < 0
             alpha1 = p;
             f2s = lambdaCRIT;
          end
    end
end



function [FamLambda,alpha] = numbering(n,alpha_chunks,A,B,C)

alpha = alpha_chunks;

z = alpha+sqrt(-1)*sqrt(1-alpha.^2);
Mfun = @(z) (eye(n)-C*z)\(A+B*z);

FamLambda = zeros(n,length(alpha));

M1 = Mfun(z(1));
[nu1,lambda1] = eig(M1);
lambda1 = diag(lambda1);
[~,tags1] = sort(real(lambda1),1,'ascend');
lambda1=lambda1(tags1);
nu1=nu1(:,tags1);

 FamLambda(:,1) = lambda1;
 

for i = 2:length(alpha)
    M2 = Mfun(z(i));
    [nu2,lambda2] = eig(M2);
    lambda2 = diag(lambda2);
    [~,tags2] = sort(real(lambda2),1,'ascend');

    nu2=nu2(:,tags2);

    lambda1=FamLambda(:,i-1);
    lambda2=lambda2(tags2);
    
    
    
    dist = (1-abs(nu1'*nu2))+sqrt( ...
    distancematrix(real(lambda1),real(lambda2)).^2+ ...
    distancematrix(imag(lambda1),imag(lambda2)).^2);

    reorder = munkres(dist)';
    FamLambda(:,i) = lambda2(reorder);
    nu2=nu2(:,reorder);
    
    % ensure the signs of each eigenvector pair were consistent if possible
    S = squeeze(real(sum(nu1.*nu2,1))) < 0;
    nu2(:,S) = -nu2(:,S);
    nu1=nu2; %update nu1
end
end


function d = distancematrix(vec1,vec2)
% simple interpoint distance matrix
[vec1,vec2] = ndgrid(vec1,vec2);
d = abs(vec1 - vec2);
end

function [assignment,cost] = munkres(costMat)
% MUNKRES   Munkres (Hungarian) Algorithm for Linear Assignment Problem. 
%
% [ASSIGN,COST] = munkres(COSTMAT) returns the optimal column indices,
% ASSIGN assigned to each row and the minimum COST based on the assignment
% problem represented by the COSTMAT, where the (i,j)th element represents the cost to assign the jth
% job to the ith worker.
%

% This is vectorized implementation of the algorithm. It is the fastest
% among all Matlab implementations of the algorithm.

% Reference:
% "Munkres' Assignment Algorithm, Modified for Rectangular Matrices", 
% http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html

% version 2.0 by Yi Cao at Cranfield University on 10th July 2008

assignment = zeros(1,size(costMat,1));
cost = 0;

costMat(costMat~=costMat)=Inf;
validMat = costMat<Inf;
validCol = any(validMat,1);
validRow = any(validMat,2);

nRows = sum(validRow);
nCols = sum(validCol);
n = max(nRows,nCols);
if ~n
    return
end

maxv=10*max(costMat(validMat));

dMat = zeros(n) + maxv;
dMat(1:nRows,1:nCols) = costMat(validRow,validCol);

%*************************************************
% Munkres' Assignment Algorithm starts here
%*************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   STEP 1: Subtract the row minimum from each row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minR = min(dMat,[],2);
minC = min(bsxfun(@minus, dMat, minR));

%**************************************************************************  
%   STEP 2: Find a zero of dMat. If there are no starred zeros in its
%           column or row start the zero. Repeat for each zero
%**************************************************************************
zP = dMat == bsxfun(@plus, minC, minR);

starZ = zeros(n,1);
while any(zP(:))
    [r,c]=find(zP,1);
    starZ(r)=c;
    zP(r,:)=false;
    zP(:,c)=false;
end

while 1
%**************************************************************************
%   STEP 3: Cover each column with a starred zero. If all the columns are
%           covered then the matching is maximum
%**************************************************************************
    if all(starZ>0)
        break
    end
    coverColumn = false(1,n);
    coverColumn(starZ(starZ>0))=true;
    coverRow = false(n,1);
    primeZ = zeros(n,1);
    [rIdx, cIdx] = find(dMat(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));
    while 1
        %**************************************************************************
        %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
        %           zero in the row containing this primed zero, Go to Step 5.  
        %           Otherwise, cover this row and uncover the column containing 
        %           the starred zero. Continue in this manner until there are no 
        %           uncovered zeros left. Save the smallest uncovered value and 
        %           Go to Step 6.
        %**************************************************************************
        cR = find(~coverRow);
        cC = find(~coverColumn);
        rIdx = cR(rIdx);
        cIdx = cC(cIdx);
        Step = 6;
        while ~isempty(cIdx)
            uZr = rIdx(1);
            uZc = cIdx(1);
            primeZ(uZr) = uZc;
            stz = starZ(uZr);
            if ~stz
                Step = 5;
                break;
            end
            coverRow(uZr) = true;
            coverColumn(stz) = false;
            z = rIdx==uZr;
            rIdx(z) = [];
            cIdx(z) = [];
            cR = find(~coverRow);
            z = dMat(~coverRow,stz) == minR(~coverRow) + minC(stz);
            rIdx = [rIdx(:);cR(z)];
            cIdx = [cIdx(:);stz(ones(sum(z),1))];
        end
        if Step == 6
            % *************************************************************************
            % STEP 6: Add the minimum uncovered value to every element of each covered
            %         row, and subtract it from every element of each uncovered column.
            %         Return to Step 4 without altering any stars, primes, or covered lines.
            %**************************************************************************
            [minval,rIdx,cIdx]=outerplus(dMat(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));            
            minC(~coverColumn) = minC(~coverColumn) + minval;
            minR(coverRow) = minR(coverRow) - minval;
        else
            break
        end
    end
    %**************************************************************************
    % STEP 5:
    %  Construct a series of alternating primed and starred zeros as
    %  follows:
    %  Let Z0 represent the uncovered primed zero found in Step 4.
    %  Let Z1 denote the starred zero in the column of Z0 (if any).
    %  Let Z2 denote the primed zero in the row of Z1 (there will always
    %  be one).  Continue until the series terminates at a primed zero
    %  that has no starred zero in its column.  Unstar each starred
    %  zero of the series, star each primed zero of the series, erase
    %  all primes and uncover every line in the matrix.  Return to Step 3.
    %**************************************************************************
    rowZ1 = find(starZ==uZc);
    starZ(uZr)=uZc;
    while rowZ1>0
        starZ(rowZ1)=0;
        uZc = primeZ(rowZ1);
        uZr = rowZ1;
        rowZ1 = find(starZ==uZc);
        starZ(uZr)=uZc;
    end
end

% Cost of assignment
rowIdx = find(validRow);
colIdx = find(validCol);
starZ = starZ(1:nRows);
vIdx = starZ <= nCols;
assignment(rowIdx(vIdx)) = colIdx(starZ(vIdx));
cost = trace(costMat(assignment>0,assignment(assignment>0)));
end

function [minval,rIdx,cIdx]=outerplus(M,x,y)
[nx,ny]=size(M);
minval=inf;
for r=1:nx
    x1=x(r);
    for c=1:ny
        M(r,c)=M(r,c)-(x1+y(c));
        if minval>M(r,c)
            minval=M(r,c);
        end
    end
end
[rIdx,cIdx]=find(M==minval);
end




