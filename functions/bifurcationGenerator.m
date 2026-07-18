% Symbolically solves for all equilibria of the model (modelEqs.m) as a
% function of rho, evaluates them over a rho grid refined near the
% bifurcation points, and classifies each as stable/unstable via the
% Jacobian eigenvalues. Used to generate the bifurcation diagrams.
%
% Inputs:
%   pars        = [mu1, mu2, r, d, kappa, lambda, xi]  (rho excluded; it
%                 is the bifurcation parameter, solved symbolically)
%   inputRhoVec = [rhoMin, rhoMax], the range to sweep
%
% Outputs (each a 1x9 cell, one entry per branch: extinct, uninfected,
% 3x cooperator-only, 3x defector-only, mixed):
%   stableCell, unStableCell: rows are [H_U,H_C,H_D,H_M,V_C,V_D,rho,
%                              eigenvalues(1:6),stability]
%   rhoVec: the (refined, sorted) rho grid used

function [stableCell,unStableCell,rhoVec] = bifurcationGenerator(pars,inputRhoVec)

syms rho Hu Hc Hd Hm Pc Pd
mu1 = pars(1);
mu2 = pars(2);
r = pars(3);
d = pars(4);
kappa = pars(5);
lambda = pars(6);
xi = pars(7);


  alpha = 1; 
  beta = 0;
  g = 1+ rho; 
 delta = rho; 

 % The model equations (see modelEqs.m), symbolic in rho
 dotHu = r.*Hu - xi.*(Hu).^2 - mu1.*Hu.*(Pc + Pd) - mu2.*Hu.*((Pc).^2 + (Pd).^2 + 2.*Pc.*Pd);
 dotHc = mu1.*Hu.*Pc + mu2.*Hu.*(Pc)^2 - d.*Hc;
 dotHd = mu1.*Hu.*Pd + mu2.*Hu.*(Pd)^2 - d .*Hd;
 dotHm = 2.*mu2.*Hu.*Pc.*Pd - d.*Hm;
 dotPc = lambda.*d.*(2.*alpha.*Hc + beta.*Hm) - d.*kappa.*Pc;
 dotPd = lambda.*d.*(g.*Hm + 2.*delta.*Hd) - d.*kappa.*Pd;

 sols = solve([dotHu,dotHc,dotHd,dotHm,dotPc,dotPd] == [0,0,0,0,0,0],[Hu,Hc,Hd,Hm,Pc,Pd]);

diffSols = {1,length(sols.Hu)};
 for i = 1:length(sols.Hu)
     diffSols{i} = [sols.Hu(i), sols.Hc(i),sols.Hd(i),sols.Hm(i),sols.Pc(i),sols.Pd(i)];
 end

 % solve() doesn't label branches, so classify each solution by
 % evaluating numerically at an arbitrary rho (0.2) and checking which
 % of H_C, H_D, H_M are nonzero
 numSols = -1*ones(length(sols.Hu),6);
for i = 1:length(sols.Hu)
     numSols(i,:) = double(subs(diffSols{i},rho,.2));
end

indexList = 1:1:length(sols.Hu);
Cindex = indexList(numSols(:,2) ~= 0 & numSols(:,3) ==0);
Dindex = indexList(numSols(:,2) == 0 & numSols(:,3) ~=0);
Mindex = indexList(numSols(:,4) ~= 0);

 conly1 = diffSols{Cindex(1)};
 conly2 = diffSols{Cindex(2)};
 conly3 = diffSols{Cindex(3)};
%
mixed1 = diffSols{Mindex};
%
donly1 = diffSols{Dindex(1)};
donly2 = diffSols{Dindex(2)};
donly3 = diffSols{Dindex(3)};

% rho values of the bifurcation points, used to refine the sampling grid
% below (collisions: mixed equilibrium meets a C-only/D-only branch;
% rhoS: saddle-node in the mixed branch; rhoTUD: uninfected/D-only
% transcritical)
indexer = double(solve(mixed1(4) == 0,rho));
 collisions = sort(indexer(indexer<1 & indexer>0 & imag(indexer) == 0 ));

 numS = (-2*mu1^3 -9*mu1*mu2*r + 2*(sqrt((mu1^2 + 3*mu2*r)^3)))*d*kappa*xi;
 denomS = 2*mu2*r^2.*(mu1^2 +4*mu2*r)*lambda;
 rhoS = numS./denomS;

 scalarThing = 1/(2*lambda);
 oldTranscritical = (kappa*d*xi)/(mu1*r);
 rhoTUD = scalarThing*oldTranscritical;

 bifVec = sort([collisions(1),collisions(2),rhoS,rhoTUD]);

 % Sample densely near each bifurcation point, sparsely elsewhere
 n=100;
  rhoVec = [linspace(inputRhoVec(1),(bifVec(1)-.001),n), ...
      linspace((bifVec(1)-.001),(bifVec(1)+.001),n), ...
      linspace((bifVec(1)+.001),bifVec(2)-.001,5*n), ...
      linspace((bifVec(2)-.001),(bifVec(2)+.001),5*n), ...
     linspace(bifVec(2)+.001,bifVec(3) - .001,5*n), ...
     linspace((bifVec(3)-.001),(bifVec(3)+.001),5*n), ...
     linspace(bifVec(3)+.001,bifVec(4)-.001,n), ...
     linspace((bifVec(4)-.001),(bifVec(4)+.001),n), ...
     linspace(bifVec(4)+.001,inputRhoVec(end),n)];

 rhoVec = sort(rhoVec);
% 

%%

%first, the c
c1Vec = zeros(size(rhoVec,2),6);
c2Vec = zeros(size(rhoVec,2),6);
c3Vec = zeros(size(rhoVec,2),6);
%next, the d
d1Vec = zeros(size(rhoVec,2),6);
d2Vec = zeros(size(rhoVec,2),6);
d3Vec = zeros(size(rhoVec,2),6);
%finally, mixed
mixedVec = zeros(size(rhoVec,2),6);

for i = 1:length(rhoVec)
c1Vec(i,:) = subs(conly1,rho,rhoVec(i));
c2Vec(i,:) = subs(conly2,rho,rhoVec(i));
c3Vec(i,:) = subs(conly3,rho,rhoVec(i));
%next, the d
d1Vec(i,:) = subs(donly1,rho,rhoVec(i));
d2Vec(i,:) = subs(donly2,rho,rhoVec(i));
d3Vec(i,:) = subs(donly3,rho,rhoVec(i));
%finally, mixed
mixedVec(i,:) = subs(mixed1,rho,rhoVec(i));
end

%and now append on rhoVec
plotVec(:,:,1) = [c1Vec,rhoVec'];
plotVec(:,:,2) = [c2Vec,rhoVec'];
plotVec(:,:,3) = [c3Vec,rhoVec'];
plotVec(:,:,4) = [d1Vec,rhoVec'];
plotVec(:,:,5) = [d2Vec,rhoVec'];
plotVec(:,:,6) = [d3Vec,rhoVec'];
plotVec(:,:,7) = [mixedVec,rhoVec'];


%% and now we'll filter out the imaginary parts

parfor i =1:size(plotVec,3)
      realSol{i}= plotVec(sum(imag(plotVec(:,:,i)),2) == 0,:,i);
end
%% The 3 defector-only roots (realSol{4:6}) aren't consistently ordered
% by the solver across rho, so branches 4/5 get relabeled below by
% tracking which of the two is larger at each rho step.
sv = [size(realSol{4},1),size(realSol{5},1),size(realSol{6},1)];
si = 4:1:6;
sl= si(sv == min(unique(sv)));
twoRoots = [realSol{sl(1)}(:,1),realSol{sl(2)}(:,1)];
maxIndex = zeros(size(twoRoots,1),1);
minIndex = zeros(size(twoRoots,1),1);
for j = 1:size(twoRoots,1)
[~,maxIndex(j)] = max(twoRoots(j,:));
[~,minIndex(j)] = min(twoRoots(j,:));
end

solHigh = zeros(size(twoRoots,1),7);
solLow = zeros(size(twoRoots,1),7);

for i = 1:size(twoRoots,1)
solHigh(i,:) = realSol{maxIndex(i)+3}(i,:);
solLow(i,:) = realSol{minIndex(i)+3}(i,:);
end 

realSolSort = realSol;
realSolSort{sl(1)} = solHigh;
realSolSort{sl(2)} = solLow;

%% Filter out solutions with any negative component (not biologically valid)
for i =1:size(realSolSort,2)
      realPosSol{i}= realSolSort{i}(sum(realSolSort{i}(:,:) <0,2) == 0,:);
end

extinctSol = [zeros(length(rhoVec),6),rhoVec'];
uninfectSol = double([[r/xi,0,0,0,0,0].*ones(length(rhoVec),6),rhoVec']);

rpCell = {extinctSol,uninfectSol,realPosSol{1},realPosSol{2},realPosSol{3},realPosSol{4},...
    realPosSol{5},realPosSol{6},realPosSol{7}};
% E,U,C1,C2,C3,D1,D2,D3,M

% Stability: evaluate the Jacobian's eigenvalues at each equilibrium row.
% Columns end up as: 1-6 equilibria, 7 rho, 8-13 eigenvalues, 14 stability
% (1 = unstable, 2 = stable, i.e. all eigenvalues have negative real part)
symVec = [Hu,Hc,Hd,Hm,Pc,Pd];
jacMat = jacobian([dotHu,dotHc,dotHd,dotHm,dotPc,dotPd],symVec);
symVecplus = [Hu,Hc,Hd,Hm,Pc,Pd,rho];

parfor i = 1:size(rpCell,2)
    for k = 1:size(rpCell{i},1);
        rpCell{i}(k,8:13) = double(eig(subs(jacMat,symVecplus,rpCell{i}(k,1:7))));
        rpCell{i}(k,14) = (sum(real(rpCell{i}(k,8:13)) > 0,2) == 0)+1;
    end
end

% Separate into stable and unstable
for i =1:size(rpCell,2)
    if size(rpCell{i},1) ==0
stableCell{i} = rpCell{i};
unStableCell{i} = rpCell{i};
    else
stableCell{i} = rpCell{i}(rpCell{i}(:,14) == 2,:);
unStableCell{i} = rpCell{i}(rpCell{i}(:,14) == 1,:);
    end
end

end %end entire function 
