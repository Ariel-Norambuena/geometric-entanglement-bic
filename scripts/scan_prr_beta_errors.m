function summary=scan_prr_beta_errors(runMode)
%SCAN_PRR_BETA_ERRORS Physical modulation-depth errors after ideal loading.
% A fractional beta error shifts zeros of J0. Errors are switched on at the
% start of the passage, not during the local loading pulse. The latter would
% require a separate multi-excitation laboratory-frame validation.
if nargin<1, runMode='paper'; end
runMode=validatestring(runMode,{'paper','dev'},mfilename,'runMode');
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','matlab'));
if strcmp(runMode,'paper')
    out=fullfile(root,'data','prr_revision');
else
    out=fullfile(root,'data','development','prr_revision');
end
if ~exist(out,'dir'), mkdir(out); end
if strcmp(runMode,'paper'), nc=804; ng=21; else, nc=204; ng=5; end
m=giant_atom_static_model(struct('Nc',nc));
ctrl=floquet_controls(0.9,20,8);
tt=linspace(0,20,4001);
[b1,b2]=ctrl.beta(tt);
interp1=griddedInterpolant(tt,b1,'pchip');
interp2=griddedInterpolant(tt,b2,'pchip');
epsilons=linspace(-0.1,0.1,ng);
result=zeros(ng*ng,8); states=zeros(m.dimension,ng*ng); count=0;
opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.05);
for iy=1:ng
    for ix=1:ng
        e1=epsilons(ix); e2=epsilons(iy);
        c=ctrl; c.u=@(t) beta_u(t,e1,e2,interp1,interp2);
        p0=zeros(m.dimension,1); p0(1)=1;
        cd=struct('counterdiabaticScale',1);
        [~,yp]=ode113(@(t,p) floquet_effective_rhs(t,p,m,c,cd),[0 10 20],p0,opts);
        pend=transpose(yp(end,:));
        [~,yh]=ode113(@(t,p) floquet_effective_rhs(t,p,m,c), ...
            linspace(20,120,501),pend,opts);
        ob=observables_single_excitation(transpose(yh));
        count=count+1;
        result(count,:)=[e1 e2 ob.bellPlusFidelity(1) ob.bellPlusFidelity(end) ...
            min(ob.bellPlusFidelity) ob.concurrence(end) ...
            max(abs(ob.norm-1)) ob.photonicPopulation(end)];
        states(:,count)=transpose(yh(end,:));
    end
    fprintf('Physical beta errors: row %d/%d\n',iy,ng);
end
vars={'epsilon_beta1','epsilon_beta2','F_preparation','F_after_hold', ...
    'minimum_F_hold','C_after_hold','norm_error','photonic_population_after_hold'};
writetable(array2table(result,'VariableNames',vars),fullfile(out,'beta_error_maps.csv'));
% Check interpolation against the existing analytic controls and test grid
% convergence at the nominal point and four corners.
checkPoints=[0 0;-0.1 -0.1;-0.1 0.1;0.1 -0.1;0.1 0.1];
checks=[];
for ncheck=[804 2004 4004]
    mm=giant_atom_static_model(struct('Nc',ncheck));
    for ic=1:size(checkPoints,1)
        ep=checkPoints(ic,:); c=ctrl;
        c.u=@(t) beta_u(t,ep(1),ep(2),interp1,interp2);
        p0=zeros(mm.dimension,1); p0(1)=1;
        [~,yy]=ode113(@(t,p) floquet_effective_rhs(t,p,mm,c, ...
            struct('counterdiabaticScale',1)),[0 10 20],p0,opts);
        [~,hh]=ode113(@(t,p) floquet_effective_rhs(t,p,mm,c), ...
            [20 70 120],transpose(yy(end,:)),opts);
        checks=[checks; ncheck ep abs(yy(end,1)+yy(end,2))^2/2 ...
            abs(hh(end,1)+hh(end,2))^2/2]; %#ok<AGROW>
    end
end
writetable(array2table(checks,'VariableNames',{'Nc','epsilon_beta1','epsilon_beta2', ...
    'F_preparation','F_after_hold'}),fullfile(out,'beta_map_grid_checks.csv'));
summary=struct('Nc',nc,'gridSize',ng,'T',20,'holdTime',100,'u0',0.9, ...
    'errorOnset','after ideal local loading','RelTol',1e-9,'AbsTol',1e-11, ...
    'minimumPreparationFidelity',min(result(:,3)), ...
    'minimumHoldFidelity',min(result(:,4)),'maximumNormError',max(result(:,7)));
save(fullfile(out,'beta_error_maps.mat'),'result','states','summary','checks','-v7');
fid=fopen(fullfile(out,'beta_error_summary.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));
end

function [u1,u2]=beta_u(t,e1,e2,f1,f2)
t=min(20,max(0,t));
u1=besselj(0,(1+e1)*f1(t));
u2=besselj(0,(1+e2)*f2(t));
end
