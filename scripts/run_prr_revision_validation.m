function report = run_prr_revision_validation(runMode)
%RUN_PRR_REVISION_VALIDATION Dressed-state, hold, and matched-resource checks.
% All propagations use the existing full-zone finite-mode Hamiltonian.
% The loading pulse is treated separately in the archived two-stage data.
% New comparisons start from |eg,0>, with matched geometry, g, T and the
% same atomic exchange pulse. No Lindblad decay is added to the explicit guide.

if nargin < 1, runMode = 'paper'; end
runMode = validatestring(runMode,{'paper','dev'},mfilename,'runMode');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','matlab'));
paper = strcmp(runMode,'paper');
if paper
    out = fullfile(root,'data','prr_revision');
else
    out = fullfile(root,'data','development','prr_revision');
end
if ~exist(out,'dir'), mkdir(out); end
if paper, nc = 2004; else, nc = 204; end
u0 = 0.9; T = 20; holdTime = 100;
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.04);
params = struct('Nc',nc,'g',0.1,'xi',1,'wc',0,'Omega0',0, ...
    'n1',6,'n2',6,'x1',0,'Dx',2);
model = giant_atom_static_model(params);
report = struct('params',params,'u0',u0,'T',T,'holdTime',holdTime, ...
    'RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.04, ...
    'initialState','eg with vacuum guide','model','central-sideband, lossless');

% Full state trajectories are retained for reproducible BIC and spatial checks.
modes = {'atomic','dressed','fixed','attach','isolated','none'};
nominal = struct();
for im = 1:numel(modes)
    name = modes{im};
    nominal.(name) = simulate(model,u0,T,holdTime,name,opts,true);
    fprintf('%s: F_Bell=%.9f, F_BIC=%.9f, F_Bell(hold)=%.9f\n', ...
        name,nominal.(name).F(401),nominal.(name).FB(401),nominal.(name).F(end));
end
save_nominal(out,nominal,report);

% Reverse the coupling path and the sign of the atomic CD after the hold.
ctrl = floquet_controls(u0,T,8);
rev = ctrl;
rev.u = @(t) ctrl.u(T-t);
rev.thetaDot = @(t) -ctrl.thetaDot(T-t);
cd = struct('counterdiabaticScale',1);
[tr,yr] = ode113(@(t,p) floquet_effective_rhs(t,p,model,rev,cd), ...
    linspace(0,T,401),nominal.atomic.finalState,opts);
retrieval = struct('time',tr+T+holdTime,'psi',transpose(yr));
retrieval.egPopulation = abs(yr(:,1)).^2;
report.retrievedEgPopulation = retrieval.egPopulation(end);
save(fullfile(out,'retrieval.mat'),'retrieval','-v7');

gridSizes = [804 2004 4004];
if ~paper, gridSizes = [204 404]; end
conv = zeros(numel(gridSizes),7);
for ii=1:numel(gridSizes)
    p=params; p.Nc=gridSizes(ii); mm=giant_atom_static_model(p);
    r=simulate(mm,u0,T,holdTime,'atomic',opts,false);
    [b,~,~]=compact_bic_branch(mm,u0,pi/4);
    c= floquet_controls(u0,T,8);
    residual=norm(floquet_effective_rhs(T,b,mm,c,struct('counterdiabaticScale',0)));
    conv(ii,:)=[p.Nc r.prepF r.prepFB r.F(end) min(r.F) residual r.normError];
end
writetable(array2table(conv,'VariableNames',{'Nc','F_Bell','F_BIC', ...
    'F_Bell_hold','minimum_F_Bell','BIC_eigen_residual','norm_error'}), ...
    fullfile(out,'grid_convergence.csv'));
report.gridConvergence=conv;

% Same wave-vector detuning and hold for all preparations. At nonzero
% detuning the nominal compact state is a target, not an exact BIC.
if paper, deltas=linspace(-0.1,0.1,21); else, deltas=[-0.1 0 0.1]; end
matchedModes={'atomic','fixed','attach'};
matched=zeros(numel(deltas)*numel(matchedModes),7); rows=0;
matchedPrepared=complex(zeros(model.dimension,size(matched,1)));
matchedFinal=matchedPrepared;
for ii=1:numel(deltas)
    mm=model;
    omegaK=-2*cos((pi/2)*(1+deltas(ii)));
    mm.wk=mm.wk-omegaK;
    for im=1:numel(matchedModes)
        r=simulate(mm,u0,T,holdTime,matchedModes{im},opts,false);
        rows=rows+1;
        matched(rows,:)=[deltas(ii) im r.prepF r.F(end) min(r.holdF) ...
            r.C(end) r.normError];
        matchedPrepared(:,rows)=r.preparedState;
        matchedFinal(:,rows)=r.finalState;
    end
    fprintf('Matched detuning %2d/%2d complete\n',ii,numel(deltas));
end
writetable(array2table(matched,'VariableNames',{'deltaK_fraction','method', ...
    'F_preparation','F_after_hold','minimum_F_hold','C_after_hold','norm_error'}), ...
    fullfile(out,'matched_k_benchmark.csv'));
report.matchedMethods=matchedModes;
save(fullfile(out,'matched_final_states.mat'),'matchedPrepared','matchedFinal', ...
    'matched','report','-v7');

durations=[5 10 15 20 30 40 60];
durationModes={'atomic','dressed','fixed','attach','none'};
durationData=zeros(numel(durations)*numel(durationModes),6); rows=0;
durationPrepared=complex(zeros(model.dimension,size(durationData,1)));
durationFinal=durationPrepared;
for ii=1:numel(durations)
    for im=1:numel(durationModes)
        r=simulate(model,u0,durations(ii),holdTime,durationModes{im},opts,false);
        rows=rows+1;
        durationData(rows,:)=[durations(ii) im r.prepF r.prepFB r.F(end) ...
            15*pi/(32*durations(ii))];
        durationPrepared(:,rows)=r.preparedState;
        durationFinal(:,rows)=r.finalState;
    end
    fprintf('Duration %g complete\n',durations(ii));
end
writetable(array2table(durationData,'VariableNames',{'T','method','F_preparation', ...
    'F_BIC_preparation','F_after_hold','peak_atomic_CD'}),fullfile(out,'duration_benchmark.csv'));
report.durationMethods=durationModes;
save(fullfile(out,'duration_final_states.mat'),'durationPrepared','durationFinal', ...
    'durationData','report','-v7');

% Carrier convergence is tested only for passage from |eg>, with all
% laboratory micromotion retained and removed before comparing the states.
frequencies=[4 6 8 12 16];
freqData=zeros(numel(frequencies),5);
for ii=1:numel(frequencies)
    c=floquet_controls(u0,T,frequencies(ii));
    p0=zeros(model.dimension,1); p0(1)=1;
    labOpts=odeset(opts,'MaxStep',min(0.01,2*pi/frequencies(ii)/50));
    [tl,yl]=ode113(@(t,p) floquet_lab_rhs(t,p,model,c,cd),linspace(0,T,401),p0,labOpts);
    slow=micromotion_to_rotating_frame(tl.',transpose(yl),c);
    eff=nominal.atomic.psi(:,1:401);
    f=abs(slow(1,:)+slow(2,:)).^2/2;
    err=max(abs(f-nominal.atomic.F(1:401)));
    infid=max(abs(1-abs(sum(conj(eff).*slow,1)).^2));
    freqData(ii,:)=[frequencies(ii) f(end) err infid max(abs(sum(abs(slow).^2,1)-1))];
    fprintf('Carrier nu/xi=%g: max Bell error %.3g\n',frequencies(ii),err);
end
writetable(array2table(freqData,'VariableNames',{'nu_over_xi','F_Bell_final', ...
    'max_Bell_error','max_state_infidelity','norm_error'}),fullfile(out,'carrier_convergence.csv'));
report.carrierConvergence=freqData;
report.atomicPreparationF=nominal.atomic.prepF;
report.atomicPreparationFB=nominal.atomic.prepFB;
report.atomicHoldMinF=min(nominal.atomic.holdF);
report.analyticZ=1/(1+(0.1*u0)^2/2);
report.residualBound=atan(0.1*u0/sqrt(2))^2;
report.maxNominalEigenResidual=max(nominal.atomic.eigenResidual);
report.dressedMaximumTrackingError=max(abs(1-nominal.dressed.FB));
save(fullfile(out,'revision_summary.mat'),'report','-v7');
write_json(fullfile(out,'revision_summary.json'),report);
end

function r=simulate(model,u0,T,holdTime,mode,opts,keepStates)
ctrl=floquet_controls(u0,T,8);
p0=zeros(model.dimension,1); p0(1)=1;
tp=linspace(0,T,401);
[t,y]=ode113(@(t,p) passage_rhs(t,p,model,ctrl,mode),tp,p0,opts);
ph=transpose(y); pre=ph(:,end);
holdCtrl=ctrl; holdCtrl.thetaDot=@(t) 0;
if strcmp(mode,'isolated'), holdCtrl.u=@(t) zero_couplings(t); end
if holdTime>0
    [th,yh]=ode113(@(t,p) floquet_effective_rhs(t,p,model,holdCtrl), ...
        linspace(T,T+holdTime,501),pre,opts);
    ph=[ph transpose(yh(2:end,:))];
    t=[t;th(2:end)];
end
obs=observables_single_excitation(ph);
theta=ctrl.theta(t.');
FB=zeros(1,numel(t)); er=zeros(size(FB)); res=zeros(size(FB));
for jj=1:numel(t)
    [b,~,info]=compact_bic_branch(model,u0,theta(jj));
    FB(jj)=abs(b'*ph(:,jj))^2;
    if keepStates
        er(jj)=norm(floquet_effective_rhs(t(jj),b,model,ctrl));
        res(jj)=abs(info.etaPrime*ctrl.thetaDot(t(jj)));
    end
end
% The inverse transform is unitary; sites 0:n+2 define the bound region.
realspace=fft(ph(3:end,:),[],1)/sqrt(model.Nc);
outside=sum(abs(realspace(model.n1+4:end,:)).^2,1);
r=struct('time',t.','F',obs.bellPlusFidelity,'C',obs.concurrence, ...
    'FB',FB,'Patom',obs.atomicPopulation,'Pphoton',obs.photonicPopulation, ...
    'Poutside',outside,'prepF',obs.bellPlusFidelity(401),'prepFB',FB(401), ...
    'holdF',obs.bellPlusFidelity(401:end),'normError',max(abs(obs.norm-1)), ...
    'preparedState',pre,'finalState',ph(:,end),'eigenResidual',er,'trackingResidual',res);
if keepStates, r.psi=ph; end
end

function dp=passage_rhs(t,p,model,ctrl,mode)
c=ctrl;
scale=1;
if strcmp(mode,'fixed'), c.u=@(t) fixed_couplings(ctrl.u0,t); end
if ismember(mode,{'attach','isolated'}), c.u=@(t) zero_couplings(t); end
if strcmp(mode,'none'), scale=0; end
dp=floquet_effective_rhs(t,p,model,c,struct('counterdiabaticScale',scale));
if strcmp(mode,'dressed')
    [~,~,info]=compact_bic_branch(model,ctrl.u0,ctrl.theta(t));
    etaDot=info.etaPrime*ctrl.thetaDot(t);
    % -i H_aux acting on psi, with H_aux=i etaDot(|p><D|-|D><p|).
    dp=dp+etaDot*(info.p*(info.D'*p)-info.D*(info.p'*p));
end
end

function [a,b]=fixed_couplings(u0,t)
a=u0/sqrt(2)+0*t; b=a;
end
function [a,b]=zero_couplings(t)
a=0*t; b=a;
end
function write_json(path,value)
fid=fopen(path,'w'); assert(fid>=0);
closer=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end

function save_nominal(out,nominal,report)
names=fieldnames(nominal);
for ii=1:numel(names)
    name=names{ii};
    trajectory=nominal.(name);
    save(fullfile(out,['nominal_' name '_states.mat']),'trajectory','report','-v7');
    nominal.(name)=rmfield(nominal.(name),'psi');
end
save(fullfile(out,'nominal_trajectories.mat'),'nominal','report','-v7');
end
