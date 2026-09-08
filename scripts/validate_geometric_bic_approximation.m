function report=validate_geometric_bic_approximation(runMode)
%VALIDATE_GEOMETRIC_BIC_APPROXIMATION Test the recovered analytical conditions.
% Separate exact conditional-atomic identities, weak photonic dressing, and
% the frozen-envelope approximation in a time-dependent memory equation.
if nargin<1, runMode='paper'; end
runMode=validatestring(runMode,{'paper','dev'},mfilename,'runMode');
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','matlab'));
paper=strcmp(runMode,'paper');
if paper, nc=2004; out=fullfile(root,'data','prr_revision','geometric_appendix');
else, nc=204; out=fullfile(root,'data','development','geometric_appendix'); end
if ~exist(out,'dir'), mkdir(out); end
u0=0.9;

% Independent real-space and momentum-space reconstructions of the eigenstate.
rows=[];
for n2=[6 10]
    m=giant_atom_static_model(struct('Nc',nc,'n1',6,'n2',n2));
    for lambda=linspace(0,2,41)
        u=u0*[lambda*n2 m.n1]/hypot(lambda*n2,m.n1);
        [B,info]=geometric_bic_state(m,u);
        ctrl=floquet_controls(); ctrl.u=@(t) fixed_u(t,u);
        residual=norm(floquet_effective_rhs(0,B,m,ctrl));
        numerator=u(1)*conj(m.g1k)*info.a(1)+u(2)*conj(m.g2k)*info.a(2);
        denom=-m.wk; regular=abs(denom)>1e-12;
        phi=zeros(m.Nc,1); phi(regular)=numerator(regular)./denom(regular);
        recovered=[info.a;phi]; recovered=recovered/norm(recovered);
        ca=recovered(1:2); z=sum(abs(ca).^2);
        rows=[rows; n2 lambda u z 2*abs(prod(ca))/z 2*abs(prod(ca)) ...
            residual norm(recovered-B)]; %#ok<AGROW>
    end
end
writetable(array2table(rows,'VariableNames',{'n2','lambda','u1','u2','Z', ...
    'C_conditional','C_unconditional','eigen_residual','reconstruction_error'}), ...
    fullfile(out,'amplitude_rule.csv'));

% Root-resolved radiation zeros: derivative-selected versus opposite channel.
m=giant_atom_static_model(struct('Nc',nc)); u=[u0 u0]/sqrt(2);
delta=logspace(-5,-1,101).'; radiation=zeros(numel(delta),5);
for ii=1:numel(delta)
    values=zeros(1,4); idx=0;
    for rootSign=[1 -1]
        k=rootSign*pi/2+delta(ii);
        v1=u(1)*(1+exp(1i*k*m.n1));
        v2=u(2)*exp(1i*k*m.Dx)*(1+exp(1i*k*m.n2));
        for sign=[1 -1]
            idx=idx+1; values(idx)=abs((v1+sign*v2)/(sqrt(2)*u0))^2;
        end
    end
    radiation(ii,:)=[delta(ii) values];
end
writetable(array2table(radiation,'VariableNames',{'delta_k','selected_plus_root', ...
    'other_plus_root','selected_minus_root','other_minus_root'}), ...
    fullfile(out,'radiation_zero.csv'));

dressing=[];
for n2=[6 10]
    for gg=logspace(-2,log10(0.4),41)
        m=giant_atom_static_model(struct('Nc',nc,'n1',6,'n2',n2,'g',gg));
        u=u0*[n2 m.n1]/hypot(n2,m.n1);
        [B,info]=geometric_bic_state(m,u);
        aVac=[info.a;zeros(nc,1)];
        dressing=[dressing; n2 gg info.Z 1-abs(aVac'*B)^2]; %#ok<AGROW>
    end
end
writetable(array2table(dressing,'VariableNames',{'n2','g_over_xi','Z','vacuum_state_infidelity'}), ...
    fullfile(out,'dressing_error.csv'));

% The frozen memory equation is integrated without assigning it a physical
% norm. Compare its atomic amplitudes directly with the unitary reference.
if paper, couplings=[0.01 0.02 0.05 0.1 0.2 0.4]; durations=[20 80 160];
else, couplings=[0.05 0.1]; durations=[20 80]; end
opts=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.04);
memory=[]; traces=struct([]); count=0;
for T=durations
    for gg=couplings
        count=count+1;
        m=giant_atom_static_model(struct('Nc',nc,'g',gg));
        ctrl=floquet_controls(u0,T,8);
        p0=zeros(m.dimension,1); p0(1)=1;
        y0=zeros(2+2*nc,1); y0(1)=1;
        tt=linspace(0,T,801);
        [t,p]=ode113(@(t,p) floquet_effective_rhs(t,p,m,ctrl, ...
            struct('counterdiabaticScale',1)),tt,p0,opts);
        [~,y]=ode113(@(t,y) frozen_envelope_rhs(t,y,m,ctrl,1),tt,y0,opts);
        err=vecnorm(p(:,1:2)-y(:,1:2),2,2);
        exactF=abs(p(:,1)+p(:,2)).^2/2;
        approxF=abs(y(:,1)+y(:,2)).^2/2;
        memory=[memory; gg T max(err) err(end) max(abs(exactF-approxF)) ...
            max(abs(sum(abs(p).^2,2)-1)) max(sum(abs(y(:,1:2)).^2,2))]; %#ok<AGROW>
        traces(count).g=gg;traces(count).T=T;traces(count).time=t;
        traces(count).exactAtomic=p(:,1:2);traces(count).frozenAtomic=y(:,1:2);
        fprintf('Frozen envelope g=%.3f T=%g: max atomic error %.6g\n',gg,T,max(err));
    end
end
writetable(array2table(memory,'VariableNames',{'g_over_xi','T','max_atomic_error', ...
    'final_atomic_error','max_Bell_overlap_difference','exact_norm_error', ...
    'maximum_frozen_atomic_norm_squared'}),fullfile(out,'frozen_memory.csv'));
report=struct('Nc',nc,'u0',u0,'RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.04, ...
    'maximumEigenResidual',max(rows(:,8)),'maximumReconstructionError',max(rows(:,9)), ...
    'maximumConditionalConcurrenceError',max(abs(rows(:,6)-2*rows(:,2)./(1+rows(:,2).^2))), ...
    'scope','band-centre static states and vacuum-input CD passage in the effective model');
save(fullfile(out,'geometric_appendix.mat'),'rows','radiation','dressing','memory','traces','report','-v7');
fid=fopen(fullfile(out,'summary.json'),'w');assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
end

function [u1,u2]=fixed_u(t,u)
u1=u(1)+0*t;u2=u(2)+0*t;
end
