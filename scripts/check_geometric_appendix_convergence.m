function check_geometric_appendix_convergence()
%CHECK_GEOMETRIC_APPENDIX_CONVERGENCE Finite-grid and solver checks for memory freezing.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src','matlab'));
out=fullfile(root,'data','prr_revision','geometric_appendix');
assert(exist(out,'dir')==7,'Run validate_geometric_bic_approximation first.');
rows=[];
settings=[804 1e-10 1e-12 0.04;2004 1e-10 1e-12 0.04; ...
    4004 1e-10 1e-12 0.04;2004 1e-11 1e-13 0.02];
for s=1:size(settings,1)
    nc=settings(s,1);
    opts=odeset('RelTol',settings(s,2),'AbsTol',settings(s,3),'MaxStep',settings(s,4));
    for T=[20 160]
        for gg=[0.1 0.4]
            m=giant_atom_static_model(struct('Nc',nc,'g',gg));
            ctrl=floquet_controls(0.9,T,8);
            p0=zeros(m.dimension,1);p0(1)=1;
            y0=zeros(2+2*nc,1);y0(1)=1;
            tt=linspace(0,T,801);
            [~,p]=ode113(@(t,p) floquet_effective_rhs(t,p,m,ctrl, ...
                struct('counterdiabaticScale',1)),tt,p0,opts);
            [~,y]=ode113(@(t,y) frozen_envelope_rhs(t,y,m,ctrl,1),tt,y0,opts);
            err=vecnorm(p(:,1:2)-y(:,1:2),2,2);
            rows=[rows; settings(s,:) gg T max(err) err(end)]; %#ok<AGROW>
        end
    end
    fprintf('Geometric appendix convergence setting %d/4 complete.\n',s);
end
writetable(array2table(rows,'VariableNames',{'Nc','RelTol','AbsTol','MaxStep', ...
    'g_over_xi','T','max_atomic_error','final_atomic_error'}),fullfile(out,'convergence.csv'));
end
