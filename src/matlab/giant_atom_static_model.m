function model = giant_atom_static_model(params)
%GIANT_ATOM_STATIC_MODEL Build the full-Brillouin-zone giant-atom model.
%
% The Floquet project uses the full-zone convention requested in the
% execution brief:
%   k_m = -pi + 2*pi*m/Nc,  m = 0,...,Nc-1.
%
% The single-excitation basis is:
%   |e,g,0>, |g,e,0>, |g,g,1_k0>, ..., |g,g,1_k(Nc-1)>.

arguments
    params struct = struct()
end

params = with_default(params, 'xi', 1.0);
params = with_default(params, 'wc', 0.0);
params = with_default(params, 'g', 0.1 * params.xi);
params = with_default(params, 'Nc', 4 * 501);
params = with_default(params, 'n1', 6);
params = with_default(params, 'n2', 6);
params = with_default(params, 'x1', 0);
params = with_default(params, 'Dx', 2);
params = with_default(params, 'Omega0', params.wc);

Nc = params.Nc;
k = -pi + 2*pi*(0:Nc-1).' / Nc;
wk = params.wc - 2 * params.xi * cos(k);

x1 = params.x1;
x2 = params.x1 + params.Dx;

g1k = giant_atom_coupling(k, params.g, Nc, x1, params.n1);
g2k = giant_atom_coupling(k, params.g, Nc, x2, params.n2);

model = params;
model.k = k;
model.wk = wk;
model.x2 = x2;
model.g1k = g1k;
model.g2k = g2k;
model.dimension = Nc + 2;
model.kConvention = "full Brillouin zone [-pi, pi)";
end

function params = with_default(params, name, value)
if ~isfield(params, name)
    params.(name) = value;
end
end

function gik = giant_atom_coupling(k, g, Nc, atomPosition, ni)
formFactor = cos(k * ni / 2);
phase = exp(-1i * k * (atomPosition + ni/2));
gik = (2 * g / sqrt(Nc)) .* formFactor .* phase;
end
