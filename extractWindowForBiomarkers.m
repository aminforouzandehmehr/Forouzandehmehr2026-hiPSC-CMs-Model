function window = extractWindowForBiomarkers(V_ode,i_tot, Istim ,t,Cai,Nai,Frc,Lsrc,svl, wspan)



flag = 0;
[pks_MDP, idx_MDP] = findpeaks(-V_ode);

idx = find(pks_MDP<0);
pks_MDP(idx)=[];
idx_MDP(idx) = [];

if length(pks_MDP)<=10 
   disp('Error')
   flag =1;
   window=[];
   return
end
numAPsBiomarkers = 9;
N = length(V_ode);
window.Time   = t(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end))*1000;% ms';
window.V_ode  = V_ode(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));
window.i_tot  = i_tot(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));
window.Istim  = Istim(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));
window.Cai    = Cai(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));
window.Nai    = Nai(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));
window.Frc    = Frc(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));
window.Lsrc   = Lsrc(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));
window.svl    = svl(idx_MDP(end-numAPsBiomarkers)-wspan:idx_MDP(end));

