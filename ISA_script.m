clc
close all
clearvars

g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11000,25000,47000,53000,79000,90000,105000];
T_0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p_0_row=[101330,22633,2488.7,120.45,58.323,1.0095,0.10444];
a_0_row=[-.0065,.003,-.0045,.004];

N_layers=length(h_G0_row)-1;
N_layer=100;

h_G_rows=nan(N_layers,N_layer);
for n=1:N_layers
    h_G_rows(n,:)=linspace(h_G0_row(n),h_G0_row(n+1),N_layer);
end

T_rows=nan(N_layers,N_layer);
p_rows=nan(N_layers,N_layer);

n_grad_layers_vec=1:2:7;
T_rows(n_grad_layers_vec,:)=T_0_row(n_grad_layers_vec).'+a_0_row.'.*(h_G_rows(n_grad_layers_vec,:)-h_G0_row(n_grad_layers_vec).');
p_rows(n_grad_layers_vec,:)=p_0_row(n_grad_layers_vec).'.*(T_rows(n_grad_layers_vec,:)./T_0_row(n_grad_layers_vec).').^(-g_0./a_0_row.'./R);

n_iso_layers_vec=2:2:7;
T_rows(n_iso_layers_vec,:)=repmat(T_0_row(n_iso_layers_vec).',1,N_layer);
p_rows(n_iso_layers_vec,:)=p_0_row(n_iso_layers_vec).'.*exp(-g_0.*(h_G_rows(n_iso_layers_vec,:)-h_G0_row(n_iso_layers_vec).')./R./T_0_row(n_iso_layers_vec).');

h_G_row=reshape(h_G_rows.',1,[]);
T_row=reshape(T_rows.',1,[]);
p_row=reshape(p_rows.',1,[]);

h_row=r.*h_G_row./(r+h_G_row);
rho_row=p_row./R./T_row;
a_row=sqrt(gamma.*R.*T_row);

plot(h_row./1e3,h_G_row./1e3)
xlabel('h (km)')
ylabel('h_G (km)')

figure
tiledlayout(1,2)
nexttile
plot(T_row-273,h_G_row./1e3)
xlabel('T (C)')
ylabel('h_G (km)')

nexttile
plot(a_row./1e3.*3600,h_G_row)
xlabel('a (km/h)')
ylabel('h_G (km)')

figure
tiledlayout(1,2)
nexttile
plot(p_row./1e5,h_G_row./1e3)
xlabel('p (bar)')
ylabel('h_G (km)')

nexttile
plot(rho_row,h_G_row)
xlabel('rho (kg/m^3)')
ylabel('h_G (km)')
