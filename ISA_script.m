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

N_layer=100;
h_G_row=nan(1,7*N_layer);
T_row=nan(1,7*N_layer);
p_row=nan(1,7*N_layer);
for n=1:4
    n_grad=2*n-1;
    n_vec=2*(n-1)*N_layer+(1:N_layer);

    h_G_row(n_vec)=linspace(h_G0_row(n_grad),h_G0_row(n_grad+1),N_layer);
    T_row(n_vec)=T_0_row(n_grad)+a_0_row(n).*(h_G_row(n_vec)-h_G0_row(n_grad));
    p_row(n_vec)=p_0_row(n_grad).*(T_row(n_vec)./T_0_row(n_grad)).^(-g_0./a_0_row(n)./R);
end

for n=1:3
    n_grad=2*n-1;
    n_iso=2*n;
    n_vec=n_grad*N_layer+(1:N_layer);

    h_G_row(n_vec)=linspace(h_G0_row(n_iso),h_G0_row(n_iso+1),N_layer);
    T_row(n_vec)=repmat(T_0_row(n_iso),1,N_layer);
    p_row(n_vec)=p_0_row(n_iso).*exp(-g_0.*(h_G_row(n_vec)-h_G0_row(n_iso))./R./T_0_row(n_iso));
end

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
