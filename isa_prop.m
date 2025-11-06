function [h_arr,T_arr,p_arr,rho_arr,a_arr]=isa_prop(h_G_arr)

g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11000,25000,47000,53000,79000,90000,105000];
T_0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p_0_row=[101330,22633,2488.7,120.45,58.323,1.0095,0.10444];
a_0_row=[-.0065,.003,-.0045,.004];

if any(h_G_arr<h_G0_row(1) | h_G_arr>h_G0_row(end))
    error("h_G_arr elements must be within ["+h_G0_row(1)+','+h_G0_row(1)+']!');
end

sz_vec=size(h_G_arr);
T_arr=nan(sz_vec);
p_arr=nan(sz_vec);

N_layers=length(h_G0_row)-1;
for n_layer=1:N_layers
    ind_vec=find(h_G_arr>h_G0_row(n_layer) & h_G_arr<=h_G0_row(n_layer+1));
    if ~isempty(ind_vec)
        if mod(n_layer,2)~=0  %n_layer is odd ==> gradient layer
            a_0=a_0_row((n_layer+1)/2);
            T_arr(ind_vec)=T_0_row(n_layer)+a_0.*(h_G_arr(ind_vec)-h_G0_row(n_layer));
            p_arr(ind_vec)=p_0_row(n_layer).*(T_arr(ind_vec)./T_0_row(n_layer)).^(-g_0./a_0./R);
        else    %n_layer is even ==> isothermal layer
            T_arr(ind_vec)=T_0_row(n_layer);
            p_arr(ind_vec)=p_0_row(n_layer).*exp(-g_0.*(h_G_arr(ind_vec)-h_G0_row(n_layer))./R./T_0_row(n_layer));
        end        
    end
end

h_arr=r.*h_G_arr./(r+h_G_arr);
rho_arr=p_arr./R./T_arr;
a_arr=sqrt(gamma.*R.*T_arr);