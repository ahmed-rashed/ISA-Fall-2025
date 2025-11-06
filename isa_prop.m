function [h_arr,T_arr,p_arr,rho_arr,a_arr]=isa_prop(h_G_arr)

g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11000,25000,47000,53000,79000,90000,105000];
T_0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p_0_row=[101330,22633,2488.7,120.45,58.323,1.0095,0.10444];
a_0_row=[-.0065,.003,-.0045,.004];

sz_vec=size(h_G_arr);
T_arr=nan(sz_vec);
p_arr=nan(sz_vec);
n=1;
for h_G=h_G_arr(:).'
    if h_G<h_G0_row(1)
        error("h_G_arr elements must be higher than "+h_G0_row(1)+'!');
    elseif h_G<=h_G0_row(2) % 1st layer
        T_arr(n)=T_0_row(1)+a_0_row(1).*(h_G-h_G0_row(1));
        p_arr(n)=p_0_row(1).*(T_arr(n)./T_0_row(1)).^(-g_0./a_0_row(1)./R);
    elseif h_G<=h_G0_row(3) % 2nd layer
        T_arr(n)=T_0_row(2);
        p_arr(n)=p_0_row(2).*exp(-g_0.*(h_G-h_G0_row(2))./R./T_0_row(2));
    elseif h_G<=h_G0_row(4) % 3rd layer
        T_arr(n)=T_0_row(3)+a_0_row(2).*(h_G-h_G0_row(3));
        p_arr(n)=p_0_row(3).*(T_arr(n)./T_0_row(3)).^(-g_0./a_0_row(2)./R);
    elseif h_G<=h_G0_row(5) % 4th layer
        T_arr(n)=T_0_row(4);
        p_arr(n)=p_0_row(4).*exp(-g_0.*(h_G-h_G0_row(4))./R./T_0_row(4));
    elseif h_G<=h_G0_row(6) % 5th layer
        T_arr(n)=T_0_row(5)+a_0_row(3).*(h_G-h_G0_row(5));
        p_arr(n)=p_0_row(5).*(T_arr(n)./T_0_row(5)).^(-g_0./a_0_row(3)./R);
    elseif h_G<=h_G0_row(7) % 6th layer
        T_arr(n)=T_0_row(6);
        p_arr(n)=p_0_row(6).*exp(-g_0.*(h_G-h_G0_row(6))./R./T_0_row(6));
    elseif h_G<=h_G0_row(8) % 7th layer
        T_arr(n)=T_0_row(7)+a_0_row(4).*(h_G-h_G0_row(7));
        p_arr(n)=p_0_row(7).*(T_arr(n)./T_0_row(7)).^(-g_0./a_0_row(4)./R);
    else
        error("h_G_arr elements must be less than "+h_G0_row(8)+'!');
    end
    
    n=n+1;
end

h_arr=r.*h_G_arr./(r+h_G_arr);
rho_arr=p_arr./R./T_arr;
a_arr=sqrt(gamma.*R.*T_arr);