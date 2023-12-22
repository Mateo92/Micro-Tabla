function c0=init_vel_KE(KE,H_ext,rho_ext,Di_int,H_int,rho_int,Di_v)
r_int=Di_int/2;
r_w=Di_v/2;
H_T=H_ext+H_int;
if r_int<=r_w
    c0= sqrt(KE/(pi*((r_w^2/2)*(rho_int*H_T))));
end

if r_int>r_w
    c0= sqrt(KE/(pi*((r_int^2/2)*(rho_int*H_T)+((r_w-r_int)^2/2)*(H_ext*rho_ext))));
end