function [r_dc_out]=eclipse(r_dc,env,rsun)

R_mod=zeros(length(r_dc(:,1)),1);
rsun_=-rsun;
for i=1:length(r_dc(:,1))
    R_mod(i,1)=norm(r_dc(i,1:3));
end
beta_ecl=asin(env.Rp ./R_mod);
ecl=zeros(length(r_dc(:,1)),1);
for i=1:length(r_dc(:,1))
    ecl(i,1)=acos(dot(rsun_,r_dc(i,1:3))/(norm(rsun_)*norm(r_dc(i,1:3))));
end
index=zeros(length(r_dc(:,1)),1);

for i=1:length(r_dc(:,1))
    if ecl(i,1)<=beta_ecl(i,1)                                              %1 se è eclisse
        index(i,1)=1;
    else
        index(i,1)=0;
    end
end
r_dc_out=[r_dc index];
end