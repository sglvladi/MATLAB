config.q = 0.01;
config.dim = 2;

CVmodel = ConstantVelocityModel(config);

xkm1 = [0; 0; 1; 1];
pkm1 = [zeros(2,1000);ones(2,1000)];
Pkm1 = [0 1 0 1; 0 1 0 1; 0 1 0 1; 1 1 0 1];
xk = CVmodel.sys(0, xkm1);
pk = CVmodel.sys(0, pkm1, CVmodel.sys_noise(1, size(pkm1,2)));
Pk = CVmodel.sys_cov(0, Pkm1);