config.q_vel = 0.01;
config.q_head = 0.16;

CHmodel = ConstantHeadingModel(config);

xkm1 = [0; 0; sqrt(2); pi/4];
pkm1 = [0 1; 0 1; sqrt(2) sqrt(2); pi/4 pi/4];
Pkm1 = [0 1 0 1; 0 1 0 1; 0 1 0 1; 1 1 0 1];
xk = CHmodel.propagate_mean(xkm1,2);
pk = CHmodel.propagate_parts(pkm1,2);
Pk = CHmodel.propagate_cov(Pkm1,2);