config.dim = 2;
config.r = 0.1;

model = PositionalObsModel(config);

xkm1 = [0; 0; sqrt(2); pi/4];
pkm1 = [0 1; 0 1; sqrt(2) sqrt(2); pi/4 pi/4];
Pkm1 = [0 1 0 1; 0 1 0 1; 0 1 0 1; 1 1 0 1];
yk = model.transform_mean(xkm1);
pk = model.transform_parts(pkm1, false);
Pk = model.transform_cov(Pkm1);