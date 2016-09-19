function val=computeres(weaklearner,ftrvalue)

p0 = exp(weaklearner.e0*(ftrvalue-weaklearner.mu0)*(ftrvalue-weaklearner.mu0))*weaklearner.n0;
p1 = exp(weaklearner.e1*(ftrvalue-weaklearner.mu1)*(ftrvalue-weaklearner.mu1))*weaklearner.n1;

val = log(p1+0.00001)-log(p0+0.00001);