function p = baseCaseParameters()
    p = struct();
    p.vaxnum = 1;
    p.T = 0;
    p.k = 0;
    p.numshot = 1;
    p.E1h=6.6;
    p.dE12=-0.6;
    p.p=0.5;
    p.masking=0;
    p.C0 = 0.2;
    p.w1 = 0;
    p.w2 = 0.5;
    p.steric = 0;
    p.memToGCFrac = 0;
    p.outputprob = 0.1;
    p.outputpcfrac = 0.1;
    p.rho = 0.95;
    p.pSER = 0;
    p.tmax=28;
    p.first = 1;
    p.last = 2000;
    p.numfrag = 10;
end