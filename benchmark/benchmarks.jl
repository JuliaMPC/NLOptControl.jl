using PkgBenchmark, NLOptControl

# TODO make a problems group and combine these with testing

@benchgroup "psMethods"  begin
    for scheme in (:lgrExplicit,:lgrImplicit)
        for Nck in ([10,8,6],[12,10,8,6])
            BrysonDenham_EXP=[:(x2[j]),:(u1[j])]; L=1/6;
            n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN]);
            dynamics!(n,BrysonDenham_EXP)
            configure!(n;(Nck=Nck),(:integrationScheme=>scheme),(:finalTimeDV=>false),(:tf=>1.0));
            obj=integrate!(n,:(0.5*u1[j]^2));
            @NLobjective(n.mdl,Min,obj);
            @bench string(scheme), Nck optimize!($n);
        end
    end
end

@benchgroup "tmMethods"  begin
    for scheme in (:trapezoidal,:bkwEuler)
        for N in (24,36)
            BrysonDenham_EXP=[:(x2[j]),:(u1[j])]; L=1/6;
            n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN]);
            dynamics!(n,BrysonDenham_EXP)
            configure!(n;(N=N),(:integrationScheme=>scheme),(:finalTimeDV=>false),(:tf=>1.0));
            obj=integrate!(n,:(0.5*u1[j]^2));
            @NLobjective(n.mdl,Min,obj);
            @bench string(scheme), N optimize!($n);
        end
    end
end
