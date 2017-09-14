using PkgBenchmark, NLOptControl

@benchgroup "psMethods"  begin
    for scheme in (:lgrExplicit,:lgrImplicit)
        for Nck in ([10],[20],[30],[40],[50],[10,8,6],[12,10,8,6])

            # Bryson Denham, basic
            BrysonDenham_EXP=[:(x2[j]),:(u1[j])]; L=1/6;
            n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN]);
            dynamics!(n,BrysonDenham_EXP)
            configure!(n;(Nck=Nck),(:integrationScheme=>scheme),(:finalTimeDV=>false),(:tf=>1.0));
            obj=integrate!(n,:(0.5*u1[j]^2));
            @NLobjective(n.mdl,Min,obj);
            @bench "BrysonDenham_basic", string(scheme), Nck optimize!($n);

            # Bryson Denham, standard
            BrysonDenham_EXP=[:(x2[j]),:(u1[j]),:(0.5*u1[j]^2)];
            n=define(numStates=3,numControls=1,X0=[0.,1.,0.],XF=[0.,-1.,NaN],XL=[0.,-10.,-10.],XU=[1/9,10.,10.],CL=[-5000.],CU=[5000.]);
            n.s.tf_max=50.0;
            dynamics!(n,BrysonDenham_EXP)
            configure!(n;(Nck=Nck),(:integrationScheme=>scheme),(:finalTimeDV=>true),(:solverSettings=>(:name=>:KNITRO)));
            x1=n.r.x[:,1];x2=n.r.x[:,2];x3=n.r.x[:,3];u1=n.r.u[:,1];
            @NLobjective(n.mdl,Min,x3[end]);
            setvalue(n.tf, 0.5)
            for i in 1:Nck[1]; setvalue(x1[i], 0.0); setvalue(x2[i], 0.0); setvalue(x3[i], 0.0); setvalue(u1[i], 0.0); end
            @bench "BrysonDenham_standard", string(scheme), Nck optimize!($n);
        end
    end
end

@benchgroup "tmMethods"  begin
    for scheme in (:trapezoidal,:bkwEuler)
        for N in (10,20,30,40,50)

            # Bryson Denham, basic
            BrysonDenham_EXP=[:(x2[j]),:(u1[j])]; L=1/6;
            n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN]);
            dynamics!(n,BrysonDenham_EXP)
            configure!(n;(N=N),(:integrationScheme=>scheme),(:finalTimeDV=>false),(:tf=>1.0));
            obj=integrate!(n,:(0.5*u1[j]^2));
            @NLobjective(n.mdl,Min,obj);
            @bench "BrysonDenham_basic", string(scheme), N optimize!($n);

            # Bryson Denham, standard
            BrysonDenham_EXP=[:(x2[j]),:(u1[j]),:(0.5*u1[j]^2)];
            n=define(numStates=3,numControls=1,X0=[0.,1.,0.],XF=[0.,-1.,NaN],XL=[0.,-10.,-10.],XU=[1/9,10.,10.],CL=[-5000.],CU=[5000.]);
            n.s.tf_max=50.0;
            dynamics!(n,BrysonDenham_EXP)
            configure!(n;(N=N),(:integrationScheme=>scheme),(:finalTimeDV=>true),(:solverSettings=>(:name=>:KNITRO)));
            x1=n.r.x[:,1];x2=n.r.x[:,2];x3=n.r.x[:,3];u1=n.r.u[:,1];
            @NLobjective(n.mdl,Min,x3[end]);
            setvalue(n.tf, 0.5)
            for i in 1:N; setvalue(x1[i], 0.0); setvalue(x2[i], 0.0); setvalue(x3[i], 0.0); setvalue(u1[i], 0.0); end
            @bench "BrysonDenham_standard", string(scheme), N optimize!($n);

        end
    end
end
