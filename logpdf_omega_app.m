function logpdf=logpdf_omega_app(omega,sigma2_U,sigma2_A,deltas)

logpdf=sum(log(normpdf(omega(1,:),0,sqrt(deltas.*sigma2_U))))+...
       sum(log(normpdf(omega(2,:),0,sqrt(deltas.*sigma2_A))));