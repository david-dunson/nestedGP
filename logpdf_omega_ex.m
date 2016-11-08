function logpdf=logpdf_omega_ex(omega,sigma2_U,sigma2_A,deltas)

J=length(deltas);

SIGMA=zeros(3,3,J);

SIGMA(1,1,:)=deltas.^3./3*sigma2_U+deltas.^5./20*sigma2_A;
SIGMA(1,2,:)=deltas.^2./2*sigma2_U+deltas.^4./8*sigma2_A;
SIGMA(1,3,:)= deltas.^3./6*sigma2_A;
SIGMA(2,2,:)= deltas*sigma2_U+deltas.^3./3*sigma2_A;
SIGMA(2,3,:)= deltas.^2./2*sigma2_A;
SIGMA(3,3,:)= deltas*sigma2_A;
SIGMA(2,1,:)=SIGMA(1,2,:);
SIGMA(3,1,:)=SIGMA(1,3,:);
SIGMA(3,2,:)=SIGMA(2,3,:);

%slower version:
%logpdf=sum(log(mvnpdf(omega',zeros(1,3),SIGMA)));

%faster version:
tmp=zeros(1,J);
for j=1:J
    tmp(j)= mvnormpdfln(omega(:,j), [], [], SIGMA(:,:,j));
end

logpdf=sum(tmp);
