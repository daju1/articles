debugmode(true);

assume(ra > 0);
assume(rq > 0);
assume(aq > 0);
assume(Rq > 0);
assume(rho0 > 0);

zq(rq, theta_q) := rq * cos(theta_q);
za(ra, theta_a) := ra * cos(theta_a);
rcq(rq, theta_q) := rq*sin(theta_q);
rca(ra, theta_a) := ra*sin(theta_a);
R (ra, theta_a, rq, theta_q, phi_q) := sqrt(rca(ra, theta_a)^2+rcq(rq, theta_q)^2+(za(ra, theta_a)-zq(rq, theta_q))^2-2*rca(ra, theta_a)*rcq(rq, theta_q)*cos(phi_q));

I1    (ra, theta_a, rq, theta_q) := quad_qags(1/R(ra, theta_a, rq, theta_q, phi_q), phi_q, 0, 2 * %pi);
rcqa2 (ra, theta_a, rq, theta_q) := (rcq(rq,theta_q)-rca(ra,theta_a))^2+(zq(rq, theta_q)-za(ra, theta_a))^2;
module(ra, theta_a, rq, theta_q) := - 4*rcq(rq,theta_q)*rca(ra,theta_a) / rcqa2(ra, theta_a, rq, theta_q);
Iphi  (ra, theta_a, rq, theta_q) := 4*elliptic_kc(module(ra, theta_a, rq, theta_q)) / sqrt(rcqa2(ra, theta_a, rq, theta_q));

I2    (ra, theta_a, rq)          := quad_qags(sin(theta_q) * I1(ra, theta_a, rq, theta_q), theta_q, 0, %pi);

rho_q (rho0, Rq, aq, r)          := rho0 / (1+exp((rq-Rq)/aq));
I3(rho0, Rq, aq, ra, theta_a) := quad_qags(rho_q(rho0, Rq, aq, rq) * rq^2 * I2(ra, theta_a, rq), rq, 0, inf);

I4(rho0, Rq, aq, ra, theta_a) := 2 * %pi * I3(rho0, Rq, aq, ra, theta_a);
I5(rho0, Rq, aq, ra)          := quad_qags(sin(theta_a) * I4(rho0, Rq, aq, ra, theta_a), theta_a, 0, %pi);
I6(rho0, Rq, aq)              := quad_qags(rho_q(rho0, Rq, aq, ra) * ra^2 * I5(rho0, Rq, aq, ra), ra, 0, inf);
I6(1.0, 1.0, 1.0);