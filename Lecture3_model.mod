// Solving the basic new keynesian model
// Lecture 3

var pi y r i y_n y_gap a d;
varexo epsilon_v epsilon_a epsilon_d;
parameters beta sigma phi theta varphi_y varphi_pi zeta rho_a rho_d;

// Parameters
beta = 0.99;            // Discount factor
sigma = 3.5;            // Intertemporal Elasticity
phi = 1;                // Frisch Elasticity
theta = 2/3;            // Slope of the Philips Curve
varphi_y = 0.125;       // Persistence of demand shock
varphi_pi = 2.5;        // Persistence of monetary policy shock
zeta = 0.9;             // Persistence of cost-push shock
rho_a = 0.9;            // Taylor rule coefficient on inflation
rho_d = 0.5;            // Taylor rule coefficient on output gap
    
// EQ Conditions
model;
// Newy Keynesian Philips Curve
pi = beta*pi(+1)+(((1-theta)*(1-beta*theta))/theta)*(sigma+phi)*y;
// IS Curve
y = y(+1)-(1/sigma)*(i-pi(+1)-r)+d;
// Natural rate of interest
r = sigma*((1+phi)/(sigma+phi))*(a(+1)-a);
// Monetary Policy rule
i = zeta*i(-1)+(1-zeta)*(varphi_y*y+varphi_pi*pi)+epsilon_v;
// Output gap
y_gap = y - y_n;
// Output
y_n = ((1+phi)/(sigma+phi))*a;
// Output and cost-push shock
a = rho_a*a(-1)+epsilon_a;
d = rho_d*d(-1)+epsilon_d;
end;

initval;
pi = 0;
y = 0;
y_n = 0;
r = 0;
i = 0;
y_gap = 0;
a = 0;
d = 0;
epsilon_v = 0;
epsilon_a = 0;
epsilon_d = 0;
end;

// Dynare
check;
steady;

// Bayesian Estimate
estimated_params;
sigma, normal_pdf, 1.5, 0.5;
beta, beta_pdf, 0.5, 0.2;
varphi_pi, normal_pdf, 1.5, 0.2;
rho_a, beta_pdf, 0.5, 0.2;
rho_d, beta_pdf, 0.2, 0.2;

// Shock standard errors
stderr epsilon_a, inv_gamma_pdf, 0.01, 2;
stderr epsilon_d, inv_gamma_pdf, 0.01, 2;
stderr epsilon_v, inv_gamma_pdf, 0.01, 2;

end; 

stoch_simul(order=1,irf=20);

options_.hessian_type = 1;

// Observed variables
varobs y pi;

// 
estimation(datafile = nk_data,mh_nblocks=2,mh_replic=100000,mh_drop=0.5,mode_compute=6,mode_check,bayesian_irf,graph_format=pdf);

// Trace plots
generate_trace_plots([1 2]);

// Shock decomposition
shock_decomposition(nograph) y pi i;
plot_shock_decomposition(nodisplay, graph_format=fig) y pi i;

