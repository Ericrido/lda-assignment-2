libname LDA '/folders/myfolders';


/* Dataset preparation */

data lda.acu2;
set lda.acu2;
by id;
prev1 = lag(frequency);
prev2= 0;
if first.id then prev1 = 0;
if last.id then prev2 = lag2(frequency);
if frequency=0 then frequency=".";
run;

/* TRANSITIONAL MODEL */

/* First order autoregressive */
proc genmod data=lda.acu2 descending;
model frequency = group chronicity age time group*time prev1 / dist=poisson;
run;

/* Second order autoregressive */
proc genmod data=lda.acu2 descending;
model frequency = group chronicity age time group*time prev1 prev2 / dist=poisson;
run;



/* CONDITIONAL MODELS */

/* Standard GEE */

proc genmod data=lda.acu2 descending;
class id timeclass;
model frequency = age chronicity group time group*time
/ dist=poisson;
repeated subject=id / withinsubject=timeclass
type=exch covb corrw modelse;
run;

/* Linearization based model */

%glimmix(data=test, procopt=%str(method=ml empirical)
,
stmts=%str(
class idnum timeclss;
model onyresp = treatn time treatn*time / solution;
repeated timeclss / subject=idnum type=cs rcorr;
),
error=binomial,
link=logit);

proc glimmix data=lda.acu2 method=RSPL empirical;
class id;
model frequency = age chronicity group time group*time
/ dist=poisson solution;
random _residual_ / subject=id type=cs;
run;


/* GENERALIZED LINEAR MIXED MODELS */

/* PQL for initial values, then adaptive quadrature */

proc glimmix data=lda.acu2 method=rspl;
class id;
model frequency = age chronicity group time group*time
/ dist=poisson solution;
random intercept / subject=id;
run;



proc nlmixed data=lda.acu2 ;
parms beta0=2.77 beta1=-0.003 beta2=-0.003 beta3=-0.07 beta4=-0.01 beta5=-0.01;
teta=beta0 + b + beta1*age + beta2*chronicity + beta3*group + beta4*time + beta5*group*time;
lambda=exp(teta);
model frequency ~ poisson(lambda);
random b ~ normal(0,sigma**2) subject=id;
run;
