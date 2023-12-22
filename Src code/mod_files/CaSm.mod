: TITLE Motoneuron Soma channels
: Calcium channels + Calcium Dynamics - Soma


NEURON {
	SUFFIX CaSm

	NONSPECIFIC_CURRENT ikca
	NONSPECIFIC_CURRENT ican
	NONSPECIFIC_CURRENT ical

	RANGE gkcabar, gcanbar, gcalbar, eca
	RANGE gkca, gcan, gcal
	RANGE mn_inf, hn_inf, ml_inf
	RANGE tau_mn, tau_hn, tau_ml
	RANGE nexp, kd
	RANGE f, kca, alpha
	RANGE thetamn, thetahn
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)

	FARADAY	= (faraday) (coulomb)
	R	= (k-mole) (joule/degC)
}

PARAMETER {

	: Calcium N-type channels
	gcanbar = 0.072837  (mho/cm2)	
	tmn	= 15	    (ms)
	thetamn	= 22	    (mV)
	thn	= 50	    (ms)
	thetahn	= 40	    (mV)

	: Calcium L-type Channels
	gcalbar = 0.0002    (mho/cm2)	
	mlexp	= 1
	tml	= 400	    (ms)
	thetaml = 45.8	    (mV)
	kml	= -3.7	    (mV)

	: Calcium-activated Potassium Channels
	gkcabar = 0.37418   (mho/cm2)	                    
	nexp	= 2
	kd	= 0.0002

	: Calcium Dynamics
	cao	= 2	    (mM)
	caio	= .0001     (mM)
	f	= 0.01
	alpha	= 1
	kca	= 4

	: General
	vtraub	= -10	    (mV)
	celsius = 36        (degC)
	ek= -80             (mV)
}

STATE {
	mn hn ml cai
}

ASSIGNED {
	dt      (ms)
	v       (mV)
	eca	(mV)
	
	ican	(mA/cm2)
	ical	(mA/cm2)
	ikca	(mA/cm2)

	gkca	(mho/cm2)
	gcan	(mho/cm2)
	gcal	(mho/cm2)
	
	mn_inf
	hn_inf
	ml_inf
	
	tau_mn	(ms)
	tau_hn	(ms)
	tau_ml	(ms)
	
	tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	eca = ((1000 * R * (celsius + 273.15)) / (2 * FARADAY)) * log(cao/cai)

	gcan = gcanbar * mn*mn*hn
	ican = gcan * (v - eca)

	gcal = gcalbar * (ml^mlexp)
	ical = gcal * (v - eca)

	gkca = gkcabar * ( (cai^nexp) / ((cai^nexp)+kd) )
	ikca = gkca * (v - ek)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations

        evaluate_fct(v)

	mn' = (mn_inf - mn) / tau_mn
	hn' = (hn_inf - hn) / tau_hn
	ml' = (ml_inf - ml) / tau_ml

	cai' = f*(-(alpha*(ican+ical))-(kca*cai))
}

UNITSOFF
INITIAL {

	:  Q10 was assumed to be 3
	tadj = 3.0 ^ ((celsius-36)/ 10 )

	cai = caio

	evaluate_fct(v)

	mn = mn_inf
	hn = hn_inf
	ml = ml_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL v2

	v2 = v - vtraub : convert to traub convention
	
	tau_mn = tmn * tadj
	mn_inf = 1 / (1+exp((v2+thetamn)/-5))
	
	tau_hn = thn * tadj
	hn_inf = 1 / (1+exp((v2+thetahn)/5))
	
	tau_ml = tml * tadj
	ml_inf = 1 / (1+exp((v2+thetaml)/kml))
}

