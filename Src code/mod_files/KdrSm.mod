: Delayed Rectifier Potassium Channels - Soma



NEURON {
	SUFFIX KdrSm

	NONSPECIFIC_CURRENT ikdr

	RANGE gkdrbar, gkdr, ek
	RANGE n_inf
	RANGE tau_n
	RANGE theta_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gkdrbar	= 0.80048	(mho/cm2)
	tn	= 5	(ms)

	theta_n = 28
	kappa_n = -15

	ek      =-80	(mV)
	celsius = 36    (degC)

	vtraub	= -10	(mV)
}

STATE {
	n
}

ASSIGNED {
	dt      (ms)
	v       (mV)

	ikdr    (mA/cm2)
	gkdr    (mho/cm2)

	n_inf

	tau_n	(ms)

	tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	gkdr = gkdrbar * n*n*n*n
	ikdr = gkdrbar * n*n*n*n * (v - ek)
}

DERIVATIVE states {   
        
	evaluate_fct(v)

	n' = (n_inf - n) / tau_n

}

UNITSOFF
INITIAL {

	:  Q10 was assumed to be 3
	tadj = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)

	n = n_inf

}

PROCEDURE evaluate_fct(v(mV)) { LOCAL v2, v22, a, b

	: convert to traub convention
	v2 = v - vtraub 
	
	tau_n = (tn * tadj) / (Exp((v2+40)/40)+Exp(-(v2+40)/50))
	n_inf = 1 / (1+Exp((v2+theta_n)/kappa_n))

}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}
