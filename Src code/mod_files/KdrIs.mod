: Delayed Rectifier Potassium Channels - Initial Segment


NEURON {
	SUFFIX KdrIs

	NONSPECIFIC_CURRENT ik

	RANGE gkbar, gk, ek
	RANGE n_inf
	RANGE tau_n
	RANGE theta_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gkbar   = 0.16552 	  (mho/cm2)

	ek      = -80.0   (mV)

	theta_n = 28
	kappa_n = -15

	celsius = 20      (degC)

	vtraub2	= -10	  (mV)

	tn  = 5		  (ms)
}

STATE {
	n
}

ASSIGNED {
	dt      (ms)
	v       (mV)

	ik      (mA/cm2)
	gk	  (mho/cm2)

	n_inf

	tau_n	(ms)

	tadj3
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	gk = gkbar * n*n*n*n
	ik = gkbar * n*n*n*n * (v - ek)

}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
 
        evaluate_fct(v)

	n' = (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {

	: Q10 adjustment
	tadj3 = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)

	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL v2

	v2 = v - vtraub2

	tau_n = (tn * tadj3) / (exp((v2+40)/40)+exp(-(v2+40)/50))
	n_inf = 1 / (1+exp((v2+theta_n)/kappa_n))
}

UNITSON
