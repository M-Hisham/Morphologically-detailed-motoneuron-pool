: Fast Sodium Channels - Initial Segment


NEURON {
	SUFFIX NafIs

	NONSPECIFIC_CURRENT ina

	RANGE gnabar, ena, gna
	RANGE m_inf, h_inf 
	RANGE tau_m, tau_h
	RANGE amA, bmA, theta_h
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gnabar  = 1.3392  (mho/cm2)
	ena     = 50.0    (mV)

	amA = 10
	bmA = 30
	th  = 30	  (ms)

	theta_h = 55
	kappa_h = 7

	celsius = 20      (degC)
	vtraub2	= -10	  (mV)
	vtraub22 = -70	  (mV)
}

STATE {
	m h
}

ASSIGNED {
	dt      (ms)
	v       (mV)

	ina	(mA/cm2)

	gna   (mho/cm2)

	m_inf
	h_inf


	tau_m   (ms)
	tau_h	(ms)

	tadj1
	tadj2
	tadj3
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	gna = gnabar * m*m*m*h
	ina = gnabar * m*m*m*h * (v - ena)

}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
        evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h

}

UNITSOFF

INITIAL {

	: Q10 adjustment
	tadj3 = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2,v22

	v2 = v - vtraub2
	v22 = v - vtraub22


	a = 0.4*vtrap(amA-v22,5)
	b = 0.4*vtrap(v22-bmA,5)
	tau_m = 1 / (a + b) / tadj3
	m_inf = a / (a + b)

	tau_h = (th * tadj3) / (exp((v2+50)/15)+exp(-(v2+50)/16))
	h_inf = 1 / (1+exp((v2+theta_h)/kappa_h))
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}

UNITSON
