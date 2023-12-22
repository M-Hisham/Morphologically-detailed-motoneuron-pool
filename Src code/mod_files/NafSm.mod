TITLE Motoneuron Soma channels
:
: Fast Sodium Channels - Soma


NEURON {
	SUFFIX NafSm
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

	gnabar  = 0.06	(mho/cm2)

	amA	= 14	(mV)
	bmA	= 38	(mV)
	th	= 30	(ms)

	theta_h = 55
	kappa_h = 7

	ena     = 50    (mV)
	celsius = 36    (degC)

	vtraub	= -10	(mV)
	vtraub2	= -70	(mV)
}

STATE {
	m h
}

ASSIGNED {
	dt      (ms)
	v       (mV)

	ina     (mA/cm2)
	gna	  (mho/cm2)

	m_inf
	h_inf

	tau_m	(ms)
	tau_h	(ms)

	tadj
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

	:  Q10 was assumed to be 3
	tadj = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)

	m = m_inf
	h = h_inf

}

PROCEDURE evaluate_fct(v(mV)) { LOCAL v2, v22, a, b

	v2 = v - vtraub : convert to traub convention
	v22 = v - vtraub2
	
	a = 0.4*vtrap(amA-v22,5)
	b = 0.4*vtrap(v22-bmA,5)
	tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)

	tau_h = (th * tadj) / (Exp((v2+50)/15)+Exp(-(v2+50)/16))
	h_inf = 1 / (1+Exp((v2+theta_h)/kappa_h))

}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(Exp(x/y)-1)
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}
