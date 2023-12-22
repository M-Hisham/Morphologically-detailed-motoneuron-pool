: Persistent Sodium Channels - Initial Segment


NEURON {
	SUFFIX NapIs

	NONSPECIFIC_CURRENT inap

	RANGE gnapbar, ena, gnap
	RANGE mp_inf
	RANGE tau_mp
	RANGE ampB, bmpB
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gnapbar = 3.2971e-5  (mho/cm2)

	ena     = 50.0    (mV)

	celsius = 20      (degC)

	vtraub1 = 10	  (mV)

	ampA = 0.01
	ampB = 38.4
	ampC = 5
	bmpA = 0.00025
	bmpB = 42.7
	bmpC = 5
}

STATE {
	mp
}

ASSIGNED {
	dt      (ms)
	v       (mV)
	inap    (mA/cm2)
	gnap    (mho/cm2)

	mp_inf

	tau_mp	(ms)

	tadj1
	tadj2
	tadj3
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	gnap = gnapbar * mp*mp*mp
	inap = gnapbar * mp*mp*mp * (v - ena)

}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
        evaluate_fct(v)
	mp' = (mp_inf - mp) / tau_mp
}

UNITSOFF

INITIAL {

	: Q10 adjustment
	tadj1 = 2.2 ^ ((celsius-20)/ 10 )

	evaluate_fct(v)
	mp = mp_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v1

	v1 = v - vtraub1 : convert to traub convention

	a = vtrap1(v1)
	b = vtrap2(v1)
	tau_mp = 1 / (a + b) / tadj1
	mp_inf = a / (a + b)
}

FUNCTION vtrap1(x) {
	if (fabs((x+ampB)/ampC) < 1e-6) {
		vtrap1 = ampA*ampC
	}else{
		vtrap1 = (ampA*(x+ampB)) / (1 - exp(-(x+ampB)/ampC))
	}
}

FUNCTION vtrap2(x) {
	if (fabs((x+bmpB)/bmpC) < 1e-6) {
		vtrap2 = -bmpA*bmpC
	}else{
		vtrap2 = (bmpA*(-(x+bmpB))) / (1 - exp((x+bmpB)/bmpC))
	}
}

UNITSON
