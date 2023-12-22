: TITLE Motoneuron L-type Calcium channels
:
: The parameters for this current come from V. Booth et al. J Neurophysiol 78:3371-3385, 1997
: Iterative equations
:


NEURON {
	SUFFIX Llva
	NONSPECIFIC_CURRENT ica
	RANGE gcaLlvabar, eca, gcaLlva
	RANGE m_inf
	RANGE theta_m
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gcaLlvabar  = 0.0003  (mho/cm2)
	eca     = 60	      (mV)	
	tau_m	= 20	      (ms)
	theta_m = -30         (mV)
	kappa_m = -6	      (mV)
}

STATE {
	m
}

ASSIGNED {
	v		(mV)
	dt		(ms)
	ica		(mA/cm2)
	m_inf
	tadj
	gcaLlva		(mho/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcaLlva = gcaLlvabar * m
	ica = gcaLlvabar * m * (v - eca)
}

DERIVATIVE states {
	evaluate_fct(v)
	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {

:
:  Q10 was assumed to be 3 for both currents
:
:	tadj = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV)) {

	m_inf = 1 / (1 + (Exp((v - theta_m)/ kappa_m)))         : / tadj

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


