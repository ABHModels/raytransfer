/* Calculate redshift */
void redshift(long double r, long double th, long double ktkp, long double& gg)
{
	long double Omega;
	long double uet;
	long double gtt,gtp,gpp,dgttdr,dgtpdr,dgppdr;
	long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31;

	/* Calculate necessary metric and derivative of metric components */
	/* If possible, use Maple or Mathematica to optimize code */

	t1 = r * r;
	t2 = pow(t1, 0.2e1);
	t3 = r * t1;
	t4 = t1 * t2;
	t5 = a22 + t1;
	t6 = sin(th);
	t6 = pow(t6, 0.2e1);
	t7 = spin * spin;
	t8 = 2;
	t9 = (double) t8 * r * t2 - t7 * (-t6 * pow(t5, 0.2e1) + t2) - t4;
	t10 = cos(th);
	t10 = t7 * pow(t10, 0.2e1);
	t11 = (t10 + t1) * r + epsi3;
	t12 = t7 + t1;
	t13 = a13 + t3;
	t5 = t5 * t6;
	t14 = t12 * t13;
	t15 = -t5 * t7 * r + t14;
	t4 = t6 * ((double) t8 * t3 * t2 * t7 - t4 * t7 * t12) + pow(t12, 0.2e1) * pow(t13, 0.2e1);
	t13 = 0.1e1 / r;
	t16 = pow(t13, 0.2e1);
	t17 = t13 * t16;
	t18 = a13 * t17 + 0.1e1;
	t19 = a22 * t16 + 0.1e1;
	t20 = t12 * t18;
	t21 = (double) t8 * r;
	t22 = epsi3 * t13 + t1 + t10;
	t23 = -t7 * t19 * t6 + t20;
	t24 = 0.3e1 * t1;
	t10 = t10 + t24;
	t25 = (0.5e1 / 0.2e1 * t1 + 0.3e1 / 0.2e1 * t7) * r + a13;
	t24 = r * t25 - t7 * (a22 + t24) * t6 / 0.2e1;
	t15 = 0.1e1 / t15;
	t26 = pow(t15, 0.2e1);
	t27 = t11 * r * t26;
	t28 = t6 * t11 * t26;
	t29 = t11 * t13;
	t30 = t7 * a22;
	t31 = 0.3e1 / 0.5e1;
	t23 = 0.1e1 / t23;
	t19 = (t20 * t19 - t1 + t21 - t7) * pow(t23, 0.2e1);
	t20 = t6 * spin;

	gtt = t27 * t9;
	gpp = t28 * t4 * t13;
	gtp = -t20 * t19 * t22;

	dgttdr = 0.4e1 * t27 * (r * (t3 * (-0.3e1 / 0.2e1 * r + 0.5e1 / 0.2e1) - t7 * (-t5 + t1)) - t24 * t9 * t15) + t9 * t26 * (r * t10 + t11);
	dgppdr = 0.4e1 * t28 * (-t4 * t24 * t15 * t13 + t14 * t25 + (0.7e1 / 0.2e1 * (-0.4e1 / 0.7e1 * r + 0.1e1) * r - 0.3e1 / 0.2e1 * t7) * t7 * t2 * t6) + t6 * t4 * t26 * t13 * (-t29 + t10);
	dgtpdr = t20 * (0.5e1 * t29 * (a13 * (t1 * (t1 / 0.5e1 + (a22 + t7) * t31) + t30) - 0.2e1 / 0.5e1 * t3 * (-t30 + t3)) * t26 + t19 * (epsi3 * t16 - t21 + (double) t8 * t22 * t23 * ((double) t8 * (t30 * t17 * t6 + r * t18) - 0.3e1 * t12 * a13 * pow(t16, 0.2e1))));

	Omega  = (-dgtpdr + sqrt(dgtpdr*dgtpdr - dgttdr*dgppdr))/dgppdr; //angular velocity

	uet = sqrt(-gtt - 2.*gtp*Omega - gpp*Omega*Omega); //t-component of 4-velocity

	gg = uet/(1. - ktkp*Omega); //redshift
}

/* Calculate specific energy at ISCO */
long double specific_energy(long double r)
{
    long double Omega, se;
    long double gtt,gtp,gpp,dgttdr,dgtpdr,dgppdr;
    long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;

	/* Calculate necessary metric and derivative of metric components */
	/* If possible, use Maple or Mathematica to optimize code */

	t1 = r * r;
	t2 = pow(t1, 0.2e1);
	t3 = r * t1;
	t4 = t1 * t2;
	t5 = a22 + t1;
	t6 = spin * spin;
	t7 = 2.0;
	t8 = (double) t7 * r * t2 - t6 * (-pow(t5, 0.2e1) + t2) - t4;
	t9 = epsi3 + t3;
	t10 = t6 + t1;
	t11 = a13 + t3;
	t12 = t10 * t11;
	t13 = t6 * r * t5 - t12;
	t4 = pow(t10, 0.2e1) * pow(t11, 0.2e1) + (double) t7 * t3 * t2 * t6 - t4 * t6 * t10;
	t11 = 0.1e1 / r;
	t14 = pow(t11, 0.2e1);
	t15 = t11 * t14;
	t16 = a13 * t15 + 0.1e1;
	t17 = a22 * t14 + 0.1e1;
	t18 = t10 * t16;
	t19 = (double) t7 * r;
	t20 = epsi3 * t11 + t1;
	t21 = -t17 * t6 + t18;
	t22 = (0.5e1 / 0.2e1 * t1 + 0.3e1 / 0.2e1 * t6) * r + a13;
	t23 = r * t22 - t6 * (a22 + 0.3e1 * t1) / 0.2e1;
	t13 = 0.1e1 / t13;
	t24 = pow(t13, 0.2e1);
	t25 = t9 * r * t24;
	t26 = t9 * t24;
	t27 = t6 * a22;
	t28 = 0.3e1 / 0.5e1;
	t21 = 0.1e1 / t21;
	t17 = (t18 * t17 - t1 + t19 - t6) * pow(t21, 0.2e1);

	gtt = t25 * t8;
	gpp = t26 * t4 * t11;
	gtp = -t17 * spin * t20;

	dgttdr = 0.4e1 * t25 * (r * (t3 * (-0.3e1 / 0.2e1 * r + 0.5e1 / 0.2e1) - t6 * (t1 - t5)) + t23 * t8 * t13) + t8 * t24 * (0.3e1 * t3 + t9);
	dgppdr = 0.4e1 * t26 * (t4 * t23 * t13 * t11 + t12 * t22 + (0.7e1 / 0.2e1 * (-0.4e1 / 0.7e1 * r + 0.1e1) * r - 0.3e1 / 0.2e1 * t6) * t6 * t2) + t4 * t24 * (-t14 * t9 + 0.3e1 * r);
	dgtpdr = spin * (t17 * (epsi3 * t14 - t19 + (double) t7 * t20 * t21 * ((double) t7 * (r * t16 + t27 * t15) - 0.3e1 * t10 * a13 * pow(t14, 0.2e1))) + 0.5e1 * t26 * (a13 * (t1 * (t1 / 0.5e1 + (a22 + t6) * t28) + t27) - 0.2e1 / 0.5e1 * t3 * (-t27 + t3)) * t11);

    Omega  = (-dgtpdr + sqrt(dgtpdr*dgtpdr - dgttdr*dgppdr))/dgppdr; //angular velocity

    se = -1.0*(gtt + Omega*gtp)/sqrt(-1.0*gtt - 2.0*Omega*gtp - Omega*Omega*gpp); //specific energy

    return se;
}

/* Calculate cosine of the emission angle modulo the redshift factor */
long double emis_angle(long double r, long double th, long double kr, long double kth)
{
	long double Zr, Zth, k, cs1, ss1, ss3, r3, angle;
	long double gurth[2];

	cs1 = cos(th);
	ss1 = sin(th);
	ss3 = ss1*ss1*ss1;
	r3 = r*r*r;

	k = 3.0 / eta * Mdl;
	Zr = 0.5 * k * sqrt(isco / (r3*ss1)) - cs1;
	Zth = 0.5 * k * cs1 * sqrt(isco / (r*ss3)) + r*ss1;

	uppermetric(r,th,gurth);

	angle = 1.0 / sqrt(gurth[0]*Zr*Zr + gurth[1]*Zth*Zth) * (Zr*kr+Zth*kth);

	if (angle < 0.0)
		angle *= -1.0;

	return angle;
}
