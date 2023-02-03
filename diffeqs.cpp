void diffeqs(long double b, long double vars[], long double diffs[])
{
	long double r, th;
	long double kt, kphi;
	long double kt2, kr2, kth2, kp2, ktp, krth;
	long double ch_rtt, ch_rtp, ch_rrr, ch_rrth, ch_rthth, ch_rpp;
	long double ch_thtt, ch_thtp, ch_thrr, ch_thrth, ch_ththth, ch_thpp;
	long double gtt, gtp, gpp, gurr, guthth;
	long double dgttdr, dgttdth, dgtpdr, dgtpdth, dgrrdr, dgrrdth, dgththdr, dgththdth, dgppdr, dgppdth;
	long double hgurr, hguthth;
	long double denom;
	long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52;

	r = vars[0];
	th = vars[1];

	/* Calculate necessary metric components and derivatives of metric components */
	/* Use code generation in Maple or Mathematica to optimize, if possible */

	t1 = r * r;
	t2 = pow(t1, 0.2e1);
	t3 = r * t1;
	t4 = t1 * t2;
	t5 = r * t2;
	t6 = a22 + t1;
	t7 = sin(th);
	t8 = spin * spin;
	t9 = pow(t6, 0.2e1);
	t10 = pow(t7, 0.2e1);
	t11 = 2;
	t12 = t5 * (double) t11 - t8 * (-t10 * t9 + t2) - t4;
	t13 = cos(th);
	t14 = t8 * pow(t13, 0.2e1);
	t15 = (t14 + t1) * r + epsi3;
	t16 = t8 + t1;
	t17 = a13 + t3;
	t18 = t10 * t8;
	t19 = t18 * r;
	t20 = t19 * t6;
	t21 = t16 * t17;
	t22 = -t21 + t20;
	t17 = t10 * ((double) t11 * t3 * t2 * t8 - t4 * t8 * t16) + pow(t16, 0.2e1) * pow(t17, 0.2e1);
	t23 = 0.1e1 / r;
	t24 = pow(t23, 0.2e1);
	t25 = t23 * t24;
	t26 = a13 * t25 + 0.1e1;
	t27 = a22 * t24 + 0.1e1;
	t28 = t16 * t26;
	t29 = (double) t11 * r;
	t30 = epsi3 * t23 + t1 + t14;
	t31 = -t18 * t27 + t28;
	t32 = t16 - t29;
	t33 = a52 * t24 + 0.1e1;
	t34 = 0.3e1 * t1;
	t35 = t14 + t34;
	t36 = (0.5e1 / 0.2e1 * t1 + 0.3e1 / 0.2e1 * t8) * r + a13;
	t34 = r * t36 - t18 * (a22 + t34) / 0.2e1;
	t22 = 0.1e1 / t22;
	t37 = pow(t22, 0.2e1);
	t38 = t22 * t37;
	t39 = t15 * r * t37;
	t40 = a52 * t8;
	t41 = a52 + t8;
	t42 = a52 + t1;
	t43 = -epsi3 * t24 + t29;
	t44 = t10 * t15 * t37;
	t45 = t15 * t23;
	t46 = a22 + t8;
	t47 = t8 * a22;
	t31 = 0.1e1 / t31;
	t27 = (t28 * t27 - t1 + t29 - t8) * pow(t31, 0.2e1);
	t28 = t10 * spin;
	t5 = -(double) t11 * t7 * t17 * t37 * t13 * (t18 - t45) - 0.4e1 * t7 * t10 * t15 * t37 * t8 * t13 * (t5 * t32 / 0.2e1 + t17 * t22 * t6);
	t29 = t8 * a13;
	t48 = 0.1e1 / t30;
	t49 = 0.1e1 / t33;
	t50 = -0.1e1 / t32;
	t51 = 0.1e1 / t32;
	t42 = pow(t42, -0.2e1);
	t52 = (double) t11 * t8 * t13 * t7;
	t40 = (-(double) t11 * (t2 * (-t41 + r) - t14 * ((t1 * (-r + 0.1e1) - a52) * r + t40)) * r + t3 * (0.4e1 * epsi3 + 0.4e1 * t40 + (-0.6e1 * a52 - 0.3e1 * epsi3) * r) + epsi3 * (-t1 * t41 + t40)) * pow(t50, 0.2e1) * t42;
	t4 = (double) t11 * t13 * (t21 * t15 + t19 * (a22 * epsi3 - ((-(a22 - t8) * r - epsi3 + a13) * r - t14 * t6) * r - t29 + t20)) * ((double) t11 * t4 + t1 * (a13 * t46 + ((a22 * r + a13) * r + t47) * r) + t29 * a22) * t7 * spin * t38;

	gtt = t39 * t12;
	gpp = t44 * t17 * t23;
	gtp = -t28 * t27 * t30;

	gurr = t48 * t32 * t33;
	guthth = t48;

	dgttdr = 0.4e1 * t39 * (r * (t3 * (-0.3e1 / 0.2e1 * r + 0.5e1 / 0.2e1) - t8 * (-t10 * t6 + t1)) + t34 * t12 * t22) + t12 * t37 * (r * t35 + t15);
	dgrrdr = t40;
	dgththdr = t43;
	dgppdr = 0.4e1 * t44 * (t21 * t36 + 0.7e1 / 0.2e1 * t18 * ((-0.4e1 / 0.7e1 * r + 0.1e1) * r - 0.3e1 / 0.7e1 * t8) * t2 + t17 * t34 * t22 * t23) + t10 * t17 * t37 * t23 * (t35 - t45);
	dgtpdr = t28 * (t27 * ((double) t11 * t30 * t31 * ((double) t11 * (t18 * a22 * t25 + r * t26) - 0.3e1 * t16 * a13 * pow(t24, 0.2e1)) - t43) + 0.5e1 * t45 * (a13 * (t1 * (0.3e1 / 0.5e1 * t46 + t1 / 0.5e1) + t47) - 0.2e1 / 0.5e1 * t3 * (-t47 + t3)) * t37);

	dgttdth = (double) t11 * t8 * r * t13 * t7 * t37 * (-r * t12 + t15 * t9) - 0.4e1 * t12 * t1 * t15 * t38 * t8 * t6 * t7 * t13;
	dgrrdth = -t52 * t51 * t49;
	dgththdth = -t52;
	dgppdth = t5;
	dgtpdth = t4;

	hgurr = 0.5*gurr;
	hguthth = 0.5*guthth;

	/* Christoffel symbols */
	ch_rtt = -hgurr*dgttdr;
	ch_rtp = -hgurr*dgtpdr;
	ch_rrr = hgurr*dgrrdr;
	ch_rrth = hgurr*dgrrdth;
	ch_rthth = -hgurr*dgththdr;
	ch_rpp = -hgurr*dgppdr;
	ch_thtt = -hguthth*dgttdth;
	ch_thtp = -hguthth*dgtpdth;
	ch_thrr = -hguthth*dgrrdth;
	ch_thrth = hguthth*dgththdr;
	ch_ththth = hguthth*dgththdth;
	ch_thpp = -hguthth*dgppdth;

	denom = (gtt*gpp-gtp*gtp);

	/* t and phi photon 4-momentum */
	kt = -(gpp+b*gtp)/denom;
	kphi = (gtp+b*gtt)/denom;

	diffs[0] = vars[3];
	diffs[1] = vars[4];
	diffs[2] = kphi;

	kt2 = kt*kt;
	kr2 = vars[3]*vars[3];
	kth2 = vars[4]*vars[4];
	kp2 = kphi*kphi;
	ktp = kt*kphi;
	krth = vars[3]*vars[4];

	/* 2nd order diff eqs for r and theta */
	diffs[3] = -(ch_rtt*kt2+ch_rrr*kr2+ch_rthth*kth2+ch_rpp*kp2+2.0*(ch_rtp*ktp+ch_rrth*krth));
	diffs[4] = -(ch_thtt*kt2+ch_thrr*kr2+ch_ththth*kth2+ch_thpp*kp2+2.0*(ch_thtp*ktp+ch_thrth*krth));
}
