void raytrace(long double xscr, long double yscr, long double traced[4], long double rdisk)
{
  long double dscr, xscr2, yscr2;
  long double r, th, phi, rau, thau, phiau, r0, th0, phi0;
  long double kr, kth, kt0, kr0, kth0, kphi0;
  long double r02, s0, s02;
  long double fact1, fact2, fact3, omega;
  long double h;

  long double height;

  long double cosem;
  long double b;

  long double g[4][4];
  long double diffs[5], vars[5], vars_temp[5], vars_4th[5], vars_5th[5], k1[5], k2[5], k3[5], k4[5], k5[5], k6[5];
  long double rgap, rmid, thmid;
  long double gfactor;
  long double err, errmin, errmax;

  int errcheck, crosscheck = 0, acccheck = 0, blockcheck = 0;
  int stop_integration = 0;
  int i;

  /* ----- Set computational parameters ----- */
  dscr = 1.0e+8;    /* distance of the observer */
  errmin = 1.0e-9;  /* error bounds for RK45 */
  errmax = 1.0e-7;
  long double cross_tol = 1.0e-8;   /* sought accuracy at disk crossing */

  // set disk height
  if(rdisk >= isco)
    height = 3.0*Mdl*(1.0-sqrt(isco/rdisk))/eta;
  else
    height = 0.0;

  h = -1.0;    /* initial step size */

  /* ----- compute photon initial conditions ----- */
  xscr2 = xscr*xscr;
  yscr2 = yscr*yscr;

  fact1 = yscr*sin(inc) + dscr*cos(inc);
  fact2 = dscr*sin(inc) - yscr*cos(inc);

  r02 = xscr2 + yscr2 + dscr*dscr;

  //Initial r, theta, and phi
  r0 = sqrt(r02);
  th0 = acos(fact1/r0);
  phi0 = atan2(xscr,fact2);

  s0  = sin(th0);
  s02 = s0*s0;

  //Initial r, theta, and phi momentum
  kr0 = dscr/r0;
  kth0 = -(cos(inc) - dscr*fact1/r02)/sqrt(r02-fact1*fact1);
  kphi0 = -xscr*sin(inc)/(xscr2+fact2*fact2);

  metric(r0, th0, g);

  //Energy used to scale affine parameter
  fact3 = sqrt(g[0][3]*g[0][3]*kphi0*kphi0-g[0][0]*(g[1][1]*kr0*kr0+g[2][2]*kth0*kth0+g[3][3]*kphi0*kphi0));

  //Initial t momentum
  kt0 = -(g[0][3]*kphi0+fact3)/g[0][0];

  //Impact parameter, b=L/E
  b = -(g[3][3]*kphi0+g[0][3]*kt0)/(g[0][0]*kt0+g[0][3]*kphi0);

  //Scale by Energy
  kr0 /= fact3;
  kth0 /= fact3;

  /* ----- RK45 ----- */

  /* set initial values */
  r = r0;
  th = th0;
  phi = phi0;

  kr = kr0;
  kth = kth0;

  omega = r02*s02*kphi0/kt0;  /* disk angular velocity */

  do {

    vars[0] = r;
    vars[1] = th;
    vars[2] = phi;
    vars[3] = kr;
    vars[4] = kth;

    do {

      errcheck = 0;

      /* ----- compute RK1 ----- */

      diffeqs(b, vars, diffs);
      for(i = 0; i <= 4; i++)
      {
        k1[i] = h*diffs[i];
        vars_temp[i] = vars[i] + a1*k1[i];
      }

      /* ----- compute RK2 ----- */

      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k2[i] = h*diffs[i];
        vars_temp[i] = vars[i] + b1*k1[i] + b2*k2[i];
      }

      /* ----- compute RK3 ----- */

      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k3[i] = h*diffs[i];
        vars_temp[i] = vars[i] + c1*k1[i] + c2*k2[i] + c3*k3[i];
      }

      /* ----- compute RK4 ----- */

      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k4[i] = h*diffs[i];
        vars_temp[i] = vars[i] + d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i];
      }

      /* ----- compute RK5 ----- */

      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k5[i] = h*diffs[i];
        vars_temp[i] = vars[i] + e1*k1[i] + e2*k2[i] + e3*k3[i] + e4*k4[i] + e5*k5[i];
      }

      /* ----- compute RK6 ----- */

      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
        k6[i] = h*diffs[i];

      /* ----- local error ----- */

      for(i=0; i<= 4; i++)
      {
        vars_4th[i] = vars[i] + f1*k1[i] + f2*k2[i] + f3*k3[i] + f4*k4[i] + f5*k5[i];
        vars_5th[i] = vars[i] + g1*k1[i] + g2*k2[i] + g3*k3[i] + g4*k4[i] + g5*k5[i] + g6*k6[i];

        err = fabs((vars_4th[i]-vars_5th[i])/max(vars_4th[i], vars[i]));

		if(err > errmax && crosscheck == 0)  /* accuracy not achieved and photon hasn't crossed disk */
          errcheck = 1;
        else if(err < errmin && errcheck != 1 && crosscheck == 0)  /* accuracy better than wanted, but photon hasn't crossed disk */
          errcheck = -1;
      }

      if(errcheck == 1)  /* accuracy not achieved, lower step size */
        h/=2.0;
      else if(errcheck == -1)  /* accuracy better than wanted, but photon hasn't crossed disk, increase step size */
        h*=2.0;

    } while (errcheck == 1);

    /* ----- cross disk/horizon check ----- */

    rau = r;
    thau = th;
    phiau = phi;
    r = vars_4th[0];
    th = vars_4th[1];
    phi = vars_4th[2];

    if ( r*cos(th) < height )  /* check if photon has crossed disk */
    {
        /* check if accuracy achieved; if so, move on to setting final values */
        if(sqrt(r*r+rau*rau-2.0*r*rau*(cos(th)*cos(thau)+sin(th)*sin(thau)*cos(phi-phiau))) <= cross_tol)
        {
            acccheck = 1;
            //printf("Crossed disk near desired radius and within error tolerance. Exit.\n");
        }
        /* otherwise, photon has crossed disk, but has not achieved accuracy; continue to zoom in on disk */
        else
        {
            crosscheck = 1;
            //printf("Crossed disk near desired radius but not within error tolerance. Continue.\n");
        }
        if(acccheck == 1)  /* set final values */
        {
            kr = vars_4th[3];
            kth = vars_4th[4];

            /* calculate average/midpoint values */
            rmid = 0.5*(r+rau);
            thmid = 0.5*(th+thau);

            //printf("%Le %Le\n", rmid, thmid);

            if (rmid*sin(thmid) >= isco-0.001 && rmid*sin(thmid) < 1.05*dscr)
            {
                stop_integration = 1; /* the photon hits the disk */
                break;
            }
            else
            {
                stop_integration = 2;   /* the photon misses the disk or other error */
                break;
            }
        }
        else if (crosscheck == 1)  /* did not achieve accuracy; go back a step, and decrease step size */
        {
            r = rau;
            th = thau;
            phi = phiau;
            h /= 2.0;
        }
    }
    else if(r <= 1.+sqrt(1.-spin2)+0.001)
    {
        stop_integration = 2;   /* photon crosses event horizon */
        //printf("Photon crossed horizon\n");
        break;
    }
    else if (h > -1.0e-20)
    {
        stop_integration = 3;   /* photon is stuck */
        //printf("Photon is stuck\n");
        break;
    }
    else  /* not done, take a step */
    {
        kr = vars_4th[3];
        kth = vars_4th[4];
    }
  } while (stop_integration == 0);

  /* ----- Calculate redshift, cosem, and return values ----- */

  if (stop_integration == 1)  /* photon hit disk, no issues */
  {
    redshift(rmid, thmid, omega, gfactor);
    cosem = gfactor * emis_angle(rmid, thmid, kr, kth);
  }
  else  /* photon crossed horizon, missed disk, or other issue */
  {
    rmid = 0.0;
    gfactor = 0.0;
  }

  traced[0] = rmid*sin(thmid);
  traced[1] = cosem;
  traced[2] = gfactor;
}
