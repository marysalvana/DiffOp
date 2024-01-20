#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

double deg2rad (double deg);
double rad2deg (double rad);

void calculatedistance (double * y, double * x, double * dist, int *n) {
  int nn = *n;
  printf("%d \n\n\n", nn);
  int i, j;
  for (i = 0; i < nn; i++) {
    for (j = 0; j < nn ; j++) {
      dist[nn * j + i] = pow(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2), 0.5);
      //printf("%f %f \n", x[i], x[j]);
    }
  }
}

double deg2rad (double deg) {
  return (deg * M_PI / 180);
}

double rad2deg (double rad) {
  return (rad * 180 / M_PI);
}

void DistanceEarth (double * y, double * x, double * dist, int *n) {
  int nn = *n;
  printf("%d \n\n\n", nn);
  int i, j;

  double lat1r, lon1r, lat2r, lon2r, u, v;
  double earthRadiusKm = 6371;

  for (i = 0; i < nn; i++) {
    for (j = 0; j < nn ; j++) {

      lat1r = deg2rad(y[i]);
      lon1r = deg2rad(x[i]);
      lat2r = deg2rad(y[j]);
      lon2r = deg2rad(x[j]);

      u = sin((lat2r - lat1r)/2);
      v = sin((lon2r - lon1r)/2);

      dist[nn * j + i] = 2.0 * earthRadiusKm * pow(pow(u, 2) + cos(lat1r) * cos(lat2r) * pow(v, 2), 0.5);
      //printf("%f %f \n", u, v);
    }
  }
}

double DistanceEarth22 (double lat1r, double lon1r, double lat2r, double lon2r,
                      double radius_of_sphere) {
  double u, v;

  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);

  //printf("%f \n", lat1r);

  return (2.0 * radius_of_sphere * pow(pow(u, 2) + cos(lat1r) * cos(lat2r) *
          pow(v, 2), 0.5));
}

double H (double lat1r, double lon1r, double z1, double lat2r, double lon2r,
          double z2, double scale_h, double scale_v, double radius_of_sphere) {
  return (pow(scale_h, 2) * pow(DistanceEarth22(lat1r, lon1r, lat2r, lon2r, radius_of_sphere), 2) +
          pow(scale_v, 2) * pow(z1 - z2, 2));
}

double H1 (double lat1r, double lon1r, double lat2r, double lon2r,
           double scale_h, double radius_of_sphere) {

  double L, l, con;

  L = lat1r - lat2r;
  l = lon1r - lon2r;

  con = 4 * pow(scale_h, 2) * pow(radius_of_sphere, 2);

  return (con * (sin(L / 2) * cos(L / 2) - sin(lat1r) * cos(lat2r) * pow(sin(l / 2), 2)));

}

double H3 (double lat1r, double lon1r, double lat2r, double lon2r,
           double scale_h, double radius_of_sphere) {

  double l, con;

  l = lon1r - lon2r;

  con = 4 * pow(scale_h, 2) * pow(radius_of_sphere, 2);

  return (con * sin(l / 2) * cos(l / 2));
}

double H33 (double lat1r, double lon1r, double lat2r, double lon2r,
            double scale_h, double radius_of_sphere) {

  double l, con;

  l = lon1r - lon2r;

  con = 2 * pow(scale_h, 2) * pow(radius_of_sphere, 2);

  return (con * (pow(cos(l / 2), 2) - pow(sin(l / 2), 2)));
}

double H12 (double lat1r, double lon1r, double lat2r, double lon2r,
            double scale_h, double radius_of_sphere) {

  double L, l, con;

  L = lat1r - lat2r;
  l = lon1r - lon2r;

  con = 4 * pow(scale_h, 2) * pow(radius_of_sphere, 2);

  return (con * (-pow(cos(L / 2), 2) / 2 + pow(sin(L / 2), 2) / 2 + sin(lat1r) * sin(lat2r) * pow(sin(l / 2), 2)));
}

double H13 (double lat1r, double lon1r, double lat2r, double lon2r,
            double scale_h, double radius_of_sphere) {

  double l, con;

  l = lon1r - lon2r;

  con = 4 * pow(scale_h, 2) * pow(radius_of_sphere, 2);

  return (-con * sin(lat1r) * sin(l / 2) * cos(l / 2));
}

double H4 (double z1, double z2, double scale_v) {
  return (2 * pow(scale_v, 2) * (z1 - z2));
}

double H44 (double scale_v) {
  return (2 * pow(scale_v, 2));
}

double C1 (double lat1r, double lon1r, double z1, double lat2r, double lon2r, double z2,
           double a1, double b1, double c1, double d1,
           double a2, double b2, double c2, double d2,
           double scale_h, double scale_v) {

  double radius_of_sphere = 6371;
  double h, h1, h2, h3, h4;

  h = H(lat1r, lon1r, z1, lat2r, lon2r, z2, scale_h, scale_v, radius_of_sphere);
  h1 = H1(lat1r, lon1r, lat2r, lon2r, scale_h, radius_of_sphere);
  h2 = H1(lat2r, lon1r, lat1r, lon2r, scale_h, radius_of_sphere);
  h3 = H3(lat1r, lon1r, lat2r, lon2r, scale_h, radius_of_sphere);
  h4 = H4(z1, z2, scale_v);

  return(0.25 * (a1 * a2 * h1 * h2 - b1 * b2 * pow(h3, 2) - c1 * c2 * pow(h4, 2)
                   - a1 * b2 * h1 * h3 + a2 * b1 * h2 * h3
                   - a1 * c2 * h1 * h4 + a2 * c1 * h2 * h4
                   - b1 * c2 * h3 * h4 - b2 * c1 * h3 * h4)
                   + h * d1 * d2);

}

double C2 (double lat1r, double lon1r, double z1, double lat2r, double lon2r, double z2,
           double a1, double b1, double c1, double d1,
           double a2, double b2, double c2, double d2,
           double scale_h, double scale_v, double nu) {

  double radius_of_sphere = 6371;
  double h12, h13, h23, h33, h44;

  h12 = H12(lat1r, lon1r, lat2r, lon2r, scale_h, radius_of_sphere);
  h13 = H13(lat1r, lon1r, lat2r, lon2r, scale_h, radius_of_sphere);
  h23 = H13(lat2r, lon1r, lat1r, lon2r, scale_h, radius_of_sphere);
  h33 = H33(lat1r, lon1r, lat2r, lon2r, scale_h, radius_of_sphere);
  h44 = H44(scale_v);

  return(-0.5 * (a1 * a2 * h12 - b1 * b2 * h33 - c1 * c2 * h44
                   - a1 * b2 * h13 + a2 * b1 * h23) + 2 * nu * d1 * d2);

}

double uni_differential (double x1, double y1, double z1, double x2, double y2, double z2,
                         double a1, double b1, double c1, double d1,
                         double a2, double b2, double c2, double d2,
                         double scale_h, double scale_v, double nu, double sigma_square) {

  double lat1r, lon1r, lat2r, lon2r;
  double con, expr, f, f_prime, C1_val, C2_val, cov_val;
  double radius_of_sphere = 6371;

  con = pow(2, nu - 1) * tgamma(nu);
  con = 1.0 / con;
  con = sigma_square * con;

  lat1r = deg2rad(y1);
  lon1r = deg2rad(x1);
  lat2r = deg2rad(y2);
  lon2r = deg2rad(x2);

  //printf("%f %f \n", y1, lat1r);

  expr = pow(H(lat1r, lon1r, z1, lat2r, lon2r, z2, scale_h, scale_v, radius_of_sphere), 0.5);

  C1_val = C1(lat1r, lon1r, z1, lat2r, lon2r, z2,
              a1, b1, c1, d1, a2, b2, c2, d2,
              scale_h, scale_v);

  C2_val = C2(lat1r, lon1r, z1, lat2r, lon2r, z2,
              a1, b1, c1, d1, a2, b2, c2, d2,
              scale_h, scale_v, nu);

  //printf("%f \n", expr);

  if(expr == 0){
    cov_val = con * (C1_val + C2_val) + sigma_square * d1 * d2;
  }else{
    f = pow(expr, nu - 1) * gsl_sf_bessel_Knu(nu - 1, expr);
    f_prime = pow(expr, nu - 2) * gsl_sf_bessel_Knu(nu - 2, expr);
    cov_val = con * (C1_val * f_prime + C2_val * f + d1 * d2 * (pow(expr, 2) * f_prime + 2 * (nu - 1) * f));
  }

  return(cov_val);
}

void CovBiDifferential (double * x, double * y, double * z,
                       double *a1, double *b1, double * c1, double *d1,
                       double *a2, double *b2, double * c2, double *d2,
                       double *scale_h, double *scale_v, double *nu,
                       double *sigsq1, double *sigsq2, double *beta,
                       double * dist, int *n) {
  int nn = *n, i, j, p = 2;

  double A1 = *a1, B1 = *b1, D1 = *d1, A2 = *a2, B2 = *b2, D2 = *d2;
  double NU = *nu, SIGSQ1 = *sigsq1, SIGSQ2 = *sigsq2, BETA = *beta, SCALE_H = *scale_h, SCALE_V = *scale_v;

  for (i = 0; i < nn; i++) {
    for (j = 0; j < nn ; j++) {

      dist[p * nn * i + j] = uni_differential(x[i], y[i], z[i], x[j], y[j], z[j],
                                          A1, B1, c1[i], D1, A1, B1, c1[j], D2,
                                          SCALE_H, SCALE_V, NU, SIGSQ1);

      dist[p * nn * i + nn + j] = dist[p * nn * (nn + j) + i] = uni_differential(x[i], y[i], z[i], x[j], y[j], z[j],
                                          A2, B2, c2[i], D2, A1, B1, c1[j], D1,
                                          SCALE_H, SCALE_V, NU, BETA * pow(SIGSQ1 * SIGSQ1, 0.5));

      dist[p * nn * (nn + i) + nn + j] = uni_differential(x[i], y[i], z[i], x[j], y[j], z[j],
                     A2, B2, c2[i], D2, A2, B2, c2[j], D2,
                     SCALE_H, SCALE_V, NU, SIGSQ2);
    }
  }
}



