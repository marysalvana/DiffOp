

void calculatedistance (double * y, double * x, double * dist, int *n);
void DistanceEarth (double * y, double * x, double * dist, int *n);
double DistanceEarth22 (double lat1r, double lon1r, double lat2r, double lon2r,
                        double radius_of_sphere);
void CovBiDifferential (double * x, double * y, double * z,
                       double *a1, double *b1, double * c1, double *d1,
                       double *a2, double *b2, double * c2, double *d2,
                       double *scale_h, double *scale_v, double *nu,
                       double *sigsq1, double *sigsq2, double *beta,
                       double * dist, int *n);

double H (double lat1r, double lon1r, double z1, double lat2r, double lon2r,
          double z2, double scale_h, double scale_v, double radius_of_sphere);
double H1 (double lat1r, double lon1r, double lat2r, double lon2r,
           double scale_h, double radius_of_sphere);
double H3 (double lat1r, double lon1r, double lat2r, double lon2r,
           double scale_h, double radius_of_sphere);
double H33 (double lat1r, double lon1r, double lat2r, double lon2r,
            double scale_h, double radius_of_sphere);
double H12 (double lat1r, double lon1r, double lat2r, double lon2r,
            double scale_h, double radius_of_sphere);
double H13 (double lat1r, double lon1r, double lat2r, double lon2r,
            double scale_h, double radius_of_sphere);
double H4 (double z1, double z2, double scale_v);
double H44 (double scale_v);
double C1 (double lat1r, double lon1r, double z1, double lat2r, double lon2r, double z2,
           double a1, double b1, double c1, double d1,
           double a2, double b2, double c2, double d2,
           double scale_h, double scale_v);
double C2 (double lat1r, double lon1r, double z1, double lat2r, double lon2r, double z2,
           double a1, double b1, double c1, double d1,
           double a2, double b2, double c2, double d2,
           double scale_h, double scale_v, double nu);

double uni_differential (double x1, double y1, double z1, double x2, double y2, double z2,
                         double a1, double b1, double c1, double d1,
                         double a2, double b2, double c2, double d2,
                         double scale_h, double scale_v, double nu, double sigma_square);
