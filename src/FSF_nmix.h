
extern "C" {

  void draw_indicators_logistic(int *r, double *z, double *lambda, int *n, 
				double *w, double *s, int *J);
  
  void draw_indicators_generic(int *r, double *res, int *n, 
			       double *w, double *m, double *s, int *J);

}
