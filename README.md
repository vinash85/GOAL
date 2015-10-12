GOAL
====

Genetic Omics Association resolve Linkage disequilibrium

#intstall from source instruction 
GOAL requires library GSL and ranmsampl to be installed from source.

1. ransampl can be installed from http://sourceforge.net/projects/ransampl/.
2. gsl libraries can be installed from http://www.gnu.org/software/gsl/

Make sure that ransampl and gsl libraries are both included in your include and library path. 
This can be done by modifying PKG_LIBS and PKG_CXXFLAGS in  src/Makevars an example is included (commented out )

To speed up the computation GOAL also uses openmp (http://openmp.org/wp/), however this is optional.
