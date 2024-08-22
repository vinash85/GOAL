<!--GOAL
====

Genetic Omics Association resolve Linkage disequilibrium

==== 
-->

# GOAL: Genetic Omics Association resolve Linkage Disequilibrium



## Step 1: Clone the Repository

First, clone the GOAL repository from GitHub:

```bash
git clone https://github.com/vinash85/GOAL.git
cd GOAL
```


## Step 2: Install Dependencies

### 2.1. Install GSL

First, install the GNU Scientific Library (GSL) using your system’s package manager:

#### For Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install libgsl-dev
```
or
```bash
wget https://mirror.team-cymru.com/gnu/gsl/gsl-latest.tar.gz -O gsl-latest.tar.gz
tar -zxvf gsl-latest.tar.gz
cd gsl-*/
./configure
make
sudo make install
```


### 2.2. Install ransampl

Download and install the `ransampl` library:

```bash
wget https://sourceforge.net/projects/ransampl/files/latest/download -O ransampl-latest.tar.gz
tar -zxvf ransampl-1.0.tar.gz
cd ransampl-*/
./configure
make
sudo make install
```

### 2.3. Install LAPACK and BLAS

Ensure that LAPACK and BLAS libraries are installed on your system:

#### For Ubuntu/Debian:

```bash
sudo apt-get install liblapack-dev libblas-dev
```



### 2.4. Install Required R Packages

Install the necessary R packages:

```r
install.packages(c("BayesLogit", "Rcpp", "RcppArmadillo", "RcppProgress"))
```

If `BayesLogit` is not available on CRAN, install it from GitHub:

```r
install.packages("devtools")
devtools::install_github("jwindle/BayesLogit")
```

## Step 3: Modify the `Makevars` File

Before compiling GOAL, update the `Makevars` file to include the paths to GSL and ransampl.

1. Navigate to the `src/` directory of the GOAL source code.
2. Open `Makevars` and modify the following lines:

   ```makefile
   PKG_CXXFLAGS = -fopenmp -I/usr/local/include -DUSE_R -DNDEBUG -DDISABLE_SINGLE -DNTHROW -DDISABLE_FIO
   PKG_LIBS = -L/usr/local/lib -lgsl -lgslcblas -lm -lransampl -fopenmp -lgomp -llapack -lblas -lgfortran -lm -lquadmath -L/usr/lib/R/lib -lR
   ```

   **Note**: If your libraries are installed in different locations, adjust the paths accordingly.

## Step 4: Install the GOAL Package

1. Navigate to the root directory of the GOAL package source code.
2. Run the installation command:

   ```bash
   sudo R CMD INSTALL .
   ```

   This will compile and install GOAL into the system-wide R library.

## Step 5: Verify the Installation

Open an R session and load the GOAL package to ensure it was installed correctly:

```r
library(GOAL)
```

If the package loads without errors, the installation was successful.

## Troubleshooting

- **Missing Dependencies**: Ensure all dependencies are installed and paths are correctly set.
- **Linker Errors**: Double-check the `Makevars` file for correct paths to GSL, ransampl, LAPACK, and BLAS.
- **Package Not Found**: Verify that the package is installed in the correct R library path using `.libPaths()`.

## Contact

For further assistance, please contact  asahu@salud.unm.edu

<!--#intstall from source instruction 
#GOAL requires library GSL and ranmsampl to be installed from source.

#1. ransampl can be installed from http://sourceforge.net/projects/ransampl/.
#2. gsl libraries can be installed from http://www.gnu.org/software/gsl/

#Make sure that ransampl and gsl libraries are both included in your include and library path. 
#This can be done by modifying PKG_LIBS and PKG_CXXFLAGS in  src/Makevars an example is included (commented out )-->

To speed up the computation GOAL also uses openmp (http://openmp.org/wp/), however this is optional.
