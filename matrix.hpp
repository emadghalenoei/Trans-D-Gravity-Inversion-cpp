//
//  matrix.h
//
//  Created by Emad Ghalenoei in 2020
//  Copyright (c) 2020 Emad Ghalenoei. All rights reserved.


#ifndef _Project_matrix_
#define _Project_matrix_

#include <stdio.h>
#include <fstream> // for file access
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <random>
#include <chrono>
#include <ctime>
#include <stdarg.h>
#include <time.h>

#define PI 3.14159265

using namespace std;

using std::vector;
using std::tuple;
class Matrix
{
private:
    unsigned m_rowSize;
    unsigned m_colSize;
    vector<vector<double> > m_matrix;
public:
    Matrix(unsigned, unsigned);
    Matrix(const char *);
    Matrix(const Matrix &);
    ~Matrix();

    // Matrix Operations
    Matrix operator+(Matrix &);
    Matrix operator-(Matrix &);
    Matrix operator*(Matrix &);
    Matrix transpose();

    // Scalar Operations
    Matrix operator+(double);
    Matrix operator-(double);
    Matrix operator*(double);
    Matrix operator/(double);
    Matrix operator^(double);

    // Aesthetic Methods
    double& operator()(const unsigned &, const unsigned &);
    void print() const;
    unsigned getRows() const;
    unsigned getCols() const;
    unsigned getLength() const;


    // Deflation
    Matrix deflation(Matrix &, double&);

    Matrix MatExtraction(const int &,const int &,const int &,const int &);
    void print_element(int, int) const;
    Matrix vectorize();
    Matrix reshape( int,  int );
    Matrix absolute();
    Matrix rounding();
    Matrix log_fun();
    Matrix toeplitz();
    Matrix exp_mat();
    Matrix cos_mat();
    //Matrix cholesky(const char * );
    Matrix inv_cholesky();
    Matrix inv_LU();
    Matrix one();
    Matrix damping();
    Matrix diag();
    Matrix eye();
    Matrix DeleteRow(const unsigned &);
    Matrix DeleteCol(const unsigned &);
    double det_cholesky();
    double log_det_cholesky();
    double sum();
    double Min();
    double Max();
    double prod();
    static Matrix linespace(double, double, int );
    static Matrix colon(const double &,const double&,const double &);
    static Matrix logspace(double, double, int );
    void write2txtfile(const char *) const;
    //static tuple<Matrix, Matrix> meshgrid(Matrix &, Matrix &, Matrix &, Matrix &,const int &,const int &);
    static void meshgrid(Matrix &, Matrix &, Matrix &, Matrix &,const int &,const int &);
    static Matrix normrnd(Matrix &, Matrix &);
    static Matrix cauchy_dist(Matrix &, Matrix &);
    static double cauchy_dist(double &, double &);
    static Matrix filter(Matrix &,Matrix &,Matrix & );
    static Matrix RandomDouble(unsigned, unsigned);
    static double rand01();
    static Matrix dot( Matrix &, Matrix &);
    static Matrix dot_division(Matrix &,Matrix &);
    static double Log_Likelihood_AR(Matrix &,Matrix &,Matrix & );
    static double Log_Likelihood_full_C(Matrix &,Matrix &);
    static double Log_Likelihood_initial(Matrix &,Matrix &, double &);
    static Matrix join_in_right(const initializer_list<Matrix> &);
    static Matrix join_in_bottom(const initializer_list<Matrix> &);
    static Matrix XZ_model(Matrix &,Matrix &, Matrix &, const int &, const int &);
    static Matrix FM_fast(Matrix &,Matrix &, const double &);
    static Matrix e2C(double &,  Matrix &, const int);
    tuple<Matrix, Matrix> LUdecomposition();
    tuple<Matrix, Matrix> cholesky();
    static void MCMC_0(Matrix &,Matrix &,const double,const double,const double,const double,const int,const int, Matrix &,Matrix &,const int, Matrix &, Matrix &, const double &, int &, int &, double &, int &,double & );
    static int select_step(double &, double &, int &, int & );
    static void birth_0(Matrix &,Matrix &,const double,const double,const double,const double,const int,const int,Matrix &,Matrix &,Matrix &,Matrix &,double &,const double &,double &,double &);
    static void death_0(Matrix &,Matrix &,const int,const int,Matrix &,Matrix &,Matrix &,Matrix &,double &,const double &,double &,double &);
    static void move_0(Matrix &,Matrix &,const double,const double,const double,const double,const int,const int,Matrix &,Matrix &,Matrix &,Matrix &,double &,const double &,double &,double &);
    static void cov_perturb_0(Matrix &,Matrix &,const int, const int,Matrix &,Matrix &,const int,Matrix &,Matrix &,double &,const double &,double &,double &);
//    static tuple<Matrix, double> birth_0(Matrix &,Matrix &,const double,const double,const double,const double,const int,const int,Matrix &,Matrix &,Matrix &,double &,const double &,Matrix &,double &);
//    static tuple<Matrix, double> death_0(Matrix &,Matrix &,const int,const int,Matrix &,Matrix &,Matrix &,double &,const double &,Matrix &,double &);
//    static tuple<Matrix, double> move_0(Matrix &,Matrix &,const double,const double,const double,const double,const int,const int,Matrix &,Matrix &,Matrix &,double &,const double &,Matrix &,double &);
//    static tuple<double, double> cov_perturb_0(Matrix &,Matrix &,const int , const int ,Matrix &,Matrix &,const int ,Matrix &,Matrix &,double &,const double &,double &,double &);









};
#endif /* defined(_Project_matrix_) */




// Constructor for Any Matrix
Matrix::Matrix(unsigned rowSize, unsigned colSize)
{
    m_rowSize = rowSize;
    m_colSize = colSize;
    m_matrix.resize(rowSize);
    for (unsigned i = 0; i < m_matrix.size(); ++i)
    {
        m_matrix[i].resize(colSize);
    }
}

// Constructor for Given Matrix
Matrix::Matrix(const char * fileName)
{
    ifstream file_A(fileName); // input file stream to open the file A.txt

    // Task 1
    // Keeps track of the Column and row sizes
    int colSize = 0;
    int rowSize = 0;

    // read it as a vector
    string line_A;
    int idx = 0;
    double element_A;
    vector <double> vector_A;

    if (file_A.is_open() && file_A.good())
    {
        // cout << "File A.txt is open. \n";
        while (getline(file_A, line_A))
        {
            rowSize += 1;
            stringstream stream_A(line_A);
            colSize = 0;
            while (1)
            {
                stream_A >> element_A;
                if (!stream_A)
                    break;
                colSize += 1;
                vector_A.push_back  (element_A);
                idx += 1;
            }
        }
    }
    else
    {
        cout << " WTF! failed to open. \n";
    }

    int j;
    idx = 0;
    m_matrix.resize(rowSize);
    for (unsigned i = 0; i < m_matrix.size(); ++i)
    {
        m_matrix[i].resize(colSize);
    }
    for (int i = 0; i < rowSize; ++i)
    {
        for (j = 0; j < colSize; ++j)
        {
            this->m_matrix[i][j] = vector_A[idx];
            idx++;
        }
    }
    m_colSize = colSize;
    m_rowSize = rowSize;

//    delete [] vector_A; // Tying up loose ends
}

// Copy Constructor
Matrix::Matrix(const Matrix &B)
{
    this->m_colSize = B.getCols();
    this->m_rowSize = B.getRows();
    this->m_matrix = B.m_matrix;

}

Matrix::~Matrix()
{

}



// Addition of Two Matrices
Matrix Matrix::operator+(Matrix &B)
{
    Matrix sum(m_rowSize, m_colSize);
    unsigned i,j;
    for (i = 0; i < m_rowSize; ++i)
    {
        for (j = 0; j < m_colSize; ++j)
        {
            sum(i,j) = this->m_matrix[i][j] + B(i,j);
        }
    }
    return sum;
}

// Subtraction of Two Matrices
Matrix Matrix::operator-(Matrix & B)
{
    Matrix diff(m_rowSize, m_colSize);
    unsigned i,j;
    for (i = 0; i < m_rowSize; ++i)
    {
        for (j = 0; j < m_colSize; ++j)
        {
            diff(i,j) = this->m_matrix[i][j] - B(i,j);
        }
    }

    return diff;
}

// Multiplication of Two Matrices
Matrix Matrix::operator*(Matrix & B)
{
    Matrix multip(m_rowSize,B.getCols());
    if(m_colSize == B.getRows())
    {
        unsigned i,j,k;
        double temp = 0.0;
        for (i = 0; i < m_rowSize; ++i)
        {
            for (j = 0; j < B.getCols(); ++j)
            {
                temp = 0.0;
                for (k = 0; k < m_colSize; ++k)
                {
                    temp += m_matrix[i][k] * B(k,j);
                }
                multip(i,j) = temp;
                //cout << multip(i,j) << " ";
            }
            //cout << endl;
        }
        return multip;
    }
    else
    {
        return "Error";
    }
}





// Scalar Addition
Matrix Matrix::operator+(double scalar)
{
    Matrix result(m_rowSize,m_colSize);
    unsigned i,j;
    for (i = 0; i < m_rowSize; ++i)
    {
        for (j = 0; j < m_colSize; ++j)
        {
            result(i,j) = this->m_matrix[i][j] + scalar;
        }
    }
    return result;
}

// Scalar Subraction
Matrix Matrix::operator-(double scalar)
{
    Matrix result(m_rowSize,m_colSize);
    unsigned i,j;
    for (i = 0; i < m_rowSize; ++i)
    {
        for (j = 0; j < m_colSize; ++j)
        {
            result(i,j) = this->m_matrix[i][j] - scalar;
        }
    }
    return result;
}

// Scalar Multiplication
Matrix Matrix::operator*(double scalar)
{
    Matrix result(m_rowSize,m_colSize);
    unsigned i,j;
    for (i = 0; i < m_rowSize; ++i)
    {
        for (j = 0; j < m_colSize; ++j)
        {
            result(i,j) = this->m_matrix[i][j] * scalar;
        }
    }
    return result;
}

// Scalar Division
Matrix Matrix::operator/(double scalar)
{
    Matrix result(m_rowSize,m_colSize);
    unsigned i,j;
    for (i = 0; i < m_rowSize; ++i)
    {
        for (j = 0; j < m_colSize; ++j)
        {
            result(i,j) = this->m_matrix[i][j] / scalar;
        }
    }
    return result;
}

// Scalar Power
Matrix Matrix::operator^(double scalar)
{
    Matrix result(m_rowSize,m_colSize);
    unsigned i,j;
    for (i = 0; i < m_rowSize; ++i)
    {
        for (j = 0; j < m_colSize; ++j)
        {
            result(i,j) = pow(this->m_matrix[i][j], scalar);
        }
    }
    return result;
}

// Returns value of given location when asked in the form A(x,y)
double& Matrix::operator()(const unsigned &rowNo, const unsigned & colNo)
{
    return this->m_matrix[rowNo][colNo];
}

// No brainer - returns row #
unsigned Matrix::getRows() const
{
    return this->m_rowSize;
}

// returns col #
unsigned Matrix::getCols() const
{
    return this->m_colSize;
}

// returns total num of matrix elements #
unsigned Matrix::getLength() const
{
    return this-> m_rowSize * m_colSize;
}

// Take any given matrices transpose and returns another matrix
Matrix Matrix::transpose()
{
    Matrix Transpose(m_colSize,m_rowSize);
    for (unsigned i = 0; i < m_colSize; ++i)
    {
        for (unsigned j = 0; j < m_rowSize; ++j)
        {
            Transpose(i,j) = this->m_matrix[j][i];
        }
    }
    return Transpose;
}

// vectorize a matrix
Matrix Matrix::vectorize()
{
    Matrix result(m_rowSize*m_colSize,1);
    for (unsigned j = 0; j < m_colSize; ++j)
    {
        for (unsigned i = 0; i < m_rowSize; ++i)
        {
            result(j*m_rowSize+i,0) = this->m_matrix[i][j];
        }
    }
    return result;
}

// absolute of a matrix
Matrix Matrix::absolute()
{
    Matrix result(m_rowSize,m_colSize);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            result(i,j) = fabs(this->m_matrix[i][j]);
        }
    }
    return result;
}


// Extract Rows and Columns of a Matrix
Matrix Matrix::MatExtraction(const int &r1,const int &r2,const int &c1,const int &c2)
{

    Matrix MatExtraction(r2-r1+1,c2-c1+1);
    for (int i = r1-1; i < r2; ++i)
    {
        for (int j = c1-1; j < c2; ++j)
        {
            MatExtraction(i-r1+1,j-c1+1) = this->m_matrix[i][j];
        }
    }
    return MatExtraction;
}

// Prints the matrix beautifully
void Matrix::print() const
{
    cout << "Matrix: " << endl;
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            cout << "[" << m_matrix[i][j] << "] ";
        }
        cout << endl;
    }
}

// Prints element of a matrix beautifully
void Matrix::print_element(int i, int j) const
{
    cout << "Matrix[" <<i<<"]"<<"["<<j<<"]: "<< m_matrix[i-1][j-1] << endl;
}

// linspace is similar to the colon operator, “:”, but gives direct control over the number of points and always includes the endpoints

// Generate linearly spaced vector
Matrix Matrix::linespace(double x1, double xn, int n)
{
    double dx = (xn-x1)/(n-1);
    Matrix A(n,1);
    for (int i=0; i<n; ++i)
    {
        A(i,0)=x1+(i*dx);
    }
    return A;
}

// linspace is similar to the colon operator, “:”, but gives direct control over the number of points and always includes the endpoints
Matrix Matrix::colon(const double &x1, const double &dx, const double &xn)
{

    vector <double> vector_A;

    int idx =0;
    while (x1+(idx*dx)<=xn)
    {

        vector_A.push_back  (x1+(idx*dx));
        idx += 1;
    }
    Matrix A (idx,1);
    int n=0;
    for(int i=0; i<idx; ++i)
    {
        A(i,0) = vector_A[n];
        n++;
    }
    return A;
}

// Generate logarithmically spaced vector
Matrix Matrix::logspace(double x1, double xn, int n)
{
    double dx = (xn-x1)/(n-1);

    Matrix A(n,1);
    for (int i=0; i<n; ++i)
    {
        A(i,0)=pow(10,(x1+(i*dx)));
    }
    return A;
}

// Normal Distribution
Matrix Matrix::normrnd(Matrix &mean, Matrix &sigma)
{
    std::random_device rd{};
    static std::mt19937 generator(std::time(0));

    Matrix result(mean.getRows(),mean.getCols());
    for (unsigned i=0; i<mean.getRows(); ++i)
    {
        for (unsigned j=0; j<mean.getCols(); ++j)
        {
            std::normal_distribution<double> distribution (mean(i,j),sigma(i,j));
            result(i,j) = distribution(generator);
        }
    }
    return result;
}


// cauchy distribution for inputted matrix
Matrix Matrix::cauchy_dist(Matrix &mean, Matrix &sigma)
{
    std::random_device rd{};
    static std::mt19937 generator(std::time(0));

    Matrix result(mean.getRows(),mean.getCols());
    for (unsigned i=0; i<mean.getRows(); ++i)
    {
        for (unsigned j=0; j<mean.getCols(); ++j)
        {
            std::cauchy_distribution<double> distribution (mean(i,j),sigma(i,j));
            result(i,j) = distribution(generator);
        }
    }
    return result;
}

// cauchy distribution for inputted double
//double Matrix::cauchy_dist(const unsigned  &mean, const unsigned  &sigma)
double Matrix::cauchy_dist( double  &mean,  double  &sigma)
{
//    std::random_device rd{};
//    static std::mt19937 generator(std::time(0));
//    std::cauchy_distribution<double> distribution (mean,sigma);
//    double result = distribution(generator);
//    return result;

    double r = Matrix::rand01();
    double x = mean+sigma*tan(PI*(r-0.5));
    return x;
}


void  Matrix::write2txtfile(const char* filename_X) const
{
    // open the file (for writing because it is an ofstream, meaning "output file stream")
    std::ofstream output(filename_X);

    for (unsigned i=0; i<m_rowSize; ++i)
    {
        for (unsigned j=0; j<m_colSize; ++j)
        {
            output << m_matrix[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output << "\n";
    }
// no need to close the file explicitly because it will be closed automatically by destructor of ofstream.
// However you can do that if you like using output.close()

}


Matrix Matrix::join_in_right(const initializer_list<Matrix> &A)
{

    int r;
    int c;
    vector <double> vector_A;
    int idx=0;
    for (auto mat: A)
    {
        r = mat.getRows();
        c = mat.getCols();
        for (int j=0; j<c ; ++j)
        {
            for (int i=0; i<r ; ++i)
            {
                vector_A.push_back  (mat(i,j));
                idx += 1;
            }
        }
    }

    c=idx/r;
    idx=0;
    Matrix output (r,c);
    for (int j=0; j<c ; ++j)
    {
        for (int i=0; i<r ; ++i)
        {
            output(i,j) = vector_A[idx];
            idx++;
        }
    }
    return output;
}

Matrix Matrix::join_in_bottom(const initializer_list<Matrix> &A)
{

    int r;
    int c;
    vector <double> vector_A;
    int idx=0;
    for (auto mat: A)
    {
        r = mat.getRows();
        c = mat.getCols();
        for (int i=0; i<r ; ++i)
        {
            for (int j=0; j<c ; ++j)
            {
                vector_A.push_back  (mat(i,j));
                idx += 1;
            }
        }
    }

    r=idx/c;
    idx=0;
    Matrix output (r,c);
    for (int i=0; i<r ; ++i)
    {
        for (int j=0; j<c ; ++j)
        {
            output(i,j) = vector_A[idx];
            idx++;
        }
    }
    return output;
}


//tuple<Matrix, Matrix> Matrix::meshgrid(Matrix &x, Matrix &z, Matrix &X, Matrix &Z,const int &NX,const int &NZ)
void Matrix::meshgrid(Matrix &x, Matrix &z, Matrix &X, Matrix &Z,const int &NX,const int &NZ)
{
    // Matrix X(NZ,NX);
    //  Matrix Z(NZ,NX);
    for (int i=0; i<NZ; ++i)
    {
        for (int j=0; j<NX; ++j)
        {
//            X[i][j]=x[j][0];
//            Z[i][j]=z[i][0];
            X(i,j)=x(j,0);
            Z(i,j)=z(i,0);
        }
    }
    // return make_tuple(X,Z);
}

// generate double random number between 0 and 1
Matrix Matrix::RandomDouble(unsigned rowSize, unsigned colSize)
{
    std::random_device rd{};
    static std::mt19937 generator(std::time(0));
    std::uniform_real_distribution<double> distribution(0, 1);
    Matrix r(rowSize,colSize);
    for (unsigned i=0; i<rowSize; ++i)
    {
        for (unsigned j=0; j<colSize; ++j)
        {
            r(i,j) = distribution(generator);

        }
    }

    return r;
}

// generate double random number between 0 and 1
double Matrix::rand01()
{

    double maxx = 100;
    double minn = 0;
    const unsigned range = maxx - minn + 1 ;
    return (rand() % range + minn)/maxx ;

//    std::random_device rd{};
//    static std::mt19937 generator(std::time(0));
//    std::uniform_real_distribution<double> distribution(0, 1);
//    Matrix r(rowSize,colSize);
//    for (unsigned i=0; i<rowSize; ++i)
//    {
//        for (unsigned j=0; j<colSize; ++j)
//        {
//            r(i,j) = distribution(generator);
//
//        }
//    }
//
//    return r;

}

// round number to closest int
Matrix Matrix::rounding()
{
    Matrix result(m_rowSize,m_colSize);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            result(i,j) = round(this->m_matrix[i][j]);
        }
    }
    return result;
}

Matrix Matrix::e2C(double & ec, Matrix &dg_obs_abs, const int n)
{
    // int n = dg_obs.getRows();
    // dg_hat=abs(dg_obs);
    // dg_hat(abs(dg_obs)<0.1)=0.1;
    Matrix T (n,n);
    T=T.eye();
    // Matrix sigma = dg_obs.absolute() * ec;  // Magnitude Scaling
    Matrix sigma = dg_obs_abs * ec;  // Magnitude Scaling
    Matrix T3 =sigma.transpose();
    Matrix T4 = sigma * T3;
    Matrix C = Matrix::dot(T,T4);
    return C;
}



// damping a matrix to make it positive definite
Matrix Matrix::damping()
{

    Matrix C (m_rowSize, m_colSize); // store result
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            C(i,j) = this->m_matrix[i][j];
        }
    }

    Matrix ind = Matrix::colon(0,1,m_rowSize-1);
    Matrix Dij=ind.toeplitz();
    Matrix T1 = Dij*(PI/(2*(m_rowSize-1)));
    Matrix T2 = T1.cos_mat();
    Matrix damp = T2 ^ 1.2 ;
    C = dot(C,damp);
    return C;
}


// diag of a matrix
Matrix Matrix::diag()
{
    Matrix result(m_rowSize,1);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {

        result(i,0) = this->m_matrix[i][i];

    }
    return result;
}



// Lower triangular matrix from cholesky decomposition
//Matrix Matrix::cholesky(const char * fileName)
tuple<Matrix, Matrix> Matrix::cholesky()
{

    const double EPSILON = 1.0E-12;
    if (m_rowSize == 0 || m_colSize == 0)
        cout <<"Inverse called on matrix with zero rows and/or columns"<<endl;

    if (m_rowSize != m_colSize)
        cout << "Non-square matrix in Cholesky Decomposition" << endl;

    // counters (better to declare in the loop decalaration)
//    int i, j, k;

    // temp object to hold the result
    Matrix A (m_rowSize, m_colSize); // store result
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j <= i; ++j)
        {
            A(i,j) = this->m_matrix[i][j];
        }
    }

    // check for positive definiteness
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        if (A(i,i) < 0.0)
            cout <<"Matrix not positive definite in Cholesky Decomposition" << endl;

        if (fabs (A(i,i)) < EPSILON)
            cout <<"Singular matrix in Cholesky Decomposition" <<endl;
    }

    // Perform Choleski decomposition
    for (unsigned j = 0; j < m_rowSize; ++j)
    {
        for(unsigned k = 0; k < j; ++k)
            A(j,j) -= A(j,k) * A(j,k);

        if (A(j,j) < 0.0)
            cout << "Square root of negative number in Cholesky Decomposition" << endl;

        A(j,j) = sqrt (A(j,j));

        for(unsigned i = j + 1; i < m_rowSize; ++i)
        {
            for (unsigned k = 0; k < j; ++k)
                A(i,j) -= A(i,k) * A(j,k);

            if (fabs(A(j,j)) < EPSILON)
                cout <<"Division by zero in Cholesky Decomposition" << endl;

            A(i,j) /= A(j,j);
        }
    }

    Matrix L (m_rowSize, m_colSize);
    Matrix U (m_rowSize, m_colSize);
    L=A;
    U=L.transpose();
    return make_tuple(L,U);

//    if (fileName == "Lower")
//    {
//        return A;
//
//    }
//    else if (fileName == "Upper")
//    {
//        return A.transpose();
//    }
//    else
//    {
//        cout <<"wrong string input" << endl;
//    }

}


// As Cholesky decomposition is twice as fast as LU-decomposition which is used to calculate general matrix determinants,
// it is recommended to use Cholesky decomposition when calculating a determinant of a symmetric positive definite matrix.
// In order to calculate the determinant of a symmetric matrix which is not positive definite algorithm on the basis of LDLT-decomposition can be used.


// matrix determinant from cholesky decomposition
double Matrix::det_cholesky()
{
    Matrix A (m_rowSize, m_colSize); // store result
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            A(i,j) = this->m_matrix[i][j];
        }
    }

    //Matrix L=A.cholesky("Lower");
    Matrix L (m_rowSize, m_colSize);
    Matrix U (m_rowSize, m_colSize);
    tie(L,U)=A.cholesky();
    Matrix DL = L.diag();
    double Pr = DL.prod();
    double det = pow(Pr, 2);
    return det;
}

// matrix log-determinant from cholesky decomposition
double Matrix::log_det_cholesky()
{
    Matrix A (m_rowSize, m_colSize); // store result
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            A(i,j) = this->m_matrix[i][j];
        }
    }

    Matrix L (m_rowSize, m_colSize);
    Matrix U (m_rowSize, m_colSize);
    tie(L,U)=A.cholesky();
    Matrix diagL = L.diag();
    Matrix Log_diagL = diagL.log_fun();
    double log_det = Log_diagL.sum();
    log_det = log_det* 2;
    return log_det;
}


// Cholesky-based matrix inversion has several benefits over LU-based one.
// First, instead of two factors (L and U) we now have only one triangular factor to invert. Less factors = less work.
// Second, there is no more row permutation matrix P. Row permutations are essential for the stability of LU decomposition,
// but Cholesky factorization is perfectly stable even without them (as long as your matrix is strictly positive definite).
// These permutations add significant overhead to algorithm running time because of two reasons:
// a) they trash CPU cache with non-local memory accesses, and b) they worsen parallelism potential by adding one more synchronization point to the algorithm.



// inverse of matrix by cholesky decomposition
Matrix Matrix::inv_cholesky()
{
    Matrix A (m_rowSize, m_colSize); // store result
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            A(i,j) = this->m_matrix[i][j];
        }
    }

    //A=A.cholesky("Lower");
    Matrix L (m_rowSize, m_colSize);
    Matrix U (m_rowSize, m_colSize);
    tie(L,U)=A.cholesky();
//    int i, j, k;
    // inversion of lower triangular matrix
    for (unsigned j = 0; j < m_rowSize; ++j)
    {
        L(j,j) = 1.0 / L(j,j);

        for (unsigned i = j + 1; i < m_rowSize; ++i)
        {
            L(i,j) *= -L(j,j) / L(i,i);

            for (unsigned k = j + 1; k < i; ++k)
                L(i,j) -= L(i,k) * L(k,j) / L(i,i);
        }
    }

    // construction of lower triangular inverse matrix
    for (unsigned j = 0; j < m_rowSize; ++j)
    {
        for (unsigned i = j; i < m_rowSize; ++i)
        {
            L(i,j) *= L(i,i);

            for (unsigned k = i + 1; k < m_rowSize; ++k)
                L(i,j) += L(k,i) * L(k,j);
        }
    }

    // fill upper diagonal
    for (unsigned i = 1; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < i; ++j)
            L(j,i) = L(i,j);
    }
    return L; // return L copy



}

// LU decomposition can be used for matrix inversion and determinent exactly like what Matlab does!
// However, Cholesky decomposion is also more useful to do these.
tuple<Matrix, Matrix> Matrix::LUdecomposition()
{
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    Matrix lower (m_rowSize, m_colSize); // store result
    Matrix upper (m_rowSize, m_colSize); // store result
    int n=m_rowSize;
    int i, j, k;
    for ( i = 0; i < n; ++i)
    {

        // Upper Triangular
        for ( k = i; k < n; ++k)
        {

            // Summation of L(i, j) * U(j, k)
            int sum = 0;
            for ( j = 0; j < i; ++j)
                sum += (lower(i,j) * upper(j,k));

            // Evaluating U(i, k)
            upper(i,k) = this->m_matrix[i][k] - sum;
        }

        // Lower Triangular
        for ( k = i; k < n; ++k)
        {
            if (i == k)
                lower(i,i) = 1; // Diagonal as 1
            else
            {

                // Summation of L(k, j) * U(j, i)
                int sum = 0;
                for ( j = 0; j < i; ++j)
                    sum += (lower(k,j) * upper(j,i));

                // Evaluating L(k, i)
                lower(k,i) = (this->m_matrix[k][i] - sum) / upper(i,i);
            }
        }
    }

    return make_tuple(lower,upper);

}
// inverse of matrix by LU decomposition (The following code actually performs the matrix inversion with same result from Matlab)
Matrix Matrix::inv_LU()
{

    Matrix A (m_rowSize, m_colSize); // store result
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            A(i,j) = this->m_matrix[i][j];
        }
    }

    int r=m_rowSize;
    int c=m_colSize;
    if (r != c)
    {
        cout << "Only Square Matrices, please" <<endl;
    }
    Matrix b (m_rowSize, m_colSize); // store result
    b=b.eye();
    for (int j=0; j<r ; ++j)
    {
        for (int i=j; i<r; ++i)
        {
            if (A(i,j)!=0)
            {
                for (int k=0; k<r; ++k)
                {
                    double s = A(j,k);
                    A(j,k) = A(i,k);
                    A(i,k) = s;
                    s = b(j,k);
                    b(j,k) = b(i,k);
                    b(i,k) = s;
                }
                double t = 1/A(j,j);
                for (int k=0; k<r ; ++k)
                {
                    A(j,k) = t * A(j,k);
                    b(j,k) = t * b(j,k);
                }
                for (int L=0; L<r; ++L)
                {
                    if (L != j)
                    {
                        t= A(L,j)* -1;
                        for (int k=0; k<r; ++k)
                        {
                            A(L,k) = A(L,k) + t * A(j,k);
                            b(L,k) = b(L,k) + t * b(j,k);

                        }
                    }
                }
            }
            break;
        }
//         if (A(i,j)==0)
//         {
//             cout << "Warning: Singular Matrix" <<endl;
//         }
    }



    return b;


}



// filter function

Matrix Matrix::filter(Matrix & b,Matrix & a,Matrix & x)
{
    Matrix a1 =a/a(0,0);
    Matrix b1 = b/a(0,0);
    unsigned sx = x.getLength();
    unsigned rowx = x.getRows();
    unsigned colx = x.getCols();

    Matrix output (rowx,colx);
    output(0,0) = b1(0,0)*x(0,0);
    for (unsigned i = 1; i < sx; ++i)
    {
        output(i,0) = 0;
        for (unsigned j = 0; j <= i; ++j)
        {
            unsigned k = i-j;
            if (j > 0)
            {
                if ((k < b1.getLength()) && (j < x.getLength()))
                {
                    output(i,0) += b1(k,0)*x(j,0);
                }
                if ((k < output.getLength()) && (j < a1.getLength()))
                {
                    output(i,0) -= a1(j,0)*output(k,0);
                }
            }
            else
            {
                if ((k < b1.getLength()) && (j < x.getLength()))
                {
                    output(i,0) += (b1(k,0)*x(j,0));
                }
            }
        }
    }
    return output;
}

// log of a matrix
Matrix Matrix::log_fun()
{
    Matrix result(m_rowSize,m_colSize);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            result(i,j) = log(this->m_matrix[i][j]);
        }
    }
    return result;
}

// sum of all elements of a matrix
double Matrix::sum()
{
    double result=0;
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            result = result + this->m_matrix[i][j];
        }
    }
    return result;
}

// Min of a matrix
double Matrix::Min()
{
    double result=m_matrix[0][0];
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            if (result > this->m_matrix[i][j])
                result = this->m_matrix[i][j];
        }
    }
    return result;
}

// Max of a matrix
double Matrix::Max()
{
    double result=m_matrix[0][0];
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            if (result < this->m_matrix[i][j])
                result = this->m_matrix[i][j];
        }
    }
    return result;
}

// Knn Search
Matrix Matrix::XZ_model(Matrix & Bc,Matrix & X, Matrix & Z, const int &NX, const int &NZ)
{
    Matrix IDX (NZ,NX);
    double dist;
    double d;
    for (int i=0; i<NZ; ++i)
    {
        for (int j=0; j< NX; ++j)
        {
            dist=sqrt(pow(Bc(0,0)-X(i,j),2) + pow(Bc(0,1)-Z(i,j),2));

            for (unsigned k=0; k<Bc.getRows(); ++k)
            {
                d = sqrt(pow(Bc(k,0)-X(i,j),2) + pow(Bc(k,1)-Z(i,j),2));
                if ( d <= dist)
                {
                    dist = d;
                    IDX(i,j)=Bc(k,2);
                }
            }
        }
    }
    return IDX;
}
//Matrix Matrix::XZ_model(Matrix & Bc,Matrix & X, Matrix & Z)
//{
//    Matrix xv = X.vectorize();
//    Matrix zv = Z.vectorize();
//    Matrix B = Bc.MatExtraction(1,Bc.getRows(),1,2);
//    Matrix IDX (xv.getRows(),xv.getCols());
//    double dist;
//    double d;
//    for (unsigned i=0; i<xv.getRows(); ++i)
//    {
//        dist=sqrt(pow(B(0,0)-xv(i,0),2) + pow(B(0,1)-zv(i,0),2));
//        for (unsigned j=0; j<B.getRows(); ++j)
//        {
//            d = sqrt(pow(B(j,0)-xv(i,0),2) + pow(B(j,1)-zv(i,0),2));
//            if ( d <= dist)
//            {
//                dist = d;
//                IDX(i,0)=Bc(j,2);
//            }
//        }
//    }
//
//    Matrix cc = IDX.reshape(X.getRows(),X.getCols());
//    return cc;
//}

// reshape
Matrix Matrix::reshape(int r, int c)
{
    Matrix result(r,c);
    int cnt= 0;

    for (unsigned j = 0; j < m_colSize; ++j)
    {
        for (unsigned i = 0; i < m_rowSize; ++i)
        {

            result(cnt % r,cnt / r) = this->m_matrix[i][j];
            cnt ++;
        }

    }
    return result;
}

// reshape
Matrix Matrix::FM_fast(Matrix & cc,Matrix & kernel, const double &density)
{
// density must be in kg/m3
    Matrix c3=cc.vectorize();
    c3=c3*density;
    Matrix dg_pre=kernel*c3;
    dg_pre=dg_pre*100000; // dg is mGal
    return dg_pre;
}

// exp function for an inputed matrix
Matrix Matrix::exp_mat()
{
    Matrix result(m_rowSize,m_colSize);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for  (unsigned j = 0; j < m_colSize; ++j)
        {
            result(i,j)=exp(this->m_matrix[i][j]);

        }

    }
    return result;
}

// cos function for an inputed matrix
Matrix Matrix::cos_mat()
{
    Matrix result(m_rowSize,m_colSize);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for  (unsigned j = 0; j < m_colSize; ++j)
        {
            result(i,j)=cos(this->m_matrix[i][j]);

        }

    }
    return result;
}


// one matrix
Matrix Matrix::one()
{
    Matrix result(m_rowSize,m_colSize);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for  (unsigned j = 0; j < m_colSize; ++j)
        {
            result(i,j)=1;

        }

    }
    return result;
}

// identical matrix
Matrix Matrix::eye()
{
    Matrix result(m_rowSize,m_colSize);
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for  (unsigned j = 0; j < m_colSize; ++j)
        {
            result(i,i)=1;

        }

    }
    return result;
}

// toeplitz
Matrix Matrix::toeplitz()
{
    Matrix result(m_rowSize,m_rowSize);
    int n = m_rowSize;
    for (int i = 0; i < n; ++i)
    {
        for  (int j = 0; j < n; ++j)
        {
            if(i>=j)
            {
                result(i,j) = this->m_matrix[(n+i-j)%n][0];
            }
            else
            {
                result(i,j) = this->m_matrix[(n-i+j)%n][0];
            }

        }

    }
    return result;
}


// dot product of two matrices
Matrix Matrix::dot( Matrix & a, Matrix  &b)
{
    Matrix result(a.getRows(), a.getCols());
    for (unsigned i = 0; i < a.getRows(); ++i)
    {
        for (unsigned j = 0; j < a.getCols(); ++j)
        {
            result(i,j) = a(i,j) * b(i,j);
        }
    }
    return result;
}

// Product of array elements (gives a double answer)
double Matrix::prod()
{
    double result = 1;
    for (unsigned i = 0; i < m_rowSize; ++i)
    {
        for (unsigned j = 0; j < m_colSize; ++j)
        {
            result = result * this->m_matrix[i][j];
        }
    }
    return result;
}

// dot_division of two matrices
Matrix Matrix::dot_division(Matrix & a,Matrix & b)
{
    Matrix result(a.getRows(), a.getCols());
    for (unsigned i = 0; i < a.getRows(); ++i)
    {
        for (unsigned j = 0; j < a.getCols(); ++j)
        {
            result(i,j) = a(i,j) / b(i,j);
        }
    }
    return result;
}

// Compute LogL for AR inversion
double Matrix::Log_Likelihood_AR(Matrix & r,Matrix & dg_obs,Matrix & AR_parameters)
{
    double e = AR_parameters(0,0);
    Matrix sigma = dg_obs.absolute() * e;
    int row = AR_parameters.getRows();
    Matrix b (row,1);
    b(0,0)=0;
    for (int i=1; i<row; ++i)
    {
        b(i,0)=AR_parameters(i,0);
    }
    Matrix a (1,1);
    a(0,0)=1;
    Matrix da=Matrix::filter(b,a,r);
    double LogL;
    Matrix t1 = sigma.log_fun();
    Matrix t2 = Matrix::dot_division(r,sigma);
    Matrix t3 = t2 ^ 2;
    Matrix t4 = dot_division(da,sigma);
    Matrix t5 = t4 ^ 2;
    Matrix t6 = r * 2;
    Matrix t7 = dot(t6,da);
    Matrix t8 = sigma^2;
    Matrix t9 = dot_division(t7,t8);

    LogL = (-1 * t1.sum()) -(0.5 * t3.sum()) -(0.5 * t5.sum())+(0.5 * t9.sum());

    return LogL;

}

// Compute LogL for a full Cov Matrix
double Matrix::Log_Likelihood_full_C(Matrix & r,Matrix & C)
{
    Matrix rt= r.transpose();
    Matrix Cinv = C.inv_cholesky();
    //Matrix Cinv = C.inv_LU();
    Matrix rtCinv = rt * Cinv;
    double LogL = (rtCinv * r)(0,0) * (-0.5);
    return LogL;

}

// Compute LogL for a full Cov Matrix at initial inversion
double Matrix::Log_Likelihood_initial(Matrix & r,Matrix & dg_obs_abs,double & ec)
{
    Matrix sigma = dg_obs_abs * ec;
    Matrix LogS = sigma.log_fun();
    double sum1 = LogS.sum();
    sum1 *=-1;
    Matrix rs = Matrix::dot_division(r,sigma);
    rs =rs^2;
    double sum2 = rs.sum();
    sum2 *= -0.5;
    return sum1 + sum2;

    //double LogL1 = C.log_det_cholesky() * (-0.5);
    //double LogL2 = Matrix::Log_Likelihood_full_C(r,Cinv);
    //return logdetC + LogL2;
}

Matrix Matrix::deflation(Matrix &X, double &eigenvalue)
{
    // Deflation formula exactly applied
    double denominator = eigenvalue / (X.transpose() * X)(0,0);
    Matrix Xtrans = X.transpose();
    Matrix RHS = (X * Xtrans);
    Matrix RHS2 = RHS * denominator;
    Matrix A2 = *this - RHS2;
    return A2;
}

// Initial MCMC
void Matrix::MCMC_0( Matrix & X,  Matrix & Z,const double x_min,const double x_max,const double z_min,const double z_max,const int NX, const int NZ,  Matrix & dg_obs,Matrix & dg_obs_abs, const int num_obs,  Matrix & kernel, Matrix & tab, const  double & density, int & Nmin, int & Nmax, double & beta, int & k, double & TrueLogL)
{
    double rstep = Matrix::rand01();
    Matrix Bc = tab.MatExtraction(1,1,4,3+(tab(0,1)*3));
    Bc = Bc.transpose();
    Bc = Bc.reshape(tab(0,1),3);
    double ec=tab(0,2);
    //Matrix C = Matrix::e2C(ec,dg_obs_abs,num_obs);
    int step= Matrix::select_step(rstep,tab(0,1),Nmin,Nmax);
    double LogLc=tab(0,0);
    switch(step)
    {
    case 91 :
        Matrix::birth_0(X,Z,x_min,x_max,z_min,z_max,NX,NZ,dg_obs,dg_obs_abs,kernel,Bc,LogLc,density,ec,beta);
        break;
    case 92 :
        Matrix::death_0(X,Z,NX,NZ,dg_obs,dg_obs_abs,kernel,Bc,LogLc,density,ec,beta);
        break;
    default:
        // Matrix::move_0(X,Z,x_min,x_max,z_min,z_max,NX,NZ,dg_obs,dg_obs_abs,kernel,Bc,LogLc,density,ec,beta);
        break;
    }

    Matrix::cov_perturb_0(X,Z,NX,NZ,dg_obs,dg_obs_abs,num_obs,kernel,Bc,LogLc,density,ec,beta);
    Matrix F_TAB (1,3);
    F_TAB (0,0)=LogLc ;
    F_TAB (0,1)=Bc.getRows() ;
    F_TAB (0,2)=ec ;
    Matrix zeros(1,(Nmax-Bc.getRows())*3);
    tab = Matrix::join_in_right({F_TAB,Bc.vectorize().transpose(),zeros});
    //if ((k%100) ==0)
    // {
    //cout << "T: " << beta << ", Iteration: " << k << ", LogL: " << LogLc << ", NumofNodes: " << Bc.getRows() << ", e: " << ec <<", TrueLogL: "<< TrueLogL<< endl;
    // printf ("T: %f  \n", beta, ", Iteration: %f  \n", k);
    printf("T:  %f , Iteration:  %i , LogL:  %f , NumofNodes: %i , e: %f , TrueLogL: %f \n",beta,k,LogLc,Bc.getRows(),ec,TrueLogL);
    // }
    //return tab;
}

// birth for initial inversion
void Matrix::birth_0(Matrix &X,Matrix &Z,const double x_min,const double x_max,const double z_min,const double z_max,const int NX, const int NZ, Matrix &dg_obs,Matrix &dg_obs_abs,Matrix &kernel,Matrix &Bc,double &LogLc,const double &density,double &ec,double &beta)
{

    double xrp =  Matrix::rand01()*(x_max-x_min) + x_min;
    double zrp =  Matrix::rand01()*(z_max-z_min) + z_min;
    Matrix PP (1,3);
    PP(0,0)=xrp;
    PP(0,1)=zrp;
    PP(0,2)=round(Matrix::rand01());
    Matrix Bp = Matrix::join_in_bottom({Bc,PP});
    Matrix cp = Matrix::XZ_model(Bp, X, Z,NX,NZ);
    Matrix dp = Matrix::FM_fast(cp,kernel,density);
    Matrix res = dg_obs - dp;
    double LogLp = Matrix::Log_Likelihood_initial(res,dg_obs_abs,ec);
    double p = exp((LogLp - LogLc)/beta);
    if (Matrix::rand01() <= p)
    {
        Bc = Bp;
        LogLc = LogLp ;
    }

    // return make_tuple(Bc,LogLc);
}

// death for initial inversion
void Matrix::death_0(Matrix &X,Matrix &Z,const int NX, const int NZ,Matrix &dg_obs,Matrix &dg_obs_abs,Matrix &kernel,Matrix &Bc,double &LogLc,const double &density,double &ec,double &beta)
{

    double r =  Matrix::rand01();
    int N = Bc.getRows();
    //double idouble = (r * (N-1)) +1 ;
    int i = round((r * (N-1)) +1);
    Matrix Bp = Bc.DeleteRow(i);
    Matrix cp = Matrix::XZ_model(Bp, X, Z,NX,NZ);
    Matrix dp = Matrix::FM_fast(cp,kernel,density);
    Matrix res = dg_obs - dp;
    double LogLp = Matrix::Log_Likelihood_initial(res,dg_obs_abs,ec);
    double p = exp((LogLp - LogLc)/beta);
    if (Matrix::rand01() <= p)
    {
        Bc = Bp;
        LogLc = LogLp ;
    }

    //return make_tuple(Bc,LogLc);
}

// death for initial inversion
void Matrix::move_0(Matrix &X,Matrix &Z,const double x_min,const double x_max,const double z_min,const double z_max,const int NX, const int NZ,Matrix &dg_obs,Matrix &dg_obs_abs,Matrix &kernel,Matrix &Bc,double &LogLc,const double &density,double &ec,double &beta)
{

    double dx = abs(X(0,1)-X(0,0));
    double dz = abs(Z(1,0)-Z(0,0));
    Matrix Bp (Bc.getRows(),Bc.getCols());
    for (unsigned i=0; i<Bc.getRows(); ++i)
    {
        for (unsigned k=0; k<Bc.getCols(); ++k)
        {
            Bp = Bc;
            switch (k)
            {
            case 0 :

                Bp(i,k)=Matrix::cauchy_dist(Bc(i,k),dx);
                while (Bp(i,k)<x_min || Bp(i,k)>x_max)
                {
                    Bp(i,k)=Matrix::cauchy_dist(Bc(i,k),dx);
                }

//                 teta =  Matrix::RandomDouble(1,1)(0,0)*(0-2*PI) + 0;
//                 radius =  Matrix::RandomDouble(1,1)(0,0)*(0-2*dx) + 0;
//                Bp(i,k)=Bc(i,k) + radius * cos(teta);
//                while (Bp(i,k)<x_min || Bp(i,k)>x_max)
//                {
//                     teta =  Matrix::RandomDouble(1,1)(0,0)*(0-2*PI) + 0;
//                     radius =  Matrix::RandomDouble(1,1)(0,0)*(0-2*dx) + 0;
//                    Bp(i,k)=Bc(i,k) + radius * cos(teta);
//                }


                break;
            case 1 :


                Bp(i,k)=Matrix::cauchy_dist(Bc(i,k),dz);
                while (Bp(i,k)<z_min || Bp(i,k)>z_max)
                {
                    Bp(i,k)=Matrix::cauchy_dist(Bc(i,k),dz);
                }


//                 teta =  Matrix::RandomDouble(1,1)(0,0)*(0-2*PI) + 0;
//                 radius =  Matrix::RandomDouble(1,1)(0,0)*(0-2*dz) + 0;
//                Bp(i,k)=Bc(i,k) + radius * sin(teta);
//                while (Bp(i,k)<z_min || Bp(i,k)>z_max)
//                {
//                     teta =  Matrix::RandomDouble(1,1)(0,0)*(0-2*PI) + 0;
//                     radius =  Matrix::RandomDouble(1,1)(0,0)*(0-2*dz) + 0;
//                    Bp(i,k)=Bc(i,k) + radius * sin(teta);
//                }


                break;
            default :
                Bp(i,k)=1-Bc(i,k);
                break;

            }

            Matrix cp = Matrix::XZ_model(Bp, X, Z, NX, NZ);
            Matrix dp = Matrix::FM_fast(cp,kernel,density);
            Matrix res = dg_obs - dp;
            double LogLp = Matrix::Log_Likelihood_initial(res,dg_obs_abs,ec);
            double p = exp((LogLp - LogLc)/beta);
            if (Matrix::rand01() <= p)
            {
                Bc = Bp;
                LogLc = LogLp ;
            }
        }
    }
    // return make_tuple(Bc,LogLc);
}



// cov purt for initial inversion
void Matrix::cov_perturb_0(Matrix &X,Matrix &Z,const int NX, const int NZ,Matrix &dg_obs,Matrix & dg_obs_abs, const int num_obs,Matrix &kernel,Matrix &Bc,double &LogLc,const double &density,double &ec,double &beta)
{

    Matrix cc = Matrix::XZ_model(Bc, X, Z,NX,NZ);
    Matrix dc = Matrix::FM_fast(cc,kernel,density);
    Matrix res = dg_obs - dc;
    double deltaE = 0.001;
    double ep =Matrix::cauchy_dist(ec,deltaE);
    while (ep < 0.00001 || ep>1)
    {
        ep =Matrix::cauchy_dist(ec,deltaE);
    }
    // Matrix C = Matrix::e2C(ep,dg_obs_abs,num_obs);
    double LogLp = Matrix::Log_Likelihood_initial(res,dg_obs_abs,ep);
    double p = exp((LogLp - LogLc)/beta);
    if (Matrix::rand01() <= p)
    {
        LogLc = LogLp ;
        ec = ep;
    }

    //return make_tuple(LogLc,ec);
}


// delete a specific row from a given matrix
Matrix Matrix::DeleteRow(const unsigned &r)
{
    Matrix result (m_rowSize-1,m_colSize);

    for (unsigned i=0; i<m_rowSize-1; ++i)
    {
        for (unsigned j =0; j<m_colSize; ++j)
        {
            if (i+1 >= r)
            {
                result(i,j)=this->m_matrix[i+1][j];
            }
            else
            {
                result(i,j)=this->m_matrix[i][j];
            }
        }
    }
    return result;
}


// delete a specific col from a given matrix
Matrix Matrix::DeleteCol(const unsigned &c)
{
    Matrix result (m_rowSize,m_colSize-1);

    for (unsigned i=0; i<m_rowSize; ++i)
    {
        for (unsigned j =0; j<m_colSize-1; ++j)
        {
            if (j+1 >= c)
            {
                result(i,j)=this->m_matrix[i][j+1];
            }
            else
            {
                result(i,j)=this->m_matrix[i][j];
            }
        }
    }
    return result;
}



// Select Step for Birth, Death, Move
int Matrix::select_step(double &r, double &n, int & Nmin, int & Nmax)
{
    int step ;
    if (n==Nmin) //only birth and move
    {
        if (r <= 0.3)
        {
            step = 91; // birth
        }
        else
        {
            step = 93; //move
        }
    }
    else if (n==Nmax)
    {
        if (r <= 0.3)
        {
            step = 92; // death
        }
        else
        {
            step = 93; //move
        }
    }
    else
    {
        if (r <= 0.3)
        {
            step = 91; // birth
        }
        else if ( r > 0.3 && r <=0.6)
        {
            step = 92; //death
        }
        else
        {
            step =93;
        }

    }
    return step;

}

