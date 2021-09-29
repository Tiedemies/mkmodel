#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::matrix<double> Matrix;

void ComputeHouseHolder(Matrix& Qv, Matrix v&)
{
    BOOST_ASSERT(Qv.size1() == v.size1() && Qv.size2() == v.size1());
    const int& n = v.size1();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            Qv(i,j) = -2*v(i,0)*v(j,0);
        }
    }
    for (int i = 0; i < n; ++i)
    {
        Qv(i,i) += 1;
    }
}

void Rescale(Matrix &x)
{
    double norm = Norm(x);
    for (int i = 0; i < x.size1();++i)
    {
        x(i,0) /= norm;
    }
}

Matrix Vmadd(const Matrix& a, const Matrix& b, double s)
{
    if(a.size1 != b.size1())
    {
        throw std::runtime_error("mismatched dimenstions");
    }
    Matrix c(a.size1(),0);
    for (int i = 0; i < a.size1();++i)
    {
        c(i,0) = a(i,0) + s*b(i,0);
    }
    return c;
}

double Norm(const Matrix& x) 
{
    double res = 0.0;
    for(int i = 0; i < x.size1(); ++i)
    {
        res += x(i,0)*x(i,0);
    }
    return res; 
}


Matrix Minor(const Matrix& M, int d)
{
    Matrix togo(M.size1(), M.size2());
    for (int i = 0; i < d; ++i)
    {
        togo(i,i) = 1.0;
    }
    for (int i = d; i < M.size1(); ++i)
    {
        for (int j = d; j < M.size2(); ++j)
        {
            togo(i,j) = M(i,j);
        }
    }
    return togo;
}

Matrix ExtractColumn(Matrix& M, int d)
{
    BOOST_ASSERT(0 <= d && d < M.size1());
    Matrix togo(M.size1(),1);
    for (int i = 0; i < M.size1(); ++i)
    {
        togo(i,0) = M(i,d);
    }
    return togo; 
}


// QR decomposition
std::pair<Matrix,Matrix> QRdecomp(Matrix& A)
{
    int m = A.size1();
    int n = A.size2(); 
    Matrix Z = A; 
    Matrix Q(m,m, 0.0);
    for (int i = 0; i < m; ++i){Q(i,i) = 1.0;}
    Matrix R(m,n);
    for (int k = 0; k < n && k < m - 1; ++k)
    {
        Matrix Qv(m,m);
        Matrix e(m,1, 0.0);
        Matrix z1 = Minor(Z,k);
        Matrix x = ExtractColumn(z1,k);
        double a = Norm(x);
        if (A(k,k) > 0)
        {
            a = -a;
        }
        e(k,0) = 1.0; 
        e = Vmadd(x,e,a);
        Rescale(e);
        ComputeHouseHolder(Qv,e);
        Z = boost::numeric::ublas::prod(Qv, z1)
    }
}

Matrix Regress(const Matrix& X, const Matrix& y)
{
    //OLS fit for y = Xb, return b
    
}
#endif