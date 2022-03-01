#ifndef SPECTRA_MY_GEN_MAT_PROD_H
#define SPECTRA_MY_GEN_MAT_PROD_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/KroneckerProduct>



namespace Spectra {


    template <typename Scalar_, typename TypeA = Eigen::Sparse, typename TypeB = Eigen::Sparse,
            int FlagsA = Eigen::ColMajor, int FlagsB = Eigen::ColMajor,
            typename StorageIndexA = int, typename StorageIndexB = int>
    class MySparseGenMatProd
    {
    public:
        using Scalar = Scalar_;

    private:
        using Index = Eigen::Index;
        using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        using MapConstVec = Eigen::Map<const Vector>;
        using MapVec = Eigen::Map<Vector>;
        using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using MapMat = Eigen::Map<Matrix>;
        using SparseMatrixA = Eigen::SparseMatrix<Scalar, FlagsA, StorageIndexA>;
        using SparseMatrixB = Eigen::SparseMatrix<Scalar, FlagsB, StorageIndexB>;

        using ConstGenericSparseMatrixA = const Eigen::Ref<const SparseMatrixA>;
        using ConstGenericSparseMatrixB = const Eigen::Ref<const SparseMatrixB>;

        using FacTypeA = Eigen::SparseLU<SparseMatrixA>;
        using FacTypeB = Eigen::SparseLU<SparseMatrixB>;

        FacTypeA a_solver;
        FacTypeB b_solver;
        ConstGenericSparseMatrixA m_mat;
        ConstGenericSparseMatrixA a_mat;
        ConstGenericSparseMatrixB b_mat;
        SparseMatrixA p;
    public:
        template <typename Derived>
        MySparseGenMatProd(
            const Eigen::SparseMatrixBase<Derived>& mat,
            const Eigen::SparseMatrixBase<Derived>& amat, 
            const Eigen::SparseMatrixBase<Derived>& bmat) 
            : m_mat(mat), a_mat(amat), b_mat(bmat)
        {
            // auto I = Eigen::SparseMatrix<double>(amat.rows(), amat.rows());
            // auto J = Eigen::SparseMatrix<double>(bmat.rows(), bmat.rows());
            // I.setIdentity();
            // J.setIdentity();
            p = Eigen::kroneckerProduct(amat, bmat);
            // auto bb = Eigen::kroneckerProduct(I, bmat);
            a_solver.compute(p);
            //b_solver.compute(bmat);
        }

        Index rows() const { return m_mat.rows(); }
        Index cols() const { return m_mat.cols(); }

        // y_out = B^-1 A * x_in
        // B y_out = A * x_in
        // (A*I) (J*B) y_out = H * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in, m_mat.cols());
            MapVec y(y_out, m_mat.rows());
            y.noalias() = m_mat * x;
            y.noalias() = a_solver.solve(y);

            //MapMat map(y_out, a_mat.cols(), b_mat.rows()); 
            //map = a_solver.transpose().solve(map.transpose());
            //map = b_solver.solve(map);
        }

        Matrix operator*(const Eigen::Ref<const Matrix>& mat_in) const
        {
            Matrix y = m_mat * mat_in;
            y = a_solver.solve(y);
            
            return y;
        }

        Scalar operator()(Index i, Index j) const
        {
            return m_mat.coef(i, j);
        }
    };


    template <typename Scalar_, typename TypeA = Eigen::Sparse, typename TypeB = Eigen::Sparse, typename TypeC = Eigen::Sparse,
            int FlagsA = Eigen::ColMajor, int FlagsB = Eigen::ColMajor, int FlagsC = Eigen::ColMajor,
            typename StorageIndexA = int, typename StorageIndexB = int, typename StorageIndexC = int>
    class ExpSparseGenMatProd
    {
    public:
        using Scalar = Scalar_;

    private:
        using Index = Eigen::Index;
        using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        using MapConstVec = Eigen::Map<const Vector>;
        using MapVec = Eigen::Map<Vector>;
        using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using MapMat = Eigen::Map<Matrix>;
        using SparseMatrixA = Eigen::SparseMatrix<Scalar, FlagsA, StorageIndexA>;
        using SparseMatrixB = Eigen::SparseMatrix<Scalar, FlagsB, StorageIndexB>;
        using SparseMatrixC = Eigen::SparseMatrix<Scalar, FlagsC, StorageIndexC>;

        using ConstGenericSparseMatrixA = const Eigen::Ref<const SparseMatrixA>;
        using ConstGenericSparseMatrixB = const Eigen::Ref<const SparseMatrixB>;
        using ConstGenericSparseMatrixC = const Eigen::Ref<const SparseMatrixC>;


        using FacTypeC = Eigen::SparseLU<Eigen::SparseMatrix<double>>;

        FacTypeC solver;
        ConstGenericSparseMatrixA a_mat;
        ConstGenericSparseMatrixB b_mat;
        ConstGenericSparseMatrixC m_mat;

        SparseMatrixC s;
        Eigen::SparseMatrix<double> f;

        double delta;
    public:
        template <typename Derived>
        ExpSparseGenMatProd(
            const Eigen::SparseMatrixBase<Derived>& mat,
            const Eigen::SparseMatrixBase<Derived>& amat, 
            const Eigen::SparseMatrixBase<Derived>& bmat, double dt) 
            : m_mat(mat), a_mat(amat), b_mat(bmat)
        {
            // auto I = Eigen::SparseMatrix<double>(amat.rows(), amat.rows());
            // auto J = Eigen::SparseMatrix<double>(bmat.rows(), bmat.rows());
            // I.setIdentity();
            // J.setIdentity();
            delta = dt;

            s = Eigen::kroneckerProduct(a_mat, b_mat);
            Eigen::SparseMatrix<double> lu = s + delta*0.5d*m_mat;
            f = (s - delta*0.5d*m_mat);
            // auto bb = Eigen::kroneckerProduct(I, bmat);
            solver.compute(lu);
            //b_solver.compute(bmat);
        }

        Index rows() const { return m_mat.rows(); }
        Index cols() const { return m_mat.cols(); }

        // y_out = B^-1 A * x_in
        // B y_out = A * x_in
        // (A*I) (J*B) y_out = H * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in, m_mat.cols());
            MapVec y(y_out, m_mat.rows());
            y.noalias() =  f * x;
            y.noalias() = solver.solve(y);

            //MapMat map(y_out, a_mat.cols(), b_mat.rows()); 
            //map = a_solver.transpose().solve(map.transpose());
            //map = b_solver.solve(map);
        }

        Matrix operator*(const Eigen::Ref<const Matrix>& mat_in) const
        {
            Matrix y = f * mat_in;
            y = solver.solve(y);
            
            return y;
        }

        Scalar operator()(Index i, Index j) const
        {
            return m_mat.coef(i, j);
        }
    };

}  // namespace Spectra

#endif  // SPECTRA_DENSE_GEN_MAT_PROD_H