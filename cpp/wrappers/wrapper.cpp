#define EIGEN_USE_MKL_ALL

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <array>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "hamiltonian2d.h"
#include "hamiltonian3d.h"
#include "bornoppenheimer.h"
#include "s32.h"
#include <omp.h> 


namespace py = pybind11;


PYBIND11_MODULE(hamiltonian, handle) {
    handle.doc() = "This is my awesome module!";

    py::class_<Hamiltonian2D>(
        handle, "PyHamiltonian2D"
    )
    .def(py::init<py::EigenDRef<Eigen::VectorXd>, int, int, std::array<double, 2>, double>())
    .def("get_eigenvalues", &Hamiltonian2D::getEigenvalues)
    .def("get_eigenvectors", &Hamiltonian2D::getEigenvectors)
    .def("get_state", &Hamiltonian2D::getEigenfunction)
    .def("get_spectrum", &Hamiltonian2D::getTheSpectrum)
    .def_property_readonly("pMatr", &Hamiltonian2D::getPMatr)
    .def_property_readonly("nMatr", &Hamiltonian2D::getNMatr)
    .def_property_readonly("d2Matr",&Hamiltonian2D::getD2Matr)
    .def_property_readonly("potential",&Hamiltonian2D::getPotential)
    ;

    py::class_<CHermiteBC>(
        handle, "PyHermiteBC"
    )
    .def(py::init<py::EigenDRef<Eigen::VectorXd>, int, int>())
    .def("fBSplineBC", &CHermiteBC::fBSplineBC)
    .def("d1BSplineBC", &CHermiteBC::d1BSplineBC)
    .def("d2BSplineBC", &CHermiteBC::d2BSplineBC)
    .def("getSpaceDim", &CHermiteBC::getSpaceDim)
    .def("getSNMatr", &CHermiteBC::generateSNMatr)
    .def("getSPMatr", &CHermiteBC::generateSPMatr)
    .def("locate", &CHermiteBC::locate)
    .def_property_readonly("collocPoints", &CHermiteBC::getCollocPoints)
    .def_property_readonly("leftPoints", &CHermiteBC::getLeftPoints)
    .def_property_readonly("midPoints", &CHermiteBC::getMidPoints)
    .def_property_readonly("rightPoints", &CHermiteBC::getRightPoints)
    ;

    py::class_<Hamiltonian3D>(
         handle, "PyHamiltonian3D"
    )
    .def(py::init<py::EigenDRef<Eigen::VectorXd>, py::EigenDRef<Eigen::VectorXd>,
                                std::array<int, 4>,
                                std::array<double, 3>,
                                double,
                                int>())
    .def("get_eigenvalues", &Hamiltonian3D::getEigenvalues)
    .def("get_eigenvectors", &Hamiltonian3D::getEigenvectors)
    .def("get_spectrum", &Hamiltonian3D::getTheSpectrum)
    .def_property_readonly("pMatr", &Hamiltonian3D::getPMatr)
    .def_property_readonly("nMatr", &Hamiltonian3D::getNMatr)
    .def_property_readonly("h", &Hamiltonian3D::getHamiltonian)
    .def("init_LU", &Hamiltonian3D::initHamiltonianLU)
    .def("init_impulse", &Hamiltonian3D::initImpulse)
    .def("init_absorption", &Hamiltonian3D::initAbsorption)
    .def("init_scaling", &Hamiltonian3D::initScaling)
    .def("get_state", &Hamiltonian3D::getEigenfunction)
    .def("evolutionStep", &Hamiltonian3D::evolutionStep)
    ;

    py::class_<BornOppenheimer2D>(
        handle, "PyBornOppenHeimer2D"
    )
    .def(py::init<py::EigenDRef<Eigen::VectorXd>, int, int, std::array<double, 3>, double, double>())
    .def("get_eigenvalues", &BornOppenheimer2D::getEigenvalues)
    .def("get_eigenvectors", &BornOppenheimer2D::getEigenvectors)
    .def("get_state", &BornOppenheimer2D::getEigenfunction)
    .def("get_spectrum", &BornOppenheimer2D::getTheSpectrum)
    .def_property_readonly("pMatr", &BornOppenheimer2D::getPMatr)
    .def_property_readonly("nMatr", &BornOppenheimer2D::getNMatr)
    ;

    handle.def("get_max_threads", &omp_get_max_threads, "Returns max number of threads");
    handle.def("set_num_threads", &omp_set_num_threads, "Set number of threads");
}
