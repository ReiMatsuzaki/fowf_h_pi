#include <boost/python.hpp>
#include "fowf.hpp"

// wrapper function for DiscreteFunc
DiscreteFunc discretize_func(boost::python::list& xs, const Grid& g) {
  int num = boost::python::len(xs);
  std::vector<CD> tmp(num);
  for(int i = 0; i < num; i++) 
    tmp[i] = boost::python::extract<CD>(xs[i]);
  DiscreteFunc d(tmp, g);
  return d;
}

BOOST_PYTHON_MODULE(fowf) {
  boost::python::def( "add", add);
  boost::python::class_<Grid>("Grid", boost::python::init<double, int>())
    .def("h", &Grid::h)
    .def("num", &Grid::num)
    .def("xs_i", &Grid::xs_i);
  boost::python::class_<HydrogenPI>("HydrogenPI", boost::python::init<double>())
    .def("ene", &HydrogenPI::ene)
    .def("l0", &HydrogenPI::l0)
    .def("l1", &HydrogenPI::l1);
  
  boost::python::class_<DiscreteFunc>("DiscreteFunc", boost::python::init<Grid>())
    .def("xs_i", &DiscreteFunc::xs_i)
    .def("ys_i", &DiscreteFunc::ys_i).
    def("size", &DiscreteFunc::size);
  boost::python::def( "discretize_func", discretize_func);

  boost::python::class_<DrivSys>("DrivSys", boost::python::init<Grid, HydrogenPI>())
    .def("solve", &DrivSys::Solve);
}
