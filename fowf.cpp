#include "fowf.hpp"

int add(int a, int b) {
  return a + b;
}

// ================ Grid ===================
Grid::Grid(double _h, int _num) : h_(_h), xs_(_num) {
  for(int i = 0; i < _num; i++)
    xs_(i) = (i + 1) * h_;
}
double Grid::h() const { return h_; }
int Grid::num() const { return xs_.rows(); }
double Grid::xs_i(int i) const { return xs_(i); }
const VD& Grid::xs() const { return xs_; }

// ================ HydrogenPI ============
HydrogenPI::HydrogenPI(double _ene): ene_(_ene), l0_(0), l1_(1) {}
double HydrogenPI::ene() const { return ene_; }
int HydrogenPI::l0() const { return l0_; }
int HydrogenPI::l1() const { return l1_; }
double HydrogenPI::Vx(double x) const {
  int l = l1_;
  return double(l * (l + 1)) / (x * x * 2.0) - 1.0 / x;
}

// ================ DiscreteFunc ==========
DiscreteFunc::DiscreteFunc(const Grid& g) : g_(g), ys_(g.num()) {}
DiscreteFunc::DiscreteFunc(const vCD& _ys, const Grid& g) :  g_(g) {
  ys_.resize(_ys.size());
  for(int i = 0; i < _ys.size(); i++) {
    ys_(i) = _ys[i];
  }
}
DiscreteFunc::DiscreteFunc(const VC& _ys, const Grid& g) : ys_(_ys), g_(g) {}
int DiscreteFunc::size() const { return g_.num(); }
double DiscreteFunc::xs_i(int i ) const { return g_.xs_i(i); }
CD DiscreteFunc::ys_i(int i) const { return ys_(i); }
const Grid& DiscreteFunc::grid() const { return g_; }
const VC& DiscreteFunc::ys() const { return ys_; }


// ================ Sparse matrix Utils ========
void AddLaplacian(const Grid& g, const HydrogenPI& h_pi, MC& mat) {
  
  int num = g.num();
  MC tmp(num, num);
  std::vector<T> v;
  
  v.push_back(T(0, 0, -2.0));
  v.push_back(T(0, 1, +1.0));

  for(int i = 1; i < num-1; i++) {
    v.push_back(T(i, i-1, +1.0 ));
    v.push_back(T(i, i  , -2.0 ));
    v.push_back(T(i, i+1, +1.0  ));
  }

  double h = g.h();
  double k = sqrt(2.0 * h_pi.ene());
  v.push_back(T(num-1,num-2,  2.0));
  v.push_back(T(num-1,num-1, CD(-2.0, 2.0*k*h)));

  tmp.setFromTriplets(v.begin(), v.end());
  
  mat += tmp;

}
void AddPotEne(const Grid& g, const HydrogenPI& h_pi, MC& mat) {

  int num = g.num();
  MC tmp(num, num);
  std::vector<T> v;

  double h = g.h();
  
  for(int i = 0; i < num; i++) {
    double x = g.xs_i(i);
    v.push_back(T(i, i,  -2 * h * h * (h_pi.Vx(x) - h_pi.ene())));
  }
  
  tmp.setFromTriplets(v.begin(), v.end());
  
  mat += tmp;


}
void AddDriv1skpLength(const Grid& g, const HydrogenPI& h_pi, VC& vec) {

  double h = g.h();
  for(int i = 0; i < g.num(); i++) {
    double x = g.xs_i(i);
    double y = -2.0 * h * h * (2.0 * x * x * exp(-x));
    vec(i) += y;
  }

}


// ================ DrivSys ================
DrivSys::DrivSys(const Grid& g, const HydrogenPI& h) : d_mat(g.num(), g.num()), m_vec(g.num()), g_(g) {

  //  d_mat.resize(g.num(), g.num());
  m_vec = VC::Zero(g.num());

  AddLaplacian(g, h, d_mat);
  AddPotEne(g, h, d_mat);
  AddDriv1skpLength(g, h, m_vec);

}
DiscreteFunc DrivSys::Solve() const {

  Eigen::SparseLU<MC> solver;
  solver.analyzePattern(d_mat);
  solver.factorize(d_mat);
  VC ys = solver.solve(m_vec);

  DiscreteFunc res(ys, g_);
  return res;
}
