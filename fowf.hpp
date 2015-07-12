#ifndef FOWF_HPP
#define FOWF_HPP

#include <vector>
#include <Eigen/Sparse>


typedef std::complex<double> CD;
typedef std::vector<CD> vCD;
typedef Eigen::SparseMatrix<CD> MC;
typedef Eigen::VectorXcd VC;
typedef Eigen::VectorXd VD;
typedef Eigen::Triplet<CD> T;

int add(int, int);

class Grid {
private:
  double h_;
  VD xs_;
public:
  Grid(double, int);
  double h() const;
  int num() const;
  double xs_i(int i) const;
  const VD& xs() const;
};

class HydrogenPI {
private:
  double ene_;
  int l0_;
  int l1_;
public:
  HydrogenPI(double);
  double ene() const;
  int l0() const;
  int l1() const;
  double Vx(double x) const;
};

// add matrix
// (T+V-E)\psi = \mu\phi
// (Laplacian -2V + 2E) = -2\mu\phi
void AddLaplacian(const Grid&, MC&);
void AddPotential(const Grid&, const HydrogenPI&, MC&);
void AddDriv1skpLength(const Grid&, const HydrogenPI&, VC&);

class DiscreteFunc {
private:
  Grid g_;
  VC ys_;
public:
  DiscreteFunc(const Grid&);
  DiscreteFunc(const std::vector<CD>& _xs, const Grid&);
  DiscreteFunc(const VC& _xs, const Grid&);
  double xs_i(int i) const;
  CD ys_i(int i) const;
  int size() const;
  const Grid& grid() const;
  const VC& ys() const;
};

class DrivSys {
private:
  MC d_mat;
  VC m_vec;
  Grid g_;
public:
  DrivSys(const Grid&, const HydrogenPI&);
  DrivSys(const MC&, const VC&, const Grid&);
  DiscreteFunc Solve() const;
};


#endif
