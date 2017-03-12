#include <assert.h>
#include <limits>
#include <vector>
#include <iostream> // for debugging
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "find_root.h"
#include "distributeND.h"

using namespace std;
using namespace Go;

namespace {

// ----------------------------------------------------------------------------
// approxmate total length of 1D curve
double total_length_1D(const pair<double, double>& bounds,
		       const DistFun& dfun,
		       const unsigned int NUM_SAMPLES = 20)
// ----------------------------------------------------------------------------    
{
  double res = 0;
  double p1 = bounds.first;
  const double interval = bounds.second - bounds.first;
  
  for (unsigned int i = 0; i != NUM_SAMPLES; ++i) {
    double p2 = bounds.first + ((i+1) * interval) / double(NUM_SAMPLES);
    res += dfun(&p1, &p2, nullptr, nullptr);
    p1 = p2;
  }
  return res;
};

// ----------------------------------------------------------------------------
  double dist_energy_1D(const double dist, const double radius, double* const grad,
			double* const jac, double eps=0)
// ----------------------------------------------------------------------------  
{
  const double dist_eps = dist + eps;
  const double dist_eps_3 = dist_eps * dist_eps * dist_eps;
  const double dist_eps_4 = dist_eps * dist_eps_3;
  const double res_sqrt = (dist > radius) ? 0 : (radius - dist) / dist_eps;

  if (grad) 
    *grad = (dist > radius) ? 0 : (-2 * (radius+eps) * (radius - dist) / dist_eps_3);

  if (jac)
    *jac = (dist > radius) ? 0 :  2 * (radius+eps) / dist_eps_3 +
                                  6 * (radius - dist) * (radius + eps) / dist_eps_4;

  return res_sqrt * res_sqrt;
}

  
// ----------------------------------------------------------------------------
class Functor1D
// ----------------------------------------------------------------------------
{
public:
  Functor1D(const unsigned int dim,
	    const pair<double, double>& bounds,
	    const DistFun& dfun,
	    const double radius,
	    const double buffer) :
    dim_(dim), bounds_(bounds), dfun_(dfun), radius_(radius),
    min_par_(dim, bounds.first + buffer),
    max_par_(dim, bounds.second - buffer), small_fac_(0.35) {}
  // the smaller small_fac_ is set, the more strongly points repel each other as they get very
  // close (but then the derivatives get very steep too)
  
  double operator()(const double* arg) const;
  void grad(const double* arg, double* grad) const;
  void jac_full(const double* arg, double* jac) const; // full (not sparse) jacobian
  
  double minPar(int n) const {return min_par_[n];}
  double maxPar(int n) const {return max_par_[n];}
  unsigned int dim() const {return dim_;}
private:
  const unsigned int dim_;
  const pair<double, double> bounds_;
  const DistFun& dfun_;
  const double radius_;
  const vector<double> min_par_;
  const vector<double> max_par_;
  const double small_fac_;
};

// ----------------------------------------------------------------------------
void minimise_newton(Functor1D, vector<double>& points);
// ----------------------------------------------------------------------------

  
}; // end anonymous namespace

namespace Go {

  
// ============================================================================
vector<double> distribute_points_1D(const double* const pts,
				    const unsigned int num_points,
				    const pair<double, double>& bounds,
				    const DistFun& dfun)
// ============================================================================  
{
  assert(bounds.first < bounds.second);

  // sort points, and cull away those outside specified bounds
  vector<double> cur_points;
  cur_points.reserve(num_points);
  copy_if(pts, pts + num_points, back_inserter(cur_points),
	  [bounds](double x) {return (x > bounds.first) && (x < bounds.second);});
  sort(cur_points.begin(), cur_points.end());
  const unsigned int N = (unsigned int) cur_points.size(); // number of internal points left
  assert(N>0);
  
  // computing estimated average distance between each point when equidistantly placed
  const double average_distance = total_length_1D(bounds, dfun) / (N + 1);

  cout << "Distance is: " << average_distance << endl;
  
  // constructing and minimizing objective
  const double BUFFER_FAC = 4.0;
  const double buffer = min(min((bounds.second - bounds.first) / (BUFFER_FAC * N),
				cur_points[0] - bounds.first),
  			    bounds.second - cur_points.back())/2;

  const double DIST_FAC = 2;

  Functor1D obj_fun(N, bounds, dfun, average_distance*DIST_FAC, buffer*0);
  minimise_newton(obj_fun, cur_points);
    
  // FunctionMinimizer<Functor1D> fmini(N, obj_fun, &cur_points[0]);
  // minimise_conjugated_gradient(fmini);
  // copy(fmini.getPar(), fmini.getPar() + N, cur_points.begin()); // keeping the new parameters

  return cur_points;
};

// ============================================================================
vector<double> distribute_points_2D(const double* const pts, // u1, v1, u2, v2, ...
				    const unsigned int num_points,
				    const vector<vector<double>>& bounds,
				    const DistFun& dfun)
// ============================================================================  
{
  // sort points, cull away those outside bounds

  // compute estimated average distance

  // construct objective function

  // minimize using Newton

  // copy and return result
  return vector<double>();
};
  
}; // end namespace Go

namespace {

// ----------------------------------------------------------------------------  
double Functor1D::operator()(const double* arg) const
// ----------------------------------------------------------------------------
{
  double result = 0;
  // @@ brute force implementation
  for (int i = -1; i != (int)dim_+1 ; ++i) {
    const double p1 = (i == -1) ? bounds_.first :
                      (i == (int)dim_) ? bounds_.second: arg[i];
    for (int j = i+1; j != (int)dim_+1; ++j) {
      const double p2 = (j==(int)dim_) ? bounds_.second : arg[j];

      if (i != -1 || j != (int)dim_) {
	const double d = dfun_(&p1, &p2, nullptr, nullptr);
	result += dist_energy_1D(d, radius_, nullptr, nullptr, small_fac_ * radius_);
      }
    }
  }
  
  return result;
}
  
// ----------------------------------------------------------------------------
void Functor1D::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  // @@ brute force implementation
  double pder[2];
  double dder;
  fill(grad, grad + dim_, 0);
  for (int i = -1; i != (int)dim_+1; ++i) {
    const double p1 = (i == -1) ? bounds_.first :
                      (i == (int)dim_) ? bounds_.second: arg[i];
    for (int j = i+1; j != (int)dim_+1; ++j) {
      const double p2 = (j==(int)dim_) ? bounds_.second : arg[j];

      if (i != -1 || j != (int) dim_) {
	dist_energy_1D(dfun_(&p1, &p2, pder, nullptr),
		       radius_, &dder, nullptr, small_fac_ * radius_);
	if ((i >= 0) && (i < (int)dim_))
	  grad[i] += dder * pder[0];
	if (j < (int)dim_)
	  grad[j] += dder * pder[1];
      }
    }
  }
};
  
// ----------------------------------------------------------------------------
void Functor1D::jac_full(const double* arg, double* jac) const
// ----------------------------------------------------------------------------
{
  double pder[2];
  double pder2[3];
  double dder, dder2;
  fill(jac, jac + dim_ * dim_, 0);
  for (int i = - 1; i != (int)dim_+1; ++i) {
    const double p1 = (i==-1)        ? bounds_.first :
                      (i==(int)dim_) ? bounds_.second : arg[i];
    for (int j = i+1; j != (int)dim_+1; ++j) {
      const double p2 = (j==(int)dim_) ? bounds_.second : arg[j];

      if (i != -1 || j != (int) dim_) {
	const double e = dist_energy_1D(dfun_(&p1, &p2, pder, pder2),
					radius_, &dder, &dder2, small_fac_ * radius_);
	if (e > 0) {
	  // nonzero contribution from this point pair
	  const bool movable_i = ( (i>=0) && (i < (int)dim_));
	  const bool movable_j = (j < (int)dim_);
	    
	  if (movable_i)  // d2u
	    jac[i*dim_ + i] += dder2 * pder[0] * pder[0] + dder * pder2[0];
	  
	  if (movable_j) //d2v
	    jac[j*dim_ + j] += dder2 * pder[1] * pder[1] + dder * pder2[1];
	  
	  if (movable_i && movable_j) {
	    const double cross_entry = dder2 * pder[0] * pder[1] + dder * pder2[2];
	    jac[i*dim_ + j] += cross_entry;
	    jac[j*dim_ + i] += cross_entry;
	  }
	}
      }
    }
  }
}
  
// ----------------------------------------------------------------------------
void minimise_newton(Functor1D functor, vector<double>& points)
// ----------------------------------------------------------------------------
{
  auto rn_to_rn_fun = [&functor] (const double* const pts, unsigned int dim) {
    ValAndJac result {vector<double>(dim, 0), vector<double>(dim*dim, 0)};
    functor.grad(pts, &result.value[0]);
    functor.jac_full(pts, &result.jacobian[0]);
    return result;
  };
  
  auto update_fun = [&functor] (double* val, const double* const update) {
    double max_step = 1;
    for (unsigned int i = 0; i != functor.dim(); ++i) {
      if ((update[i] > 0)  && val[i] + update[i] > functor.maxPar(i))
	max_step = min(max_step, (functor.maxPar(i) - val[i]) / update[i]);
      else if ((update[i] < 0) && val[i] + update[i] < functor.minPar(i))
	max_step = min(max_step, (val[i] - functor.minPar(i)) / (-update[i]));
    }
    
    if (max_step < 1)
      max_step = max_step * 0.9; // avoid sending parameters all away to the boundary

    cout << "max_step " << max_step << endl;
    
    for (unsigned int i = 0; i != functor.dim(); ++i)
      val[i] += max_step * update[i];
  };

  auto update_fun_2 = [&functor] (double* val, const double* const update) {

    vector<int> outside_ixs;
    for (unsigned int i = 0; i != functor.dim(); ++i) {
      val[i] += update[i];
      if ((val[i] < functor.minPar(i)) || (val[i] > functor.maxPar(i))) {
	outside_ixs.push_back(i);
      }
    }
    // move outside points back in again, where there is room
    vector<double> tmp(val, val+functor.dim());

    for (int i = outside_ixs.size()-1; i >= 0; --i) 
      tmp.erase(tmp.begin() + outside_ixs[i]);
    sort(tmp.begin(), tmp.end());

    for (int i = 0; i != outside_ixs.size(); ++i) {
      const int ix = outside_ixs[i];
      vector<double> diff;
      for (int i = 1; i != tmp.size(); ++i) {
	diff.push_back(tmp[i] - tmp[i-1]);
      }
      int it = (max_element(diff.begin(), diff.end()) - diff.begin());
      double new_val = (tmp[it] + tmp[it+1])/2;
      tmp.push_back(new_val);
    }
    assert(tmp.size() == functor.dim());
    
  };
  
  
  const double TOL = 1e-12;
  const unsigned int MAX_ITER = 80;
  points = find_root(rn_to_rn_fun,
		     &points[0],
		     functor.dim(),
		     update_fun, TOL, MAX_ITER);
}
  
}; // end anonymous namespace 


