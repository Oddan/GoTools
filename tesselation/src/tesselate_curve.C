#include "tesselate_curve.h"
#include "make_curve_spacing_fun.h"
#include "find_root.h"

#include <fstream>
#include "GoTools/geometry/PointCloud.h"

using namespace std;
using namespace Go;

namespace {

  vector<double> define_parvec(double minparam, double maxparam, unsigned int num_intparams) {
    vector<double> result(num_intparams+2);
    result.front() = minparam;
    result.back()  = maxparam;
    for (size_t i = 1; i != num_intparams+1; ++i)
      result[i] = (maxparam-minparam)/(num_intparams + 1) * double(i);
    return result;
  }

  double norm2(const double* const x, unsigned int num) {
    double res = 0;
    for (unsigned int i = 0; i != num; ++i) {
      res = res + x[i] * x[i];
    }
    return res;
  }
    
  
  // @@ Debug purposes
  void store_points_and_curve(const ParamCurve& pc,
			      double* t,
			      unsigned int num_points,
			      const string& filename)
  {

    ofstream os(filename);
    pc.writeStandardHeader(os);  pc.write(os);
    vector<double> coords;
    for (double* pt = t; pt != t+num_points; ++pt) {
      Point cur_point;
      pc.point(cur_point, *pt);
      coords.push_back(cur_point[0]);
      coords.push_back(cur_point[1]);
      coords.push_back(cur_point[2]);
    }
    PointCloud<3> pcloud(coords.begin(), (int)coords.size()/3);
    pcloud.writeStandardHeader(os);
    pcloud.write(os);
    os.close();
  }

  
}; // end anonymous namespace

namespace Go {
  
vector<double> tesselate_curve(const ParamCurve& pc, unsigned int num_internal_points)
{
  const double TOL = 1e-9;
  const int ITER = 200;
  const vector<double> p = define_parvec(pc.startparam(), pc.endparam(), num_internal_points);
  const auto fun = make_curve_spacing_fun(pc);

  // ------------------------- specify updater function -------------------------

  UpdateFun updater = [&] (double* x, const double* const dx) {

    const double MAX_STEP = (pc.endparam() - pc.startparam())/20;
    const double dx_max = max( *max_element(dx, dx + num_internal_points),
  			      -*min_element(dx, dx + num_internal_points));
    //cout << dx_max << endl;
    const double fac = (dx_max < MAX_STEP) ? 1 : MAX_STEP/dx_max;
    
    // trying to update function using jacobian.  If this fails, use the fact
    // that the linear system is really a gradient, and employ steepest descent.
    const double orig_funval = get<0>(fun)(x, (int)num_internal_points);
    vector<double> x_tmp(x, x + num_internal_points);

    for (size_t i = 0; i != (size_t)num_internal_points; ++i) {
      x_tmp[i] += fac * dx[i];
      x_tmp[i] = max(x_tmp[i], pc.startparam());
      x_tmp[i] = min(x_tmp[i], pc.endparam());
    }
    for (size_t i = 0; i != (size_t)num_internal_points-1; ++i) {
      x_tmp[i] = min(x_tmp[i], x_tmp[i+1]);
    }
    double new_funval = get<0>(fun)(&x_tmp[0], (int)num_internal_points);
    if (new_funval < orig_funval) {
      copy(x_tmp.begin(), x_tmp.end(), x);
      return;
    }

    // if we got here, the use of the jacobian did not decrease our function
    // value.  Let us try a steepest descent instead.
    cout << "Resort to SD" << endl;
    const auto res = get<1>(fun)(x, num_internal_points);

    // minimize along direction given by -
    directional_minimum(fun, x, &res.value[0], num_internal_points);
    
    // cout << "orig funval: " << orig_funval << endl;
    // cout << "jac funval: " << new_funval << endl;
    // const auto res = get<1>(fun)(x, num_internal_points);
    // double scaling = sqrt(norm2(dx, num_internal_points) /
    // 			  norm2(&(res.value[0]), num_internal_points));
    // const int MAX_ITER = 100;
    
    // for (int i = 0; i != MAX_ITER; ++i) {
    //   cout << i << endl;
    //   copy(x, x+num_internal_points, x_tmp.begin());

    //   for (size_t i = 0; i != (size_t)num_internal_points; ++i) {
    // 	x_tmp[i] -= scaling * res.value[i];
    // 	x_tmp[i] = max(x_tmp[i], pc.startparam());
    // 	x_tmp[i] = min(x_tmp[i], pc.endparam());
    //   }
    //   for (size_t i = 0; i != (size_t)num_internal_points-1; ++i) {
    // 	x_tmp[i] = min(x_tmp[i], x_tmp[i+1]);
    //   }
    //   new_funval = get<0>(fun)(&x_tmp[0], (int)num_internal_points);
    //   cout << "updated funval: " << new_funval << endl;
    //   if (new_funval < orig_funval) {

    // 	const auto new_grad =get<1>(fun)(&x_tmp[0], num_internal_points);
    // 	double sprod = 0;
    // 	for (int i = 0; i != num_internal_points; ++i) {
    // 	  sprod = new_grad.value[i] * res.value[i];
    // 	}
    // 	if (sprod < 0) {
    // 	  // we can probably do better by further reducing the stepsize
    // 	  copy(x, x+num_internal_points, x_tmp.begin());
    // 	  for (size_t i = 0; i != (size_t)num_internal_points; ++i) {
    // 	    x_tmp[i] -= scaling/2 * res.value[i];
    // 	    x_tmp[i] = max(x_tmp[i], pc.startparam());
    // 	    x_tmp[i] = min(x_tmp[i], pc.endparam());
    // 	  }
    // 	  for (size_t i = 0; i != (size_t)num_internal_points-1; ++i) {
    // 	    x_tmp[i] = min(x_tmp[i], x_tmp[i+1]);
    // 	  }
    // 	}
    // 	break;
    //   }
    //   scaling /= 1.5;
    //   //cout << scaling << endl;
    }
    // copy whatever we found to the original vector.
  //copy(x_tmp.begin(), x_tmp.end(), x);
  };

  // UpdateFun updater = [&] (double* x, const double* const dx) {

  //   double min_spacing = (pc.endparam() - pc.startparam());
  //   for (int i = 0; i != num_internal_points-1; ++i) 
  //     min_spacing = min(min_spacing, x[i+1] - x[i]);

  //   const double dx_max = max( *max_element(dx, dx + num_internal_points),
  // 			      -*min_element(dx, dx + num_internal_points));

  //   double fac = (dx_max > min_spacing) ? min_spacing / dx_max : 1;
  //   //fac = fac * 0.8;
    
  //   for (size_t i = 0; i != (size_t)num_internal_points; ++i) {
  //     x[i] += fac * dx[i];
  //     x[i] = max(x[i], pc.startparam());
  //     x[i] = min(x[i], pc.endparam());
  //   }
  //   store_points_and_curve(pc, x, num_internal_points, "krull.g2");
  // };

  // UpdateFun updater = [&] (double* x, const double* const dx) {

  //   const double MAX_STEP = (pc.endparam() - pc.startparam())/20;
  //   const double dx_max = max( *max_element(dx, dx + num_internal_points),
  // 			      -*min_element(dx, dx + num_internal_points));
  //   //cout << dx_max << endl;
  //   const double fac = (dx_max < MAX_STEP) ? 1 : MAX_STEP/dx_max;
    
  //   for (size_t i = 0; i != (size_t)num_internal_points; ++i) {
  //     x[i] += fac * dx[i];
  //     x[i] = max(x[i], pc.startparam());
  //     x[i] = min(x[i], pc.endparam());
  //   }
  //   for (size_t i = 0; i != (size_t)num_internal_points-1; ++i) {
  //     x[i] = min(x[i], x[i+1]);
  //   }
  // };

  // find optimal set of points
  return find_root(get<1>(fun), &p[1], num_internal_points, updater, TOL, ITER); 
}
  
};
