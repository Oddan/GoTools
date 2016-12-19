#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"

#include "GoTools/lrsplines2D/SSurfTraceIsocontours.h"

using namespace std;
using namespace Go;

namespace {
  CurveVec compute_levelset(const SplineSurface& ss,
			    SISLSurf* ss_sisl, // SISL surface corresponding to ss
			    const double isoval,
			    const bool include_3D_curves,
			    const bool use_sisl_marching);
}; // end anonymous namespace 

namespace Go {


// ============================================================================
vector<CurveVec> SSurfTraceIsocontours(const SplineSurface& ss,
				       const vector<double>& isovals,
				       bool include_3D_curves,
				       bool use_sisl_marching)
// ============================================================================
{
  // assert(!use_sisl_marching); // @@ Unimplemented for now
  
  // Compute topology for each requested level-set.  We use SISL for this
  SISLSurf* sislsurf1D = GoSurf2SISL(ss, false);

  // Defining function tracing out the level set for a specified isovalue
  const function<CurveVec(double)> comp_lset = [&] (double ival)
    {return compute_levelset(ss, sislsurf1D, ival, include_3D_curves, use_sisl_marching);};

  // Computing all level-set curves for all isovalues ("transforming" each
  // isovalue into its corresponding level-set)
  const vector<CurveVec> result = apply_transform(isovals, comp_lset);

  // Cleaning up after use of SISL objects
  freeSurf(sislsurf1D);

  // Returning result
  return result;
}

} // end namespace Go;


namespace {

// ----------------------------------------------------------------------------  
pair<SISLIntcurve**, int> compute_topology(SISLSurf* ss_sisl, double isoval)
// ----------------------------------------------------------------------------
{
  // This function does the equivalent of SISL sh1851, but assumes that the
  //provided surface is 1D (bivariate spline function), noncyclic and k-regular.

  
  // Finding intersections, using SISL function sh1761.
  int jpt = 0; // number of single intersection point
  int jcrv = 0; // number of intersection curves
  int jsurf = 0; 
  double* gpar = SISL_NULL; // parameter values of the single intersection points
  double* spar = SISL_NULL; // dummy array
  int *pretop=SISL_NULL;
  SISLIntcurve** wcurve; // array containing descriptions of the intersection curves
  SISLIntsurf** wsurf; 

  auto qp  = newPoint(&isoval, 1, 1);
  auto qo1 = newObject(SISLSURFACE);
  qo1->s1 = ss_sisl;
  qo1->o1 = qo1;
  auto qo2 = newObject(SISLPOINT);
  qo2->p1 = qp;
  
  SISLIntdat* qintdat = SISL_NULL; // intersection result
  auto freeall = [&] (bool skip_wcurves = false) {
    if (gpar)                         free(gpar);
    if (spar)                         free(spar);
    if (pretop)                       free(pretop);
    if (qintdat)                      freeIntdat(qintdat);
    if ((bool)wcurve & (jcrv > 0) & !skip_wcurves) freeIntcrvlist(wcurve, jcrv);
    for (int i = 0; i < jsurf; ++i)   freeIntsurf(wsurf[i]);
    if ((bool)wsurf & (jsurf > 0))          free(wsurf);
  };
  
  auto cleanup_and_throw = [&freeall] (string str) {
    freeall();
    throw runtime_error(str.c_str());
  };
  
  const double epsge = 1e-6;
  int kstat = 0;
  
  // find intersections
  sh1761(qo1, qo2, epsge, &qintdat, &kstat);
  if (kstat < 0) cleanup_and_throw("SISL error in sh1761.");

  const int kdeg = 1;
    
  // express intersections on output format (SISL)
  if (qintdat)
    hp_s1880 (qo1, qo1, kdeg, 2, 0, qintdat, &jpt, &gpar, &spar, &pretop,
	      &jcrv, &wcurve, &jsurf, &wsurf, &kstat);
  
  freeall(true);
  
  return pair<SISLIntcurve**, int>(wcurve, jcrv);
}

// ----------------------------------------------------------------------------  
pair<CurvePtr, CurvePtr> trace_isoval_curve(const SplineSurface& ss,
					      SISLIntcurve* ic, double ival)
// ----------------------------------------------------------------------------
{
  return pair<CurvePtr, CurvePtr>(); // @@ dummy
}
  
// ----------------------------------------------------------------------------
CurveVec compute_levelset(const SplineSurface& ss,
			  SISLSurf* ss_sisl, // SISL surface corresponding to ss
			  const double isoval,
			  const bool include_3D_curves,
			  const bool use_sisl_marching)
// ----------------------------------------------------------------------------
{
  assert(!use_sisl_marching); // @@ not yet implemented
  
  // compute topology
  pair<SISLIntcurve**, int> cv_topo_pts = compute_topology(ss_sisl, isoval);

  // marching out curves
  CurveVec result(size_t(cv_topo_pts.second));
  transform(cv_topo_pts.first, cv_topo_pts.first + cv_topo_pts.second, back_inserter(result), 
	    [&] (SISLIntcurve* ic) {return trace_isoval_curve(ss, ic, isoval);});

  // cleaning up
  if (cv_topo_pts.second > 0)
    freeIntcrvlist(cv_topo_pts.first, cv_topo_pts.second);

  return result;

}

  //using CurveVec = std::vector<std::pair<CurvePtr, CurvePtr>>;
  
}; // end anonymous namespace 
