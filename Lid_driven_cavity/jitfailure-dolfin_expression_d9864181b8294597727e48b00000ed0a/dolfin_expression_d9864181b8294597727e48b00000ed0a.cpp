
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_d9864181b8294597727e48b00000ed0a : public Expression
  {
     public:
       double nu_0;
double nu_1;
double sig_tight;
double vort_x;
double vort_y;
double sig_center;


       dolfin_expression_d9864181b8294597727e48b00000ed0a()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = nu_0+(nu_1-nu_0)/(1+exp(sig_tight*((sqrt(pow(x[0]-vort_x,2)+pow(x[1]-vort_y,2)))-sig_center))));

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "nu_0") { nu_0 = _value; return; }          if (name == "nu_1") { nu_1 = _value; return; }          if (name == "sig_tight") { sig_tight = _value; return; }          if (name == "vort_x") { vort_x = _value; return; }          if (name == "vort_y") { vort_y = _value; return; }          if (name == "sig_center") { sig_center = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "nu_0") return nu_0;          if (name == "nu_1") return nu_1;          if (name == "sig_tight") return sig_tight;          if (name == "vort_x") return vort_x;          if (name == "vort_y") return vort_y;          if (name == "sig_center") return sig_center;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {

       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {

       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_d9864181b8294597727e48b00000ed0a()
{
  return new dolfin::dolfin_expression_d9864181b8294597727e48b00000ed0a;
}

