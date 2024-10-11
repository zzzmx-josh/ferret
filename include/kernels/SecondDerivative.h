#pragma once

#include "ADKernel.h"

class SecondDerivative : public ADKernel
{
public:
  static InputParameters validParams();
  SecondDerivative(const InputParameters &parameters);

protected:
  virtual ADReal computeQpResidual() override;
  // virtual Real computeQpJacobian() override;
  // virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  // const VariableValue & _variable_to_take_the_derivative_of; // Coupled variable for v_x
  // const unsigned int _ii, _jj, _kk;
  // const unsigned int _ii;
  // const Real _len_scale;
  const MaterialProperty<Real> & _lambda11;
  const MaterialProperty<Real> & _lambda110;
  const VariableGradient & _grad_v;
  /// Coupled variable number
  // unsigned int _v_var;
  /// Coupled variable
  // const GenericVariableValue<is_ad> & _v;
  /// Multiplier for the coupled force term
};
