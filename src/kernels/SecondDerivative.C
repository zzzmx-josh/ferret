/*
 * File: SecondDerivative.C
 * File Created: 9th October 2024 10:20:17 pm
 * Author: zzzmx_josh (zmxzmx1997@outlook.com)
 * -----
 * Last Modified: 9th October 2024 10:20:17 pm
 * Modified By: zzzmx_josh (zmxzmx1997@outlook.com>)
 * -----
 */

#include "SecondDerivative.h"

class SecondDerivative;

registerMooseObject("FerretApp", SecondDerivative);

InputParameters SecondDerivative::validParams()
{
  InputParameters params = ADKernel::validParams();
  // params.addRequiredParam<unsigned int>("component", "The component (0: x, 1: y, 2: z) to apply this ADKernel to.");
  params.addRequiredCoupledVar("variable_to_take_the_derivative_of", "The coupled velocity variable v.");
  return params;
}

SecondDerivative:: SecondDerivative(const InputParameters &parameters)
    : ADKernel(parameters),
      // _v(coupledValue("variable_to_take_the_derivative_of")),
      // _component(getParam<unsigned int>("component")),
      // _ii(getParam<unsigned int>("component")),
      _lambda11(getMaterialProperty<Real>("lambda11")),
      _lambda110(getMaterialProperty<Real>("lambda11_lambda110")),
      _grad_v(coupledGradient("variable_to_take_the_derivative_of"))
{
}

ADReal  SecondDerivative::computeQpResidual()
{
  return -_grad_v[_qp] * _grad_test[_i][_qp] * _lambda11[_qp] * _lambda110[_qp];
}

// Real  SecondDerivative::computeQpJacobian()
// {
//   return 0;
// }

// Real  SecondDerivative::computeQpOffDiagJacobian(unsigned int jvar)
// {
//   return 0;
// }