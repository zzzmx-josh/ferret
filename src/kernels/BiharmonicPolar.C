/*
 * File: BiharmonicPolar.C
 * File Created: 9th October 2024 10:19:03 pm
 * Author: zzzmx_josh (zmxzmx1997@outlook.com)
 * -----
 * Last Modified: 9th October 2024 10:19:04 pm
 * Modified By: zzzmx_josh (zmxzmx1997@outlook.com>)
 * -----
 */


#include "BiharmonicPolar.h"

registerMooseObject("FerretApp", BiharmonicPolar);

InputParameters
BiharmonicPolar::validParams()
{
  InputParameters params = Kernel::validParams();
  return params;
}

BiharmonicPolar::BiharmonicPolar(const InputParameters & parameters)
  : Kernel(parameters),
   _second_u(second()), 
   _second_phi(secondPhi()),
  _second_test(secondTest()),
  _lambda110(getMaterialProperty<Real>("lambda11_lambda110")),
   _lambda11(getMaterialProperty<Real>("lambda11"))
{
}

Real
BiharmonicPolar::computeQpResidual()
{

  Real trace = 0.0;

  for (unsigned int i = 0; i < _mesh.dimension(); ++i)
  {
      trace += _second_u[_qp](i, i) * _second_test[_i][_qp](i, i) ; 
  }

  return trace * _lambda110[_qp] *  _lambda11[_qp];
}

Real
BiharmonicPolar::computeQpJacobian()
{

  Real trace = 0.0;

  for (unsigned int i = 0; i < _mesh.dimension(); ++i)
  {
      trace += _second_phi[_j][_qp](i, i) * _second_test[_i][_qp](i, i);  
  }

  return trace  * _lambda110[_qp] * _lambda11[_qp];

}

