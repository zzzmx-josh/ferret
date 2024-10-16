/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_FE_phi.h"

registerMooseObject("FerretApp",PMP_FE_phi);

InputParameters PMP_FE_phi::validParams()
{
    InputParameters params = Kernel::validParams();
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");
    return params;
}

PMP_FE_phi::PMP_FE_phi(const InputParameters & parameters)
    :Kernel(parameters),
    _npolar(coupledComponents("polarizations")),
    _P_var(_npolar),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _K(getMaterialProperty<RankTwoTensor>("K")),
    _D(getMaterialProperty<RealVectorValue>("D")),
    _deps0dP(getMaterialProperty<RankThreeTensor>("deps0dP"))
{
    for (unsigned int i = 0; i < _npolar; ++i)
    {  _P_var[i] = coupled("polarizations", i); }
}

Real PMP_FE_phi::computeQpResidual()
{
    Real temp = 0.0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
      temp = temp + _D[_qp](i)*_grad_test[_i][_qp](i);
    }
    return temp;
}

Real PMP_FE_phi::computeQpJacobian()
{
    Real temp = 0.0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            temp = temp - _K[_qp](i,j)*_grad_phi[_j][_qp](j)*_grad_test[_i][_qp](i);
        }
    }
    return temp;
}

Real PMP_FE_phi::computeQpOffDiagJacobian(unsigned int jvar)
{
    // K_phi_P
    RankTwoTensor I = RankTwoTensor(RankTwoTensor::initIdentity);
    Real temp;
    for(unsigned int jj = 0; jj < _npolar; ++jj)
    {
        temp = 0.0;
        if(jvar==_P_var[jj])
        {
            for(unsigned int i = 0; i < _npolar; ++i)
            {
                temp = temp + I(i,jj)*_phi[_j][_qp]*_grad_test[_i][_qp](i);
            }
            return temp;
        }
    }

    return 0;
}
