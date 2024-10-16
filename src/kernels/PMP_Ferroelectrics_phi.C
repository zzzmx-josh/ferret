/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 1.1          */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#include "PMP_Ferroelectrics_phi.h"

registerMooseObject("FerretApp",PMP_Ferroelectrics_phi);

InputParameters PMP_Ferroelectrics_phi::validParams()
{
    InputParameters params = Kernel::validParams();
    params.addClassDescription("Kernel for Gauss's law");
    params.addRequiredCoupledVar("varphi","Coupled variable");
    params.addRequiredCoupledVar("displacements", "Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");
    return params;
}

PMP_Ferroelectrics_phi::PMP_Ferroelectrics_phi(const InputParameters & parameters)
    :Kernel(parameters),
    _ndisp(coupledComponents("displacements")),
    _u_var(_ndisp),
    _varphi_var(coupled("varphi")),
    _npolar(coupledComponents("polarizations")),
    _P_var(_npolar),
    _D(getMaterialProperty<RealVectorValue>("D")),
    _dDddu(getMaterialProperty<RankThreeTensor>("dDddu")),
    _dDddphi(getMaterialProperty<RankTwoTensor>("dDddphi")),
    _dDdP(getMaterialProperty<RankTwoTensor>("dDdP"))
{
    for (unsigned int i = 0; i < _ndisp; ++i)
    { _u_var[i] = coupled("displacements", i); }
    for (unsigned int i = 0; i < _npolar; ++i)
    { _P_var[i] = coupled("polarizations", i); }
}

Real PMP_Ferroelectrics_phi::computeQpResidual()
{
    Real temp = 0.0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
      temp = temp + _D[_qp](i)*_grad_test[_i][_qp](i);
    }
    return temp;
}

Real PMP_Ferroelectrics_phi::computeQpJacobian()
{
    Real temp = 0.0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            temp = temp + _dDddphi[_qp](i,j)*_grad_test[_i][_qp](i)*_grad_phi[_j][_qp](j);
        }
    }
    return temp;
}

Real PMP_Ferroelectrics_phi::computeQpOffDiagJacobian(unsigned int jvar)
{
    // K_phi_u
    Real temp;
    for(unsigned int jj = 0; jj < _ndisp; ++jj)
    {
        temp = 0.0;
        if(jvar==_u_var[jj])
        {
            for(unsigned int i = 0; i < _ndisp; ++i)
            {
                for(unsigned int k = 0; k < _ndisp; ++k)
                {
                    temp = temp + _dDddu[_qp](i,jj,k)*_grad_test[_i][_qp](i)*_grad_phi[_j][_qp](k);
                }
            }
            return temp;
        }
    }
    // K_phi_P
    for(unsigned int jj = 0; jj < _npolar; ++jj)
    {
        temp = 0.0;
        if(jvar==_P_var[jj])
        {
            for(unsigned int i = 0; i < _npolar; ++i)
            {
                temp = temp + _dDdP[_qp](i,jj)*_grad_test[_i][_qp](i)*_phi[_j][_qp];
            }
            return temp;
        }
    }

    return 0;
}
