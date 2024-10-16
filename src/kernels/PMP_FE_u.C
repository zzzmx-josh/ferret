/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_FE_u.h"

registerMooseObject("FerretApp",PMP_FE_u);

InputParameters PMP_FE_u::validParams()
{
    InputParameters params=Kernel::validParams();
    params.addRequiredParam<unsigned int>("component","The component of displacements");
    params.addRequiredCoupledVar("varphi","Coupled variable");
    params.addRequiredCoupledVar("displacements", "Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    return params;
}

PMP_FE_u::PMP_FE_u(const InputParameters &parameters)
    :Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _u_var(_ndisp),
    _varphi_var(coupled("varphi")),
    _npolar(coupledComponents("polarizations")),
    _P_var(_npolar),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _E(getMaterialProperty<RealVectorValue>("E")),
    _deps0dP(getMaterialProperty<RankThreeTensor>("deps0dP")),
    _f(getMaterialProperty<RankFourTensor>("f")),
    _GradientP(getMaterialProperty<RankTwoTensor>("GradientP"))
{
    for (unsigned int i = 0; i < _ndisp; ++i)
    { _u_var[i] = coupled("displacements", i); }
    for (unsigned int i = 0; i < _npolar; ++i)
    { _P_var[i] = coupled("polarizations", i); }
}

Real PMP_FE_u::computeQpResidual()
{
    Real temp1 = 0.0;
    Real temp2 = 0.0;
    for(unsigned int j = 0; j < _ndisp; ++j)
    {
        temp1 = temp1 + _stress[_qp](_component,j)*_grad_test[_i][_qp](j);
        for(unsigned int k = 0; k < _npolar; ++k)
        {
            for(unsigned int l = 0; l < _npolar; ++l)
            {
                temp2 = temp2 + 0.5*(_f[_qp](_component,j,k,l)+_f[_qp](j,_component,k,l))*_GradientP[_qp](k,l)*_grad_test[_i][_qp](j);
            }
        }
    }
    return temp1+temp2;
}

Real PMP_FE_u::computeQpJacobian()
{
    Real temp = 0.0;
    for(unsigned int j = 0; j < _ndisp; ++j)
    {
        for (unsigned int l = 0; l < _ndisp; ++l)
        {
            temp = temp + _elasticity_tensor[_qp](_component,j,_component,l)*_grad_phi[_j][_qp](l)*_grad_test[_i][_qp](j);
        }
    }
    return temp;
}

Real PMP_FE_u::computeQpOffDiagJacobian(unsigned int jvar)
{
    // K_u_u
    Real temp1;
    for(unsigned int ii = 0; ii < _ndisp; ++ii)
    {
        temp1 = 0.0;
        if(jvar==_u_var[ii])
        {
            for(unsigned int j = 0; j < _ndisp; ++j)
            {
                for (unsigned int l = 0; l < _ndisp; ++l)
                {
                    temp1 = temp1 + _elasticity_tensor[_qp](_component,j,ii,l)*_grad_phi[_j][_qp](l)*_grad_test[_i][_qp](j);
                }
            }
            return temp1;
        }
    }
    // K_u_P
    Real temp2;
    for(unsigned int jj = 0; jj < _npolar; ++jj)
    {
        temp1 = 0.0;
        temp2 = 0.0;
        if(jvar==_P_var[jj])
        {
            for(unsigned int j = 0; j < _npolar; ++j)
            {
                for(unsigned int l = 0; l < _npolar; ++l)
                {
                    temp2 = temp2 + 0.5*(_f[_qp](_component,j,jj,l)+_f[_qp](j,_component,jj,l))*_grad_phi[_j][_qp](l)*_grad_test[_i][_qp](j);
                    for(unsigned int k = 0; k < _npolar; ++k)
                    {
                        temp1 = temp1 - _elasticity_tensor[_qp](_component,j,k,l)*_deps0dP[_qp](k,l,jj)*_phi[_j][_qp]*_grad_test[_i][_qp](j);
                    }
                }
            }
            return temp1+temp2;
        }
    }

    return 0.0;
}
