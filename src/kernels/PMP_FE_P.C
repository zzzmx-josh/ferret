/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_FE_P.h"
registerMooseObject("FerretApp",PMP_FE_P);

InputParameters PMP_FE_P::validParams()
{
    InputParameters params = TimeKernel::validParams();
    params.addRequiredParam<Real>("inv_alpha","Inverse of mobility");
    params.addRequiredParam<unsigned int>("component","The component of polarizations");
    params.addRequiredCoupledVar("displacements", "Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    params.addCoupledVar("varphi","Coupled varphi");
    return params;
}

PMP_FE_P::PMP_FE_P(const InputParameters & parameters)
    :TimeKernel(parameters),
    _inv_alpha(getParam<Real>("inv_alpha")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _u_var(_ndisp),
    _varphi_var(coupled("varphi")),
    _npolar(coupledComponents("polarizations")),
    _P_var(_npolar),
    _K(getMaterialProperty<RankTwoTensor>("K")),
    _g(getMaterialProperty<RankFourTensor>("g")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _E(getMaterialProperty<RealVectorValue>("E")),
    _GradientP(getMaterialProperty<RankTwoTensor>("GradientP")),
    _GradientG(getMaterialProperty<RankTwoTensor>("GradientG")),
    _GradientF(getMaterialProperty<RankTwoTensor>("GradientF")),
    _deps0dP(getMaterialProperty<RankThreeTensor>("deps0dP")),
    _d2eps0dPdP(getMaterialProperty<RankFourTensor>("d2eps0dPdP")),
    _dpsidP(getMaterialProperty<RealVectorValue>("dpsidP")),
    _d2psidPdP(getMaterialProperty<RankTwoTensor>("d2psidPdP")),
    _f(getMaterialProperty<RankFourTensor>("f"))
{
    for (unsigned int i = 0; i < _ndisp; ++i)
    { _u_var[i] = coupled("displacements", i); }
    for (unsigned int i = 0; i < _npolar; ++i)
    { _P_var[i] = coupled("polarizations", i); }
}

Real PMP_FE_P::computeQpResidual()
{
    I = RankTwoTensor(RankTwoTensor::initIdentity);

    Real temp1 = 0.0;
    Real temp2 = 0.0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
        temp2 = temp2 + _GradientG[_qp](_component,i)*_grad_test[_i][_qp](i)+_GradientF[_qp](_component,i)*_grad_test[_i][_qp](i);
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            temp1 = temp1 - _deps0dP[_qp](i,j,_component)*_stress[_qp](i,j)*_test[_i][_qp];
        }
    }
    return temp1-_E[_qp](_component)*_test[_i][_qp]+_dpsidP[_qp](_component)*_test[_i][_qp]
          +temp2+_inv_alpha*_test[_i][_qp]*_u_dot[_qp];
}

Real PMP_FE_P::computeQpJacobian()
{
    // K_P_P
    Real temp1 = 0.0;
    Real temp2 = 0.0;
    Real temp3 = 0.0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            temp1 = temp1 - _d2eps0dPdP[_qp](i,j,_component,_component)*_stress[_qp](i,j);
            temp3 = temp3 + _g[_qp](_component,i,_component,j)*_grad_test[_i][_qp](i)*_grad_phi[_j][_qp](j);
            for(unsigned int k = 0; k < _npolar; ++k)
            {
                for(unsigned int l = 0; l < _npolar; ++l)
                {
                     temp2 = temp2 + _deps0dP[_qp](i,j,_component)*_elasticity_tensor[_qp](i,j,k,l)*_deps0dP[_qp](k,l,_component);
                }
            }
        }
    }
    return ( temp1+temp2+_d2psidPdP[_qp](_component,_component) )*_phi[_j][_qp]*_test[_i][_qp]
          +temp3+_inv_alpha*_test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp];
}

Real PMP_FE_P::computeQpOffDiagJacobian(unsigned int jvar)
{
    I = RankTwoTensor(RankTwoTensor::initIdentity);
    
    // K_P_u
    Real temp1;
    Real temp2;
    for(unsigned int ii = 0; ii < _ndisp; ++ii)
    {
        temp1 = 0.0;
        temp2 = 0.0;
        if(jvar==_u_var[ii])
        {
            for(unsigned int j = 0; j < _ndisp; ++j)
            {
                for(unsigned int l = 0; l < _ndisp; ++l)
                {
                    temp2 = temp2 + 0.5*(_f[_qp](ii,l,_component,j)+_f[_qp](ii,l,j,_component))*_grad_phi[_j][_qp](l)*_grad_test[_i][_qp](j);
                    for(unsigned int k = 0; k < _ndisp; ++k)
                    {
                        temp1 = temp1 - _elasticity_tensor[_qp](k,l,ii,j)*_deps0dP[_qp](k,l,_component)*_grad_phi[_j][_qp](j)*_test[_i][_qp];
                    }
                }
            }
            return temp1+temp2;
        }
    }
    // K_P_varphi
    if(jvar==_varphi_var)
    {
        temp1 = 0.0;
        for(unsigned int i = 0; i < _npolar; ++i)
        {
            temp1 = temp1 + I(_component,i)*_grad_phi[_j][_qp](i)*_test[_i][_qp];
        }
        return temp1;
    }
    // K_P_P
    Real temp3;
    for(unsigned int jj = 0; jj < _npolar; ++jj)
    {
        temp1 = 0.0;
        temp2 = 0.0;
        temp3 = 0.0;
        if(jvar==_P_var[jj])
        {
            for(unsigned int i = 0; i < _npolar; ++i)
            {
                for(unsigned int j = 0; j < _npolar; ++j)
                {
                    temp1 = temp1 - _d2eps0dPdP[_qp](i,j,_component,jj)*_stress[_qp](i,j);
                    temp3 = temp3 + _g[_qp](_component,i,jj,j)*_grad_test[_i][_qp](i)*_grad_phi[_j][_qp](j);
                    for(unsigned int l = 0; l < _npolar; ++l)
                    {
                        for(unsigned int k = 0; k < _npolar; ++k)
                        {
                             temp2 = temp2 + _deps0dP[_qp](i,j,_component)*_elasticity_tensor[_qp](i,j,k,l)*_deps0dP[_qp](k,l,jj);
                        }
                    }
                }
            }
            return ( temp1+temp2+_d2psidPdP[_qp](_component,jj) )*_phi[_j][_qp]*_test[_i][_qp]+temp3;
        }
    }

    return 0;
}
