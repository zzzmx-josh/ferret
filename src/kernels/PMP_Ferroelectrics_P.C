/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 1.1          */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#include "PMP_Ferroelectrics_P.h"
registerMooseObject("FerretApp",PMP_Ferroelectrics_P);

InputParameters PMP_Ferroelectrics_P::validParams()
{
    InputParameters params = TimeKernel::validParams();
    params.addClassDescription("Kernel for Allen-Cahn equation");
    params.addRequiredParam<Real>("inv_alpha","Inverse of mobility");
    params.addRequiredParam<unsigned int>("component","The component of polarizations");
    params.addRequiredCoupledVar("displacements", "Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");
    params.addCoupledVar("varphi","Coupled varphi");
    return params;
}

PMP_Ferroelectrics_P::PMP_Ferroelectrics_P(const InputParameters & parameters)
    :TimeKernel(parameters),
    _inv_alpha(getParam<Real>("inv_alpha")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _u_var(_ndisp),
    _varphi_var(coupled("varphi")),
    _npolar(coupledComponents("polarizations")),
    _P_var(_npolar),
    _second_u(second()),
    _second_test(secondTest()),
    _second_phi(secondPhi()),
    _dfdP(getMaterialProperty<RealVectorValue>("dfdP")),
    _dfddP(getMaterialProperty<RankTwoTensor>("dfddP")),
    _d2fdPdeps(getMaterialProperty<RankThreeTensor>("d2fdPdeps")),
    _d2fdPddphi(getMaterialProperty<RankTwoTensor>("d2fdPddphi")),
    _d2fdPdP(getMaterialProperty<RankTwoTensor>("d2fdPdP")),
    _d2fddPddP(getMaterialProperty<RankFourTensor>("d2fddPddP")),
    _flexo_status(getMaterialProperty<Real>("flexo_status")),
    _d2fddPdeps(getMaterialProperty<RankFourTensor>("d2fddPdeps")),
    _d2fdPddeps(getMaterialProperty<RankFourTensor>("d2fdPddeps"))
{
    for (unsigned int i = 0; i < _ndisp; ++i)
    { _u_var[i] = coupled("displacements", i); }
    for (unsigned int i = 0; i < _npolar; ++i)
    { _P_var[i] = coupled("polarizations", i); }
}

Real PMP_Ferroelectrics_P::computeQpResidual()
{
    Real temp = 0.0;
    for(unsigned int j = 0; j < _npolar; ++j)
    {
        temp = temp + _dfddP[_qp](_component,j)*_grad_test[_i][_qp](j);
    }
    return _dfdP[_qp](_component)*_test[_i][_qp]+temp+_inv_alpha*_u_dot[_qp]*_test[_i][_qp];
}

Real PMP_Ferroelectrics_P::computeQpJacobian()
{
    // K_P_P
    Real temp = _d2fdPdP[_qp](_component,_component)*_test[_i][_qp]*_phi[_j][_qp];
    for(unsigned int j = 0; j < _npolar; ++j)
    {
        for(unsigned int l = 0; l < _npolar; ++l)
        {
            temp = temp + _d2fddPddP[_qp](_component,j,_component,l)*_grad_test[_i][_qp](j)*_grad_phi[_j][_qp](l);
        }
    }
    return temp+_inv_alpha*_du_dot_du[_qp]*_test[_i][_qp]*_phi[_j][_qp];
}

Real PMP_Ferroelectrics_P::computeQpOffDiagJacobian(unsigned int jvar)
{
    // K_P_u
    Real temp;
    for(unsigned int jj = 0; jj < _ndisp; ++jj)
    {
        temp = 0.0;
        if(jvar==_u_var[jj])
        {
            for(unsigned int j = 0; j < _ndisp; ++j)
            {
                temp = temp + _d2fdPdeps[_qp](_component,jj,j)*_test[_i][_qp]*_grad_phi[_j][_qp](j);
                if (_flexo_status[_qp] == 1)
                {
                    for(unsigned int l = 0; l < _ndisp; ++l)
                    {
                        // Flexoelectricity, contribution to dfdP
                        temp = temp + _d2fdPddeps[_qp](_component,jj,l,j)*_test[_i][_qp]*_second_phi[_j][_qp](l,j);
                        // Flexoelectricity, contribution to dfddP
                        temp = temp + _d2fddPdeps[_qp](_component,j,jj,l)*_grad_test[_i][_qp](j)*_grad_phi[_j][_qp](l);
                    }
                }
            }
            return temp;
        }
    }
    // K_P_varphi
    if(jvar==_varphi_var)
    {
        temp = 0.0;
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            temp = temp + _d2fdPddphi[_qp](_component,j)*_test[_i][_qp]*_grad_phi[_j][_qp](j);
        }
        return temp;
    }
    // K_P_P
    for(unsigned int jj = 0; jj < _npolar; ++jj)
    {
        temp = 0.0;
        if(jvar==_P_var[jj])
        {
            temp = _d2fdPdP[_qp](_component,jj)*_test[_i][_qp]*_phi[_j][_qp];
            for(unsigned int j = 0; j < _npolar; ++j)
            {
                for(unsigned int l = 0; l < _npolar; ++l)
                {
                    temp = temp + _d2fddPddP[_qp](_component,j,jj,l)*_grad_test[_i][_qp](j)*_grad_phi[_j][_qp](l);
                }
            }
            return temp;
        }
    }

    return 0;
}
