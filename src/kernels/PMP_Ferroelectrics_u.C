/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 1.1          */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#include "PMP_Ferroelectrics_u.h"

registerMooseObject("FerretApp",PMP_Ferroelectrics_u);

InputParameters PMP_Ferroelectrics_u::validParams()
{
    InputParameters params=Kernel::validParams();
    params.addClassDescription("Kernel for stress equilibrium without body force");
    params.addRequiredParam<unsigned int>("component","The component of displacements");
    params.addRequiredCoupledVar("varphi","Coupled variable");
    params.addRequiredCoupledVar("displacements", "Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    return params;
}

PMP_Ferroelectrics_u::PMP_Ferroelectrics_u(const InputParameters &parameters)
    :Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _u_var(_ndisp),
    _varphi_var(coupled("varphi")),
    _npolar(coupledComponents("polarizations")),
    _P_var(_npolar),
    _second_u(second()),
    _second_test(secondTest()),
    _second_phi(secondPhi()),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _dstressdeps(getMaterialProperty<RankFourTensor>("dstressdeps")),
    _dstressddphi(getMaterialProperty<RankThreeTensor>("dstressddphi")),
    _dstressdP(getMaterialProperty<RankThreeTensor>("dstressdP")),
    _dstressddP(getMaterialProperty<RankFourTensor>("dstressddP")),
    _flexo_status(getMaterialProperty<Real>("flexo_status")),
    _t(getMaterialProperty<RankThreeTensor>("t")),
    _dtdP(getMaterialProperty<RankFourTensor>("dtdP"))
{
    for (unsigned int i = 0; i < _ndisp; ++i)
    { _u_var[i] = coupled("displacements", i); }
    for (unsigned int i = 0; i < _npolar; ++i)
    { _P_var[i] = coupled("polarizations", i); }
}

Real PMP_Ferroelectrics_u::computeQpResidual()
{
    Real temp = 0.0;
    for(unsigned int j = 0; j < _ndisp; ++j)
    {
        temp = temp + _stress[_qp](_component,j)*_grad_test[_i][_qp](j);
        // Flexoelectricity, contribution from the method with  strain gradient
        if (_flexo_status[_qp] == 1)
        {
            for(unsigned int l = 0; l < _npolar; ++l)
            {
                temp = temp + _t[_qp](_component,j,l)*_second_test[_i][_qp](j,l);
            }
        }
    }

    return temp;
}

Real PMP_Ferroelectrics_u::computeQpJacobian()
{
    Real temp = 0.0;
    for(unsigned int j = 0; j < _ndisp; ++j)
    {
        for (unsigned int l = 0; l < _ndisp; ++l)
        {
            temp = temp + _dstressdeps[_qp](_component,j,_component,l)*_grad_test[_i][_qp](j)*_grad_phi[_j][_qp](l);
        }
    }
    return temp;
}

Real PMP_Ferroelectrics_u::computeQpOffDiagJacobian(unsigned int jvar)
{
    // K_u_u
    Real temp;
    for(unsigned int jj = 0; jj < _ndisp; ++jj)
    {
        temp = 0.0;
        if(jvar==_u_var[jj])
        {
            for(unsigned int j = 0; j < _ndisp; ++j)
            {
                for (unsigned int l = 0; l < _ndisp; ++l)
                {
                    temp = temp + _dstressdeps[_qp](_component,j,jj,l)*_grad_test[_i][_qp](j)*_grad_phi[_j][_qp](l);
                }
            }
            return temp;
        }
    }
    // K_u_varphi
    temp = 0.0;
    if(jvar==_varphi_var)
    {
        for(unsigned int j = 0; j < _ndisp; ++j)
        {
            for(unsigned int k = 0; k < _npolar; ++k)
            {
                temp = temp + _dstressddphi[_qp](_component,j,k)*_grad_test[_i][_qp](j)*_grad_phi[_j][_qp](k);
            }
        }
        return temp;
    }
    // K_u_P
    for(unsigned int jj = 0; jj < _npolar; ++jj)
    {
        temp = 0.0;
        if(jvar==_P_var[jj])
        {
            for(unsigned int j = 0; j < _ndisp; ++j)
            {
                temp = temp + _dstressdP[_qp](_component,j,jj)*_grad_test[_i][_qp](j)*_phi[_j][_qp];
                if (_flexo_status[_qp] == 1)
                {
                    for(unsigned int l = 0; l < _npolar; ++l)
                    {
                        // Only one of the following is non-zero for different methods
                        // Flexoelectricity, contribution from the method with reduced order
                        temp = temp + _dstressddP[_qp](_component,j,jj,l)*_grad_test[_i][_qp](j)*_grad_phi[_j][_qp](l);
                        // Flexoelectricity, contribution from the method with strain gradient
                        temp = temp + _dtdP[_qp](_component,j,l,jj)*_second_test[_i][_qp](j,l)*_phi[_j][_qp];
                    }
                }
            }
            return temp;
        }
    }

    return 0.0;
}
