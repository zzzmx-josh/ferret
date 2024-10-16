/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_FE_ComputeStress.h"

registerMooseObject("FerretApp",PMP_FE_ComputeStress);

InputParameters PMP_FE_ComputeStress::validParams()
{
    InputParameters params = ComputeStressBase::validParams();
    params.addRequiredCoupledVar("displacements","Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    return params;
}

PMP_FE_ComputeStress::PMP_FE_ComputeStress(const InputParameters &parameters)
    :ComputeStressBase(parameters),
    _ndisp(coupledComponents("displacements")),
    _npolar(coupledComponents("polarizations")),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _GradientG(declareProperty<RankTwoTensor>("GradientG")),
    _GradientF(declareProperty<RankTwoTensor>("GradientF")),
    _f_elastic(declareProperty<Real>("f_elastic")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _g(getMaterialProperty<RankFourTensor>("g")),
    _f(getMaterialProperty<RankFourTensor>("f")),
    _GradientP(getMaterialProperty<RankTwoTensor>("GradientP"))
{}

void PMP_FE_ComputeStress::computeQpStress()
{
    _stress[_qp].zero();
    for(unsigned int i = 0; i < _ndisp; ++i)
    {
        for(unsigned int j = 0; j < _ndisp; ++j)
        {
            for(unsigned int k = 0; k < _ndisp; ++k)
            {
                for(unsigned int l = 0; l < _ndisp; ++l)
                {
                    _stress[_qp](i,j) = _stress[_qp](i,j) + _elasticity_tensor[_qp](i,j,k,l)*_mechanical_strain[_qp](k,l);
                }
            }
        }
    }

    // elastic energy density
    _f_elastic[_qp]  = 0;
    for(unsigned int i = 0; i < _ndisp; ++i)
    {
        for(unsigned int j = 0; j < _ndisp; ++j)
        {
            _f_elastic[_qp]   = _f_elastic[_qp] +0.5*_stress[_qp](i,j)*_mechanical_strain[_qp](i,j);
        }
    }

    // Micro-force tensor
    _GradientG[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                for (unsigned int l = 0; l < _npolar; ++l)
                {
                    _GradientG[_qp](i,j) = _GradientG[_qp](i,j) + _g[_qp](i,j,k,l)*_GradientP[_qp](k,l);
                }
            }
        }
    }

    // Flexoelectricity
    _GradientF[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                for (unsigned int l = 0; l < _npolar; ++l)
                {
                    _GradientF[_qp](i,j) = _GradientF[_qp](i,j) + _f[_qp](k,l,i,j)*_total_strain[_qp](k,l);
                }
            }
        }
    }

    // Jac=dSigma/dStrain
    // for linear elastic problem, Jac=Cijkl
    _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
