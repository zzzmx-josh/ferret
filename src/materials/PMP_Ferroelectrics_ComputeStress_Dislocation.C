/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 2.0.0        */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#include "PMP_Ferroelectrics_ComputeStress_Dislocation.h"

registerMooseObject("FerretApp",PMP_Ferroelectrics_ComputeStress_Dislocation);

InputParameters PMP_Ferroelectrics_ComputeStress_Dislocation::validParams()
{
    InputParameters params = ComputeStressBase::validParams();
    params.addClassDescription("Compute stress field and derivatives of the free energy density with dislocation induced eigenstress");
    params.addRequiredCoupledVar("displacements","Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    return params;
}

PMP_Ferroelectrics_ComputeStress_Dislocation::PMP_Ferroelectrics_ComputeStress_Dislocation(const InputParameters &parameters)
    :ComputeStressBase(parameters),
    _ndisp(coupledComponents("displacements")),
    _npolar(coupledComponents("polarizations")),
    _f_elastic(declareProperty<Real>("f_elastic")),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _dstressdeps(declareProperty<RankFourTensor>("dstressdeps")),
    _dstressddphi(declareProperty<RankThreeTensor>("dstressddphi")),
    _dstressdP(declareProperty<RankThreeTensor>("dstressdP")),
    _dDddu(declareProperty<RankThreeTensor>("dDddu")),
    _dDddphi(declareProperty<RankTwoTensor>("dDddphi")),
    _dDdP(declareProperty<RankTwoTensor>("dDdP")),
    _dfdP(declareProperty<RealVectorValue>("dfdP")),
    _dfddP(declareProperty<RankTwoTensor>("dfddP")),
    _d2fdPdeps(declareProperty<RankThreeTensor>("d2fdPdeps")),
    _d2fdPddphi(declareProperty<RankTwoTensor>("d2fdPddphi")),
    _d2fdPdP(declareProperty<RankTwoTensor>("d2fdPdP")),
    _d2fddPddP(declareProperty<RankFourTensor>("d2fddPddP")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _E(getMaterialProperty<RealVectorValue>("E")),
    _K(getMaterialProperty<RankTwoTensor>("K")),
    _g(getMaterialProperty<RankFourTensor>("g")),
    _flexo(getMaterialProperty<RankFourTensor>("flexo")),
    _GradientP(getMaterialProperty<RankTwoTensor>("GradientP")),
    _dpsidP(getMaterialProperty<RealVectorValue>("dpsidP")),
    _d2psidPdP(getMaterialProperty<RankTwoTensor>("d2psidPdP")),
    _deps0dP(getMaterialProperty<RankThreeTensor>("deps0dP")),
    _d2eps0dPdP(getMaterialProperty<RankFourTensor>("d2eps0dPdP")),
    _stress_D(getMaterialProperty<RankTwoTensor>("stress_D")),
    _beta_eD(getMaterialProperty<RankTwoTensor>("beta_eD"))
{}

void PMP_Ferroelectrics_ComputeStress_Dislocation::computeQpStress()
{
    // Dislocation
    // Displacement gradient from dislocation part
    // Mechanical strain from dislocation part
    RankTwoTensor mechanical_strain_D;
    for(unsigned int i = 0; i < _ndisp; ++i)
    {
        for(unsigned int j = 0; j < _ndisp; ++j)
        {
            mechanical_strain_D(i,j) = 0.5*(_beta_eD[_qp](i,j) + _beta_eD[_qp](j,i));
        }
    }

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

    // add eigen (additional) stress from dislocations
    _stress[_qp] = _stress[_qp] + _stress_D[_qp];

    // elastic energy density
    _f_elastic[_qp]  = 0;
    for(unsigned int i = 0; i < _ndisp; ++i)
    {
        for(unsigned int j = 0; j < _ndisp; ++j)
        {
            _f_elastic[_qp]   = _f_elastic[_qp] + 0.5*_stress[_qp](i,j)*(_mechanical_strain[_qp](i,j)+mechanical_strain_D(i,j));
        }
    }

    // Calculate derivatives of the free energy density for kernels
    RankTwoTensor I = RankTwoTensor(RankTwoTensor::initIdentity);
    // for u_i
    _dstressdeps[_qp] = _elasticity_tensor[_qp];
    _dstressddphi[_qp].zero();
    _dstressdP[_qp].zero();
    for (unsigned int i = 0; i < _ndisp; ++i)
    {
        for (unsigned int j = 0; j < _ndisp; ++j)
        {
            for (unsigned int m = 0; m < _npolar; ++m)
            {
                for (unsigned int k = 0; k < _ndisp; ++k)
                {
                    for (unsigned int l = 0; l < _ndisp; ++l)
                    {
                       _dstressdP[_qp](i,j,m) = _dstressdP[_qp](i,j,m) - _elasticity_tensor[_qp](i,j,k,l)*_deps0dP[_qp](k,l,m);
                    }
                }
            }
        }
    }

    // for varphi
    _dDddu[_qp].zero();
    _dDddphi[_qp] = -_K[_qp];
    _dDdP[_qp] = I;

    // for P_i
    _dfdP[_qp] = _dpsidP[_qp] -_E[_qp];
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                _dfdP[_qp](k) = _dfdP[_qp](k) -_deps0dP[_qp](i,j,k)*_stress[_qp](i,j);
            }
        }
    }

    _dfddP[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                for (unsigned int l = 0; l < _npolar; ++l)
                {
                    _dfddP[_qp](i,j) = _dfddP[_qp](i,j) + _g[_qp](i,j,k,l)*_GradientP[_qp](k,l);
                }
            }
        }
    }

    _d2fdPdeps[_qp].zero();
    for (unsigned int m = 0; m < _npolar; ++m)
    {
        for (unsigned int n = 0; n < _ndisp; ++n)
        {
            for (unsigned int l = 0; l < _ndisp; ++l)
            {
                for (unsigned int i = 0; i < _ndisp; ++i)
                {
                    for (unsigned int j = 0; j < _ndisp; ++j)
                    {
                        _d2fdPdeps[_qp](m,n,l) = _d2fdPdeps[_qp](m,n,l) - _deps0dP[_qp](i,j,m)*_elasticity_tensor[_qp](i,j,n,l);
                    }
                }
            }
        }
    }

    _d2fdPddphi[_qp] = I;

    _d2fdPdP[_qp] = _d2psidPdP[_qp];
    for (unsigned int m = 0; m < _npolar; ++m)
    {
        for (unsigned int n = 0; n < _npolar; ++n)
        {
            for (unsigned int i = 0; i < _npolar; ++i)
            {
                for (unsigned int j = 0; j < _npolar; ++j)
                {
                    _d2fdPdP[_qp](m,n) = _d2fdPdP[_qp](m,n) - _d2eps0dPdP[_qp](i,j,m,n)*_stress[_qp](i,j);
                    for (unsigned int k = 0; k < _npolar; ++k)
                    {
                        for (unsigned int l = 0; l < _npolar; ++l)
                        {
                            _d2fdPdP[_qp](m,n) = _d2fdPdP[_qp](m,n) + _deps0dP[_qp](i,j,m)*_elasticity_tensor[_qp](i,j,k,l)*_deps0dP[_qp](k,l,n);
                        }
                    }
                }
            }
        }
    }

    _d2fddPddP[_qp] = _g[_qp];

    // Jac=dSigma/dStrain
    // for linear elastic problem, Jac=Cijkl
    _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
