/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 2.0.0        */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#include "PMP_Ferroelectrics_SP_ComputeStress_Dislocation.h"

registerMooseObject("FerretApp",PMP_Ferroelectrics_SP_ComputeStress_Dislocation);

InputParameters PMP_Ferroelectrics_SP_ComputeStress_Dislocation::validParams()
{
    InputParameters params = ComputeStressBase::validParams();
    params.addClassDescription("Compute stress field and derivatives of the free energy density");
    params.addRequiredCoupledVar("displacements","Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    return params;
}

PMP_Ferroelectrics_SP_ComputeStress_Dislocation::PMP_Ferroelectrics_SP_ComputeStress_Dislocation(const InputParameters &parameters)
    :ComputeStressBase(parameters),
    _ndisp(coupledComponents("displacements")),
    _npolar(coupledComponents("polarizations")),
    _f_elastic(declareProperty<Real>("f_elastic")),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _dstressdeps(declareProperty<RankFourTensor>("dstressdeps")),
    _dstressddphi(declareProperty<RankThreeTensor>("dstressddphi")),
    _dstressdP(declareProperty<RankThreeTensor>("dstressdP")),
    _D(declareProperty<RealVectorValue>("D")),
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
    _e(getMaterialProperty<RankThreeTensor>("e")),
    _E(getMaterialProperty<RealVectorValue>("E")),
    _K(getMaterialProperty<RankTwoTensor>("K")),
    _g(getMaterialProperty<RankFourTensor>("g")),
    _flexo(getMaterialProperty<RankFourTensor>("flexo")),
    _dedP(getMaterialProperty<RankFourTensor>("dedP")),
    _d2NPdPdP(getMaterialProperty<RankTwoTensor>("d2NPdPdP")),
    _e_const(getMaterialProperty<RankThreeTensor>("e_const")),
    _vec_P(getMaterialProperty<RealVectorValue>("vec_P")),
    _GradientP(getMaterialProperty<RankTwoTensor>("GradientP")),
    _dpsidP(getMaterialProperty<RealVectorValue>("dpsidP")),
    _d2psidPdP(getMaterialProperty<RankTwoTensor>("d2psidPdP")),
    _deps0dP(getMaterialProperty<RankThreeTensor>("deps0dP")),
    _d2eps0dPdP(getMaterialProperty<RankFourTensor>("d2eps0dPdP")),
    _stress_D(getMaterialProperty<RankTwoTensor>("stress_D")),
    _beta_eD(getMaterialProperty<RankTwoTensor>("beta_eD"))
{}

void PMP_Ferroelectrics_SP_ComputeStress_Dislocation::computeQpStress()
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
                _stress[_qp](i,j) = _stress[_qp](i,j) - _e[_qp](k,i,j)*_E[_qp](k);
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
            _f_elastic[_qp] = _f_elastic[_qp] + 0.5*_stress[_qp](i,j)*(_mechanical_strain[_qp](i,j)+mechanical_strain_D(i,j));
            for(unsigned int k = 0; k < _ndisp; ++k)
            {
                _f_elastic[_qp] = _f_elastic[_qp] -_e[_qp](i,j,k)*(_mechanical_strain[_qp](i,j)+mechanical_strain_D(i,j))*_E[_qp](k);
            }
        }
    }

    // Electric displacement
    _D[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        _D[_qp](i) = _D[_qp](i)+_vec_P[_qp](i);
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            _D[_qp](i) = _D[_qp](i)+_K[_qp](i,j)*_E[_qp](j);
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                _D[_qp](i) = _D[_qp](i)+_e[_qp](i,j,k)*_mechanical_strain[_qp](j,k);
            }
        }
    }

    // Calculate derivatives of the free energy density for kernels
    RankTwoTensor I = RankTwoTensor(RankTwoTensor::initIdentity);
    // for u_i
    _dstressdeps[_qp]  = _elasticity_tensor[_qp];
    _dstressddphi[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                _dstressddphi[_qp](i,j,k) = _e[_qp](k,i,j);
            }
        }
    }
    _dstressdP[_qp].zero();
    for (unsigned int i = 0; i < _ndisp; ++i)
    {
        for (unsigned int j = 0; j < _ndisp; ++j)
        {
            for (unsigned int m = 0; m < _npolar; ++m)
            {
                for (unsigned int k = 0; k < _ndisp; ++k)
                {
                    _dstressdP[_qp](i,j,m) = _dstressdP[_qp](i,j,m) -_dedP[_qp](k,i,j,m)*_E[_qp](k);
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
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                _dDddu[_qp](i,j,k) = _e[_qp](i,j,k);
            }
        }
    }
    _dDddphi[_qp] = -_K[_qp];
    _dDdP[_qp]    = I;
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int m = 0; m < _npolar; ++m)
        {
            for (unsigned int j = 0; j < _npolar; ++j)
            {
                for (unsigned int k = 0; k < _npolar; ++k)
                {
                    _dDdP[_qp](i,m) = _dDdP[_qp](i,m) -_dedP[_qp](i,j,k,m)*_mechanical_strain[_qp](j,k)-_e[_qp](i,j,k)*_deps0dP[_qp](j,k,m);
                }
            }
        }
    }

    // for P_i
    _dfdP[_qp] = _dpsidP[_qp] -_E[_qp];
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int m = 0; m < _npolar; ++m)
            {
                _dfdP[_qp](m) = _dfdP[_qp](m) -_deps0dP[_qp](i,j,m)*_stress[_qp](i,j);
                for (unsigned int k = 0; k < _npolar; ++k)
                {
                    _dfdP[_qp](m) = _dfdP[_qp](m) - _dedP[_qp](k,i,j,m)*_mechanical_strain[_qp](i,j)*_E[_qp](k);
                }
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
                    _d2fdPdeps[_qp](m,n,l) = _d2fdPdeps[_qp](m,n,l) - _dedP[_qp](i,n,l,m)*_E[_qp](i);
                    for (unsigned int j = 0; j < _ndisp; ++j)
                    {
                        _d2fdPdeps[_qp](m,n,l) = _d2fdPdeps[_qp](m,n,l) - _deps0dP[_qp](i,j,m)*_elasticity_tensor[_qp](i,j,n,l);
                    }
                }
            }
        }
    }

    _d2fdPddphi[_qp] = I;
    for (unsigned int m = 0; m < _npolar; ++m)
    {
        for (unsigned int n = 0; n < _npolar; ++n)
        {
            for (unsigned int i = 0; i < _npolar; ++i)
            {
                for (unsigned int j = 0; j < _npolar; ++j)
                {
                    _d2fdPdP[_qp](m,n) = _d2fdPdP[_qp](m,n) + _deps0dP[_qp](n,i,j)*_e[_qp](n,i,j)
                                                            + _dedP[_qp](n,i,j,m)*_mechanical_strain[_qp](i,j);
                }
            }
        }
    }

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
                        _d2fdPdP[_qp](m,n) = _d2fdPdP[_qp](m,n) + _deps0dP[_qp](i,j,m)*_dedP[_qp](k,i,j,n)*_E[_qp](k)
                                                                + _deps0dP[_qp](i,j,n)*_dedP[_qp](k,i,j,m)*_E[_qp](k)
                                                                + _d2NPdPdP[_qp](m,n)*_e_const[_qp](k,i,j)*_mechanical_strain[_qp](i,j)*_E[_qp](k);
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
