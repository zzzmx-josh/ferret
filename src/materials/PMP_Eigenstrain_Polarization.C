/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_Eigenstrain_Polarization.h"

registerMooseObject("FerretApp",PMP_Eigenstrain_Polarization);

//add dislocation width, replace eta with f(eta) for calculating Eigenstrain
InputParameters PMP_Eigenstrain_Polarization::validParams()
{
    InputParameters params=ComputeEigenstrainBase::validParams();
    params.addClassDescription("Computes an Eigenstrain");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    params.addCoupledVar("P1","polarization component 1");
    params.addCoupledVar("P2","polarization component 2");
    params.addCoupledVar("P3","polarization component 3");
    return params;
}
PMP_Eigenstrain_Polarization::PMP_Eigenstrain_Polarization(const InputParameters &parameters)
    :ComputeEigenstrainBase(parameters),
    _couple_P1(coupledValue("P1")),
    _couple_P2(coupledValue("P2")),
    _couple_P3(_mesh.dimension() >= 3 ? coupledValue("P3") : _zero),
    _npolar(coupledComponents("polarizations")),
    _epsilon0(declareProperty<RankTwoTensor>("epsilon0")),
    _deps0dP(declareProperty<RankThreeTensor>("deps0dP")),
    _d2eps0dPdP(declareProperty<RankFourTensor>("d2eps0dPdP")),
    _Q(getMaterialProperty<RankFourTensor>("Q"))
{}

void
PMP_Eigenstrain_Polarization::computeQpEigenstrain()
{
    RankTwoTensor I = RankTwoTensor(RankTwoTensor::initIdentity);
    
    RealVectorValue vec_P;
    vec_P(0) = _couple_P1[_qp];
    vec_P(1) = _couple_P2[_qp];
    vec_P(2) = _couple_P3[_qp];

    // Polarization eigenstrain
    _epsilon0[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                for (unsigned int l = 0; l < _npolar; ++l)
                {
                    _epsilon0[_qp](i,j) = _epsilon0[_qp](i,j) + _Q[_qp](i,j,k,l)*vec_P(k)*vec_P(l);
                }
            }
        }
    }

    _eigenstrain[_qp] = _epsilon0[_qp];

    // derivatives
    _deps0dP[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                for (unsigned int l = 0; l < _npolar; ++l)
                {
                    for (unsigned int p = 0; p < _npolar; ++p)
                    {
                        _deps0dP[_qp](i,j,p) = _deps0dP[_qp](i,j,p)+_Q[_qp](i,j,k,l)*(I(k,p)*vec_P(l)+vec_P(k)*I(l,p));
                    }
                }
            }
        }
    }


    _d2eps0dPdP[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            for (unsigned int k = 0; k < _npolar; ++k)
            {
                for (unsigned int l = 0; l < _npolar; ++l)
                {
                    for (unsigned int p = 0; p < _npolar; ++p)
                    {
                        for (unsigned int q = 0; q < _npolar; ++q)
                        {
                            _d2eps0dPdP[_qp](i,j,p,q) = _d2eps0dPdP[_qp](i,j,p,q)+_Q[_qp](i,j,k,l)*(I(k,p)*I(l,q)+I(k,q)*I(l,p));
                        }
                    }
                }
            }
        }
    }
}
