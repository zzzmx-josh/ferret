/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_Eigenstrain_SP_Polarization.h"

registerMooseObject("FerretApp",PMP_Eigenstrain_SP_Polarization);

//add dislocation width, replace eta with f(eta) for calculating Eigenstrain
InputParameters PMP_Eigenstrain_SP_Polarization::validParams()
{
    InputParameters params=ComputeEigenstrainBase::validParams();
    params.addClassDescription("Computes an Eigenstrain");
    params.addRequiredParam<Real>("P0","equilibrium polarization");
    params.addRequiredParam<Real>("eps0","strain for fully poled state");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    params.addCoupledVar("P1","polarization component 1");
    params.addCoupledVar("P2","polarization component 2");
    params.addCoupledVar("P3","polarization component 3");
    return params;
}
PMP_Eigenstrain_SP_Polarization::PMP_Eigenstrain_SP_Polarization(const InputParameters &parameters)
    :ComputeEigenstrainBase(parameters),
    _P0(getParam<Real>("P0")),
    _eps0(getParam<Real>("eps0")),
    _couple_P1(coupledValue("P1")),
    _couple_P2(coupledValue("P2")),
    _couple_P3(_mesh.dimension() >= 3 ? coupledValue("P3") : _zero),
    _npolar(coupledComponents("polarizations")),
    _epsilon0(declareProperty<RankTwoTensor>("epsilon0")),
    _deps0dP(declareProperty<RankThreeTensor>("deps0dP")),
    _d2eps0dPdP(declareProperty<RankFourTensor>("d2eps0dPdP"))
{}

void
PMP_Eigenstrain_SP_Polarization::computeQpEigenstrain()
{
    RankTwoTensor I = RankTwoTensor(RankTwoTensor::initIdentity);

    RealVectorValue vec_P;
    vec_P(0) = _couple_P1[_qp];
    vec_P(1) = _couple_P2[_qp];
    vec_P(2) = _couple_P3[_qp];

    Real NormP = sqrt(_couple_P1[_qp]*_couple_P1[_qp]+_couple_P2[_qp]*_couple_P2[_qp]+_couple_P3[_qp]*_couple_P3[_qp]);
    if(NormP < 1e-8)
    {
        NormP = 1e-8;
    }

    RealVectorValue dNormPdP = vec_P/NormP;
    RankTwoTensor d2NormPdPdP;
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            d2NormPdPdP(i,j) = I(i,j)/NormP-vec_P(i)*vec_P(j)/(NormP*NormP*NormP);
        }
    }

    // Polarization eigenstrain
    RankTwoTensor PP;
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        for (unsigned int j = 0; j < _npolar; ++j)
        {
          PP(i,j) = vec_P(i)*vec_P(j);
          _epsilon0[_qp](i,j) = 1.5*_eps0/_P0*( PP(i,j)/NormP-1/3*I(i,j)*NormP );
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
                _deps0dP[_qp](i,j,k) = 1.5*_eps0/_P0*( PP(i,j)*(-1/NormP/NormP*dNormPdP(k))
                                                      +(I(i,k)*vec_P(j)+vec_P(i)*I(j,k))/NormP
                                                      -1/3*I(i,j)*dNormPdP(k) );
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
                    _d2eps0dPdP[_qp](i,j,k,l) = 1.5*_eps0/_P0*( (I(i,l)*vec_P(j)+vec_P(i)*I(j,l))*(-1/NormP/NormP*dNormPdP(k))
                                               +PP(i,j)*(2/NormP/NormP/NormP*dNormPdP(l)*dNormPdP(k)-1/NormP/NormP*d2NormPdPdP(k,l))
                                               +(I(i,k)*I(j,l)+I(i,l)*I(j,k))/NormP
                                               +(I(i,k)*vec_P(j)+vec_P(i)*I(j,k))*(-1/NormP/NormP*dNormPdP(l))
                                               -1/3*I(i,j)*d2NormPdPdP(k,l) );
                }
            }
        }
    }
}
