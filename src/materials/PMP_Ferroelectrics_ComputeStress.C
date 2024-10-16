/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou                    ver 2.0.0        */
/* Sichuan University                          2023-07-09       */
/*==============================================================*/

#include "PMP_Ferroelectrics_ComputeStress.h"

registerMooseObject("FerretApp",PMP_Ferroelectrics_ComputeStress);

InputParameters PMP_Ferroelectrics_ComputeStress::validParams()
{
    InputParameters params = ComputeStressBase::validParams();
    params.addClassDescription("Compute stress field and derivatives of the free energy density");
    params.addParam<Real>("strain_gradient",1,"Method to deal with strain gradient. 1 for interpolated strain gradient in elements");
    params.addRequiredCoupledVar("displacements","Coupled displacements");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    params.addCoupledVar("u1",0,"displacement component 1");
    params.addCoupledVar("u2",0,"displacement component 2");
    params.addCoupledVar("u3",0,"displacement component 3");

    return params;
}

PMP_Ferroelectrics_ComputeStress::PMP_Ferroelectrics_ComputeStress(const InputParameters &parameters)
    :ComputeStressBase(parameters),
    _strain_gradient(getParam<Real>("strain_gradient")),
    _ndisp(coupledComponents("displacements")),
    _npolar(coupledComponents("polarizations")),
    _second_grad_couple_u1(coupledSecond("u1")),
    _second_grad_couple_u2(coupledSecond("u2")),
    _second_grad_couple_u3(coupledSecond("u3")),
    _f_elastic(declareProperty<Real>("f_elastic")),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _dstressdeps(declareProperty<RankFourTensor>("dstressdeps")),
    _dstressddphi(declareProperty<RankThreeTensor>("dstressddphi")),
    _dstressdP(declareProperty<RankThreeTensor>("dstressdP")),
    _dstressddP(declareProperty<RankFourTensor>("dstressddP")),
    _dDddu(declareProperty<RankThreeTensor>("dDddu")),
    _dDddphi(declareProperty<RankTwoTensor>("dDddphi")),
    _dDdP(declareProperty<RankTwoTensor>("dDdP")),
    _dfdP(declareProperty<RealVectorValue>("dfdP")),
    _dfddP(declareProperty<RankTwoTensor>("dfddP")),
    _d2fdPdeps(declareProperty<RankThreeTensor>("d2fdPdeps")),
    _d2fdPddphi(declareProperty<RankTwoTensor>("d2fdPddphi")),
    _d2fdPdP(declareProperty<RankTwoTensor>("d2fdPdP")),
    _d2fddPddP(declareProperty<RankFourTensor>("d2fddPddP")),
    //---Contributed from Flexoelectricity---
    _deps(declareProperty<RankThreeTensor>("deps")),
    _d2fddPdeps(declareProperty<RankFourTensor>("d2fddPdeps")),
    _d2fdPddeps(declareProperty<RankFourTensor>("d2fdPddeps")),
    _t(declareProperty<RankThreeTensor>("t")),
    _dtdP(declareProperty<RankFourTensor>("dtdP")),
    // Get material properties
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _E(getMaterialProperty<RealVectorValue>("E")),
    _K(getMaterialProperty<RankTwoTensor>("K")),
    _g(getMaterialProperty<RankFourTensor>("g")),
    _flexo(getMaterialProperty<RankFourTensor>("flexo")),
    _flexo_status(getMaterialProperty<Real>("flexo_status")),
    _vec_P(getMaterialProperty<RealVectorValue>("vec_P")),
    _GradientP(getMaterialProperty<RankTwoTensor>("GradientP")),
    _dpsidP(getMaterialProperty<RealVectorValue>("dpsidP")),
    _d2psidPdP(getMaterialProperty<RankTwoTensor>("d2psidPdP")),
    _deps0dP(getMaterialProperty<RankThreeTensor>("deps0dP")),
    _d2eps0dPdP(getMaterialProperty<RankFourTensor>("d2eps0dPdP"))
{}

void PMP_Ferroelectrics_ComputeStress::computeQpStress()
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
            _f_elastic[_qp]   = _f_elastic[_qp] + 0.5*_stress[_qp](i,j)*_mechanical_strain[_qp](i,j);
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

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Flexoelectricity
    if (_flexo_status[_qp]==1)
    {
        RankThreeTensor _grad_second_u;
        _grad_second_u.zero();
        for (unsigned int i = 0; i < _ndisp; ++i)
        {
            for (unsigned int j = 0; j < _ndisp; ++j)
            {
                _grad_second_u(0,i,j) = _second_grad_couple_u1[_qp](i,j);
                _grad_second_u(1,i,j) = _second_grad_couple_u2[_qp](i,j);
                _grad_second_u(2,i,j) = _second_grad_couple_u3[_qp](i,j);
            }
        }
        // Strain gradient
        _deps[_qp].zero();
        for (unsigned int i = 0; i < _ndisp; ++i)
        {
            for (unsigned int j = 0; j < _ndisp; ++j)
            {
                for (unsigned int k = 0; k < _ndisp; ++k)
                {
                    _deps[_qp](i,j,k) = 0.5*( _grad_second_u(i,j,k) + _grad_second_u(j,i,k) );
                }

            }
        }
        //-----------contribution to residual and stiffness matrix--------------
        // ----------for u_i------------
        _dstressddP[_qp].zero();    // from flexoelectricty
        _t[_qp].zero();             // double stress tensor
        _dtdP[_qp].zero();
        if (_strain_gradient == 0)
        {
            // residual constructed with integration by parts to reduce order of derivatives in strain gradient
            for (unsigned int i = 0; i < _npolar; ++i)
            {
                for (unsigned int j = 0; j < _npolar; ++j)
                {
                    for (unsigned int k = 0; k < _npolar; ++k)
                    {
                        for (unsigned int l = 0; l < _npolar; ++l)
                        {
                            _stress[_qp](i,j) = _stress[_qp](i,j) + _flexo[_qp](i,j,k,l)*_GradientP[_qp](k,l);
                        }
                    }
                }
            }
            _dstressddP[_qp] = _flexo[_qp];
        }
        else
        {
            // with strain gradient, for IGA
            for (unsigned int i = 0; i < _npolar; ++i)
            {
                for (unsigned int j = 0; j < _npolar; ++j)
                {
                    for (unsigned int k = 0; k < _npolar; ++k)
                    {
                        for (unsigned int l = 0; l < _npolar; ++l)
                        {
                            _stress[_qp](i,j) = _stress[_qp](i,j) + 0.5*_flexo[_qp](i,j,k,l)*_GradientP[_qp](k,l);
                            // double stress tensor
                            _t[_qp](i,j,l) = _t[_qp](i,j,l) - 0.5*_flexo[_qp](i,j,k,l)*_vec_P[_qp](k);
                            _dtdP[_qp](i,j,l,k) = -0.5*_flexo[_qp](i,j,k,l);  // pay attention to the sequence of index k 
                        }
                    }
                }
            }
            _dstressddP[_qp] = 0.5*_flexo[_qp];

        }

        // ---------for P_i--------------
        _d2fddPdeps[_qp].zero();   // flexoelectricity: method with reduced order
        _d2fdPddeps[_qp].zero();      // flexoelectricity: method with interpolated strain gradient
        if (_strain_gradient == 0)
        {
            // residual constructed with integration by parts to reduce order of derivatives in strain gradient
            for (unsigned int i = 0; i < _npolar; ++i)
            {
                for (unsigned int j = 0; j < _npolar; ++j)
                {
                    for (unsigned int k = 0; k < _npolar; ++k)
                    {
                        for (unsigned int l = 0; l < _npolar; ++l)
                        {
                            _dfddP[_qp](i,j) = _dfddP[_qp](i,j) + _flexo[_qp](k,l,i,j)*_total_strain[_qp](k,l);
                            _d2fddPdeps[_qp](i,j,k,l) = _flexo[_qp](k,l,i,j);
                        }
                    }
                }
            }
        }
        else
        {
            // redidual constructed with interpolated strain gradient within the element
            for(unsigned int i = 0; i < _ndisp; ++i)
            {
               for(unsigned int j = 0; j < _ndisp; ++j)
                 {
                    for(unsigned int k = 0; k < _ndisp; ++k)
                    {
                        for (unsigned int l = 0; l < _npolar; ++l)
                        {
                            _dfdP[_qp](i)    = _dfdP[_qp](i) - 0.5*_flexo[_qp](k,l,i,j)*_deps[_qp](k,l,j);
                            _dfddP[_qp](i,j) = _dfddP[_qp](i,j) + 0.5*_flexo[_qp](k,l,i,j)*_total_strain[_qp](k,l);
                            // contribution to stiffness matrix
                            _d2fdPddeps[_qp](i,k,l,j) = -0.5*_flexo[_qp](k,l,i,j);  // pay attention to the sequence of index
                            _d2fddPdeps[_qp](i,j,k,l) = 0.5*_flexo[_qp](k,l,i,j);
                        }
                    }
                }
            }
        }

    } // end of the contribution from flexoelectricity

    // Jac=dSigma/dStrain
    // for linear elastic problem, Jac=Cijkl
    _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
