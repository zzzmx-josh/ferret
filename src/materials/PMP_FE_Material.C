/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_FE_Material.h"

registerMooseObject("FerretApp",PMP_FE_Material);

InputParameters PMP_FE_Material::validParams()
{
    InputParameters params = Material::validParams();
    params.addParam<Real>("a0",0,"coefficients of potential, important for configurational force");
    params.addRequiredParam<Real>("a1","coefficients of potential");
    params.addRequiredParam<Real>("a11","coefficients of potential");
    params.addRequiredParam<Real>("a12","coefficients of potential");
    params.addRequiredParam<Real>("a111","coefficients of potential");
    params.addRequiredParam<Real>("a112","coefficients of potential");
    params.addRequiredParam<Real>("a123","coefficients of potential");
    params.addRequiredParam<Real>("a1111","coefficients of potential");
    params.addRequiredParam<Real>("a1112","coefficients of potential");
    params.addRequiredParam<Real>("a1122","coefficients of potential");
    params.addRequiredParam<Real>("a1123","coefficients of potential");
    params.addParam<Real>("theta",0,"coordinate rotation");
    params.addRequiredParam<std::vector<Real>>("Permittivity", "Permittivity");
    params.addRequiredParam<std::vector<Real>>("ElecotrictConstants", "Electrostrictive coefficient ");
    params.addRequiredParam<std::vector<Real>>("FlexoConstants", "Flexoelectricity coefficient ");
    params.addRequiredParam<std::vector<Real>>("InterfaceConstants", "Interface coefficient ");
    params.addCoupledVar("E00",0,"Applied external electric load in x direction, saved in auxvariable");
    params.addCoupledVar("E01",0,"Applied external electric load in y direction, saved in auxvariable");
    params.addCoupledVar("E02",0,"Applied external electric load in z direction, saved in auxvariable");
    params.addRequiredCoupledVar("polarizations","Coupled polarizations");

    params.addCoupledVar("varphi","electric potential");
    params.addCoupledVar("P1","polarization component 1");
    params.addCoupledVar("P2","polarization component 2");
    params.addCoupledVar("P3","polarization component 3");
    return params;
}

PMP_FE_Material::PMP_FE_Material(const InputParameters & parameters)
    :Material(parameters),
    _a0(getParam<Real>("a0")),
    _a1(getParam<Real>("a1")),
    _a11(getParam<Real>("a11")),
    _a12(getParam<Real>("a12")),
    _a111(getParam<Real>("a111")),
    _a112(getParam<Real>("a112")),
    _a123(getParam<Real>("a123")),
    _a1111(getParam<Real>("a1111")),
    _a1112(getParam<Real>("a1112")),
    _a1122(getParam<Real>("a1122")),
    _a1123(getParam<Real>("a1123")),
    _theta(getParam<Real>("theta")),
    _Permittivity(getParam<std::vector<Real>>("Permittivity")),
    _ElecotrictConstants(getParam<std::vector<Real>>("ElecotrictConstants")),
    _FlexoConstants(getParam<std::vector<Real>>("FlexoConstants")),
    _InterfaceConstants(getParam<std::vector<Real>>("InterfaceConstants")),
    _npolar(coupledComponents("polarizations")),
    _couple_P1(coupledValue("P1")),
    _couple_P2(coupledValue("P2")),
    _couple_P3(_mesh.dimension() >= 3 ? coupledValue("P3") : _zero),
    _grad_couple_varphi(coupledGradient("varphi")),
    _grad_couple_P1(coupledGradient("P1")),
    _grad_couple_P2(coupledGradient("P2")),
    _grad_couple_P3(coupledGradient("P3")),
    _E00(coupledValue("E00")),
    _E01(coupledValue("E01")),
    _E02(coupledValue("E02")),
    _K(declareProperty<RankTwoTensor>("K")),
    _Q(declareProperty<RankFourTensor>("Q")),
    _g(declareProperty<RankFourTensor>("g")),
    _f(declareProperty<RankFourTensor>("f")),
    _E(declareProperty<RealVectorValue>("E")),
    _D(declareProperty<RealVectorValue>("D")),
    _GradientP(declareProperty<RankTwoTensor>("GradientP")),
    _dpsidP(declareProperty<RealVectorValue>("dpsidP")),
    _d2psidPdP(declareProperty<RankTwoTensor>("d2psidPdP")),
    _f_bulk(declareProperty<Real>("f_bulk")),
    _f_electric(declareProperty<Real>("f_electric")),
    _f_gradient(declareProperty<Real>("f_gradient"))
{
    _KCoe.zero();
    _KCoe(0,0)=_Permittivity[0];
    _KCoe(1,1)=_Permittivity[1];
    _KCoe(2,2)=_Permittivity[2];

    // Electrostrictive constants
    _QCoe(0,0,0,0) = _ElecotrictConstants[0];  // Q11
    _QCoe(1,1,1,1) = _ElecotrictConstants[0];  // Q22
    _QCoe(2,2,2,2) = _ElecotrictConstants[0];  // Q33
    _QCoe(0,0,1,1) = _ElecotrictConstants[1];  // Q12
    _QCoe(1,1,0,0) = _ElecotrictConstants[1];  // Q21
    _QCoe(0,0,2,2) = _ElecotrictConstants[1];  // Q13
    _QCoe(2,2,0,0) = _ElecotrictConstants[1];  // Q31
    _QCoe(1,1,2,2) = _ElecotrictConstants[1];  // Q23
    _QCoe(2,2,1,1) = _ElecotrictConstants[1];  // Q32
    _QCoe(0,1,0,1) = 0.5*_ElecotrictConstants[2];  // Q44
    _QCoe(0,1,1,0) = 0.5*_ElecotrictConstants[2];  // Q44
    _QCoe(1,0,1,0) = 0.5*_ElecotrictConstants[2];  // Q44
    _QCoe(1,0,0,1) = 0.5*_ElecotrictConstants[2];  // Q44
    _QCoe(0,2,0,2) = 0.5*_ElecotrictConstants[2];  // Q55
    _QCoe(0,2,2,0) = 0.5*_ElecotrictConstants[2];  // Q55
    _QCoe(2,0,0,2) = 0.5*_ElecotrictConstants[2];  // Q55
    _QCoe(2,0,2,0) = 0.5*_ElecotrictConstants[2];  // Q55
    _QCoe(1,2,1,2) = 0.5*_ElecotrictConstants[2];  // Q66
    _QCoe(1,2,2,1) = 0.5*_ElecotrictConstants[2];  // Q66
    _QCoe(2,1,1,2) = 0.5*_ElecotrictConstants[2];  // Q66
    _QCoe(2,1,2,1) = 0.5*_ElecotrictConstants[2];  // Q66

    // Interface constants
    _gCoe.zero();
    // g(i,j,k,l) != g(j,i,k,l)
    _gCoe(0,0,0,0) = _InterfaceConstants[0];  // g11
    _gCoe(1,1,1,1) = _InterfaceConstants[0];  // g22
    _gCoe(2,2,2,2) = _InterfaceConstants[0];  // g33
    _gCoe(0,0,1,1) = _InterfaceConstants[1];  // g12
    _gCoe(1,1,0,0) = _InterfaceConstants[1];  // g21
    _gCoe(0,0,2,2) = _InterfaceConstants[1];  // g13
    _gCoe(2,2,0,0) = _InterfaceConstants[1];  // g31
    _gCoe(1,1,2,2) = _InterfaceConstants[1];  // g23
    _gCoe(2,2,1,1) = _InterfaceConstants[1];  // g32
    _gCoe(0,1,0,1) = _InterfaceConstants[2];  // g44
    _gCoe(1,0,1,0) = _InterfaceConstants[2];  // g44
    _gCoe(0,2,0,2) = _InterfaceConstants[2];  // g55
    _gCoe(2,0,2,0) = _InterfaceConstants[2];  // g55
    _gCoe(1,2,1,2) = _InterfaceConstants[2];  // g66
    _gCoe(2,1,2,1) = _InterfaceConstants[2];  // g66

    // Flexoelectric constants
    _FlexoCoe.zero();
    // f=[f11   f12   0
    //    f12   f11   0
    //    0     0     f44];
    _FlexoCoe(0,0,0,0) = _FlexoConstants[0];  // f11
    _FlexoCoe(1,1,1,1) = _FlexoConstants[0];  // f22
    _FlexoCoe(2,2,2,2) = _FlexoConstants[0];  // f33
    _FlexoCoe(0,0,1,1) = _FlexoConstants[1];  // f12
    _FlexoCoe(1,1,0,0) = _FlexoConstants[1];  // f21
    _FlexoCoe(0,0,2,2) = _FlexoConstants[1];  // f13
    _FlexoCoe(2,2,0,0) = _FlexoConstants[1];  // f31
    _FlexoCoe(1,1,2,2) = _FlexoConstants[1];  // f23
    _FlexoCoe(2,2,1,1) = _FlexoConstants[1];  // f32
    _FlexoCoe(0,1,0,1) = 0.5*_FlexoConstants[2];  // f44
    _FlexoCoe(0,1,1,0) = 0.5*_FlexoConstants[2];  // f44
    _FlexoCoe(1,0,1,0) = 0.5*_FlexoConstants[2];  // f44
    _FlexoCoe(1,0,0,1) = 0.5*_FlexoConstants[2];  // f44
    _FlexoCoe(0,2,0,2) = 0.5*_FlexoConstants[2];  // f55
    _FlexoCoe(0,2,2,0) = 0.5*_FlexoConstants[2];  // f55
    _FlexoCoe(2,0,0,2) = 0.5*_FlexoConstants[2];  // f55
    _FlexoCoe(2,0,2,0) = 0.5*_FlexoConstants[2];  // f55
    _FlexoCoe(1,2,1,2) = 0.5*_FlexoConstants[2];  // f66
    _FlexoCoe(1,2,2,1) = 0.5*_FlexoConstants[2];  // f66
    _FlexoCoe(2,1,1,2) = 0.5*_FlexoConstants[2];  // f66
    _FlexoCoe(2,1,2,1) = 0.5*_FlexoConstants[2];  // f66

}

void PMP_FE_Material::computeQpProperties()
{
    // permitivity and electrostrictive constants
    _K[_qp] = _KCoe;
    // Electrostrictive
    _Q[_qp] = _QCoe;
    // interface energy coefficient
    _g[_qp] = _gCoe;
    // Flexoelectricity constants
    _f[_qp] = _FlexoCoe;

    RealVectorValue vec_P;
    vec_P(0) = _couple_P1[_qp];
    vec_P(1) = _couple_P2[_qp];
    vec_P(2) = _couple_P3[_qp];

    // Electric field
    RealVectorValue E0;
    E0(0) = _E00[_qp];                       // applied electric field
    E0(1) = _E01[_qp];
    E0(2) = _E02[_qp];
    _E[_qp] = -_grad_couple_varphi[_qp]+E0;  // total electric field

    // Electric displacement
    _D[_qp].zero();
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        _D[_qp](i) = _D[_qp](i)+vec_P(i);
        for (unsigned int j = 0; j < _npolar; ++j)
        {
            _D[_qp](i) = _D[_qp](i)+_K[_qp](i,j)*_E[_qp](j);
        }
    }

    // electric energy density
    _f_electric[_qp] = 0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
        _f_electric[_qp] = _f_electric[_qp] -vec_P(i)*_E[_qp](i);
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            _f_electric[_qp]  = _f_electric[_qp] -0.5*_K[_qp](i,j)*_E[_qp](i)*_E[_qp](j);
        }
    }

    // Polarization gradient and Gradient energy density
    for (unsigned int i = 0; i < _npolar; ++i)
    {
        _GradientP[_qp](0,i) = _grad_couple_P1[_qp](i);
        _GradientP[_qp](1,i) = _grad_couple_P2[_qp](i);
        _GradientP[_qp](2,i) = _grad_couple_P3[_qp](i);
    }

    _f_gradient[_qp] = 0;
    for(unsigned int i = 0; i < _npolar; ++i)
    {
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            for(unsigned int k = 0; k < _npolar; ++k)
            {
                for(unsigned int l = 0; l < _npolar; ++l)
                {
                    _f_gradient[_qp] = _f_gradient[_qp] +0.5*_g[_qp](i,j,k,l)*_GradientP[_qp](i,j)*_GradientP[_qp](k,l);
                }
            }
        }
    }

    //  Coordinate transformation
    RankTwoTensor R;
    R.zero();
    R(0,0) = cos(_theta/180*3.141592653);
    R(0,1) = -sin(_theta/180*3.141592653);
    R(1,0) = sin(_theta/180*3.141592653);
    R(1,1) = cos(_theta/180*3.141592653);
    R(2,2) = 1;

    Real p01 = _couple_P1[_qp];
    Real p02 = _couple_P2[_qp];
    Real p03 = _couple_P3[_qp];

    Real P1 = p01*cos(_theta/180*3.141592653) - p02*sin(_theta/180*3.141592653);
    Real P2 = p01*sin(_theta/180*3.141592653) + p02*cos(_theta/180*3.141592653);
    Real P3 = p03;

    // Landau energy
    Real P12 = P1*P1;
    Real P13 = P12*P1;
    Real P14 = P13*P1;
    Real P15 = P14*P1;
    Real P16 = P15*P1;
    Real P17 = P16*P1;
    Real P18 = P17*P1;

    Real P22 = P2*P2;
    Real P23 = P22*P2;
    Real P24 = P23*P2;
    Real P25 = P24*P2;
    Real P26 = P25*P2;
    Real P27 = P26*P2;
    Real P28 = P27*P2;

    Real P32 = P3*P3;
    Real P33 = P32*P3;
    Real P34 = P33*P3;
    Real P35 = P34*P3;
    Real P36 = P35*P3;
    Real P37 = P36*P3;
    Real P38 = P37*P3;

    _f_bulk[_qp]  = _a0 + _a1*(P12+P22+P32)+_a11*(P14+P24+P34)+_a12*(P12*P22+P22*P32+P12*P32)
                   +_a111*(P16+P26+P36)+_a112*(P12*(P24+P34)+P22*(P14+P34)+P32*(P14+P24))
                   +_a123*(P12*P22*P32)+_a1111*(P18+P28+P38)
                   +_a1112*(P16*(P22+P32)+P26*(P12+P32)+P36*(P12+P22))
                   +_a1122*(P14*P24+P24*P34+P14*P34)
                   +_a1123*(P14*P22*P32+P24*P32*P12+P34*P12*P22);

    RealVectorValue dpdP;
    RankTwoTensor d2pdPdP;
    dpdP.zero();
    d2pdPdP.zero();

    dpdP(0) = _a112*(2*P1*(P24 + P34) + 4*P13*P22 + 4*P13*P32) + 2*P1*_a1
                    + 4*P13*_a11 + 6*P15*_a111 + 8*P17*_a1111 + _a1122*(4*P13*P24 + 4*P13*P34)
                    + _a1123*(4*P13*P22*P32 + 2*P1*P24*P32 + 2*P1*P22*P34)
                    + _a1112*(2*P1*P26 + 2*P1*P36 + 6*P15*(P22 + P32))
                    + _a12*(2*P1*P22 + 2*P1*P32) + 2*P1*P22*P32*_a123;
    dpdP(1) = _a112*(2*P2*(P14 + P34) + 4*P12*P23 + 4*P23*P32) + 2*P2*_a1
                    + 4*P23*_a11 + 6*P25*_a111 + 8*P27*_a1111 + _a1122*(4*P14*P23 + 4*P23*P34)
                    + _a1123*(2*P14*P2*P32 + 4*P12*P23*P32 + 2*P12*P2*P34)
                    + _a1112*(2*P16*P2 + 2*P2*P36 + 6*P25*(P12 + P32))
                    + _a12*(2*P2*P12 + 2*P2*P32) + 2*P12*P2*P32*_a123;
    dpdP(2) = _a112*(2*P3*(P14 + P24) + 4*P12*P33 + 4*P22*P33) + 2*P3*_a1
                    + 4*P33*_a11 + 6*P35*_a111 + 8*P37*_a1111 + _a1122*(4*P14*P33 + 4*P24*P33)
                    + _a1123*(2*P14*P22*P3 + 2*P12*P24*P3 + 4*P12*P22*P33)
                    + _a1112*(2*P16*P3 + 2*P26*P3 + 6*P35*(P12 + P22))
                    + _a12*(2*P3*P12 + 2*P3*P22) + 2*P12*P22*P3*_a123;

    d2pdPdP(0,0) = 2*_a1 + _a1123*(12*P12*P22*P32 + 2*P24*P32 + 2*P22*P34)
                         + _a1112*(2*P26 + 2*P36 + 30*P14*(P22 + P32))
                         + _a12*(2*P22 + 2*P32) + 12*P12*_a11 + 30*P14*_a111
                         + 56*P16*_a1111 + _a1122*(12*P12*P24 + 12*P12*P34)
                         + _a112*(12*P12*P22 + 12*P12*P32 + 2*P24 + 2*P34) + 2*P22*P32*_a123;
    d2pdPdP(0,1) = _a112*(8*P13*P2 + 8*P1*P23) + _a1112*(12*P15*P2 + 12*P1*P25)
                         + _a1123*(8*P13*P2*P32 + 8*P1*P23*P32 + 4*P1*P2*P34)
                         + 16*P13*P23*_a1122 + 4*P1*P2*_a12 + 4*P1*P2*P32*_a123;
    d2pdPdP(0,2) = _a112*(8*P13*P3 + 8*P1*P33) + _a1112*(12*P15*P3 + 12*P1*P35)
                         + _a1123*(8*P13*P22*P3 + 4*P1*P24*P3 + 8*P1*P22*P33)
                         + 16*P13*P33*_a1122 + 4*P1*P3*_a12 + 4*P1*P22*P3*_a123;

    d2pdPdP(1,0) = _a112*(8*P13*P2 + 8*P1*P23) + _a1112*(12*P15*P2 + 12*P1*P25)
                         + _a1123*(8*P13*P2*P32 + 8*P1*P23*P32 + 4*P1*P2*P34)
                         + 16*P13*P23*_a1122 + 4*P1*P2*_a12 + 4*P1*P2*P32*_a123;
    d2pdPdP(1,1) = 2*_a1 + _a1123*(2*P14*P32 + 12*P12*P22*P32 + 2*P12*P34)
                         + _a1112*(2*P16 + 2*P36 + 30*P24*(P12 + P32))
                         + _a12*(2*P12 + 2*P32) + 12*P22*_a11 + 30*P24*_a111
                         + 56*P26*_a1111 + _a1122*(12*P14*P22 + 12*P22*P34)
                         + _a112*(2*P14 + 12*P12*P22 + 12*P22*P32 + 2*P34) + 2*P12*P32*_a123;
    d2pdPdP(1,2) = _a112*(8*P23*P3 + 8*P2*P33) + _a1112*(12*P25*P3 + 12*P2*P35)
                         + _a1123*(4*P14*P2*P3 + 8*P12*P23*P3 + 8*P12*P2*P33)
                         + 16*P23*P33*_a1122 + 4*P2*P3*_a12 + 4*P12*P2*P3*_a123;

    d2pdPdP(2,0) = _a112*(8*P13*P3 + 8*P1*P33) + _a1112*(12*P15*P3 + 12*P1*P35)
                         + _a1123*(8*P13*P22*P3 + 4*P1*P24*P3 + 8*P1*P22*P33)
                         + 16*P13*P33*_a1122 + 4*P1*P3*_a12 + 4*P1*P22*P3*_a123;
    d2pdPdP(2,1) = _a112*(8*P23*P3 + 8*P2*P33) + _a1112*(12*P25*P3 + 12*P2*P35)
                         + _a1123*(4*P14*P2*P3 + 8*P12*P23*P3 + 8*P12*P2*P33)
                         + 16*P23*P33*_a1122 + 4*P2*P3*_a12 + 4*P12*P2*P3*_a123;
    d2pdPdP(2,2) = 2*_a1 + _a1123*(2*P14*P22 + 2*P12*P24 + 12*P12*P22*P32)
                         + _a1112*(2*P16 + 2*P26 + 30*P34*(P12 + P22))
                         + _a12*(2*P12 + 2*P22) + 12*P32*_a11 + 30*P34*_a111
                         + 56*P36*_a1111 + _a1122*(12*P14*P32 + 12*P24*P32)
                         + _a112*(2*P14 + 12*P12*P32 + 2*P24 + 12*P22*P32) + 2*P12*P22*_a123;

    _dpsidP[_qp].zero();
    _d2psidPdP[_qp].zero();
    for(unsigned int i = 0; i < _npolar; ++i)
    {
        for(unsigned int j = 0; j < _npolar; ++j)
        {
            _dpsidP[_qp](i) = _dpsidP[_qp](i) + R(j,i)*dpdP(j);
            for(unsigned int k = 0; k < _npolar; ++k)
            {
                for(unsigned int l = 0; l < _npolar; ++l)
                {
                    _d2psidPdP[_qp](i,j) = _d2psidPdP[_qp](i,j) + R(k,i)*d2pdPdP(k,l)*R(l,j);
                }
            }
        }
    }

}
