/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_Dislocation_AnalyticalSolution.h"

registerMooseObject("FerretApp",PMP_Dislocation_AnalyticalSolution);

InputParameters PMP_Dislocation_AnalyticalSolution::validParams()
{
    InputParameters params = Material::validParams();
    params.addRequiredParam<std::vector<Real>>("b_n","bx by bz nx ny nz");
    params.addRequiredParam<std::vector<Real>>("r0","x0 y0 z0");
    params.addRequiredParam<Real>("h","distribution of dislocation width");
    params.addParam<Real>("angle",0,"coordinate rotation");
    params.addRequiredCoupledVar("displacements", "Coupled displacements");

    return params;
}

PMP_Dislocation_AnalyticalSolution::PMP_Dislocation_AnalyticalSolution(const InputParameters &parameters)
    :Material(parameters),
    _b_n(getParam<std::vector<Real> >("b_n")),
    _r0(getParam<std::vector<Real> >("r0")),
    _h(getParam<Real>("h")),
    _angle(getParam<Real>("angle")),
    _ndisp(coupledComponents("displacements")),
    _stress_D(declareProperty<RankTwoTensor>("stress_D")),
    _beta_eD(declareProperty<RankTwoTensor>("beta_eD")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor"))
{}

void PMP_Dislocation_AnalyticalSolution::computeQpProperties()
{
    // Singular solution of dislocations in anisotropic media from Love (1982)
    RankTwoTensor stress, GradientDisp;
    stress.zero();
    GradientDisp.zero();

    unsigned int num = _r0.size()/3;

    // constants in plane strain problem or 3D
    Real c11 = _elasticity_tensor[_qp](0,0,0,0);
    Real c22 = _elasticity_tensor[_qp](1,1,1,1);
    Real c12 = _elasticity_tensor[_qp](0,0,1,1);
    Real c66 = _elasticity_tensor[_qp](0,1,0,1);
    Real cb11 = pow(_elasticity_tensor[_qp](0,0,0,0)*_elasticity_tensor[_qp](1,1,1,1),0.5);
    Real M = (cb11+c12)*pow((cb11-c12)/(c22*c66)/(cb11+c12+2*c66),0.5);
    Real v = c12/c11;

    Real pi = 3.1415927;
    Real r4 = 0;
    Real x = 0;
    Real y = 0;
    Real z = 0;
    Real r = 0;
    Real x0 = 0;
    Real y0 = 0;
    Real z0 = 0;
    RealVectorValue vb,n,r0,e2;
    RealVectorValue coord;
    RealVectorValue ban;
    RankTwoTensor R;

    for (unsigned int i = 0; i < num; ++i)
    {
        // Burgers vector and location
        for(unsigned int j = 0; j < 3; ++j)
        {
            vb(j) = _b_n[i*6+j];
            n(j)  = _b_n[i*6+j+3];
            r0(j) = _r0[i*3+j];
        }
        ban = vb+n;

        if(ban(2)==0)
        {
            x  = _q_point[_qp](0)-r0(0);  // coordinate transformation 1 reverse x coordinate
            y  = _q_point[_qp](1)-r0(1);
            r  = pow(x*x+y*y,0.5);
            r4 = pow(x*x+cb11/c22*y*y,2)+(cb11+c12)*(cb11-c12-2*c66)*x*x*y*y/c22/c66;

            // Stress in local coordinate system
            // cut off radius r = 0.4 nm
            if(r>0.8)
            {
                stress(0,0) = stress(0,0) + M*vb(0)/(2*pi*r4*c22)*( ((cb11-c12)*(cb11+c12+2*c66)-cb11*c66)*x*x*y+cb11*cb11*c66/c22*y*y*y)
                                          + M*vb(1)*c66/(2*pi*r4)*(cb11*x*y*y/c22-x*x*x);
                stress(1,1) = stress(1,1) + M*vb(0)*c66/(2*pi*r4)*(-x*x*y+cb11*y*y*y/c22)
                                          - M*vb(1)/(2*pi*r4*cb11)*( ((cb11-c12)*(cb11+c12+2*c66)-cb11*c66)*x*y*y+c22*c66*x*x*x);
                stress(0,1) = stress(0,1) + M*vb(0)*c66/(2*pi*r4)*(-x*x*x+cb11*x*y*y/c22)+M*vb(1)*c66/(2*pi*r4)*(-x*x*y+cb11*y*y*y/c22);
                stress(1,0) = stress(0,1);

                stress(2,2) = v*(stress(0,0)+stress(1,1));
            }
        }
        if(ban(1)==0)
        {
            x0 = _q_point[_qp](0)-r0(0);  // coordinate transformation 1 reverse x coordinate
            y0 = _q_point[_qp](1)-r0(1);
            z0 = _q_point[_qp](2)-r0(2);

            // Rotate
            R.zero();
            R(0,0) = cos(_angle/180*3.141592653);
            R(0,1) = -sin(_angle/180*3.141592653);
            R(1,0) = sin(_angle/180*3.141592653);
            R(1,1) = cos(_angle/180*3.141592653);
            R(2,2) = 1;

            x = R(0,0)*x0+R(0,1)*y0;
            y = R(1,0)*x0+R(1,1)*y0;
            z = z0;

            r  = pow(x*x+z*z,0.5);
            r4 = pow(x*x+cb11/c22*z*z,2)+(cb11+c12)*(cb11-c12-2*c66)*x*x*z*z/c22/c66;

            // Stress in local coordinate system
            // cut off radius r = 0.4 nm
            if(r>0.8)
            {
                stress(0,0) = stress(0,0) + M*vb(0)/(2*pi*r4*c22)*( ((cb11-c12)*(cb11+c12+2*c66)-cb11*c66)*x*x*z+cb11*cb11*c66/c22*z*z*z)
                                          + M*vb(2)*c66/(2*pi*r4)*(cb11*x*z*z/c22-x*x*x);
                stress(2,2) = stress(2,2) + M*vb(0)*c66/(2*pi*r4)*(-x*x*z+cb11*z*z*z/c22)
                                          - M*vb(2)/(2*pi*r4*cb11)*( ((cb11-c12)*(cb11+c12+2*c66)-cb11*c66)*x*z*z+c22*c66*x*x*x);
                stress(0,2) = stress(0,2) + M*vb(0)*c66/(2*pi*r4)*(-x*x*x+cb11*x*z*z/c22)+M*vb(2)*c66/(2*pi*r4)*(-x*x*z+cb11*z*z*z/c22);
                stress(2,0) = stress(0,2);

                stress(1,1) = v*(stress(0,0)+stress(2,2));
            }

        }

        // Displacement gradient in local coordinate system  unfinished
        GradientDisp(0,0) = GradientDisp(0,0) - 0;
        GradientDisp(0,1) = GradientDisp(0,1) + 0;
        GradientDisp(1,0) = GradientDisp(1,0) - 0;
        GradientDisp(1,1) = GradientDisp(1,1) + 0;

    }

    _stress_D[_qp] = stress;
    _beta_eD[_qp]  = GradientDisp;
}
