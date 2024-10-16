/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_Dislocation_NumericalSolution.h"

registerMooseObject("FerretApp",PMP_Dislocation_NumericalSolution);

InputParameters PMP_Dislocation_NumericalSolution::validParams()
{
    InputParameters params = Material::validParams();

    params.addCoupledVar("sigma_D_11",0, "Component of additional stress field");
    params.addCoupledVar("sigma_D_22",0, "Component of additional stress field");
    params.addCoupledVar("sigma_D_33",0, "Component of additional stress field");
    params.addCoupledVar("sigma_D_12",0, "Component of additional stress field");
    params.addCoupledVar("sigma_D_13",0, "Component of additional stress field");
    params.addCoupledVar("sigma_D_23",0, "Component of additional stress field");

    // This part is not required. It is used for Eshelby tensor
    params.addCoupledVar("beta_eD_11",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_12",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_13",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_21",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_22",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_23",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_31",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_32",0, "Component of additional mechanical strain");
    params.addCoupledVar("beta_eD_33",0, "Component of additional mechanical strain");

    return params;
}

PMP_Dislocation_NumericalSolution::PMP_Dislocation_NumericalSolution(const InputParameters &parameters)
    :Material(parameters),
     _sigma_D_11(coupledValue("sigma_D_11")),
     _sigma_D_22(coupledValue("sigma_D_22")),
     _sigma_D_33(_mesh.dimension() >= 3 ? coupledValue("sigma_D_33") : _zero),
     _sigma_D_12(coupledValue("sigma_D_12")),
     _sigma_D_13(_mesh.dimension() >= 3 ? coupledValue("sigma_D_13") : _zero),
     _sigma_D_23(_mesh.dimension() >= 3 ? coupledValue("sigma_D_23") : _zero),
     _beta_eD_11(coupledValue("beta_eD_11")),
     _beta_eD_12(coupledValue("beta_eD_12")),
     _beta_eD_13(_mesh.dimension() >= 3 ? coupledValue("beta_eD_13") : _zero),
     _beta_eD_21(coupledValue("beta_eD_21")),
     _beta_eD_22(coupledValue("beta_eD_22")),
     _beta_eD_23(_mesh.dimension() >= 3 ? coupledValue("beta_eD_23") : _zero),
     _beta_eD_31(_mesh.dimension() >= 3 ? coupledValue("beta_eD_31") : _zero),
     _beta_eD_32(_mesh.dimension() >= 3 ? coupledValue("beta_eD_32") : _zero),
     _beta_eD_33(_mesh.dimension() >= 3 ? coupledValue("beta_eD_33") : _zero),
     _stress_D(declareProperty<RankTwoTensor>("stress_D")),
     _beta_eD(declareProperty<RankTwoTensor>("beta_eD"))
{}

void PMP_Dislocation_NumericalSolution::computeQpProperties()
{
    // Non-singular solution of dislocations from Cai Wei (2006)
    // Numerical solution
    _stress_D[_qp](0,0) = _sigma_D_11[_qp];
    _stress_D[_qp](0,1) = _sigma_D_12[_qp];
    _stress_D[_qp](0,2) = _sigma_D_13[_qp];
    _stress_D[_qp](1,0) = _sigma_D_12[_qp];
    _stress_D[_qp](1,1) = _sigma_D_22[_qp];
    _stress_D[_qp](1,2) = _sigma_D_23[_qp];
    _stress_D[_qp](2,0) = _sigma_D_13[_qp];
    _stress_D[_qp](2,1) = _sigma_D_23[_qp];
    _stress_D[_qp](2,2) = _sigma_D_33[_qp];

    _beta_eD[_qp](0,0) = _beta_eD_11[_qp];
    _beta_eD[_qp](0,1) = _beta_eD_12[_qp];
    _beta_eD[_qp](0,2) = _beta_eD_13[_qp];
    _beta_eD[_qp](1,0) = _beta_eD_21[_qp];
    _beta_eD[_qp](1,1) = _beta_eD_22[_qp];
    _beta_eD[_qp](1,2) = _beta_eD_23[_qp];
    _beta_eD[_qp](2,0) = _beta_eD_31[_qp];
    _beta_eD[_qp](2,1) = _beta_eD_32[_qp];
    _beta_eD[_qp](2,2) = _beta_eD_33[_qp];
}
