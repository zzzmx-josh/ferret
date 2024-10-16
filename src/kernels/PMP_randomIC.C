/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#include "PMP_randomIC.h"

registerMooseObject("FerretApp",PMP_randomIC);

InputParameters PMP_randomIC::validParams()
{
    InputParameters params=InitialCondition::validParams();
    params.addRequiredParam<Real>("factor","factor");

    return params;
}

PMP_randomIC::PMP_randomIC(const InputParameters &parameters)
    :InitialCondition(parameters),
    _factor(getParam<Real>("factor"))
{}


Real PMP_randomIC::value(const Point & /*p*/)
{
    return _factor*pow(-1,rand()%2)*(rand()%1001/(float)1000);
}
