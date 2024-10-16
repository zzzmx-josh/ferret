/****************************************************************/
/*                Panda's Multi-Physics (PMP)                   */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/*             TU-Darmstadt & Sichuan University                */
/****************************************************************/
/*==============================================================*/
/* Created by Xiandong Zhou        ver1.0.0      2022-07-23     */
/*==============================================================*/

#pragma once

#include "InitialCondition.h"

class PMP_randomIC;

class PMP_randomIC:public InitialCondition
{
public:
    static InputParameters validParams();
    PMP_randomIC(const InputParameters & parameters);

    virtual Real value(const Point & /*p*/);

protected:
    const Real &_factor;
};
