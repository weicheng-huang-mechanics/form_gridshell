#include "contactForce.h"

contactForce::contactForce(elasticRod &m_rod, timeStepper &m_stepper, double m_stiffTimes)
{
	rod = &m_rod;
    stepper = &m_stepper;

    stiffness = rod->EI * m_stiffTimes;

    Id3<<1,0,0,
         0,1,0,
         0,0,1;
}

contactForce::~contactForce()
{
	;
}

void contactForce::computeFc(MatrixXd m_contactNodes, VectorXd m_contactRefLen)
{
    
    for (int i = 0; i < rod->contact_local_num; i++)
    {
        local_nv = rod->contactMapping(i, 0);

        currecntNode1(0) = rod->x(4 * local_nv + 0);
        currecntNode1(1) = rod->x(4 * local_nv + 1);
        currecntNode1(2) = rod->x(4 * local_nv + 2);

        for (int j = 0; j < m_contactNodes.rows(); j++)
        {
            if ( rod->contactMapping(i, 1) == m_contactNodes(j, 3) )
            {
                currecntNode2(0) = m_contactNodes(j, 0);
                currecntNode2(1) = m_contactNodes(j, 1);
                currecntNode2(2) = m_contactNodes(j, 2);

                refLen = m_contactRefLen(rod->contactMapping(i, 1));
            }
        }

        tangent = currecntNode1 - currecntNode2;

        edgeLen = tangent.norm();

        tangent = tangent / tangent.norm(); 

        epsX = edgeLen / refLen - 1.0;
        Fc = stiffness * tangent * epsX;

        for (int k = 0; k < 3; k++)
        {
            currentDof = 4 * local_nv + k;

            stepper->addForce(currentDof,  Fc(k));
        }
    }

}

void contactForce::computeJc(MatrixXd m_contactNodes, VectorXd m_contactRefLen)
{
    for (int i = 0; i < rod->contact_local_num; i++)
    {
        local_nv = rod->contactMapping(i, 0);

        currecntNode1(0) = rod->x(4 * local_nv + 0);
        currecntNode1(1) = rod->x(4 * local_nv + 1);
        currecntNode1(2) = rod->x(4 * local_nv + 2);

        for (int j = 0; j < m_contactNodes.rows(); j++)
        {
            if ( rod->contactMapping(i, 1) == m_contactNodes(j, 3) )
            {
                currecntNode2(0) = m_contactNodes(j, 0);
                currecntNode2(1) = m_contactNodes(j, 1);
                currecntNode2(2) = m_contactNodes(j, 2);

                refLen = m_contactRefLen(rod->contactMapping(i, 1));
            }
        }

        tangent = currecntNode1 - currecntNode2;

        u = tangent;
        v = u.transpose();

        edgeLen = tangent.norm();

        Jc = - stiffness * ( (1/refLen - 1/edgeLen) * Id3 + ( 1/edgeLen ) * (u*v) / ( u.norm() * u.norm() ) );

        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                currentDof1 = 4 * local_nv + j;
                currentDof2 = 4 * local_nv + k;
                stepper->addJacobian( currentDof1, currentDof2, - Jc(k,j) );
            }
        }
    }
}