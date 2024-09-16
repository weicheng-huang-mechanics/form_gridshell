#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean

	// Physical parameters
    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	numRods = m_inputData.GetIntOpt("num-rods");			// number of rods
	deltaLen = m_inputData.GetScalarOpt("deltaLen");
	stiffTimes = m_inputData.GetScalarOpt("stiffTimes");
	speed = m_inputData.GetScalarOpt("speed");
	
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) 
		return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	time_t current_time = time(0);

	// Open an input file named after the current time
	ostringstream name;
    name << "datafiles/simDER.txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
		return;

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
		return;

	if (timeStep == Nstep)
    {
        Vector3d outputNode;
        
        for (int i = 0; i < numRods; i++)
        {
        	for (int j = 0; j < rodsVector[i]->nv; j++)
        	{
        		outputNode(0) = rodsVector[i]->x(4 * j + 0);
        		outputNode(1) = rodsVector[i]->x(4 * j + 1);
        		outputNode(2) = rodsVector[i]->x(4 * j + 2);

        		outfile << outputNode(0) << " " << outputNode(1) << " " << outputNode(2) << endl;
        	}
        }
    }
}

void world::setRodStepper()
{
	rodStart = MatrixXd::Zero(numRods, 4);

	ifstream in1("inputfile/rodStart.txt");
	for(int i = 0; i < numRods; ++i)
	{
		for(int j = 0; j < 4; ++j)
		{
			in1 >> rodStart(i,j);
		}
	}
	in1.close();

	rodEnd = MatrixXd::Zero(numRods, 4);

	ifstream in2("inputfile/rodEnd.txt");
	for(int i = 0; i < numRods; ++i)
	{
		for(int j = 0; j < 4; ++j)
		{
			in2 >> rodEnd(i,j);
		}
	}
	in2.close();

	computeContactPair();

	double dm = deltaLen * M_PI * pow(rodRadius, 2.0) * density;
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM;
	forceTol = tol * characteristicForce / dm;

	for (int i=0; i < numRods; i++)
	{
		// Set up geometry
		rodGeometry(i);	

		// Create the rod 
		rodsVector.push_back( new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
			youngM, shearM, rodLength, theta, contactMapping) );

		// Set up boundary condition
		rodBoundaryCondition(i);
	
		// setup the rod so that all the relevant variables are populated
		rodsVector[i]->setup();
		// End of rod setup

		// set up the time stepper
		stepperVector.push_back( new timeStepper(*rodsVector[i]) );
		totalForceVector.push_back( stepperVector[i]->getForce() );

		// declare the forces
		v_stretchForce.push_back( new elasticStretchingForce( *rodsVector[i], *stepperVector[i]) );
		v_bendingForce.push_back( new elasticBendingForce( *rodsVector[i], *stepperVector[i]) );
		v_twistingForce.push_back( new elasticTwistingForce( *rodsVector[i], *stepperVector[i]) );
		v_inertialForce.push_back( new inertialForce( *rodsVector[i], *stepperVector[i]) );
		v_gravityForce.push_back( new externalGravityForce( *rodsVector[i], *stepperVector[i], gVector) );
		v_contactForce.push_back( new contactForce( *rodsVector[i], *stepperVector[i], stiffTimes) );

		// Allocate every thing to prepare for the first iteration
		rodsVector[i]->updateTimeStep();
	}

	Nstep = totalTime/deltaTime;
	
	timeStep = 0;
	currentTime = 0.0;

	computeContactRefLen();
}

// Setup geometry
void world::rodGeometry(int i)
{
	xStart = rodStart(i, 0);
	yStart = rodStart(i, 1);
	xEnd = rodStart(i, 2);
	yEnd = rodStart(i, 3);

	rodLength = sqrt( (xEnd - xStart) * (xEnd - xStart) + (yEnd - yStart) * (yEnd - yStart) );

	ne_local = ceil(rodLength / deltaLen);
	nv_local = ne_local + 1;

	vertices = MatrixXd::Zero(nv_local, 3);

	direction(0) = xEnd - xStart;
	direction(1) = yEnd - yStart;
	direction(2) = 0;

	direction = direction / direction.norm();

	delta_local = rodLength / ne_local;

	for (int k = 0; k < nv_local; k++)
	{
		vertices(k, 0) = xStart + direction(0) * k * delta_local;
		vertices(k, 1) = yStart + direction(1) * k * delta_local;
		vertices(k, 2) = abs(1e-6 * rand() / (double)(RAND_MAX));
	}

	int local_contact_num;
	local_contact_num = 0;
	
	for (int k = 0; k < contactPair.rows(); k++)
	{
		if (contactPair(k, 0) == i || contactPair(k, 1) == i)
		{
			local_contact_num = local_contact_num + 1;
		}
	}

	contactMapping = MatrixXi::Zero(local_contact_num, 2);

	int ind;
	ind = 0;
	for (int k = 0; k < contactPair.rows(); k++)
	{
		if (contactPair(k, 0) == i || contactPair(k, 1) == i)
		{
			double min_dis;
			min_dis = 1e8;
			double local_dis;

			for (int m = 0; m < nv_local; m++)
			{
				local_dis = ComputeDistance( vertices(m, 0), vertices(m, 1), contactPair(k, 2), contactPair(k, 3) );

				if (local_dis < min_dis)
				{
					min_dis = local_dis;
					contactMapping(ind, 0) = m; 
				}
			}

			vertices(contactMapping(ind, 0), 0) = contactPair(k, 2) + 1e-6 * rand() / (double)(RAND_MAX);
			vertices(contactMapping(ind, 0), 1) = contactPair(k, 3) + 1e-6 * rand() / (double)(RAND_MAX);

			contactMapping(ind, 1) = k;
			ind = ind + 1;
		}
	}

    // initial theta should be zeros
    theta = VectorXd::Zero(ne_local);
}

void world::rodBoundaryCondition(int n)
{
	// for start node
	xCurrent = rodsVector[n]->getVertex(0);
	target(0) = rodEnd(n, 0);
	target(1) = rodEnd(n, 1);
	target(2) = 0;

	towards = target - xCurrent;

	if ( towards.norm() > speed * deltaTime)
	{
		xCurrent = xCurrent + towards / towards.norm() * speed * deltaTime;
		rodsVector[n]->setVertexBoundaryCondition(xCurrent,0);
	}

	// for end node
	xCurrent = rodsVector[n]->getVertex(rodsVector[n]->ne);
	target(0) = rodEnd(n, 2);
	target(1) = rodEnd(n, 3);
	target(2) = 0;

	towards = target - xCurrent;

	if ( towards.norm() > speed * deltaTime)
	{
		xCurrent = xCurrent + towards / towards.norm() * speed * deltaTime;
		rodsVector[n]->setVertexBoundaryCondition(xCurrent, rodsVector[n]->ne);
	}

}

void world::updateTimeStep()
{
	for (int i = 0 ;i < numRods; i++)
	{
		updateEachRod(i);
	}

	currentTime += deltaTime;
		
	timeStep++;

	cout << currentTime << endl;
}
	

void world::updateEachRod(int n)
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	iter = 0;

	rodBoundaryCondition(n);

	// Start with a trial solution for our solution x
	rodsVector[n]->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		rodsVector[n]->prepareForIteration();
		
		stepperVector[n]->setZero();

		// Compute the forces
		v_inertialForce[n]->computeFi();
		v_stretchForce[n]->computeFs();	
		v_bendingForce[n]->computeFb();
		v_twistingForce[n]->computeFt();

		if (currentTime < 10.0)
		{
			v_gravityForce[n]->computeFg();
		}

		prepareContact(n);

		v_contactForce[n]->computeFc(contactNodes, contactRefLen);

		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < rodsVector[n]->uncons; i++)
		{
			normf += totalForceVector[n][i] * totalForceVector[n][i];
		}

		normf = sqrt(normf);

		if (iter == 0)
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		
		if (solved == false)
		{
			// compute jacobian
			v_inertialForce[n]->computeJi();
			v_stretchForce[n]->computeJs();
			v_bendingForce[n]->computeJb();
			v_twistingForce[n]->computeJt();
			v_gravityForce[n]->computeJg();
			v_contactForce[n]->computeJc(contactNodes, contactRefLen);

			stepperVector[n]->integrator(); // Solve equations of motion
			rodsVector[n]->updateNewtonX(totalForceVector[n]); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
	
	rodsVector[n]->updateTimeStep();

	if (render) 
	{
	//	cout << "n= " << n << " t=" << currentTime << " iter=" << iter << endl;
	}
	
	if (solved == false)
	{
		cout << "Convergence error with rod " << n << endl;
		// timeStep = Nstep; // we are exiting
	}


}

int world::simulationRunning()
{
	if (timeStep<Nstep) 
		return 1;
	else 
	{
		return -1;
	}
}

int world::numPoints(int n)
{
	return rodsVector[n]->nv;
}

double world::getScaledCoordinate(int n, int i)
{
	return rodsVector[n]->x[i] * 4;
}

int world::getnumRods()
{
	return numRods;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

double world::ComputeDistance(double x1, double y1, double x2, double y2)
{
  	double ans;
  	ans = sqrt( (x1 - x2) * (x1 - x2) 
              + (y1 - y2) * (y1 - y2) );
  	return ans;
}

int world::getIntersection(double p0_x, double p0_y, double p1_x, double p1_y, 
    double p2_x, double p2_y, double p3_x, double p3_y, double &i_x, double &i_y)
{
    double s02_x, s02_y, s10_x, s10_y, s32_x, s32_y, s_numer, t_numer, denom, t;
    s10_x = p1_x - p0_x;
    s10_y = p1_y - p0_y;
    s32_x = p3_x - p2_x;
    s32_y = p3_y - p2_y;

    denom = s10_x * s32_y - s32_x * s10_y;
    if (denom == 0)
    {
      return 0;
    } 
    bool denomPositive = denom > 0;

    s02_x = p0_x - p2_x;
    s02_y = p0_y - p2_y;
    s_numer = s10_x * s02_y - s10_y * s02_x;
    if ((s_numer < 0) == denomPositive)
    {
        return 0; 
    }

    t_numer = s32_x * s02_y - s32_y * s02_x;
    if ((t_numer < 0) == denomPositive)
    {
        return 0; 
    }

    if (((s_numer > denom) == denomPositive) || ((t_numer > denom) == denomPositive))
        return 0; 

    // Collision detected
    t = t_numer / denom;

    i_x = p0_x + (t * s10_x);
    i_y = p0_y + (t * s10_y);

    return 1;
}

void world::computeContactPair()
{
	double x1Start, y1Start, x1End, y1End;
	double x2Start, y2Start, x2End, y2End;
	double contact_x, contact_y;

	contact_num = 0;

	for (int i = 0; i < numRods; i++)
	{
		x1Start = rodStart(i, 0);
		y1Start = rodStart(i, 1);
		x1End = rodStart(i, 2);
		y1End = rodStart(i, 3);

		for (int j = i+1; j < numRods; j++)
		{
			x2Start = rodStart(j, 0);
			y2Start = rodStart(j, 1);
			x2End = rodStart(j, 2);
			y2End = rodStart(j, 3);

			if (getIntersection(x1Start, y1Start, x1End, y1End, 
                          x2Start, y2Start, x2End, y2End, contact_x, contact_y) == 1)
				contact_num = contact_num + 1;
		}
	}

	contactPair = MatrixXd::Zero(contact_num, 4);

	int temp;
	temp = 0;

	for (int i = 0; i < numRods; i++)
	{
		x1Start = rodStart(i, 0);
		y1Start = rodStart(i, 1);
		x1End = rodStart(i, 2);
		y1End = rodStart(i, 3);

		for (int j = i+1; j < numRods; j++)
		{
			x2Start = rodStart(j, 0);
			y2Start = rodStart(j, 1);
			x2End = rodStart(j, 2);
			y2End = rodStart(j, 3);

			if (getIntersection(x1Start, y1Start, x1End, y1End, 
                          x2Start, y2Start, x2End, y2End, contact_x, contact_y) == 1)
			{
				contactPair(temp, 0) = i;
				contactPair(temp, 1) = j;
				contactPair(temp, 2) = contact_x;
				contactPair(temp, 3) = contact_y;

				temp = temp + 1;
			}
		}
	}
}

void world::computeContactRefLen()
{
	contactRefLen = VectorXd::Zero(contact_num);

	int local_nv_1, local_nv_2;

	Vector3d xCurrent1;
	Vector3d xCurrent2;

	for (int i = 0; i < contact_num; i++)
	{
		// find first node
		for (int j = 0; j < numRods; j++)
		{
			if ( contactPair(i, 0) == j )
			{
				for (int k = 0; k < rodsVector[j]->contact_local_num; k++)
				{
					if (rodsVector[j]->contactMapping(k, 1) == i)
					{
						local_nv_1 = rodsVector[j]->contactMapping(k, 0);
						xCurrent1(0) = rodsVector[j]->nodes(local_nv_1, 0);
						xCurrent1(1) = rodsVector[j]->nodes(local_nv_1, 1);
						xCurrent1(2) = rodsVector[j]->nodes(local_nv_1, 2);
					}
				}
			}
		}


		// find second node
		for (int j = 0; j < numRods; j++)
		{
			if ( contactPair(i, 1) == j )
			{
				for (int k = 0; k < rodsVector[j]->contact_local_num; k++)
				{
					if (rodsVector[j]->contactMapping(k, 1) == i)
					{
						local_nv_2 = rodsVector[j]->contactMapping(k, 0);
						xCurrent2(0) = rodsVector[j]->nodes(local_nv_2, 0);
						xCurrent2(1) = rodsVector[j]->nodes(local_nv_2, 1);
						xCurrent2(2) = rodsVector[j]->nodes(local_nv_2, 2);
					}
				}
			}
		}

		contactRefLen(i) = (xCurrent2 - xCurrent1).norm();
	}
}

void world::prepareContact(int n)
{
	contactNodes = MatrixXd::Zero(rodsVector[n]->contact_local_num, 4);

	int temp;
	temp = 0;

	int contactRodNum;
	int local_contact_nv;

	for (int i = 0; i < contact_num; i++)
	{
		bool ifContact;
		ifContact = false;

		if ( contactPair(i, 0) == n )
		{
			contactRodNum = contactPair(i, 1);
			ifContact = true;
		}

		if ( contactPair(i, 1) == n )
		{
			contactRodNum = contactPair(i, 0);
			ifContact = true;
		}

		if (ifContact)
		{
			for (int k = 0; k < rodsVector[contactRodNum]->contact_local_num; k++)
			{
				if (rodsVector[contactRodNum]->contactMapping(k, 1) == i)
				{
					local_contact_nv = rodsVector[contactRodNum]->contactMapping(k, 0);
				}
			}

			contactNodes(temp, 0) = rodsVector[contactRodNum]->x(4 * local_contact_nv + 0);
			contactNodes(temp, 1) = rodsVector[contactRodNum]->x(4 * local_contact_nv + 1);
			contactNodes(temp, 2) = rodsVector[contactRodNum]->x(4 * local_contact_nv + 2);
			contactNodes(temp, 3) = i;

			temp = temp + 1;
		}
	}

}
