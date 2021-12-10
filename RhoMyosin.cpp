// This program attempts to reproduce the Rho Myosin module of the 
// Rangamani et al. 2016 paper, using adaptive step-size control.
// Created on 9-12-2021 
//

#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cmath>

// Parameters of the integration algorithm 

const int		nvar = 46;			// number of variables
const double	tEnd = 300.0;		// max simulation time 
const double	dt0 = 0.01;			// initial step size 
const double	dtsav = 0.1;		// time interval between data points 
const double	tolerance = 1.0e-6;	// acceptable local error during numerical integration

// Create a list of variable names 
enum Variables
{
	Ca, CaM, CaCaM, Ng, NgCaM, CaMKII, Factin, CaMKIIFactin, Gactin, CaMKIIGactin, CaMKIIp, CaN, CaNact, I1, I1act, PP1, PP1act,
	Cdc42GEF, Cdc42GEFact, Cdc42GDP, Cdc42GTP, GAP, GAPact, WASP, WASPact, Arp23, Arp23act,
	SSH1, SSH1act, LIMK, LIMKact, Cofilin, Cofilinact,
	Fnewactin, B, Bp,
	RhoGEF, RhoGEFact, RhoGDP, RhoGTP, ROCK, ROCKact, MyoPpase, MyoPpaseact, MLC, MLCact
};

// Definition of the right-hand side of the ODE model
// t is the current time 
// vector x contains the current values of the variables 
// vector dxdt is filled with the values of the derivatives dx/dt
// Split up per module 

void rhsCaMKII(const double& t, const std::vector<double>& x, std::vector<double>& dxdt)
{
	// Write down the reaction fluxes 
	double v1 = 7.75 * pow(x[Ca], 3.0) - x[CaCaM];
	double v2 = 5.0 * x[Ng] * x[CaM] - x[NgCaM];
	double v3 = x[CaMKII] * x[Factin] - 4.0 * x[CaMKIIFactin];
	double v4 = x[CaMKII] * x[Gactin] - 4.0 * x[CaMKIIGactin];
	double v5 = (120.0 * pow(x[CaCaM], 4.0) * x[CaMKII]) / (pow(4.0, 4.0) + pow(x[CaCaM], 4.0)) + (x[CaMKIIp] * x[CaMKII]) / (10.0 + x[CaMKII]);
	double v6 = (15.0 * x[PP1act] * x[CaMKIIp]) / (3.0 + x[CaMKIIp]);
	double v7 = (127.0 * pow(x[CaCaM], 4.0) * x[CaN]) / (pow(0.34, 4.0) + pow(x[CaCaM], 4.0));
	double v8 = (0.34 * x[CaMKIIp] * x[CaN]) / (127.0 + x[CaN]);
	double v9 = (0.034 * x[CaNact] * x[I1]) / (4.97 + x[I1]);
	double v10 = (0.0688 * x[CaMKIIp] * x[I1act]) / (127.0 + x[I1act]);
	double v11 = (50.0 * x[I1act] * x[PP1]) / (80.0 + x[PP1]) + (2.0 * x[PP1act] * x[PP1]) / (80.0 + x[PP1]);
	double v12 = (0.07166 * x[CaMKIIp] * x[PP1act]) / (4.97 + x[PP1act]);

	// Determine all dxdt's
	dxdt[Ca] = -3.0 * v1;
	dxdt[CaM] = -v1 - v2;
	dxdt[CaCaM] = v1;
	dxdt[Ng] = -v2;
	dxdt[NgCaM] = v2;
	dxdt[CaMKII] = -v3 - v4 - v5 + v6;
	dxdt[Factin] = -v3;
	dxdt[CaMKIIFactin] = v3;
	dxdt[Gactin] = -v4;
	dxdt[CaMKIIGactin] = v4;
	dxdt[CaMKIIp] = v5 - v6;
	dxdt[CaN] = -v7 + v8;
	dxdt[CaNact] = v7 - v8;
	dxdt[I1] = -v9 + v10;
	dxdt[I1act] = v9 - v10;
	dxdt[PP1] = -v11 + v12;
	dxdt[PP1act] = v11 - v12;
}

void rhsArp23(const double& t, const std::vector<double>& x, std::vector<double>& dxdt)
{
	rhsCaMKII(t, x, dxdt);

	// Write down the reaction fluxes 
	double v1 = (0.01 * x[CaMKIIp] * x[Cdc42GEF]) / (1.0 + x[Cdc42GEF]);
	double v2 = (0.01 * x[PP1act] * x[Cdc42GEFact]) / (1.0 + x[Cdc42GEFact]);
	double v3 = (0.75 * x[Cdc42GEFact] * x[Cdc42GDP]) / (1.0 + x[Cdc42GDP]);
	double v4 = (0.1 * x[GAPact] * x[Cdc42GTP]) / (1.0 + x[Cdc42GTP]);
	double v5 = (0.01 * x[CaMKIIp] * x[GAP]) / (1.0 + x[GAP]);
	double v6 = (0.01 * x[PP1act] * x[GAPact]) / (1.0 + x[GAPact]);
	double v7 = 0.02 * x[Cdc42GTP] * x[WASP] - 0.001 * x[WASPact];
	double v8 = 0.1 * x[Arp23] * x[WASPact] - 0.0 * x[Arp23act];

	// Determine all dxdt's
	dxdt[Cdc42GEF] = -v1 + v2;
	dxdt[Cdc42GEFact] = v1 - v2;
	dxdt[Cdc42GDP] = -v3 + v4;
	dxdt[Cdc42GTP] = v3 - v4 - v7;
	dxdt[GAP] = -v5 + v6;
	dxdt[GAPact] = v5 - v6;
	dxdt[WASP] = -v7;
	dxdt[WASPact] = v7 - v8;
	dxdt[Arp23] = -v8;
	dxdt[Arp23act] = v8;
}

void rhsCofilin(const double& t, const std::vector<double>& x, std::vector<double>& dxdt)
{
	rhsArp23(t, x, dxdt);

	// Write down the reaction fluxes 
	double v1 = (0.34 * x[CaNact] * x[SSH1]) / (4.97 + x[SSH1]);
	double v2 = (127 * x[CaMKIIp] * x[SSH1act]) / (0.34 + x[SSH1act]);
	double v3 = (0.9 * x[ROCKact] * x[LIMK]) / (0.3 + x[LIMK]);
	double v4 = (0.34 * x[SSH1act] * x[LIMKact]) / (4.0 + x[LIMKact]);
	double v5 = (0.34 * x[SSH1act] * x[Cofilin]) / (4.0 + x[Cofilin]);
	double v6 = (0.34 * x[LIMKact] * x[Cofilinact]) / (4.0 + x[Cofilinact]);

	// Determine all dxdt's
	dxdt[SSH1] = -v1 + v2;
	dxdt[SSH1act] = v1 - v2;
	dxdt[LIMK] = -v3 + v4;
	dxdt[LIMKact] = v3 - v4;
	dxdt[Cofilin] = -v5 + v6;
	dxdt[Cofilinact] = v5 - v6;
}

void rhsActin(const double& t, const std::vector<double>& x, std::vector<double>& dxdt)
{
	rhsCofilin(t, x, dxdt);

	// Write down the severing and nucleation functions
	double fsev = (0.1 * 0.0002 * pow(x[Cofilinact], 4.0) * x[Factin]) / 0.0001;
	double fnuc = (15.3 * x[Arp23act] * x[Factin] * x[Gactin]) / (2.0 + x[Arp23act]);

	// Write down the function of membrane velocity 
	double vmb = 0.1 * x[Bp] / (x[Bp] + 10.0 * exp(50.0 / x[Bp]));

	// Write down the reaction fluxes 
	double v1 = 0.001 * x[Fnewactin];
	double v2 = fsev + 0.1 * x[Factin] + 0.01 * x[Factin];
	double v3 = fnuc;
	double v4 = 106.0 * (fsev + fnuc) - 0.04 * x[B];
	double v5 = (0.1 - vmb) * x[B] - 0.04 * x[Bp];

	// Determine all dxdt's
	dxdt[Fnewactin] = -v1;
	dxdt[Factin] += v1 - v2;
	dxdt[Gactin] += v2 - v3;
	dxdt[Arp23act] += -v3;
	dxdt[B] = v4;
	dxdt[Bp] = v5;
}

void rhs(const double& t, const std::vector<double>& x, std::vector<double>& dxdt)
{
	rhsActin(t, x, dxdt);

	// Write down the reaction fluxes 
	double v1 = (0.01 * x[CaMKIIp] * x[RhoGEF]) / (1.0 + RhoGEF);
	double v2 = (0.1 * x[PP1act] * x[RhoGEFact]) / (1.0 + x[RhoGEFact]);
	double v3 = (0.75 * x[RhoGEFact] * x[RhoGDP]) / (1.0 + x[RhoGDP]);
	double v4 = (0.1 * x[GAPact] * x[RhoGTP]) / (1.0 + x[RhoGTP]);
	double v5 = 0.02 * x[RhoGTP] * x[ROCK] - 0.001 * x[ROCKact];
	double v6 = 0.01 * x[MyoPpase] + (3.0 * x[MyoPpaseact] * x[MyoPpase]) / (16.0 + x[MyoPpase]);
	double v7 = (2.357 * x[ROCKact] * x[MyoPpaseact]) / (0.1 + x[MyoPpaseact]);
	double v8 = 0.01 * x[MLC] + (1.8 * x[ROCKact] * x[MLC]) / (2.47 + x[MLC]);
	double v9 = (1.0 * x[MyoPpaseact] * x[MLCact]) / (16.0 + x[MLCact]);

	// Determine all dxdt's
	dxdt[RhoGEF] = -v1 + v2;
	dxdt[RhoGEFact] = v1 - v2;
	dxdt[RhoGDP] = -v3 + v4;
	dxdt[RhoGTP] = v3 - v4 - v5;
	dxdt[ROCK] = -v5;
	dxdt[ROCKact] = v5;
	dxdt[MyoPpase] = -v6 + v7;
	dxdt[MyoPpaseact] = v6 - v7;
	dxdt[MLC] = -v8 + v9;
	dxdt[MLCact] = v8 - v9;
}

// ODE integration routine 

const double kdShrinkMax = 0.1;		// decrease step size by no more than this factor 
const double kdGrowMax = 1.2;		// increase step size by no more than this factor 
const double kdSafety = 0.9;		// safety factor in adaptive stepsize control 
const double kdMinH = 1.0e-6;		// minimum step size 

bool bogackiShampineStepper(double& t, std::vector<double>& x, std::vector<double>& dxdt1, double& h)
{
	// step 1 with function 

	// step 2
	std::vector<double> xtmp(nvar);
	for (int i = 0; i < nvar; ++i)
		xtmp[i] = x[i] + 0.5 * h * dxdt1[i];
	std::vector<double> dxdt2(nvar);
	rhs(t + 0.5 * h, xtmp, dxdt2);

	// step 3
	for (int i = 0; i < nvar; ++i)
		xtmp[i] = x[i] + 0.75 * h * dxdt2[i];
	std::vector<double> dxdt3(nvar);
	rhs(t + 0.75 * h, xtmp, dxdt3);

	// step 4
	for (int i = 0; i < nvar; ++i)
		xtmp[i] = x[i] + (1.0 / 9.0) * h * (2.0 * dxdt1[i] + 3.0 * dxdt2[i] + 4.0 * dxdt3[i]);
	std::vector<double> dxdt4(nvar);
	rhs(t + h, xtmp, dxdt4);

	double errMax = 0.0;
	for (int i = 0; i < nvar; ++i)
	{
		// Propagate solution by lowest order method 
		xtmp[i] = x[i] + h / 24.0 * (7.0 * dxdt1[i] + 6.0 * dxdt2[i] + 8.0 * dxdt3[i] + 3.0 * dxdt4[i]);

		// Compute error 
		double erri = fabs(h * ((5.0 / 72.0) * dxdt1[i] - (1.0 / 12.0) * dxdt2[i] - (1.0 / 9.0) * dxdt3[i] + (1.0 / 8.0) * dxdt4[i])) / tolerance;
		if (erri > errMax)
			errMax = erri;
	}

	// Adjust step size 
	const double fct = errMax > 0.0 ? kdSafety / std::pow(errMax, 1.0 / 3.0) : kdGrowMax;
	if (errMax > 1.0)
	{
		// Reduce step size and reject step 
		if (fct < kdShrinkMax)
			h *= kdShrinkMax;
		else
			h *= fct;
		if (h < kdMinH)
			throw std::runtime_error("step size underflow in bogackiShampineStepper().");
		return false;
	}
	else
	{
		// Update solution and increase step size 
		dxdt1 = dxdt4;
		x = xtmp;
		t += h;
		if (fct > kdGrowMax)
			h *= kdGrowMax;
		else
			h *= fct;
		return true;
	}
}

// Function main()

int main()
{
	try
	{
		// Open data file 
		std::ofstream ofs("data.csv");
		if (!ofs.is_open())
			throw std::runtime_error("unable to open file.\n");

		// Write to the datafile what you want to save 
		ofs << 't' << ',' << "Ca,CaM,CaCaM,Ng,NgCaM,CaMKII,Factin,CaMKIIFactin,Gactin,CaMKIIGactin,CaMKIIp,CaN,CaNact,I1,I1act,PP1,PP1act" << ','
			<< "Cdc42GEF,Cdc42GEFact,Cdc42GDP,Cdc42GTP,GAP,GAPact,WASP,WASPact,Arp23,Arp23act" << ','
			<< "SSH1, SSH1act, LIMK, LIMKact, Cofilin, Cofilinact" << ','
			<< "Fnewactin,B,Bp" << ','
			<< "RhoGEF,RhoGEFact,RhoGDP,RhoGTP,ROCK,ROCKact,MyoPpase,MyoPpaseact,MLC,MLCact" << '\n';

		// Set initial conditions 
		std::vector<double> x(nvar, 0.0);

		// Change the initial conditions of some to the ones listed in Table S8 
		x[Ca] = 1.0;				// Eventually take out 
		x[CaMKIIFactin] = 10.0;
		x[CaMKIIGactin] = 10.0;
		x[CaN] = 1.0;
		x[CaM] = 10.0;
		x[Ng] = 20.0;
		x[I1] = 1.8;
		x[PP1] = 0.27;
		x[WASP] = 1.0;
		x[Arp23] = 1.0;
		x[Cdc42GDP] = 1.0;
		x[Cdc42GEF] = 0.1;
		x[LIMK] = 2.0;
		x[SSH1] = 2.0;
		x[Cofilin] = 2.0;
		x[Bp] = 1.0;
		x[B] = 30.0;
		x[MyoPpaseact] = 0.1;
		x[RhoGEF] = 0.1;
		x[RhoGDP] = 1.0;
		x[ROCK] = 1.0;
		x[MyoPpase] = 1.1;
		x[MLC] = 5.0;
		x[GAP] = 0.1;

		// Step 1 once
		std::vector<double> dxdt1(nvar);
		rhs(0.0, x, dxdt1);

		// Start numerical integration 
		int nOK = 0, nStep = 0;
		double dtMin = dt0, dtMax = kdMinH;

		for (double t = 0.0, tsav = 0.0, dt = dt0; t < tEnd; ++nStep)
		{
			if (bogackiShampineStepper(t, x, dxdt1, dt))
				++nOK;

			if (dt < dtMin)
				dtMin = dt;
			else if (dt > dtMax)
				dtMax = dt;

			if (t > tsav)
			{
				// Write the data to the datafile
				ofs << t << ',';
				for (size_t i = 0; i < x.size() - 1; ++i)
					ofs << x[i] << ',';
				ofs << x.back() << '\n';

				// Update the new saving time
				tsav += dtsav;
			}
		}

		// Close datafile 
		ofs.close();

		// Report integration data 
		std::cout << "integration complete.\n"
			<< "number of steps: " << nStep << '\n'
			<< "proportion bad steps: " << 1.0 - nOK * 1.0 / nStep << '\n'
			<< "average step size: " << tEnd / nStep << '\n'
			<< "min step size: " << dtMin << '\n'
			<< "max step size: " << dtMax << '\n';
	}
	catch (std::exception& error)
	{
		std::cerr << "error: " << error.what();
		exit(EXIT_FAILURE);
	}

	// Terminate the program 
	return 0;
}
