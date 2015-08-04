#pragma hdrstop
#include "NLEQSolver.h"
#include "rrExecutableModel.h"
#include "rrStringUtils.h"
#include "rrUtils.h"
#include "nleq/nleq1.h"
#include "rrLogger.h"
#include "rrUtils.h"
#include "rrException.h"
#include "rrConfig.h"

#include <Poco/Mutex.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <unordered_map>

namespace rr
{

// NLEQ is an ancient Fortran77 routine that assumes that there is only
// one program which has a hard coded function in it.
// So, there is no concept of a user suplied data block, have to store
// the model in this static location -- only a single thread at a time
// may use the nleq steady state.
static ExecutableModel* callbackModel = NULL;

// mutex to ensure only one thead
using Poco::Mutex;
static Mutex mutex;

// the NLEQ callback, we use same data types as f2c here.
static void ModelFunction(int* nx, double* y, double* fval, int* pErr);

static string ErrorForStatus(int error);

static bool isError(int e)
{
    return e != 4 && e != 5;
}

static bool isWarning(int e)
{
    //case 4:     return ("Warning: Superlinear or quadratic convergence slowed down near the solution");
    //case 5:     return ("Warning: Error Tolerance reached but solution is suspect");

    return e == 4 || e == 5;
}

NLEQSolver::NLEQSolver(ExecutableModel *_model) :
    IWK(0),
    LIWK(0),
    LWRK(0),
    RWK(0),
    XScal(0),
    ierr(0),
    iopt(0),
    model(0),
    nOpts(50)
{
	addSetting("maximum_iterations", 100, "Specifies ... (int)", "(int) ...");
	addSetting("minimum_damping", 1e-4, "Specifies ... (double).", "(double) ...");
	addSetting("relative_tolerance", 1e-16, "Specifies ... (double).", "(double) ...");
	NLEQSolver::loadConfigSettings();

    model = _model;

    if(model)
    {
        setup();
    }
}

NLEQSolver::~NLEQSolver()
{
    delete[] IWK;
    delete[] RWK;
    delete[] XScal;
    delete[] iopt;
}

void NLEQSolver::loadConfigSettings()
{
	Solver::loadConfigSettings();
	// Load settings specific to solver integrator

	NLEQSolver::setValue("maximum_iterations", Config::getInt(Config::STEADYSTATE_MAXIMUM_NUM_STEPS));
	NLEQSolver::setValue("relative_tolerance", Config::getDouble(Config::STEADYSTATE_RELATIVE));
	NLEQSolver::setValue("minimum_damping", Config::getDouble(Config::STEADYSTATE_MINIMUM_DAMPING));
}

void NLEQSolver::setup()
{
    // size of state vector
    n = model->getStateVector(0);

    // Allocate space, see NLEQ docs for details
    LWRK = (n + 2 + 15)*n + 61;
    LIWK = n + 52;

    XScal = new double[n];
    // Set up default scaling factors
    for (int i = 0; i < n; i++)
    {
        XScal[i] = 1.0;
    }

    iopt = new long[nOpts];
    for (int i = 0; i < nOpts; i++)
    {
        iopt[i] = 0;
    }

    // Set for Highly nonlinear problem
    iopt[31 - 1] = 4;

    // Initialise all array elements to 0.0
    IWK = new long[LIWK];
    for (int i = 0; i < LIWK; i++)
    {
        IWK[i] = 0;
    }
	IWK[31 - 1] = getValueAsInt("maximum_iterations");


    RWK = new double[LWRK];
    for (int i = 0; i < LWRK; i++)
    {
        RWK[i] = 0.0;
    }

	RWK[22 - 1] = getValueAsDouble("minimum_damping");
}

bool NLEQSolver::isAvailable()
{
	NLEQSolver *temp = new NLEQSolver(NULL);
    if(temp)
    {
        bool val = temp->getNumberOfModelEvaluations() ? true : false;
        delete temp;
        return val;
    }

    return false;
}

std::string NLEQSolver::getSolverName() const {
	return NLEQSolver::getName();
}

std::string NLEQSolver::getName() {
	return "nleq";
}

std::string NLEQSolver::getSolverDescription() const {
	return NLEQSolver::getDescription();
}

std::string NLEQSolver::getDescription() {
	return "CVODE is a deterministic ODE solver from the SUNDIALS suite "
		"of timecourse integrators. It implements an Adams-Moulton solver "
		"for non-stiff problems and a backward differentiation formula "
		"(BDF) solver for stiff problems.";
}

std::string NLEQSolver::getSolverHint() const {
	return NLEQSolver::getHint();
}

std::string NLEQSolver::getHint() {
	return "Steady-state nonlinear systems of equations solver";
}

Solver::SolverMethod NLEQSolver::getSolverMethod() const
{
	return Solver::SolverMethod::SteadyState;
}

void NLEQSolver::setValue(string key, const Variant& val)
{
	NLEQSolver::setValue(key, val);

	/*	In addition to typically value-setting behavior, some settings require further changes
	within NLEQ. */
	if (key == "maximum_iterations")
	{
		IWK[31 - 1] = getValueAsInt("maximum_iterations");
	}
	else if (key == "minimum_damping")
	{
		RWK[22 - 1] = getValueAsDouble("minimum_damping");
	}
	// Don't have to worry about relative tolerance because it is always assigned as a temp variable with the solve() method.
}

double NLEQSolver::solve(const vector<double>& yin)
{
    // lock so only one thread can be here.
    Mutex::ScopedLock lock(mutex);

    if (yin.size() == 0)
    {
        return 0;
    }

    // Set up a dummy Jacobian, actual Jacobian is computed
    // by NLEQ using finite differences
    //    double* Jacobian = new double[1];

    ierr = 0;
	IWK[31 - 1] = getValueAsInt("maximum_iterations");

    // Set up default scaling factors
    for (int i = 0; i < n; i++)
    {
        XScal[i] = 1.0;
    }

    for (int i = 0; i < nOpts; i++)
    {
        iopt[i] = 0;
    }

    iopt[31 - 1] = 3; // Set for Highly nonlinear problem

    // Initialise all array elements to 0.0
    for (int i = 0; i < LIWK; i++)
    {
        IWK[i] = 0;
    }

	IWK[31 - 1] = getValueAsInt("maximum_iterations");
    for (int i = 0; i < LWRK; i++)
    {
        RWK[i] = 0.0;
    }

    RWK[22 - 1] = 1E-20; // Minimal allowed damping factor

    // For some reason NLEQ modifies the tolerance value, use a local copy instead
	double tmpTol = getValueAsDouble("relative_tolerance");

    // set up the thread local variables, only this thread
    // access them.
    if (callbackModel)
    {
        throw Exception("global callbackModel is set, this should never happen!");
    }

    // should be exception safe, just make sure

    try
    {
        callbackModel = model;
        vector<double> stateVector(n);
        model->getStateVector(&stateVector[0]);

        NLEQ1(  &n,
                &ModelFunction,
                NULL,
                &stateVector[0],
                XScal,
                &tmpTol,
                iopt,
                &ierr,
                &LIWK,
                IWK,
                &LWRK,
                RWK);

        // done, clear it.
        callbackModel = NULL;
    }
    catch(...)
    {
        // clear the global model and re-throw the exception.
        callbackModel = NULL;
        throw;
    }

    if (ierr == 2) // retry
    {
        for (int i = 0; i < nOpts; i++)
        {
            iopt[i] = 0;
        }

        iopt[31 - 1] = 3; // Set for Highly nonlinear problem
        iopt[0] = 1; // Try again but tell NLEQ not to reinitialize
		tmpTol = getValueAsDouble("relative_tolerance");

    }

    if(ierr > 0 )
    {
        if (isWarning(ierr)) {
            Log(Logger::LOG_WARNING) << ErrorForStatus(ierr);
        } else {
            string err = ErrorForStatus(ierr);
            Log(Logger::LOG_ERROR)<<"Error :"<<err;
            throw NLEQException(err);
        }
    }

    return computeSumsOfSquares();
}



/*     FCN(N,X,F,IFAIL) Ext    Function subroutine */
/*       N              Int    Number of vector components (input) */
/*       X(N)           Dble   Vector of unknowns (input) */
/*       F(N)           Dble   Vector of function values (output) */
/*       IFAIL          Int    FCN evaluation-failure indicator. (output) */
/*                             On input:  Has always value 0 (zero). */
/*                             On output: Indicates failure of FCN eval- */
/*                                uation, if having a value <= 2. */
/*                             If <0: NLEQ1 will be terminated with */
/*                                    error code = 82, and IFAIL stored */
/*                                    to IWK(23). */
/*                             If =1: A new trial Newton iterate will */
/*                                    computed, with the damping factor */
/*                                    reduced to it's half. */
/*                             If =2: A new trial Newton iterate will */
/*                                    computed, with the damping factor */
/*                                    reduced by a reduct. factor, which */
/*                                    must be output through F(1) by FCN, */
/*                                    and it's value must be >0 and < 1. */
/*                             Note, that if IFAIL = 1 or 2, additional */
/*                             conditions concerning the damping factor, */
/*                             e.g. the minimum damping factor or the */
/*                             bounded damping strategy may also influ- */
/*                             ence the value of the reduced damping */
/*                             factor. */

void ModelFunction(int* nx, double* y, double* fval, int* pErr)
{
    ExecutableModel* model = callbackModel;
    assert(model && "model is NULL");

    assert(*nx == model->getStateVector(0) && "incorrect state vector size");

    // sets the state vector
    model->setStateVector(y);

    model->getStateVectorRate(0, y, fval);

    if (rr::Logger::getLevel() >= Logger::LOG_DEBUG)
    {
        std::stringstream ss;

        ss << "y: [";
        for (int i = 0; i < *nx; ++i)
        {
            ss << y[i];
            if (i + 1 < *nx)
            {
                ss << ", ";
            }
        }
        ss << "]" << std::endl;

        ss << "dydt: [";
        for (int i = 0; i < *nx; ++i)
        {
            ss << fval[i];
            if (i + 1 < *nx)
            {
                ss << ", ";
            }
        }
        ss << "]" << std::endl;

        Log(Logger::LOG_DEBUG) << ss.str();
    }

    *pErr = 0;
}

void NLEQSolver::setScalingFactors(const vector<double>& sx)
{
    for (int i = 0; i < n; i++)
    {
        XScal[i] = sx[i];
    }
}

int NLEQSolver::getNumberOfNewtonIterations()
{
    return IWK[0];
}

int NLEQSolver::getNumberOfCorrectorSteps()
{
    return IWK[2];
}

int NLEQSolver::getNumberOfModelEvaluations()
{
    return IWK[3];
}

int NLEQSolver::getNumberOfJacobianEvaluations()
{
    return IWK[4];
}

int NLEQSolver::getNumberOfModelEvaluationsForJacobian()
{
    return IWK[7];
}

string ErrorForStatus(int error)
{
    switch (error)
    {
    case 1:     return ("Jacobian matrix singular in NLEQ");
    case 2:     return ("Maximum iterations exceeded");
    case 3:     return ("Damping factor has became to small to continue");
    case 4:     return ("Warning: Superlinear or quadratic convergence slowed down near the solution");
    case 5:     return ("Warning: Error Tolerance reached but solution is suspect");
    case 10:    return ("Integer or real workspace too small in NLEQ");
    case 20:    return ("Bad input to size of model parameter");
    case 21:    return ("Nonpositive value for RTOL supplied to NLEQ");
    case 22:    return ("Negative scaling value via vector XSCAL supplied");
    case 30:    return ("One or more fields specified in IOPT are invalid (NLEQ)");
    case 80:    return ("Error signalled by linear solver routine N1FACT, in NLEQ");
    case 81:    return ("Error signalled by linear solver routine N1SOLV, in NLEQ");
    case 82:    return ("Possible negative concentrations in solution (NLEQ)");
    case 83:    return ("Error signalled by user routine JAC in NLEQ");
    default:    return (format("Unknown error in NLEQ, errCode = {0}", error));
    }
}

double NLEQSolver::computeSumsOfSquares()
{
    double sum = 0;
    vector<double> rates(model->getStateVector(0));
    model->getStateVectorRate(0, 0, &rates[0]);

    for (int i = 0; i < n; i++)
    {
        sum = sum + pow(rates[i], 2.0);
    }
    return sqrt(sum);
}



}    //end of namespace


