#ifndef rrNLEQInterfaceH
#define rrNLEQInterfaceH
#include <vector>
#include "Solver.h"
#include "rrExporter.h"
#include "rrExecutableModel.h"
#include "rrSteadyStateSolver.h"
using std::vector;

namespace rr
{

/**
 * @internal
 */
class RR_DECLSPEC NLEQSolver : public Solver
{

public:
    /**
     * Creates a new Instance of NLEQ for the given Model
     */
    NLEQSolver(ExecutableModel *_model = NULL);
    ~NLEQSolver();

	void loadConfigSettings();
	void loadSBMLSettings(const std::string& filename);

	/**
	* @author WBC
	* @brief Get the name for this Solver
	* @note Delegates to @ref getName
	*/
	std::string getSolverName() const;

	/**
	* @author JKM
	* @brief Get the name for this Solver
	*/
	static std::string getName();

	/**
	* @author WBC
	* @brief Get the description for this Solver
	* @note Delegates to @ref getDescription
	*/
	std::string getSolverDescription() const;

	/**
	* @author JKM
	* @brief Get the description for this Solver
	*/
	static std::string getDescription();

	/**
	* @author WBC
	* @brief Get the hint for this Solver
	* @note Delegates to @ref getHint
	*/
	std::string getSolverHint() const;

	/**
	* @author JKM
	* @brief Get the hint for this Solver
	*/
	static std::string getHint();

	// ** Getters / Setters ************************************************

	/**
	* @author WBC, ETS, MTK
	* @brief Always deterministic for CVODE
	*/
	SolverMethod getSolverMethod() const;

	/**
	* @author WBC, ETS, MTK
	* @brief Sets the value of an Solver setting (e.g. maximum_iterations)
	*/
	void setValue(std::string setting, const Variant& value);


	// ** Solver routines
	double solve(const vector<double>& yin);

public:


private:
	ExecutableModel *model; // Model generated from the SBML. Static so we can access it from standalone function
	
	int nOpts;
    long *IWK;
    long LIWK;
    long LWRK;
    double *RWK;
    double *XScal;
    long ierr;
    long *iopt;
    
    long n;
    void setup();

    bool isAvailable();

    //int maxIterations;
    //double relativeTolerance;
    //double minDamping;


    /// <summary>
    /// Sets the Scaling Factors
    /// </summary>
    /// <param name="sx">Array of Scaling factors</param>
    void                            setScalingFactors(const vector<double>& sx);

    /// <summary>
    /// Returns the Number of Newton Iterations
    /// </summary>
    /// <returns>the Number of Newton Iterations</returns>
    int                             getNumberOfNewtonIterations();

    /// <summary>
    /// Returns the Number of Corrector steps
    /// </summary>
    /// <returns>Returns the Number of Corrector steps</returns>
    int                             getNumberOfCorrectorSteps();

    /// <summary>
    /// Returns the Number of Model Evaluations
    /// </summary>
    /// <returns>the Number of Model Evaluations</returns>
    int                             getNumberOfModelEvaluations();

    /// <summary>
    /// Returns the Number Of Jacobian Evaluations
    /// </summary>
    /// <returns>the Number Of Jacobian Evaluations</returns>
    int                             getNumberOfJacobianEvaluations();

    /// <summary>
    /// Returns the Number of Model Evaluations For Jacobian
    /// </summary>
    /// <returns>the Number of Model Evaluations For Jacobian</returns>
    int                             getNumberOfModelEvaluationsForJacobian();


    double                          computeSumsOfSquares();



};
}

#endif
