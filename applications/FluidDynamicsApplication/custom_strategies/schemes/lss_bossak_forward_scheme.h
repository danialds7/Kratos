//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_LSS_BOSSAK_FORWARD_SCHEME_H_INCLUDED)
#define KRATOS_LSS_BOSSAK_FORWARD_SCHEME_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"
#include "utilities/time_discretization.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/fluid_adjoint_slip_utilities.h"
#include "custom_utilities/fluid_lss_variable_utilities.h"
#include "custom_utilities/fluid_lss_sensitivity.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class LSSBossakForwardScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LSSBossakForwardScheme);

    using IndexType = std::size_t;

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using EquationIdVectorType = std::vector<IndexType>;

    using AssembledVectorType = std::unordered_map<IndexType, double>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit LSSBossakForwardScheme(
        AdjointResponseFunction::Pointer pResponseFunction,
        FluidLSSSensitivity::Pointer pFluidLeastSquaresShadowingSensitivity,
        FluidLSSVariableUtilities::Pointer pFluidLeastSquaresShadowingUtilities,
        const Variable<double>& rResponseDesignTotalDerivativeVariable,
        const Variable<double>& rDeltaTimeDesignTotalDerivativeVariable,
        const double BossakAlpha,
        const double DeltaTimeDialationAlpha,
        const double FinalResponseValue,
        const IndexType Dimension,
        const IndexType BlockSize,
        const IndexType EchoLevel)
        : BaseType(),
          mpResponseFunction(pResponseFunction),
          mpFluidLeastSquaresShadowingSensitivity(pFluidLeastSquaresShadowingSensitivity),
          mpFluidLeastSquaresShadowingUtilities(pFluidLeastSquaresShadowingUtilities),
          mrResponseDesignTotalDerivativeVariable(rResponseDesignTotalDerivativeVariable),
          mrDeltaTimeDesignTotalDerivativeVariable(rDeltaTimeDesignTotalDerivativeVariable),
          mDeltaTimeDialationAlpha(DeltaTimeDialationAlpha),
          mFinalResponseValue(FinalResponseValue),
          mDimension(Dimension),
          mEchoLevel(EchoLevel),
          mBossak(BossakAlpha),
          mBossakConstants(mBossak),
          mAdjointSlipUtilities(Dimension, BlockSize)
    {
        KRATOS_TRY

        // Allocate auxiliary memory.
        const int number_of_threads = ParallelUtilities::GetNumThreads();
        mTLS.resize(number_of_threads);

        KRATOS_ERROR_IF_NOT(mpFluidLeastSquaresShadowingUtilities->GetPrimalVariablePointersList().size() == BlockSize)
                << "Provided block size does not match with the number of primal variables provided in least squares shadowing utilities. [ block size = "
                << BlockSize << ", number of primal variables in least squares shadowing utilities = " << mpFluidLeastSquaresShadowingUtilities->GetPrimalVariablePointersList().size() << " ].\n";

        KRATOS_ERROR_IF(mrResponseDesignTotalDerivativeVariable == mrDeltaTimeDesignTotalDerivativeVariable)
                << "Response design total derivative variable is same as the delta time design total derivative variable. [ variable = " << mrDeltaTimeDesignTotalDerivativeVariable.Name() << " ].\n";

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Created [ Dimensionality = " << Dimension << ", BlockSize = " << BlockSize << " ].\n";

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~LSSBossakForwardScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        IndexPartition<IndexType>(rModelPart.NumberOfElements()).for_each([&](const IndexType iElement) {
            const auto& r_element = *(rModelPart.ElementsBegin() + iElement);
            mpFluidLeastSquaresShadowingUtilities->CheckVariables(r_element);
        });

        const IndexType buffer_size = rModelPart.GetBufferSize();
        KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size needs to be greater than 1 in " << rModelPart.FullName() << " [ buffer size = " << buffer_size << " ].\n";

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        KRATOS_ERROR_IF(delta_time < 0.0) << "LSSBossakForwardScheme should be run in forward time with positive delta time. [ delta_time = " << delta_time << " ].\n";

        return BaseType::Check(rModelPart);

        KRATOS_CATCH("");
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BaseType::Initialize(rModelPart);

        mAdjointSlipUtilities.Initialize(rModelPart);

        rModelPart.GetProcessInfo()[mrResponseDesignTotalDerivativeVariable] = 0.0;

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        // Calculate the delta time total design derivative
        struct TLS {
            Matrix mRotatedResidualPartialTimeDerivative;
            Matrix mResidualPartialTimeDerivative;
            Vector mAdjointSolution;
        };

        auto& r_process_info = rModelPart.GetProcessInfo();

        const double elemental_eta = block_for_each<SumReduction<double>>(rModelPart.Elements(), TLS(), [&](ModelPart::ElementType& rElement, TLS& rTLS) -> double {
            mpFluidLeastSquaresShadowingUtilities->GetAdjointValues(rTLS.mAdjointSolution, rElement, 0);
            rElement.CalculateSensitivityMatrix(TIME_STEP_SENSITIVITY, rTLS.mResidualPartialTimeDerivative, r_process_info);

            mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
                rTLS.mRotatedResidualPartialTimeDerivative, rTLS.mResidualPartialTimeDerivative, rElement.GetGeometry());

            return inner_prod(rTLS.mAdjointSolution, row(rTLS.mRotatedResidualPartialTimeDerivative, 0));
        });

        const double delta_time = r_process_info[DELTA_TIME];
        const double coeff = delta_time * delta_time * mBossak.GetGamma();
        const auto& r_primal_variable_pointers_list = mpFluidLeastSquaresShadowingUtilities->GetPrimalVariablePointersList();
        const auto& r_adjoint_first_derivative_variable_pointers_list = mpFluidLeastSquaresShadowingUtilities->GetAdjointFirstDerivativeVariablePointersList();
        const IndexType number_of_variables = r_primal_variable_pointers_list.size();

        const double nodal_eta = block_for_each<SumReduction<double>>(rModelPart.GetCommunicator().LocalMesh().Nodes(), [&](ModelPart::NodeType& rNode) -> double {
            double value = 0.0;
            for (IndexType i = 0; i < number_of_variables; ++i) {
                value += (rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i], 0) - rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i], 1)) * rNode.FastGetSolutionStepValue(*r_adjoint_first_derivative_variable_pointers_list[i]);
            }
            return value;
        });

        const double eta = elemental_eta + nodal_eta / coeff;
        r_process_info[mrDeltaTimeDesignTotalDerivativeVariable] = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(eta) / (-2.0 * mDeltaTimeDialationAlpha * mDeltaTimeDialationAlpha);

        // set bossak constants
        mBossakConstants.SetConstants(delta_time);

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Computed time dialation total design derivative is " << r_process_info[mrDeltaTimeDesignTotalDerivativeVariable] << " in " << rModelPart.FullName()
            << " and stored in " << mrDeltaTimeDesignTotalDerivativeVariable.Name() << ".\n";

        mCurrentResponseValue = mpResponseFunction->CalculateValue(rModelPart);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentElement, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];

        CalculateEntityLHSContribution(
            rLHS_Contribution, rEquationId, rCurrentElement, r_tls.mResidualFirstDerivatives, r_tls.mRotatedResidualFirstDerivatives,
            r_tls.mResidualSecondDerivatives, r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentCondition, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];

        CalculateEntityLHSContribution(
            rLHS_Contribution, rEquationId, rCurrentCondition, r_tls.mResidualFirstDerivatives, r_tls.mRotatedResidualFirstDerivatives,
            r_tls.mResidualSecondDerivatives, r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);

    }

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BaseType::FinalizeSolutionStep(rModelPart, A, Dx, b);

        // update the vdot
        auto& r_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_process_info[DELTA_TIME];
        const double coeff_1 = 1 / (delta_time * mBossak.GetGamma());
        const double coeff_2 = (mBossak.GetGamma() - 1) / mBossak.GetGamma();
        const double coeff_3 = r_process_info[mrDeltaTimeDesignTotalDerivativeVariable] / (delta_time * delta_time * mBossak.GetGamma());

        const auto& r_primal_variable_pointers_list = mpFluidLeastSquaresShadowingUtilities->GetPrimalVariablePointersList();
        const auto& r_lss_variable_pointers_list = mpFluidLeastSquaresShadowingUtilities->GetLSSVariablePointersList();
        const auto& r_lss_first_derivative_variable_pointers_list = mpFluidLeastSquaresShadowingUtilities->GetLSSFirstDerivativeVariablePointersList();
        const IndexType number_of_variables = r_primal_variable_pointers_list.size();

        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            for (IndexType i = 0; i < number_of_variables; ++i){
                double value = 0.0;

                value += coeff_1 * rNode.FastGetSolutionStepValue(*r_lss_variable_pointers_list[i]);
                value -= coeff_1 * rNode.FastGetSolutionStepValue(*r_lss_variable_pointers_list[i], 1);
                value += coeff_2 * rNode.FastGetSolutionStepValue(*r_lss_first_derivative_variable_pointers_list[i], 1);
                value -= coeff_3 * rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i]);
                value += coeff_3 * rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i], 1);

                rNode.FastGetSolutionStepValue(*r_lss_first_derivative_variable_pointers_list[i]) = value;
            }
        });

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Computed primal first derivative total design derivatives in " << rModelPart.FullName() << ".\n";

        struct TLS
        {
            Matrix mResidualFirstDerivatives;
            Matrix mResidualSecondDerivatives;
            Matrix mResidualTimeStepDerivatives;
            Vector mResponseFirstDerivatives;
            Vector mResponseSecondDerivatives;
            Vector mResponseTimeStepDerivatives;
            Vector mCurrentLSSSolution;
            Vector mCurrentLSSSecondDerivativeSolution;
            Vector mPreviousLSSSecondDerivativeSolution;
        };

        const double elemental_design_deriv_contribution = block_for_each<SumReduction<double>>(rModelPart.Elements(), TLS(), [&](ModelPart::ElementType& rElement, TLS& rTLS) -> double {
            return this->CalculateEntityResponseFunctionTotalDerivativeContributions(
                            rElement,
                            rTLS.mResidualFirstDerivatives,
                            rTLS.mResidualSecondDerivatives,
                            rTLS.mResidualTimeStepDerivatives,
                            rTLS.mResponseFirstDerivatives,
                            rTLS.mResponseSecondDerivatives,
                            rTLS.mResponseTimeStepDerivatives,
                            rTLS.mCurrentLSSSolution,
                            rTLS.mCurrentLSSSecondDerivativeSolution,
                            rTLS.mPreviousLSSSecondDerivativeSolution,
                            r_process_info);
        });

        const double condition_design_deriv_contribution = block_for_each<SumReduction<double>>(rModelPart.Conditions(), TLS(), [&](ModelPart::ConditionType& rCondition, TLS& rTLS) -> double {
            return this->CalculateEntityResponseFunctionTotalDerivativeContributions(
                            rCondition,
                            rTLS.mResidualFirstDerivatives,
                            rTLS.mResidualSecondDerivatives,
                            rTLS.mResidualTimeStepDerivatives,
                            rTLS.mResponseFirstDerivatives,
                            rTLS.mResponseSecondDerivatives,
                            rTLS.mResponseTimeStepDerivatives,
                            rTLS.mCurrentLSSSolution,
                            rTLS.mCurrentLSSSecondDerivativeSolution,
                            rTLS.mPreviousLSSSecondDerivativeSolution,
                            r_process_info);
        });

        const double current_response_design_total_derivative = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(elemental_design_deriv_contribution + condition_design_deriv_contribution)
                                                          + r_process_info[mrDeltaTimeDesignTotalDerivativeVariable] * (mCurrentResponseValue - mFinalResponseValue);

        auto& response_design_total_derivative = r_process_info[mrResponseDesignTotalDerivativeVariable];
        const double time = r_process_info[TIME];
        response_design_total_derivative = (response_design_total_derivative * (time - delta_time) + current_response_design_total_derivative) / time;

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Computed response function total design derivative is "
                << response_design_total_derivative << " in " << rModelPart.FullName() << " and stored in " << mrResponseDesignTotalDerivativeVariable.Name() <<".\n";

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "LSSBossakForwardScheme";
    }

    ///@}

private:
    ///@name Private Classes
    ///@{

    struct BossakConstants
    {
        double mC1;
        double mC2;
        double mC3;
        double mCurrentStepSecondDerivativeCoefficient;
        double mPreviousStepSecondDerivativeCoefficient;
        const TimeDiscretization::Bossak& mrBossak;

        explicit BossakConstants(const TimeDiscretization::Bossak& rBossak) : mrBossak(rBossak)
        {
            mCurrentStepSecondDerivativeCoefficient = 1 - mrBossak.GetAlphaM();
            mPreviousStepSecondDerivativeCoefficient = mrBossak.GetAlphaM();
        }

        void SetConstants(const double DeltaTime)
        {
            mC1 = 1 / (DeltaTime * mrBossak.GetGamma());
            mC2 = (mrBossak.GetGamma() - 1) / mrBossak.GetGamma();
            mC3 = 1 / (DeltaTime * DeltaTime * mrBossak.GetGamma());
        }
    };

    struct TLS
    {
        Vector mLSSValues;

        // LHS TLS
        Matrix mResidualFirstDerivatives;
        Matrix mRotatedResidualFirstDerivatives;
        Matrix mResidualSecondDerivatives;
        Matrix mRotatedResidualSecondDerivatives;

        // RHS TLS
        Vector mCurrentPrimalSolution;
        Vector mPreviousPrimalSolution;
        Vector mPreviousLSSSolution;
        Vector mPreviousLSSSecondDerivativeSolution;
        Matrix mResidualTimeStepDerivatives;
        Matrix mRotatedResidualTimeStepDerivatives;
        Vector mRotatedResidualDesignDerivatives;
    };

    ///@}
    ///@name Private Member Variables
    ///@{

    AdjointResponseFunction::Pointer mpResponseFunction;

    FluidLSSSensitivity::Pointer mpFluidLeastSquaresShadowingSensitivity;

    FluidLSSVariableUtilities::Pointer mpFluidLeastSquaresShadowingUtilities;

    const Variable<double>& mrResponseDesignTotalDerivativeVariable;

    const Variable<double>& mrDeltaTimeDesignTotalDerivativeVariable;

    const double mDeltaTimeDialationAlpha;

    const double mFinalResponseValue;

    const IndexType mDimension;

    const IndexType mEchoLevel;

    const TimeDiscretization::Bossak mBossak;

    BossakConstants mBossakConstants;

    double mCurrentResponseValue;

    FluidAdjointSlipUtilities mAdjointSlipUtilities;

    std::vector<TLS> mTLS;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void CalculateEntitySystemContributions(
        TEntityType& rEntity,
        LocalSystemMatrixType& rLHS,
        LocalSystemVectorType& rRHS,
        typename TEntityType::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];

        CalculateEntityLHSContribution<TEntityType>(
            rLHS, rEquationId, rEntity, r_tls.mResidualFirstDerivatives, r_tls.mRotatedResidualFirstDerivatives,
            r_tls.mResidualSecondDerivatives, r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);

        CalculateEntityRHSContribution<TEntityType>(
            rRHS, rEntity, r_tls.mCurrentPrimalSolution, r_tls.mPreviousPrimalSolution, r_tls.mPreviousLSSSolution,
            r_tls.mPreviousLSSSecondDerivativeSolution, r_tls.mResidualTimeStepDerivatives,
            r_tls.mRotatedResidualTimeStepDerivatives, r_tls.mRotatedResidualDesignDerivatives,
            r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);

        // Calculate system contributions in residual form.
        if (rLHS.size1() != 0) {
            rEntity.GetValuesVector(r_tls.mLSSValues);
            noalias(rRHS) -= prod(rLHS, r_tls.mLSSValues);
        }

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityLHSContribution(
        Matrix& rLHS,
        typename TEntityType::EquationIdVectorType& rEquationId,
        TEntityType& rEntity,
        Matrix& rEntityResidualFirstDerivatives,
        Matrix& rEntityRotatedResidualFirstDerivatives,
        Matrix& rEntityResidualSecondDerivatives,
        Matrix& rEntityRotatedResidualSecondDerivatives,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rEntity.EquationIdVector(rEquationId, rCurrentProcessInfo);

        rEntity.CalculateFirstDerivativesLHS(rEntityResidualFirstDerivatives, rCurrentProcessInfo);
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            rEntityRotatedResidualFirstDerivatives, rEntityResidualFirstDerivatives, rEntity.GetGeometry());

        rEntityRotatedResidualFirstDerivatives = trans(rEntityRotatedResidualFirstDerivatives);

        rEntity.CalculateSecondDerivativesLHS(rEntityResidualSecondDerivatives, rCurrentProcessInfo);
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rEntityRotatedResidualSecondDerivatives, rEntityResidualSecondDerivatives, rEntity.GetGeometry());
        rEntityRotatedResidualSecondDerivatives = trans(rEntityRotatedResidualSecondDerivatives);

        if (rLHS.size1() != rEntityRotatedResidualFirstDerivatives.size1() || rLHS.size2() != rEntityRotatedResidualFirstDerivatives.size2()) {
            rLHS.resize(rEntityRotatedResidualFirstDerivatives.size1(), rEntityRotatedResidualFirstDerivatives.size2(), false);
        }

        noalias(rLHS) = rEntityRotatedResidualFirstDerivatives + rEntityRotatedResidualSecondDerivatives * mBossakConstants.mC1 * mBossakConstants.mCurrentStepSecondDerivativeCoefficient;

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityRHSContribution(
        Vector& rRHS,
        TEntityType& rEntity,
        Vector& rCurrentPrimalSolution,
        Vector& rPreviousPrimalSolution,
        Vector& rPreviousLSSSolution,
        Vector& rPreviousLSSSecondDerivativeSolution,
        Matrix& rResidualTimeStepDerivatives,
        Matrix& rRotatedResidualTimeStepDerivatives,
        Vector& rRotatedResidualDesignDerivatives,
        const Matrix& rRotatedResidualSecondDerivatives,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const double eta = rCurrentProcessInfo[mrDeltaTimeDesignTotalDerivativeVariable];

        rEntity.CalculateSensitivityMatrix(TIME_STEP_SENSITIVITY, rResidualTimeStepDerivatives, rCurrentProcessInfo);
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rRotatedResidualTimeStepDerivatives, rResidualTimeStepDerivatives, rEntity.GetGeometry());

        rEntity.GetValuesVector(rPreviousLSSSolution, 1);
        rEntity.GetSecondDerivativesVector(rPreviousLSSSecondDerivativeSolution, 1);
        mpFluidLeastSquaresShadowingUtilities->GetPrimalValues(rCurrentPrimalSolution, rEntity, 0);
        mpFluidLeastSquaresShadowingUtilities->GetPrimalValues(rPreviousPrimalSolution, rEntity, 1);

        mpFluidLeastSquaresShadowingSensitivity->CalculateResidualSensitivity(rRotatedResidualDesignDerivatives,
            rEntity, mAdjointSlipUtilities, rCurrentProcessInfo);

        if (rRHS.size() != rCurrentPrimalSolution.size()) {
            rRHS.resize(rCurrentPrimalSolution.size());
        }

        noalias(rRHS)  = prod(rRotatedResidualSecondDerivatives, rPreviousLSSSolution) * (mBossakConstants.mC1 * mBossakConstants.mCurrentStepSecondDerivativeCoefficient);
        noalias(rRHS) -= prod(rRotatedResidualSecondDerivatives, rPreviousLSSSecondDerivativeSolution) * (
                             mBossakConstants.mC2 * mBossakConstants.mCurrentStepSecondDerivativeCoefficient
                            + mBossakConstants.mPreviousStepSecondDerivativeCoefficient);
        noalias(rRHS) += prod(rRotatedResidualSecondDerivatives, rCurrentPrimalSolution - rPreviousPrimalSolution) * (mBossakConstants.mC3 * mBossakConstants.mCurrentStepSecondDerivativeCoefficient * eta);
        noalias(rRHS) -= row(rRotatedResidualTimeStepDerivatives, 0) * eta;
        noalias(rRHS) -= rRotatedResidualDesignDerivatives;

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    double CalculateEntityResponseFunctionTotalDerivativeContributions(
        TEntityType& rEntity,
        Matrix& rResidualFirstDerivatives,
        Matrix& rResidualSecondDerivatives,
        Matrix& rResidualTimeStepDerivatives,
        Vector& rResponseFirstDerivatives,
        Vector& rResponseSecondDerivatives,
        Vector& rResponseTimeStepDerivatives,
        Vector& rCurrentLSSSolution,
        Vector& rCurrentLSSSecondDerivativeSolution,
        Vector& rPreviousLSSSecondDerivativeSolution,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        const double eta = rProcessInfo[mrDeltaTimeDesignTotalDerivativeVariable];

        double value = 0.0;

        rEntity.CalculateFirstDerivativesLHS(rResidualFirstDerivatives, rProcessInfo);
        mpResponseFunction->CalculateFirstDerivativesGradient(rEntity, rResidualFirstDerivatives, rResponseFirstDerivatives, rProcessInfo);
        mpFluidLeastSquaresShadowingUtilities->GetLSSValues(rCurrentLSSSolution, rEntity);
        value += inner_prod(rResponseFirstDerivatives, rCurrentLSSSolution);

        rEntity.CalculateSecondDerivativesLHS(rResidualSecondDerivatives, rProcessInfo);
        mpResponseFunction->CalculateSecondDerivativesGradient(rEntity, rResidualSecondDerivatives, rResponseSecondDerivatives, rProcessInfo);
        mpFluidLeastSquaresShadowingUtilities->GetLSSFirstDerivativeValues(rCurrentLSSSecondDerivativeSolution, rEntity);
        mpFluidLeastSquaresShadowingUtilities->GetLSSFirstDerivativeValues(rPreviousLSSSecondDerivativeSolution, rEntity, 1);
        value += inner_prod(rResponseSecondDerivatives, rCurrentLSSSecondDerivativeSolution) * mBossakConstants.mCurrentStepSecondDerivativeCoefficient;
        value += inner_prod(rResponseSecondDerivatives, rPreviousLSSSecondDerivativeSolution) * mBossakConstants.mPreviousStepSecondDerivativeCoefficient;

        rEntity.CalculateSensitivityMatrix(TIME_STEP_SENSITIVITY, rResidualTimeStepDerivatives, rProcessInfo);
        mpResponseFunction->CalculatePartialSensitivity(rEntity, SHAPE_SENSITIVITY, rResidualTimeStepDerivatives, rResponseTimeStepDerivatives, rProcessInfo);
        value += rResponseTimeStepDerivatives[0] * eta;

        value += mpFluidLeastSquaresShadowingSensitivity->CalculateResponseSensitivity(rEntity, *mpResponseFunction, mAdjointSlipUtilities, rProcessInfo);

        return value * rProcessInfo[DELTA_TIME];

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class LSSBossakForwardScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_LSS_BOSSAK_FORWARD_SCHEME_H_INCLUDED defined */