// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opti/bcl_opti.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace opti
} // namespace bcl
// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opti/bcl_opti_ensemble_filter.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> EnsembleFilter::s_Instance
    (
      util::Enumerated< OptimizationInterface< assemble::Ensemble< assemble::ProteinModel> > >::AddInstance( new EnsembleFilter())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EnsembleFilter::EnsembleFilter() :
      m_ScoreFunction(),
      m_KeepPercentage()
    {
    }

    //! @brief construct from members
    //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
    //! @param KEEP percentage of model to keep
    EnsembleFilter::EnsembleFilter( const score::ProteinModel &SCORE_FUNCTION, double KEEP) :
      m_ScoreFunction( SCORE_FUNCTION),
      m_KeepPercentage( KEEP)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new EnsembleFilter
    EnsembleFilter *EnsembleFilter::Clone() const
    {
      return new EnsembleFilter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &EnsembleFilter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &EnsembleFilter::GetAlias() const
    {
      static const std::string s_alias( "EnsembleFilter");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EnsembleFilter::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Filters an ensemble using a provided scoring function.");
      serializer.AddInitializer
      (
        "score function",
        "score function to evaluate the sampled protein models",
        io::Serialization::GetAgent( &m_ScoreFunction)
      );
      serializer.AddInitializer
      (
        "keep",
        "percentage of models to keep",
        io::Serialization::GetAgent( &m_KeepPercentage)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief filters the provided ensemble
    //! @param ENSEMBLE ensemble to be filtered
    void EnsembleFilter::Optimize( assemble::Ensemble< assemble::ProteinModel> &ENSEMBLE) const
    {
      // compute the score for each element of the ensemble
      typedef assemble::Ensemble< assemble::ProteinModel>::Element Element;
      storage::Vector< storage::Pair< double, Element> > elements;
      storage::Vector< double> scores;
      for( auto el_it( ENSEMBLE.Begin()), el_it_end( ENSEMBLE.End()); el_it != el_it_end; ++el_it)
      {
        const assemble::ProteinModel &model( el_it->GetElement());
        const double score( ( *m_ScoreFunction)( model));
        const storage::Pair< double, Element> element_score( score, *el_it);
        elements.PushBack( element_score);
        scores.PushBack( score);
      }

      // sort the scores and determine the cutoff
      scores.Sort( std::less< double>());
      const size_t number_models( ENSEMBLE.GetSize());
      const size_t number_keep( m_KeepPercentage * number_models);
      const double cutoff( scores( number_keep - 1));

      // only keep elements with a score below the cutoff value
      assemble::Ensemble< assemble::ProteinModel> ensemble;
      for( auto el_it( elements.Begin()), el_it_end( elements.End()); el_it != el_it_end; ++el_it)
      {
        const double score( el_it->First());
        if( score <= cutoff)
        {
          ensemble.AddElement( el_it->Second());
        }
      }
      ENSEMBLE = ensemble;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opti
} // namespace bcl
// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math.h"
#include "opti/bcl_opti_improvement_type.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    //! @brief returns ImprovementType as string
    //! @param TYPE the improvement type
    //! @return the string for the improvement type
    const std::string &GetImprovementTypeName( const ImprovementType &TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "SmallerIsBetter",         //!< smaller values indicate improvement
        "LargerIsBetter",          //!< larger values indicate improvement
        "SmallerAbsIsBetter",      //!< smaller absolute values indicate improvement (typically error)
        "SmallerEqualIsBetter",    //!< smaller or equal values indicate improvement
        "LargerEqualIsBetter",     //!< larger or equal values indicate improvement
        GetStaticClassName< ImprovementType>()
      };

      return s_descriptors[ TYPE];
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    template< typename t_DataType>
    bool TestDoesImprove( const t_DataType &TEST, const t_DataType &PRIOR, const ImprovementType &TYPE)
    {
      bool improved( false);
      switch( TYPE)
      {
        case e_SmallerIsBetter:      improved = TEST <  PRIOR; break;
        case e_LargerIsBetter:       improved = TEST >  PRIOR; break;
        case e_SmallerEqualIsBetter: improved = TEST <= PRIOR; break;
        case e_LargerEqualIsBetter:  improved = TEST >= PRIOR; break;
        case e_SmallerAbsIsBetter:   improved = math::Absolute( TEST) < math::Absolute( PRIOR); break;
        default: break;
      }
      return improved;
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const double &TEST, const double &PRIOR, const ImprovementType &TYPE)
    {
      return TestDoesImprove( TEST, PRIOR, TYPE);
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const float &TEST, const float &PRIOR, const ImprovementType &TYPE)
    {
      return TestDoesImprove( TEST, PRIOR, TYPE);
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const int &TEST, const int &PRIOR, const ImprovementType &TYPE)
    {
      return TestDoesImprove( TEST, PRIOR, TYPE);
    }

  } // namespace opti
} // namespace bcl
// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opti/bcl_opti_phase.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    //! @brief returns Phase as string
    //! @param TYPE the improvement type
    //! @return the string for the improvement type
    const std::string &GetPhaseName( const Phase &TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Start",
        "Iteration",
        "End",
        "Always",
        GetStaticClassName< Phase>()
      };

      return s_descriptors[ TYPE];
    }

    //! @brief test whether one phase logically equals another
    //! @param A, B the phases to compare
    //! @return true if A == B
    bool PhasesEqual( const Phase &A, const Phase &B)
    {
      return A == B || A == e_Always || B == e_Always;
    }

  } // namespace opti
} // namespace bcl
// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opti/bcl_opti_step_status.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    //! @brief conversion to a string from a StepStatus
    //! @param STEP_STATUS the step status to get a string for
    //! @return a string representing that step status
    const std::string &GetStepStatusName( const StepStatus &STEP_STATUS)
    {
      static const std::string s_descriptors[] =
      {
        "improved",
        "accepted",
        "rejected",
        "skipped",
        GetStaticClassName< StepStatus>()
      };
      return s_descriptors[ size_t( STEP_STATUS)];
    }

  } // namespace opti
} // namespace bcl
// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opti/bcl_opti_template_instantiations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    template class BCL_API CriterionCombine< assemble::ProteinModel, double>;

    template class BCL_API CriterionNumberIterations< assemble::ProteinModel, double>;

    template class BCL_API CriterionUnimproved< assemble::ProteinModel, double>;

    template class BCL_API EnsembleNode< assemble::ProteinModel>;

    template class BCL_API OptimizationIdentity< assemble::ProteinModel>;

    template class BCL_API OptimizationIdentity< assemble::Ensemble< assemble::ProteinModel> >;

    template class BCL_API Pipeline< assemble::ProteinModel>;

  } // namespace opti
} // namespace bcl
// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opti/bcl_opti_tracker_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TrackerBase::TrackerBase() :
      m_IterationNumber( 0),
      m_BestModelIterationNumber( 0),
      m_ImprovementType( e_SmallerIsBetter),
      m_Phase( e_Start),
      m_Counts( s_NumberStepStatus, size_t( 0)),
      m_IntervalStart( s_NumberStepStatus, size_t( 0)),
      m_IntervalEnd( s_NumberStepStatus, size_t( 0)),
      m_StatusOfLastStep()
    {
    }

    //! @brief construct from ImprovementType
    //! @param IMPROVEMENT_TYPE change type that indicates an improvement
    TrackerBase::TrackerBase( const ImprovementType &IMPROVEMENT_TYPE) :
      m_IterationNumber( 0),
      m_BestModelIterationNumber( 0),
      m_ImprovementType( IMPROVEMENT_TYPE),
      m_Phase( e_Start),
      m_Counts( s_NumberStepStatus, size_t( 0)),
      m_IntervalStart( s_NumberStepStatus, size_t( 0)),
      m_IntervalEnd( s_NumberStepStatus, size_t( 0)),
      m_StatusOfLastStep()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the number of iterations since the last improvement
    //! @return number of iterations since the last improvement
    size_t TrackerBase::GetIterationsSinceLastImprovement() const
    {
      return m_IterationNumber - m_BestModelIterationNumber;
    }

    //! @brief resets the tracker
    void TrackerBase::Reset()
    {
      // reset members
      m_IterationNumber = 0;
      m_BestModelIterationNumber = 0;
      m_Phase = e_Start;
      m_Counts.SetAllElements( 0);
      m_IntervalStart.SetAllElements( 0);
      m_IntervalEnd.SetAllElements( 0);
      m_StatusOfLastStep = e_Rejected;
    }

    //! @brief sets the phase of the iteration
    //! @PHASE the now current phase of the iteration
    void TrackerBase::SetPhase( const Phase &PHASE)
    {
      m_Phase = PHASE;
    }

    //! @brief return number of steps in a row of specified STEP_STATUS
    //! @param STEP_STATUS StepStatus of interest
    //! @return number of steps in a row of specified STEP_STATUS
    size_t TrackerBase::GetNumberStepsInARow( const StepStatus STEP_STATUS) const
    {
      // if this step status was observed in the last step
      if( m_IntervalEnd( STEP_STATUS) == m_IterationNumber)
      {
        return m_IterationNumber - m_IntervalStart( STEP_STATUS);
      }

      // otherwise it means some other step status type was observed last, therefore return 0
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TrackerBase::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_IterationNumber, ISTREAM);
      io::Serialize::Read( m_BestModelIterationNumber, ISTREAM);
      io::Serialize::Read( m_ImprovementType, ISTREAM);
      io::Serialize::Read( m_Phase, ISTREAM);
      io::Serialize::Read( m_Counts, ISTREAM);
      io::Serialize::Read( m_IntervalStart, ISTREAM);
      io::Serialize::Read( m_IntervalEnd, ISTREAM);
      io::Serialize::Read( m_StatusOfLastStep, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &TrackerBase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write memebers
      io::Serialize::Write( m_IterationNumber, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BestModelIterationNumber, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ImprovementType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Phase, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Counts, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IntervalStart, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IntervalEnd, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StatusOfLastStep, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief Update the status of this tracker
    void TrackerBase::Update( const StepStatus STEP_STATUS)
    {
      // increment counts for the provided StepStatus
      ++m_Counts( STEP_STATUS);

      if( STEP_STATUS == e_Skipped)
      {
        // otherwise ignore skipped steps
        return;
      }

      // if previous iteration was not of this step status
      // in which case the end of interval would be not equal to last iteration number
      if( m_IntervalEnd( STEP_STATUS) < m_IterationNumber - 1)
      {
        // update the begin of the interval to this iteration
        m_IntervalStart( STEP_STATUS) = m_IterationNumber;
      }

      // update the end of interval to current iteration
      m_IntervalEnd( STEP_STATUS) = m_IterationNumber;

      // update last step
      m_StatusOfLastStep = STEP_STATUS;

      ++m_IterationNumber;

      if( STEP_STATUS == e_Improved)
      {
        m_BestModelIterationNumber = m_IterationNumber;
      }
    }

  } // namespace opti

} // namespace bcl
