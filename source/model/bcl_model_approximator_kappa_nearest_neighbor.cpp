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
#include "model/bcl_model_approximator_kappa_nearest_neighbor.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_kappa_nearest_neighbor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ApproximatorKappaNearestNeighbor::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorKappaNearestNeighbor())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ApproximatorKappaNearestNeighbor::ApproximatorKappaNearestNeighbor() :
      ApproximatorBase(),
      m_MinKappa( 1),
      m_Kappa( 1),
      m_MaxKappa( 10),
      m_LastObjectiveFunctionResult( util::GetUndefined< float>())
    {
    }

    //! @brief clone function
    //! @return pointer to new ApproximatorKappaNearestNeighbor
    ApproximatorKappaNearestNeighbor *ApproximatorKappaNearestNeighbor::Clone() const
    {
      return new ApproximatorKappaNearestNeighbor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief sets the training data
    //! @param DATA training data set to be set
    void ApproximatorKappaNearestNeighbor::SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA)
    {
      ApproximatorBase::m_TrainingData = DATA;
      BCL_Assert( !DATA->IsEmpty(), "no training data loaded");

      // rescale the input features
      ApproximatorBase::m_TrainingData->GetFeatures().Rescale( KappaNearestNeighbor::s_DefaultInputRange);

      BCL_Message
      (
        util::Message::e_Standard,
        "Setting up training data with " + util::Format()( ApproximatorBase::m_TrainingData->GetSize()) + " points"
      );
    }

    //! @brief returns the current model
    //! @return current model
    util::ShPtr< Interface> ApproximatorKappaNearestNeighbor::GetCurrentModel() const
    {
      return util::ShPtr< Interface>( new KappaNearestNeighbor( this->m_TrainingData, m_Kappa));
    }

    //! @brief returns the current approximation
    //! @return ShPtr to the current argument result pair
    const util::ShPtr
    <
      storage::Pair< util::ShPtr< Interface>, float>
    > ApproximatorKappaNearestNeighbor::GetCurrentApproximation() const
    {
      // create ShPtr to the model
      util::ShPtr< Interface> sp_model( new KappaNearestNeighbor( this->m_TrainingData, m_Kappa));

      // return ShPtr to model result pair
      return util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      (
        new storage::Pair< util::ShPtr< Interface>, float>
        (
          sp_model,
          m_ObjectiveFunction->operator()( sp_model)
        )
      );
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorKappaNearestNeighbor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorKappaNearestNeighbor::GetAlias() const
    {
      static const std::string s_alias( "KappaNearestNeighbor");
      return s_alias;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorKappaNearestNeighbor::CanContinue() const
    {
      return m_Kappa < m_MaxKappa || !util::IsDefined( m_LastObjectiveFunctionResult);
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorKappaNearestNeighbor::Next()
    {
      // assert when there is no interval
      BCL_Assert( this->m_TrainingData.IsDefined() && !this->m_TrainingData->IsEmpty(), "no training data given!");

      // ShPtr to the model
      util::ShPtr< Interface> model;

      // check for first round
      if( !util::IsDefined( m_LastObjectiveFunctionResult))
      {
        m_Kappa = m_MinKappa;
        model = GetCurrentApproximation()->First();
        m_LastObjectiveFunctionResult = m_ObjectiveFunction->operator()( model);

        BCL_Message
        (
          util::Message::e_Standard,
          "Kappa: " + util::Format()( m_Kappa) + " ObjFunction: "
          + util::Format()( m_LastObjectiveFunctionResult)
        );
      }
      else if( m_Kappa < m_MaxKappa)
      {
        ++m_Kappa;
        model = GetCurrentApproximation()->First();
        m_LastObjectiveFunctionResult = m_ObjectiveFunction->operator()( model);
        BCL_Message
        (
          util::Message::e_Standard,
          "Kappa: " + util::Format()( m_Kappa) + " ObjFunction: "
          + util::Format()( m_LastObjectiveFunctionResult)
        );
      }
      else
      {
        model = GetCurrentApproximation()->First();
      }

      // track the current model
      this->Track
      (
        util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>
          (
            model,
            m_LastObjectiveFunctionResult
          )
        )
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorKappaNearestNeighbor::Read( std::istream &ISTREAM)
    {
      // call Read of base class
      ApproximatorBase::Read( ISTREAM);

      // read members
      io::Serialize::Read( m_MinKappa, ISTREAM);
      io::Serialize::Read( m_Kappa, ISTREAM);
      io::Serialize::Read( m_MaxKappa, ISTREAM);
      io::Serialize::Read( m_LastObjectiveFunctionResult, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorKappaNearestNeighbor::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // call Write of base class
      ApproximatorBase::Write( OSTREAM, INDENT) << '\n';

      // write members
      io::Serialize::Write( m_MinKappa, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Kappa, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxKappa, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LastObjectiveFunctionResult, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorKappaNearestNeighbor::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "A k-nearest-neighbor predictor; iteration optimizes k. see http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm"
      );

      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each batch step",
        io::Serialization::GetAgent( &ApproximatorBase::m_ObjectiveFunction->GetImplementation()),
        "RMSD"
      );
      parameters.AddInitializer
      (
        "min kappa",
        "minimum # of nearest neighbors for kNN selection",
        io::Serialization::GetAgentWithMin( &m_MinKappa, size_t( 1)),
        "1"
      );
      parameters.AddInitializer
      (
        "max kappa",
        "maximum # of nearest neighbors for kNN selection",
        io::Serialization::GetAgentWithMin( &m_MaxKappa, size_t( 1)),
        "10"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
