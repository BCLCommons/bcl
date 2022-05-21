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
#include "score/bcl_score_protein_model_sse.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSE::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelSSE())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelSSE::ProteinModelSSE() :
      m_ScoreSingle(),
      m_Normalize( false)
    {
    }

    //! @brief construct from a function
    //! @param SP_SINGLE_FUNCTION single function to be used
    //! @param NORMALIZE
    //! @param SCORE_TYPE score type
    //! @param READABLE_SCHEME scheme that is more human readable
    ProteinModelSSE::ProteinModelSSE
    (
      const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &SP_SINGLE_FUNCTION,
      const bool NORMALIZE,
      const ProteinModel::Type &SCORE_TYPE,
      const std::string &READABLE_SCHEME
    ) :
      m_ScoreSingle( *SP_SINGLE_FUNCTION),
      m_Normalize( NORMALIZE),
      m_ScoreType( SCORE_TYPE),
      m_ReadableScheme( READABLE_SCHEME)
    {
      if( m_ReadableScheme.empty())
      {
        m_ReadableScheme = GetScheme();
      }
    }

    //! @brief virtual copy constructor
    ProteinModelSSE *ProteinModelSSE::Clone() const
    {
      return new ProteinModelSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelSSE::GetScheme() const
    {
      return m_ScoreSingle->GetScheme();
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &ProteinModelSSE::GetAlias() const
    {
      static const std::string s_alias( "ProteinModelSSE");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSE::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Wrapping class for scoring functions evaluating individual SSEs.");
      serializer.AddInitializer
      (
        "score function",
        "score function to evaluate the individual SSEs",
        io::Serialization::GetAgent( &m_ScoreSingle)
      );
      serializer.AddInitializer
      (
        "normalize",
        "normalize by the number of evaluated SSEs",
        io::Serialization::GetAgent( &m_Normalize)
      );
      // serializer.AddInitializer
      // (
      //   "type",
      //   "type of the score function",
      //   io::Serialization::GetAgent( &m_ScoreType)
      // );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an PROTEIN_MODEL and returning a t_ResultType object
    //! @param PROTEIN_MODEL Protein Model to be used to evaluate the function
    //! @return function value of the given argument
    double ProteinModelSSE::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      BCL_MessageDbg( "score singles");
      //instantiate the scoresum
      double score( 0.0);

      size_t scored_entities( 0);

      //instantiate a util::SiPtrVector to the elements in the ProteinModel
      const util::SiPtrVector< const assemble::SSE> all_elements( PROTEIN_MODEL.GetSSEs());

      // get any associated membrane
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      //iterate over all elements
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          element_itr( all_elements.Begin()),
          element_itr_end( all_elements.End());
        element_itr != element_itr_end; ++element_itr)
      {
        //call the pairscore
        const storage::Pair< double, size_t> current_score
        (
          m_ScoreSingle->operator()
          (
            **element_itr,
            sp_membrane.IsDefined() ? *sp_membrane : biol::Membrane::GetUndefinedMembrane()
          )
        );
        score += current_score.First();
        scored_entities += current_score.Second();
      }

      // if normalize flag is given and at least 1 entity was scored
      if( m_Normalize && scored_entities != 0)
      {
        score /= double( scored_entities);
      }

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelSSE::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {
      //instantiate a util::SiPtrVector to the elements in the ProteinModel
      const util::SiPtrVector< const assemble::SSE> all_elements( PROTEIN_MODEL.GetSSEs());

      // get any associated membrane
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      //iterate over all possible pairs of elements
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          element_itr( all_elements.Begin()), element_itr_end( all_elements.End());
        element_itr != element_itr_end;
        ++element_itr
      )
      {
        //write the detailed scheme and value for the pair score
        m_ScoreSingle->WriteDetailedSchemeAndValues
        (
          **element_itr,
          sp_membrane.IsDefined() ? *sp_membrane : biol::Membrane::GetUndefinedMembrane(),
          OSTREAM
        );
      }

      //end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
