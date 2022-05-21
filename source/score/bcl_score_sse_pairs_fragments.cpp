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
#include "score/bcl_score_sse_pairs_fragments.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPairsFragments::s_Instance
    (
      util::Enumerated< math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double> >::AddInstance
      (
        new SSEPairsFragments()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPairsFragments::SSEPairsFragments() :
      m_Packer( assemble::GetSSEGeometryPackingListPickers().GetUndefined()),
      m_SPScorePair(),
      m_NormalizeOverNumberOfPairs( false)
    {
    }

    //! @brief construct from a Function Interface and a normalization flag boolean
    //! @param PACKING_PICKER Function that will be used to calculate the packings between fragments
    //! @param PACKING_FUNCTION Function that will be used to calculate scores for all packings
    //! @param NORMALIZE_OVER_NUMBER_OF_PAIRS boolean flag to determine whether to normalize over number of pairs,
    SSEPairsFragments::SSEPairsFragments
    (
      const assemble::SSEGeometryPackingListPicker &PACKING_PICKER,
      const SSEPackInterface &PACKING_FUNCTION,
      const bool NORMALIZE_OVER_NUMBER_OF_PAIRS
    ) :
      m_Packer( PACKING_PICKER),
      m_SPScorePair( PACKING_FUNCTION.Clone()),
      m_NormalizeOverNumberOfPairs( NORMALIZE_OVER_NUMBER_OF_PAIRS)
    {
      // make sure the minimal interface lengths of the packer and the scoring function match each other
      BCL_Assert
      (
        ( *m_Packer)->GetMinimalInterfaceLength() == m_SPScorePair->GetMinimalInterfaceLength(),
        "The minimal interface lengths of the packer and the scoring function do not match " +
        util::Format()( ( *m_Packer)->GetMinimalInterfaceLength()) + " for packer vs " +
        util::Format()( m_SPScorePair->GetMinimalInterfaceLength()) + " for scoring function"
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new SSEPairsFragments instance copied from this one
    SSEPairsFragments *SSEPairsFragments::Clone() const
    {
      return new SSEPairsFragments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPairsFragments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPairsFragments::GetScheme() const
    {
      return m_SPScorePair->GetScheme();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SSEPairsFragments::GetAlias() const
    {
      static const std::string s_name( "SSEPairsFragments");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEPairsFragments::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Evaluates fragment interactions between SSEs.");
      serializer.AddInitializer
      (
        "packer",
        "packer used to evaluate packing",
        io::Serialization::GetAgent( &m_Packer)
      );
      serializer.AddInitializer
      (
        "score function",
        "function used to score the fragment interactions",
        io::Serialization::GetAgent( &m_SPScorePair)
      );
      serializer.AddInitializer
      (
        "normalize",
        "normalize over the number of pairs",
        io::Serialization::GetAgent( &m_NormalizeOverNumberOfPairs)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns whether the given pair of SSEs are valid
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return whether the given pair of SSEs are valid
    bool SSEPairsFragments::AreValidSSEs( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      return m_SPScorePair->AreValidSSEs( SSE_A, SSE_B);
    }

    //! @brief f to sum up all pairwise scores for the given SSE Pair
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return summed up score for the given SSE Pair
    double SSEPairsFragments::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      //instatiate the scoresum
      double score( 0.0);

      // if both SSEs are not helix or strand
      if( !SSE_A.GetType()->IsStructured() || !SSE_B.GetType()->IsStructured())
      {
        // return 0
        return double( 0.0);
      }

      // make sure the SSEs are valid
      if( !AreValidSSEs( SSE_A, SSE_B))
      {
        return double( 0.0);
      }

      BCL_MessageDbg( "scoring fragments with " + m_SPScorePair->GetScheme());
      BCL_MessageDbg
      (
        "SSEs " + SSE_A.GetIdentification() + " and " + SSE_B.GetIdentification()
      );

      // get the list of packings
      const storage::List< assemble::SSEGeometryPacking> packing_list( ( *m_Packer)->operator ()( SSE_A, SSE_B));

      // iterate over the packings
      for
      (
        storage::List< assemble::SSEGeometryPacking>::const_iterator
          pack_itr( packing_list.Begin()), pack_itr_end( packing_list.End());
        pack_itr != pack_itr_end; ++pack_itr
      )
      {

        // call the pairscore with pair of fragments
        const double this_score( m_SPScorePair->operator()( *pack_itr));

        BCL_MessageDbg
        (
          "fragments " +
          pack_itr->GetFirstSSEGeometry()->GetIdentification() + " and " +
          pack_itr->GetSecondSSEGeometry()->GetIdentification() +
          " weight: " + util::Format()( pack_itr->GetInteractionWeight()) +
          " contact " + pack_itr->GetContactType().GetName() +
          " distance " + util::Format()( pack_itr->GetDistance()) +
          " strand pairing w " + util::Format()( pack_itr->GetStrandStrandPairingWeight()) +
          " relative position w " + util::Format()( pack_itr->GetRelativePositionWeight()) +
          " twist " + util::Format()( pack_itr->GetTwistAngle()) +
          " score: " + util::Format()( this_score)
        );

        // sum up the score
        score += this_score;
      }

      // if given SSEs has fragments and normalize flag is set
      if( m_NormalizeOverNumberOfPairs && !packing_list.IsEmpty())
      {
        // then normalize the score by the fragments_a list since it is smaller
        return score / packing_list.GetSize();
      }

      //end
      return score;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &SSEPairsFragments::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Packer, ISTREAM);
      io::Serialize::Read( m_SPScorePair, ISTREAM);
      io::Serialize::Read( m_NormalizeOverNumberOfPairs, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &SSEPairsFragments::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Packer, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SPScorePair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NormalizeOverNumberOfPairs, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @param OSTREAM std::sotream to be written to
    //! @return std::ostream which was written to
    std::ostream &SSEPairsFragments::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B,
      std::ostream &OSTREAM
    ) const
    {

      // if both SSEs are not helix or strand
      if
      (
        ( SSE_A.GetType() != biol::GetSSTypes().HELIX && SSE_A.GetType() != biol::GetSSTypes().STRAND)
        ||
        ( SSE_B.GetType() != biol::GetSSTypes().HELIX && SSE_B.GetType() != biol::GetSSTypes().STRAND)
      )
      {
        // return 0
        return OSTREAM;
      }

      OSTREAM << SSE_A.GetIdentification() << '\t'
              << SSE_B.GetIdentification() << '\t'
              << operator()( SSE_A, SSE_B) << '\n';

      // get the list of packings
      const storage::List< assemble::SSEGeometryPacking> packing_list( ( *m_Packer)->operator ()( SSE_A, SSE_B));

      // iterate over the packings
      for
      (
        storage::List< assemble::SSEGeometryPacking>::const_iterator
          pack_itr( packing_list.Begin()), pack_itr_end( packing_list.End());
        pack_itr != pack_itr_end; ++pack_itr
      )
      {
        OSTREAM << "\t\t"
                << pack_itr->GetFirstSSEGeometry()->GetIdentification() << '\t'
                << pack_itr->GetSecondSSEGeometry()->GetIdentification() << '\t';

        //write the detailed scheme and value for the pair score
        m_SPScorePair->WriteDetailedSchemeAndValues( *pack_itr, OSTREAM);
      }

      //end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
