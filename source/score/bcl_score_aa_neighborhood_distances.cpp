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
#include "score/bcl_score_aa_neighborhood_distances.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AANeighborhoodDistances::s_Instance
    (
      util::Enumerated< AANeighborhoodInterface>::AddInstance( new AANeighborhoodDistances())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AANeighborhoodDistances::AANeighborhoodDistances() :
      m_ScoreDistance()
    {
    }

    //! @brief constructor from a distance score
    //! @param SP_AA_PAIR_DISTANCE_SCORE aa pair distance score
    AANeighborhoodDistances::AANeighborhoodDistances
    (
      const AAPairDistanceInterface &SP_AA_PAIR_DISTANCE_SCORE
    ) :
      m_ScoreDistance( SP_AA_PAIR_DISTANCE_SCORE)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AANeighborhoodDistances object copied from this one
    AANeighborhoodDistances *AANeighborhoodDistances::Clone() const
    {
      return new AANeighborhoodDistances( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborhoodDistances::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AANeighborhoodDistances::GetScheme() const
    {
      return m_ScoreDistance->GetScheme();
    }

    //! @brief access to the minimal sequence separation
    //! @return minimal sequence separation for neighbors of that exposure score between amino acids in the same chain
    size_t AANeighborhoodDistances::GetMinimalSequenceSeparation() const
    {
      return m_ScoreDistance->GetMinimalSequenceSeparation();
    }

    //! @brief access to the distance cutoff
    //! @return distance cutoff above which the neighbor does not have influence on the score anymore
    double AANeighborhoodDistances::GetDistanceCutoff() const
    {
      return m_ScoreDistance->GetDistanceCutoff();
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &AANeighborhoodDistances::GetAlias() const
    {
      static const std::string s_name( "AANeighborhoodDistances");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate neighbor potential for a given amino acid and its AANeighborList
    //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
    //! @return neighbor potential for a given amino acid and its AANeighborList
    double AANeighborhoodDistances::operator()
    (
      const assemble::AANeighborList &AA_NEIGHBOR_LIST,
      const util::SiPtr< const biol::Membrane> &MEMBRANE
    ) const
    {
      // check that the given neighbor list has the proper parameters
      BCL_Assert
      (
        AA_NEIGHBOR_LIST.GetDistanceCutoff() == m_ScoreDistance->GetDistanceCutoff() &&
        AA_NEIGHBOR_LIST.GetMinimalSequenceSeparation() == m_ScoreDistance->GetMinimalSequenceSeparation(),
        "given neighbor list does not have the proper distance cutoff or sequence separation: " +
        util::Format()( AA_NEIGHBOR_LIST.GetDistanceCutoff()) + " == " + util::Format()( m_ScoreDistance->GetDistanceCutoff()) + " && " +
        util::Format()( AA_NEIGHBOR_LIST.GetMinimalSequenceSeparation()) = " != " + util::Format()( m_ScoreDistance->GetMinimalSequenceSeparation())
      );

      // store the address of the current aa
      const biol::AABase &current_aa( *AA_NEIGHBOR_LIST.GetCenterAminoAcid());

      // check that the current amino acid has a defined coordinate
      if( !current_aa.GetFirstSidechainAtom().GetCoordinates().IsDefined())
      {
        return 0.0;
      }

      double score( 0);
      size_t count( 0);

      // iterate over all neighbors
      for
      (
        assemble::AANeighborList::const_iterator itr( AA_NEIGHBOR_LIST.Begin()), itr_end( AA_NEIGHBOR_LIST.End());
        itr != itr_end;
        ++itr
      )
      {
        score += m_ScoreDistance->operator ()( current_aa, *itr->First(), itr->Second());
        ++count;
      }

      // normalize by number of neighbors
      if( count > 0)
      {
        score /= count;
      }

      // end
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AANeighborhoodDistances::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ScoreDistance, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @return returns the output stream
    std::ostream &AANeighborhoodDistances::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ScoreDistance, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AA_NEIGHBOR_LIST AANeighborList to be used for a single AABase
    //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AANeighborhoodDistances::WriteDetailedSchemeAndValues
    (
      const assemble::AANeighborList &AA_NEIGHBOR_LIST,
      const util::SiPtr< const biol::Membrane> &MEMBRANE,
      std::ostream &OSTREAM
    ) const
    {
      // check that the given neighbor list has the proper parameters
      BCL_Assert
      (
        AA_NEIGHBOR_LIST.GetDistanceCutoff() == m_ScoreDistance->GetDistanceCutoff() &&
        AA_NEIGHBOR_LIST.GetMinimalSequenceSeparation() == m_ScoreDistance->GetMinimalSequenceSeparation(),
        "given neighbor list does not have the proper distance cutoff or sequence separation: " +
        util::Format()( AA_NEIGHBOR_LIST.GetDistanceCutoff()) + " == " + util::Format()( m_ScoreDistance->GetDistanceCutoff()) + " && " +
        util::Format()( AA_NEIGHBOR_LIST.GetMinimalSequenceSeparation()) = " != " + util::Format()( m_ScoreDistance->GetMinimalSequenceSeparation())
      );

      // store the address of the current aa
      const biol::AABase &current_aa( *AA_NEIGHBOR_LIST.GetCenterAminoAcid());

      // check that the current aminoacid has a defined coordinate
      if( !current_aa.GetFirstSidechainAtom().GetCoordinates().IsDefined())
      {
        return OSTREAM;
      }

      // iterate over all neighbors
      for
      (
        assemble::AANeighborList::const_iterator itr( AA_NEIGHBOR_LIST.Begin()), itr_end( AA_NEIGHBOR_LIST.End());
        itr != itr_end;
        ++itr
      )
      {
        m_ScoreDistance->WriteDetailedSchemeAndValues( current_aa, *itr->First(), OSTREAM) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief return parameters for data members that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AANeighborhoodDistances::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Scores the neighbors by adding all aa pair distance scores");
      return serializer;
    }

  } // namespace score
} // namespace bcl
