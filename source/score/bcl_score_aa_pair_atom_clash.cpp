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
#include "score/bcl_score_aa_pair_atom_clash.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AAPairAtomClash::s_Instance
    (
      GetObjectInstances().AddInstance( new AAPairAtomClash())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AAPairAtomClash::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "atomclash");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAPairAtomClash::AAPairAtomClash() :
      m_SigmoidWidth( 1.0),
      m_DistanceCutoff( 12.0),
      m_MinimalSequenceSeparation( 0),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief constructor from a specified histogram file
    //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
    //! @param MINIMAL_SEQUENCE_SEPARATION minimal sequence separation
    //! @param SCHEME scheme to be used
    AAPairAtomClash::AAPairAtomClash
    (
      const double SIGMOID_WIDTH,
      const size_t MINIMAL_SEQUENCE_SEPARATION,
      const std::string &SCHEME
    ) :
      m_SigmoidWidth( SIGMOID_WIDTH),
      m_DistanceCutoff( 12.0),
      m_MinimalSequenceSeparation( MINIMAL_SEQUENCE_SEPARATION),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AAPairAtomClash
    AAPairAtomClash *AAPairAtomClash::Clone() const
    {
      return new AAPairAtomClash( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairAtomClash::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param ATOM_A first atom of interest
    //! @param ATOM_B second atom of interest
    //! @return clash score for the given atom pair
    double AAPairAtomClash::operator()
    (
      const biol::Atom &ATOM_A,
      const biol::Atom &ATOM_B
    ) const
    {
      // calculate the sum covalent radius
      // van der waals radius would consider N-C-CA N-CA distance as clash (2.4 < 3.39)
      const double clash_distance
      (
        1.25 *
        (
          ATOM_A.GetType()->GetElementType()->GetProperty( chemistry::ElementTypeData::e_CovalentRadius) +
          ATOM_B.GetType()->GetElementType()->GetProperty( chemistry::ElementTypeData::e_CovalentRadius)
        )
      );

//      const linal::Vector3D &coord_a( ATOM_A.GetCoordinates());
//      const linal::Vector3D &coord_b( ATOM_B.GetCoordinates());
//
//      // check if at least all individual coordinates are too close, only then calculation of actual distance is required
//      for
//      (
//        const double *ptr_a( coord_a.Begin()), *ptr_b( coord_b.Begin()), *ptr_a_end( coord_a.End());
//        ptr_a != ptr_a_end;
//        ++ptr_a, ++ptr_b
//      )
//      {
//        if( math::Absolute( *ptr_a - *ptr_b) > clash_distance)
//        {
//          return 0.0;
//        }
//      }

      // calculate the distance between these two atoms
      const double distance( biol::Distance( ATOM_A, ATOM_B));

      // calculate and return score
      return CalculateRepulsiveTerm( distance, clash_distance);
    }

    //! @brief evaluate clashes for all atoms pairs and return it
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return pair of clash score and the number of scored entities
    double AAPairAtomClash::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    ) const
    {
      // initialize score
      double score( 0.0);

      const int seq_dist( AMINO_ACID_A.GetSeqID() - AMINO_ACID_B.GetSeqID());

      // iterate over the atoms of the first amino acid
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator
          atom_itr_a( AMINO_ACID_A.GetAtoms().Begin()), atom_itr_a_end( AMINO_ACID_A.GetAtoms().End());
        atom_itr_a != atom_itr_a_end; ++atom_itr_a
      )
      {
        const biol::Atom &atom_a( **atom_itr_a);
        const biol::AtomType &type_a( atom_a.GetType());

        for
        (
          util::SiPtrVector< const biol::Atom>::const_iterator
            atom_itr_b( AMINO_ACID_B.GetAtoms().Begin()), atom_itr_b_end( AMINO_ACID_B.GetAtoms().End());
          atom_itr_b != atom_itr_b_end; ++atom_itr_b
        )
        {
          const biol::Atom &atom_b( **atom_itr_b);
          const biol::AtomType &type_b( atom_b.GetType());

          // B following A allows peptide bond
          if( seq_dist == -1 && type_a == biol::GetAtomTypes().C && type_b == biol::GetAtomTypes().N)
          {
            continue;
          }

          // A following B allows peptide bond
          if( seq_dist == 1 && type_a == biol::GetAtomTypes().N && type_b == biol::GetAtomTypes().C)
          {
            continue;
          }

          // calculate and sum up the score
          const double current_repulsive( operator()( atom_a, atom_b));
          score += current_repulsive;
        }
      }

      return score;
    }

    //! @brief evaluate clashes for all atoms pairs and return it
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @return pair of clash score and the number of scored entities
    double AAPairAtomClash::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      double DISTANCE
    ) const
    {
      // call the other operator since CB distance by itself is not enough to calculate the score
      return operator()( AMINO_ACID_A, AMINO_ACID_B);
    }

    //! @brief calculate clash score for a protein model
    //! @param MODEL the protein model of interest
    //! @return amino acid pairing potential for given protein
    double AAPairAtomClash::operator()( const assemble::ProteinModel &MODEL) const
    {
      // first, hash the structured SSE index of each
      return 0.0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AAPairAtomClash::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SigmoidWidth, ISTREAM);
      io::Serialize::Read( m_MinimalSequenceSeparation, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AAPairAtomClash::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SigmoidWidth, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinimalSequenceSeparation, OSTREAM, INDENT);
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairAtomClash::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      std::ostream &OSTREAM
    ) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief repulsive term from distance and sum of van der waals radii
    //! @param DISTANCE the actual distance between atoms of interest
    //! @param VDW_DISTANCE sum of van der waals radii
    //! @return repulsive term
    double AAPairAtomClash::CalculateRepulsiveTerm
    (
      const double DISTANCE,
      const double VDW_DISTANCE
    ) const
    {
      const double difference( VDW_DISTANCE - DISTANCE);

      // no repulsion
      if( difference <= 0.0)
      {
        return 0.0;
      }
      // full repulsion
      else if( difference >= m_SigmoidWidth)
      {
        return 1.0;
      }

      return math::WeightBetweenZeroAndPi( ( ( m_SigmoidWidth - difference) / m_SigmoidWidth) * math::g_Pi);
    }

  } // namespace score
} // namespace bcl
