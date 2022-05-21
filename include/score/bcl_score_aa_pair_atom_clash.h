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

#ifndef BCL_SCORE_AA_PAIR_ATOM_CLASH_H_
#define BCL_SCORE_AA_PAIR_ATOM_CLASH_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_pair_distance_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairAtomClash
    //! @brief Scores clashes between atoms of a given pair of amino acids
    //! @details Potential that evalutes clashes between amino acids. Unlike AAPairClash score, this score actually
    //! evaluate all atom pairs and uses Van Der Waals radii to determine clashes
    //!
    //! @see @link example_score_aa_pair_atom_clash.cpp @endlink
    //! @author karakam, alexanns
    //! @date Jun 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairAtomClash :
      public AAPairDistanceInterface
    {

    private:

    //////////
    // data //
    //////////

      //! width of the sigmoidal function to be used
      double m_SigmoidWidth;

      //! distance cutoff above which score will always be 0
      double m_DistanceCutoff;

      //! minimal sequence separation
      size_t m_MinimalSequenceSeparation;

      //! scheme to be used in outputting
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AAPairAtomClash();

      //! @brief constructor from a specified histogram file
      //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
      //! @param MINIMAL_SEQUENCE_SEPARATION minimal sequence separation
      //! @param SCHEME scheme to be used
      AAPairAtomClash
      (
        const double SIGMOID_WIDTH,
        const size_t MINIMAL_SEQUENCE_SEPARATION,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new AAPairAtomClash
      AAPairAtomClash *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get sigmoid width
      //! @return sigmoid width
      double GetSigmoidWidth() const
      {
        return m_SigmoidWidth;
      }

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation in sequence distance
      size_t GetMinimalSequenceSeparation() const
      {
        return m_MinimalSequenceSeparation;
      }

      //! @brief access to the maximal distance that is considered
      //! @return the maximal distance between two amino acids that is scored
      double GetDistanceCutoff() const
      {
        return m_DistanceCutoff;
      }

      //! @brief are amino acids of different chains considered
      //! @return return true, if amino acid pairs of different chains are considered
      bool GetConsiderDifferentChain() const
      {
        return true;
      }

      //! @brief return scheme
      //! @return scheme
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate clashes for given atom pair
      //! @param ATOM_A first atom of interest
      //! @param ATOM_B second atom of interest
      //! @return clash score for the given atom pair
      double operator()
      (
        const biol::Atom &ATOM_A,
        const biol::Atom &ATOM_B
      ) const;

      //! @brief evaluate clashes for all atoms pairs and return it
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @return pair of clash score and the number of scored entities
      double operator()
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B
      ) const;

      //! @brief evaluate clashes for all atoms pairs and return it
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @param DISTANCE distance between the amino acid pair
      //! @return pair of clash score and the number of scored entities
      double operator()
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        double DISTANCE
      ) const;

      //! @brief calculate clash score for a protein model
      //! @param MODEL the protein model of interest
      //! @return amino acid pairing potential for given protein
      double operator()( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        std::ostream &OSTREAM
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief repulsive term from distance and sum of van der waals radii
      //! @param DISTANCE the actual distance between atoms of interest
      //! @param VDW_DISTANCE sum of van der waals radii
      //! @return repulsive term
      double CalculateRepulsiveTerm
      (
        const double DISTANCE,
        const double VDW_DISTANCE
      ) const;

    }; // class AAPairAtomClash

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_AA_PAIR_ATOM_CLASH_H_
