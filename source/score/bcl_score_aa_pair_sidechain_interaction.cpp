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
#include "score/bcl_score_aa_pair_sidechain_interaction.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> AAPairSidechainInteraction::s_Instance
    (
      util::Enumerated< AAPairDistanceInterface>::AddInstance( new AAPairSidechainInteraction())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AAPairSidechainInteraction::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_identifier( "aaside");

      // end
      return s_default_identifier;

    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    AAPairSidechainInteraction::AAPairSidechainInteraction() :
        m_Scheme(),
        m_AAPairDistance()
    {
    }

    //! @brief constructor from a specified histogram file
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    AAPairSidechainInteraction::AAPairSidechainInteraction
    (
      const std::string &HISTOGRAM_FILENAME,
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME),
      m_AAPairDistance( new AAPairDistance( HISTOGRAM_FILENAME))
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairDistance object that is copied from this one
    AAPairSidechainInteraction *AAPairSidechainInteraction::Clone() const
    {
      return new AAPairSidechainInteraction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairSidechainInteraction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AAPairSidechainInteraction::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operators //
  ///////////////

    //! @brief calculate the weighted energy according to distance and AATypes for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return the weighted energy according to distance and AATypes for given amino acid pair
    double AAPairSidechainInteraction::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    ) const
    {
      // return the weighted aa pair distance score
      return
        WeightOfInteraction( AMINO_ACID_A, AMINO_ACID_B) *
        m_AAPairDistance->operator()( AMINO_ACID_A, AMINO_ACID_B);
    }

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairSidechainInteraction::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const double DISTANCE
    ) const
    {
      // return the weighted aa pair distance score
      return
        WeightOfInteraction( AMINO_ACID_A, AMINO_ACID_B) *
        m_AAPairDistance->operator()( AMINO_ACID_A, AMINO_ACID_B, DISTANCE);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the weight of interaction for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return the weight of interaction for given amino acid pair
    double AAPairSidechainInteraction::WeightOfInteraction
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    )
    {
      // calculate the angle and the corresponding weight for the first amino acid
      const double alpha_a
      (
        linal::ProjAngle
        (
          AMINO_ACID_A.GetCA().GetCoordinates(),
          AMINO_ACID_A.GetFirstSidechainAtom().GetCoordinates(),
          AMINO_ACID_B.GetCA().GetCoordinates()
        )
      );
      const double weight_a( Weight90To45Transition( alpha_a));

      // calculate the angle and the corresponding weight for the second amino acid
      const double alpha_b
      (
        linal::ProjAngle
        (
          AMINO_ACID_B.GetCA().GetCoordinates(),
          AMINO_ACID_B.GetFirstSidechainAtom().GetCoordinates(),
          AMINO_ACID_A.GetCA().GetCoordinates()
        )
      );
      const double weight_b( Weight90To45Transition( alpha_b));

      // return the multiplication of calculated weights
      return weight_a * weight_b;
    }

    //! @brief calculates a weight for the angle
    //! 0 - 45 returns 0; 45 - 90 is a cos transition to 0; larger 90 returns 0
    //! @param ANGLE_RAD angle in readians
    //! @return calculated weight for the angle
    double AAPairSidechainInteraction::Weight90To45Transition( const double ANGLE_RAD)
    {
      // if angle less than 45 return 1
      if( ANGLE_RAD <= math::g_Pi / 4)
      {
        return double( 1);
      }
      // if angle more than 90 return 1
      if( ANGLE_RAD >= math::g_Pi / 2)
      {
        return double( 0);
      }

      // otherwise ( if between 45 and 90) return a cosine transition value based on the angle
      return ( std::cos( ( ANGLE_RAD - math::g_Pi / 4) * 4) + 1) / 2;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAPairSidechainInteraction::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_AAPairDistance, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAPairSidechainInteraction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AAPairDistance, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairSidechainInteraction::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      std::ostream &OSTREAM
    ) const
    {
      //write Scheme
      OSTREAM << AMINO_ACID_A.GetSeqID() << '\t'
              << AMINO_ACID_A.GetType()->GetThreeLetterCode() << '\t'
              << AMINO_ACID_B.GetSeqID()<< '\t'
              << AMINO_ACID_B.GetType()->GetThreeLetterCode() << '\t'
              << biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B) << '\t'
              << operator()( AMINO_ACID_A, AMINO_ACID_B) << '\n';

      //end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
