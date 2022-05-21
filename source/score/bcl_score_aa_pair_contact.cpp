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
#include "score/bcl_score_aa_pair_contact.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "contact/bcl_contact_prediction_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAPairContact::s_Instance
    (
      GetObjectInstances().AddInstance( new AAPairContact())
    );

    //! lower and higher range of CB distance to be considered as a contact
    const storage::Pair< double, double> AAPairContact::s_CBDistanceRange = storage::Pair< double, double>( 4.0, 12.0);

    //! probability shift for each contact type
    const double AAPairContact::s_ProbabilityShift[ contact::Types::s_NumberValidTypes] = { 0.5, 0.5, 0.5, 0.5, 0.5};

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAPairContact::AAPairContact()
    {
    }

    //! @brief constructor from a contact type and a PredictionMap
    //! @param CONTACT_TYPE corresponding contact type
    //! @param SP_PREDICTION_MAP ShPtr to PredictionMap of interest
    AAPairContact::AAPairContact
    (
      const contact::Type &CONTACT_TYPE,
      const util::ShPtr< contact::PredictionMap> &SP_PREDICTION_MAP
    ) :
      m_ContactType( CONTACT_TYPE),
      m_PredictionMap( SP_PREDICTION_MAP)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairContact copied from this one
    AAPairContact *AAPairContact::Clone() const
    {
      return new AAPairContact( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairContact::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator to calculate the contact score for the given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return the contact score for the given amino acid pair
    double AAPairContact::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    ) const
    {
      BCL_MessageDbg( "score aa pair contact");
      // retrieve and store the CB distance between AA pair.
      const double cb_distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // if the CB distance is greater then right side of CB range, there is no contact so energy contribution is 0.
      if( cb_distance >= s_CBDistanceRange.Second())
      {
        return double( 0);
      }

      // store the iterator in the map that corresponds to these two amino acids
      const double energy
      (
        s_ProbabilityShift[ m_ContactType] -
        m_PredictionMap->GetPredictions
                         (
                           storage::VectorND< 2, util::SiPtr< const biol::AAData> >
                           (
                             ( *AMINO_ACID_A.GetData()),
                             ( *AMINO_ACID_B.GetData())
                           )
                         ).First()( m_ContactType)
      );

      // if there is no available prediction for these aapairs ( probably they are in the border or close neighbors) return 0
      if( !util::IsDefined( energy))
      {
        return double( 0);
      }

      // if the CB distance is smaller than the left side of CB range, return -prediction;
      else if( cb_distance <= s_CBDistanceRange.First())
      {
        return energy;
      }

      // if the CB distance is within range, calculate a score and return it
      return CalculateScore( cb_distance, energy);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAPairContact::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ContactType, ISTREAM);
      io::Serialize::Read( m_PredictionMap, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @return returns the output stream
    std::ostream &AAPairContact::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ContactType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PredictionMap, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairContact::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      std::ostream &OSTREAM
    ) const
    {
      //write both amino acids, cb distance, contacttype and value
      OSTREAM << AMINO_ACID_A.GetSeqID() << '\t'
              << AMINO_ACID_A.GetType()->GetThreeLetterCode() << '\t'
              << AMINO_ACID_B.GetSeqID() << '\t'
              << AMINO_ACID_B.GetType()->GetThreeLetterCode() << '\t'
              << biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B) << '\t'
              << m_ContactType << '\t' << operator()( AMINO_ACID_A, AMINO_ACID_B) << '\n';

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the score given a distance and a prediction vaue from the ANNs
    //! @param DISTANCE Distance between amino acids
    //! @param ENERGY energy value from the ANN
    //! @return the score given a distance and a prediction vaue from the ANNs
    double AAPairContact::CalculateScore( const double DISTANCE, const double ENERGY) const
    {
      // calculate the angle
      const double x
      (
        ( DISTANCE - s_CBDistanceRange.First()) *
        ( math::g_Pi) /
        ( s_CBDistanceRange.Second() - s_CBDistanceRange.First())
      );

      //  Cosine of the calculated value is multiplied with the prediction and returned
      return ENERGY * ( cos( x) + 1) / 2;
    }

  } // namespace score
} // namespace bcl
