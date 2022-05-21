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

#ifndef BCL_SCORE_AA_PAIR_CONTACT_H_
#define BCL_SCORE_AA_PAIR_CONTACT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "contact/bcl_contact.fwd.hh"

// includes from bcl - sorted alphabetically
#include "contact/bcl_contact_types.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairContact
    //! @brief This is a Function derived template class for scoring AA pair Contacts
    //!
    //! @see @link example_score_aa_pair_contact.cpp @endlink
    //! @author karakam, woetzen
    //! @date 05.06.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairContact :
      public math::BinaryFunctionInterfaceSerializable< biol::AABase, biol::AABase, double>
    {

    private:

    //////////
    // data //
    //////////

      //! lower and higher range of CB distance to be considered as a contact
      static const storage::Pair< double, double> s_CBDistanceRange;

      //! probability shift for each contact type
      static const double s_ProbabilityShift[];

      //! corresponding contact type
      contact::Type m_ContactType;

      //! ShPtr to PredictionMap to be used
      util::ShPtr< contact::PredictionMap> m_PredictionMap;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AAPairContact();

      //! @brief constructor from a contact type and a PredictionMap
      //! @param CONTACT_TYPE corresponding contact type
      //! @param SP_PREDICTION_MAP ShPtr to PredictionMap of interest
      AAPairContact
      (
        const contact::Type &CONTACT_TYPE,
        const util::ShPtr< contact::PredictionMap> &SP_PREDICTION_MAP
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new AAPairContact copied from this one
      AAPairContact *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        // initalize static scheme and return it
        static const std::string s_scheme( "aapaircontact");
        return s_scheme;
      }
    ///////////////
    // operators //
    ///////////////

      //! @brief operator to calculate the contact score for the given amino acid pair
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @return the contact score for the given amino acid pair
      double operator()( const biol::AABase &AMINO_ACID_A, const biol::AABase &AMINO_ACID_B) const;

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT number of indentations
      //! @return returns the output stream
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

      //! @brief calculates the score given a distance and a prediction vaue from the ANNs
      //! @param DISTANCE Distance between amino acids
      //! @param ENERGY energy value from the ANN
      //! @return the score given a distance and a prediction vaue from the ANNs
      double CalculateScore( const double DISTANCE, const double ENERGY) const;

    }; //class AAPairContact

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_PAIR_CONTACT_H_
