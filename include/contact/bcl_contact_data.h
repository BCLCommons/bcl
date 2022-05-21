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

#ifndef BCL_CONTACT_DATA_H_
#define BCL_CONTACT_DATA_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_contact_types.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Data
    //! @brief stores contact prediction for each contact type and an merged contact prediction
    //!
    //! @see @link example_contact_data.cpp @endlink
    //! @author heinzes1
    //! @date Aug 3, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Data :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector< double> m_Data; //!< contact prediction data: probability for each contact type
      double m_Merged; //!< probability of amino acids being in contact

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //! single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Data();

      //! @brief Clone function
      //! @return pointer to new ContactData
      Data *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the merged contact prediction
      //! @return a const reference to the prediction value
      const double &MergedPrediction() const;

      //! @brief return the merged contact prediction
      //! @return a changeable reference to the prediction value
      double &MergedPrediction();

      //! @brief returns if all contact prediction data is defined
      //! @return if data is defined
      bool IsDefined() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief swaps prediction values at positions POS_A and POS_B
      //! @param CONTACT_TYPE_A first type position to swap
      //! @param CONTACT_TYPE_B second type position to swap
      void Swap( const Type CONTACT_TYPE_A, const Type CONTACT_TYPE_B);

    ///////////////
    // operators //
    ///////////////

      //! @brief returns the prediction for contact type CONTACT_TYPE
      //! @param CONTACT_TYPE the contact type for which to return the prediction
      //! @return a const reference to the prediction value
      const double &operator[]( const Type &CONTACT_TYPE) const;

      //! @brief returns the prediction for contact type CONTACT_TYPE
      //! @param CONTACT_TYPE the contact type for which to return the prediction
      //! @return a changeable reference to the prediction value
      double &operator[]( const Type &CONTACT_TYPE);

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

    }; // class Data

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_DATA_H_
