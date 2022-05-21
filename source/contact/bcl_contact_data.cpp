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
#include "contact/bcl_contact_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    //! single instance
    const util::SiPtr< const util::ObjectInterface> Data::s_Instance( GetObjectInstances().AddInstance( new Data()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Data::Data() :
      m_Data( Types::s_NumberValidTypes, util::GetUndefined< double>()), // set all probabilities to undefined
      m_Merged( util::GetUndefined< double>())
    {
    }

    //! @brief Clone function
    //! @return pointer to new Data
    Data *Data::Clone() const
    {
      return new Data( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Data::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the merged contact prediction
    //! @return a const reference to the prediction value
    const double &Data::MergedPrediction() const
    {
      return m_Merged;
    }

    //! @brief return the merged contact prediction
    //! @return a changeable reference to the prediction value
    double &Data::MergedPrediction()
    {
      return m_Merged;
    }

    //! @brief returns if all contact prediction data is defined
    //! @return if data is defined
    bool Data::IsDefined() const
    {
      return m_Data.IsDefined() && util::IsDefined( m_Merged);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief swaps prediction values at positions CONTACT_TYPE_A and CONTACT_TYPE_B
    //! @param CONTACT_TYPE_A first type position to swap
    //! @param CONTACT_TYPE_B second type position to swap
    void Data::Swap( const Type CONTACT_TYPE_A, const Type CONTACT_TYPE_B)
    {
      std::swap( m_Data( CONTACT_TYPE_A), m_Data( CONTACT_TYPE_B));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief returns the prediction for contact type CONTACT_TYPE
    //! @param CONTACT_TYPE the contact type for which to return the prediction
    //! @return a const reference to the prediction value
    const double &Data::operator[]( const Type &CONTACT_TYPE) const
    {
      return m_Data( CONTACT_TYPE);
    }

    //! @brief returns the prediction for contact type CONTACT_TYPE
    //! @param CONTACT_TYPE the contact type for which to return the prediction
    //! @return a changeable reference to the prediction value
    double &Data::operator[]( const Type &CONTACT_TYPE)
    {
      return m_Data( CONTACT_TYPE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Data::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);
      io::Serialize::Read( m_Merged, ISTREAM);

      return ISTREAM; // return the stream
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Data::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_Merged, OSTREAM, INDENT);

      return OSTREAM; // return the stream
    }

  } // namespace contact
} // namespace bcl
