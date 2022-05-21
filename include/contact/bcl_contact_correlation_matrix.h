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

#ifndef BCL_CONTACT_CORRELATION_MATRIX_H_
#define BCL_CONTACT_CORRELATION_MATRIX_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_symmetric_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CorrelationMatrix
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_contact_correlation_matrix.cpp @endlink
    //! @author teixeipl
    //! @date Jul 30, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CorrelationMatrix :
      public storage::SymmetricMatrix< double>
    {
    private:

      //! Member value to use for replacements
      double m_Value;
//      //! Member bool stating whether values are to be replaced when returned
//      bool m_ReplaceValues;

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Constructor taking in dimension
      CorrelationMatrix( const size_t DIMENSIONS = 0, const double VALUE = 0.0, const bool REPLACE = true);

      //! @brief Clone function
      //! @return pointer to new CorrelationMatrix
      CorrelationMatrix *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief Returns the average for the symmetric matrix
      //! @return double corresponding to average correlation value of matrix excluding NaN and the diagonal
      double GetAverage() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns value for replacing NaN
      //! @return double corresponding to replacing value
      const double &GetValue() const
      {
        return m_Value;
      }

      //! @brief Sets the value for replacing NaN
      //! @param VALUE to be set
      void SetValue( const double &VALUE)
      {
        m_Value = VALUE;
      }

    ///////////////
    // operators //
    ///////////////

      //! operator( POS) return const reference to element at POS
      const double &operator()( const size_t I_POS, const size_t J_POS) const;

      //! operator( POS) return reference to changeable element at POS
      double &operator()( const size_t I_POS, const size_t J_POS);

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

    }; // class CorrelationMatrix

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_CORRELATION_MATRIX_H_
