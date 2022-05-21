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
#include "contact/bcl_contact_correlation_matrix.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CorrelationMatrix::s_Instance
    (
      GetObjectInstances().AddInstance( new CorrelationMatrix())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Constructor taking in dimension
    CorrelationMatrix::CorrelationMatrix( const size_t DIMENSIONS, const double VALUE, const bool REPLACE) :
      storage::SymmetricMatrix< double>( DIMENSIONS),
      m_Value( VALUE)
    {
      // TODO: make another class to handle this case        m_ReplaceValues( REPLACE)
    }

    //! @brief Clone function
    //! @return pointer to new CorrelationMatrix
    CorrelationMatrix *CorrelationMatrix::Clone() const
    {
      return new CorrelationMatrix( *this);
    }

  /////////////////
  // data access //
  /////////////////

    // TODO: Make example test for GetAverage

    //! @brief Returns the average for the symmetric matrix
    //! @return double corresponding to average correlation value of matrix excluding NaN and the diagonal
    double CorrelationMatrix::GetAverage() const
    {
      math::RunningAverage< double> average;
      for( size_t i( 0), size( GetSize()); i < size; ++i)
      {
        for( size_t j( i + 1); j < size; ++j)
        {
          const double value( storage::SymmetricMatrix< double>::operator()( i, j));
          if( util::IsDefined( value))
          {
            average += value;
          }
        }
      }

      return average;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CorrelationMatrix::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! Returns value but substitutes in m_Value if NaN is encountered
    //! TODO: Will this allow changing of NaN?  Should I try and stop this?
    const double &CorrelationMatrix::operator()( const size_t I_POS, const size_t J_POS) const
    {
      const double &value( storage::SymmetricMatrix< double>::operator()( I_POS, J_POS));
      return util::IsDefined( value) ? value : m_Value;
    }

    //! operator( POS) return reference to changeable element at POS
    double &CorrelationMatrix::operator()( const size_t I_POS, const size_t J_POS)
    {
      // Call const method and return
      return storage::SymmetricMatrix< double>::operator()( I_POS, J_POS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CorrelationMatrix::Read( std::istream &ISTREAM)
    {
      // read members
      storage::SymmetricMatrix< double>::Read( ISTREAM);
      io::Serialize::Read( m_Value, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CorrelationMatrix::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      storage::SymmetricMatrix< double>::Write( OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Value, OSTREAM, INDENT);              // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace contact
} // namespace bcl
