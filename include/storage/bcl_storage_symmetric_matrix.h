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

#ifndef BCL_STORAGE_SYMMETRIC_MATRIX_H_
#define BCL_STORAGE_SYMMETRIC_MATRIX_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_storage_vector.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SymmetricMatrix
    //! @brief SymmetricMatrix represents a matrix where all positions i,j are = j,i and either query returns the
    //!              same value
    //! @details
    //! a.) has a constant number of elements it holds
    //! b.) no iterators exist for this container
    //! c.) random access to elements
    //! d.) no insertion or deletion of elements
    //! e.) can return const reference or changeable reference to elements in the matrix
    //!
    //! @see @link example_storage_symmetric_matrix.cpp @endlink
    //! @author teixeipl
    //! @date Aug 16, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SymmetricMatrix :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      size_t m_Size;

      Vector< t_DataType> m_StorageVector;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SymmetricMatrix() :
        m_Size( 0),
        m_StorageVector()
      {
      }

      //! @brief constructor with size specified
      //! @param N is size_t specifying the length along both sides of the symmetric matrix
      SymmetricMatrix( const size_t N) :
        m_Size( N),
        m_StorageVector( N * ( N + 1) / 2)
      {
      }

      //! @brief Clone function
      //! @return pointer to new SymmetricMatrix
      SymmetricMatrix< t_DataType> *Clone() const
      {
        return new SymmetricMatrix< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief GetSize returns size of the container
      //! @return returns size_t size, i.e. number of elements stored
      size_t GetSize() const
      {
        return m_Size;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief fills the matrix with the given value, overwrites current matrix
      //! @param FILL_VALUE to use
      void Fill( const t_DataType FILL_VALUE)
      {
        m_StorageVector.Resize( m_Size, FILL_VALUE);
      }

      //! @brief Resets the SymmetricMatrix with default size of 0, deletes all data
      //! @param N sets the size for the newly reset SymmetricMatrix
      void Reset()
      {
        m_StorageVector.Reset();
        m_Size = 0;
      }

      //! @brief Resizes matrix and fills new vector with value provided
      //! @param NEW_SIZE to change matrix to
      //! @param FILL_VALUE to use to fill new matrix
      void Resize( const size_t NEW_SIZE, const t_DataType FILL_VALUE = t_DataType())
      {
        m_StorageVector.Resize( NEW_SIZE * ( NEW_SIZE + 1) / 2, FILL_VALUE);
        m_Size = NEW_SIZE;
      }

//! TODO: examine container / linal inheritance

    ///////////////
    // operators //
    ///////////////

      //! operator( POS) return const reference to element at POS
      const t_DataType &operator()( const size_t I_POS, const size_t J_POS) const
      {
        // Ensure that i !> j
        if( I_POS <= J_POS)
        {
          // swaps i and j in the case that i > j since only half the matrix is stored and i,j = j,i
          return m_StorageVector( ( ( ( J_POS + 1) * J_POS) / 2) + I_POS);
        }
        // else
        // Get position i,j in the matrix by getting position j + ((i * (i-1))/2 in the vector
        return m_StorageVector( ( ( ( I_POS + 1) * I_POS) / 2) + J_POS);
      }

      //! operator( POS) return reference to changeable element at POS
      t_DataType &operator()( const size_t I_POS, const size_t J_POS)
      {
        // Ensure that i !> j
        if( I_POS <= J_POS)
        {
          // swaps i and j in the case that i > j since only half the matrix is stored and i,j = j,i
          return m_StorageVector( ( ( ( J_POS + 1) * J_POS) / 2) + I_POS);
        }
        // else
        // Get position i,j in the matrix by getting position j + ((i * (i-1))/2 in the vector
        return m_StorageVector( ( ( ( I_POS + 1) * I_POS) / 2) + J_POS);
      }

      //! operator== checks whether two SymmetricMatrices are equal
      bool operator==( const SymmetricMatrix< t_DataType> &MATRIX_A) const
      {
        return m_StorageVector == MATRIX_A.m_StorageVector;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_Size, ISTREAM);
        io::Serialize::Read( m_StorageVector, ISTREAM);

        return ISTREAM; // return the stream
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Size, OSTREAM, INDENT) << "\n";
        io::Serialize::Write( m_StorageVector, OSTREAM, INDENT);

        return OSTREAM; // return the stream
      }

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class SymmetricMatrix

  } // namespace storage
} // namespace bcl

#endif // BCL_STORAGE_SYMMETRIC_MATRIX_H_
