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

#ifndef BCL_LINAL_VECTOR_CONST_REFERENCE_H_
#define BCL_LINAL_VECTOR_CONST_REFERENCE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_const_interface.h"
#include "io/bcl_io_serialize.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorConstReference
    //! @brief const reference to a particular vector range
    //!
    //! @see @link example_linal_vector_const_reference.cpp @endlink
    //! @author mendenjl
    //! @date Apr 02, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorConstReference :
      public VectorConstInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      //! length of VectorConstReference
      size_t m_Size;

      //! range of dynamically allocated memory of size m_Size
      const t_DataType *m_Data;

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
      VectorConstReference() :
        m_Size( 0),
        m_Data( NULL)
      {
      }

      //! @brief construct from dimension and pointer to data
      //! @param SIZE number of elements in DATA to reference
      //! @param DATA pointer to field of data
      VectorConstReference
      (
        const size_t SIZE,
        const t_DataType *DATA
      ) :
        m_Size( SIZE),
        m_Data( DATA)
      {
        // check if allocation was successful
        BCL_Assert
        (
          m_Data != NULL || m_Size == size_t( 0),
          "cannot reference NULL data elements of type: " + GetStaticClassName< t_DataType>()
        );
      }

      //! @brief construct from dimension and pointer to data
      //! @param SIZE number of elements in DATA to reference
      //! @param DATA pointer to field of data
      VectorConstReference
      (
        const size_t SIZE,
        t_DataType *DATA
      ) :
        m_Size( SIZE),
        m_Data( DATA)
      {
        // check if allocation was successful
        BCL_Assert
        (
          m_Data != NULL || m_Size == size_t( 0),
          "cannot reference NULL data elements of type: " + GetStaticClassName< t_DataType>()
        );
      }

      //! @brief copy constructor from VectorConstInterface
      //! @param VECTOR vector interface to be copied from
      explicit VectorConstReference( const VectorConstInterface< t_DataType> &VECTOR) :
        m_Size( VECTOR.GetSize()),
        m_Data( VECTOR.Begin())
      {
      }

      //! @brief Clone function
      //! @return pointer to new VectorConstReference< t_DataType>
      VectorConstReference< t_DataType> *Clone() const
      {
        return new VectorConstReference< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief size of VectorConstReference
      //! @return size of VectorConstReference
      size_t GetSize() const
      {
        return m_Size;
      }

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of VectorConstReference
      const t_DataType *Begin() const
      {
        return m_Data;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in VectorConstReference
      const t_DataType *End() const
      {
        return m_Data + m_Size;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief create subVectorConstReference from this VectorConstReference
      //! @param SIZE number of elements in subVectorConstReference
      //! @param POS starting position in VectorConstReference - default 0
      //! @return VectorConstReference, that is a subVectorConstReference of this
      VectorConstReference< t_DataType> CreateSubVectorConstReference
      (
        const size_t SIZE,
        const size_t POS = 0
      ) const
      {
        // check that there are enough cols and rows
        BCL_Assert( SIZE + POS <= m_Size, "this VectorConstReference is too small for subVectorConstReference");

        // if size is equal (pos have to be 0 as asserted)
        if( SIZE == m_Size)
        {
          return *this;
        }

        // create VectorConstReference
        return VectorConstReference< t_DataType>( SIZE, m_Data + POS);
      }

    ///////////////
    // operators //
    ///////////////

      //! return reference of element ( POS)
      const t_DataType &operator()( const size_t POS) const
      {
        AssertValidPosition( POS);
        return m_Data[ POS];
      }

      //! return pointer to element
      const t_DataType *operator[]( const size_t POS) const
      {
        AssertValidPosition( POS);
        return m_Data + POS;
      }

      //! @brief equal operator
      //! @param VectorConstReference source VectorConstReference
      //! @return reference to this assigned VectorConstReference
      VectorConstReference< t_DataType> &operator =( const VectorConstReference< t_DataType> &VECTOR)
      {
        m_Data = VECTOR.m_Data;
        m_Size = VECTOR.m_Size;
        // return reference to this VectorConstReference
        return *this;
      }

      //! @brief equal operator
      //! @param VectorConstReference source VectorConstReference
      //! @return reference to this assigned VectorConstReference
      VectorConstReference< t_DataType> &operator =( const VectorConstInterface< t_DataType> &VECTOR)
      {
        m_Data = VECTOR.Begin();
        m_Size = VECTOR.GetSize();
        // return reference to this VectorConstReference
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! write VectorConstReference to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const util::Format &FORMAT) const
      {
        // write size
        OSTREAM << m_Size << '\n';

        //write data
        for( const t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
        {
          OSTREAM << ' ' << FORMAT( *ptr);
        }

        // end
        return OSTREAM;
      }

    protected:

      //! write VectorConstReference to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write size
        io::Serialize::Write( m_Size, OSTREAM, INDENT) << '\n';
        io::Serialize::InsertIndent( OSTREAM, INDENT);

        //write data
        for( const t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
        {
          io::Serialize::Write( *ptr, OSTREAM) << '\t';
        }

        // end
        return OSTREAM;
      }

      //! read VectorConstReference from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        BCL_Exit( "You can't read in a vector reference because it does not have ownership!", -1);

        // end
        return ISTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! check whether position is valid
      void AssertValidPosition( const size_t POS) const
      {
        BCL_Assert
        (
          POS < m_Size,
          "cannot access element outside range! " + util::Format()( POS) + " >= " + util::Format()( m_Size)
        );
      }

    }; // template class VectorConstReference

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_CONST_REFERENCE_H_
