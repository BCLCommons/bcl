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

#ifndef BCL_MODEL_FEATURE_REFERENCE_H_
#define BCL_MODEL_FEATURE_REFERENCE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_interface.h"
#include "io/bcl_io_serialize.h"
#include "linal/bcl_linal_vector_const_interface.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically
#include <iterator>
#include <ostream>

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureReference
    //! @brief References a feature from a feature data set
    //!
    //! @tparam t_DataType can be float, float, int, complex, etc...
    //!
    //! @see @link example_model_feature_reference.cpp @endlink
    //! @author loweew, woetzen
    //! @date Sep 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class FeatureReference :
      public linal::VectorConstInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      //! pointer to the data that is referenced
      const t_DataType *m_Data;

      //! total number of elements in the reference
      size_t            m_Size;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FeatureReference() :
        m_Data( NULL),
        m_Size( 0)
      {
      }

      //! @brief constructs from size and pointer to data
      //! @param SIZE the number of elements to be referenced
      //! @param DATA the pointer to the data
      FeatureReference( const size_t SIZE, const t_DataType *DATA) :
        m_Data( DATA),
        m_Size( SIZE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new FeatureReference< t_DataType>
      FeatureReference< t_DataType> *Clone() const
      {
        return new FeatureReference< t_DataType>( *this);
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

      //! @brief get size of DataSet
      //! @return the number of data items in the DataSet
      size_t GetSize() const
      {
        return m_Size;
      }

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Matrix
      const t_DataType *Begin() const
      {
        return m_Data;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Matrix
      const t_DataType *End() const
      {
        return m_Data + m_Size;
      }

      //! @brief bool indicating whether or not feature has ownership
      //! @return bool
      bool HasOwnership() const
      {
        return false;
      }

    ///////////////
    // operators //
    ///////////////

      //! return copy of element ( POS)
      const t_DataType &operator()( const size_t POS) const
      {
        AssertValidPosition( POS);
        return m_Data[ POS];
      }

      bool operator <( const FeatureReference< t_DataType> &DATA) const
      {
        BCL_Assert( DATA.GetSize() == m_Size, "unequal size");
        for( size_t i( 0), size( m_Size); i < size; ++i)
        {
          if( m_Data[ i] != DATA( i))
          {
            return m_Data[ i] < DATA( i);
          }
        }
        return false;
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

        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Size, OSTREAM, INDENT) << '\n';
        io::Serialize::InsertIndent( OSTREAM, INDENT);

        // write elements
        std::ostream_iterator< t_DataType> itr( OSTREAM, "\t");
        std::copy( m_Data, m_Data + m_Size, itr);

        // return the stream
        return OSTREAM;
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

    private:

    }; // template class FeatureReference

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> FeatureReference< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new FeatureReference< t_DataType>())
    );

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_REFERENCE_H_
