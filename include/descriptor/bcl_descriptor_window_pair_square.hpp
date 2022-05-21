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

// include header of this class
#include "bcl_descriptor_window_pair_square.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType, typename t_ReturnType>
    WindowPairSquare< t_DataType, t_ReturnType>::WindowPairSquare() :
      m_Descriptor(),
      m_Size( 0)
    {
    }

    //! @brief constructor from a descriptor
    //! @param DESCRIPTOR descriptor to be used
    template< typename t_DataType, typename t_ReturnType>
    WindowPairSquare< t_DataType, t_ReturnType>::WindowPairSquare
    (
      const BasePair< t_DataType, t_ReturnType> &DESCRIPTOR,
      const size_t &SIZE
    ) :
      m_Descriptor( DESCRIPTOR),
      m_Size( SIZE)
    {
    }

    //! @brief constructor from a descriptor
    //! @param DESCRIPTOR descriptor to be used
    template< typename t_DataType, typename t_ReturnType>
    WindowPairSquare< t_DataType, t_ReturnType>::WindowPairSquare
    (
      const util::Implementation< Base< t_DataType, t_ReturnType> > &DESCRIPTOR,
      const size_t &SIZE
    ) :
      m_Descriptor( DESCRIPTOR),
      m_Size( SIZE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new WindowPairSquare
    template< typename t_DataType, typename t_ReturnType>
    WindowPairSquare< t_DataType, t_ReturnType> *WindowPairSquare< t_DataType, t_ReturnType>::Clone() const
    {
      return new WindowPairSquare( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType, typename t_ReturnType>
    const std::string &WindowPairSquare< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &WindowPairSquare< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_window_name( "WindowPairSquare");
      return s_window_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType, typename t_ReturnType>
    size_t WindowPairSquare< t_DataType, t_ReturnType>::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetSizeOfFeatures() * math::Sqr( 2 * m_Size + 1);
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    template< typename t_DataType, typename t_ReturnType>
    Type::Symmetry WindowPairSquare< t_DataType, t_ReturnType>::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    template< typename t_DataType, typename t_ReturnType>
    bool WindowPairSquare< t_DataType, t_ReturnType>::ConsiderRepeatedElements() const
    {
      return m_Descriptor->GetType().ConsiderRepeatedObjects();
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer WindowPairSquare< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes descriptor windows for " + GetStaticClassName< t_DataType>() +
        " returning NaNs for positions beyond the end of the sequence"
      );
      parameters.AddInitializer( "", "descriptor to use for the window", io::Serialization::GetAgent( &m_Descriptor));
      parameters.AddInitializer
      (
        "size",
        "desired size of the window (0 = just central value, 1 = central value and neighboring values for both members of the pair, ...)",
        io::Serialization::GetAgent( &m_Size)
      );

      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT_A, ELEMENT_B the elements of the sequence of interest
    //! @param STORAGE storage for the descriptor
    //! @return true, if the recalculation was completed
    template< typename t_DataType, typename t_ReturnType>
    void WindowPairSquare< t_DataType, t_ReturnType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT_A,
      const iterate::Generic< const t_DataType> &ELEMENT_B,
      linal::VectorReference< t_ReturnType> &STORAGE
    )
    {
      // size of the sub feature
      size_t sub_feature_size( m_Descriptor->GetSizeOfFeatures());

      const Type type( this->GetType());

      // get positions for both iterators
      size_t position_a( ELEMENT_A.GetPosition());
      size_t position_b( ELEMENT_B.GetPosition());
      const size_t seq_size( this->GetCurrentObject()->GetSize());

      // Rewind to the earliest part of the window, if possible
      const size_t start_position_a( position_a - std::min( position_a, m_Size));
      const size_t start_position_b( position_b - std::min( position_b, m_Size));

      // get the ending positions
      const size_t end_position_a( std::min( seq_size, position_a + m_Size + size_t( 1)));
      const size_t end_position_b( std::min( seq_size, position_b + m_Size + size_t( 1)));

      iterate::Generic< const t_DataType> itr_a( ELEMENT_A);
      iterate::Generic< const t_DataType> itr_b( ELEMENT_B);

      const size_t col_size( m_Descriptor->GetSizeOfFeatures() * ( 2 * m_Size + 1));

      // current position in the storage vector reference
      size_t storage_pos( 0);

      if( position_a - start_position_a < m_Size)
      {
        // could not rewind far enough, insert NaNs for first rows
        const size_t row_diff( m_Size + start_position_a - position_a);
        storage_pos += row_diff * col_size;
      }
      itr_a.GotoPosition( start_position_a);
      itr_b.GotoPosition( start_position_b);

      const size_t offset_start_row( start_position_b + m_Size - position_b);
      const size_t offset_end_row( position_b + m_Size + 1 - end_position_b);

      const size_t n_rows( end_position_a - start_position_a);
      const size_t n_used_cols( end_position_b - start_position_b);

      for( size_t row_number( 0); row_number < n_rows; ++row_number, ++itr_a)
      {
        // copy the iterator
        Iterator< t_DataType> itr_seq( type, itr_a, itr_b);
        storage_pos += offset_start_row;
        for( size_t col_number( 0); col_number < n_used_cols; ++col_number, ++itr_seq, storage_pos += sub_feature_size)
        {
          STORAGE.CreateSubVectorReference( sub_feature_size, storage_pos).CopyValues( ( *m_Descriptor)( itr_seq));
        }
        storage_pos += offset_end_row;
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType, typename t_ReturnType>
    iterate::Generic< Base< t_DataType, t_ReturnType> > WindowPairSquare< t_DataType, t_ReturnType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, t_ReturnType> >( &m_Descriptor, &m_Descriptor + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool WindowPairSquare< t_DataType, t_ReturnType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_Descriptor.IsDefined() && m_Descriptor->GetType().GetDimension() != size_t( 2))
      {
        ERR_STREAM << "WindowPairSquare descriptors only work with pair-wise descriptors";
        return false;
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
