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
#include "bcl_descriptor_partial.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief virtual copy constructor
    template< typename t_DataType, typename t_ReturnType>
    Partial< t_DataType, t_ReturnType> *Partial< t_DataType, t_ReturnType>::Clone() const
    {
      return new Partial( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Partial< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Partial< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_name( "Partial");
      return s_name;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
    template< typename t_DataType, typename t_ReturnType>
    void Partial< t_DataType, t_ReturnType>::SetDimensionHook()
    {
      // Check that setting the dimension is valid in this case
      if( m_Property->GetType().GetDimension() != size_t( 1))
      {
        BCL_Assert( m_MaxColumn < m_Property->GetSizeOfFeatures(), "Bad column index: " + util::Format()( m_MaxColumn));
      }
      else
      {
        const size_t feature_size( m_Property->GetSizeOfFeatures());
        BCL_Assert
        (
          m_MaxColumn < feature_size,
          "Changing dimension to " + util::Format()( this->GetDimensionSetting())
          + " invalidate the partial descriptor: max col id: " + util::Format()( m_MaxColumn)
          + " >= max feature id: " + util::Format()( feature_size)
        );
      }
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool Partial< t_DataType, t_ReturnType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      const size_t dimension( m_Property->GetType().GetDimension());
      m_MaxColumn = 0;
      if( m_Columns.IsEmpty())
      {
        // in no columns are specified, return nothing
        return true;
      }

      // find the maximum index specified
      m_MaxColumn = math::Statistics::MaximumValue( m_Columns.Begin(), m_Columns.End());

      // determine the normal feature size of the sub property
      const size_t normal_feature_size( m_Property->GetSizeOfFeatures());

      if( dimension != size_t( 1))
      {
        if( m_MaxColumn >= normal_feature_size)
        {
          ERR_STREAM << "Tried to extract parts of the descriptor that do not exist (index "
                     << m_MaxColumn << " >= " << normal_feature_size;
          return false;
        }
        return true;
      }

      // handle the case that the dimensions are already plausible given the maximum index
      if( normal_feature_size > m_MaxColumn)
      {
        return true;
      }

      const size_t min_dimension( m_MaxColumn / normal_feature_size + 1);
      size_t old_max_column( m_MaxColumn);
      // set max column to 0 so that SetDimension does not assert
      m_MaxColumn = 0;
      Base< t_DataType, t_ReturnType>::SetDimension( min_dimension);
      m_MaxColumn = old_max_column;

      // check that the dimension was really increased
      if( m_Property->GetSizeOfFeatures() <= m_MaxColumn)
      {
        ERR_STREAM << "Tried to extract parts of the descriptor that do not exist (index "
                   << m_MaxColumn << " >= " << m_Property->GetSizeOfFeatures();
        return false;
      }
      return true;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType, typename t_ReturnType>
    size_t Partial< t_DataType, t_ReturnType>::GetNormalSizeOfFeatures() const
    {
      return m_Columns.GetSize();
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_DataType, typename t_ReturnType>
    size_t Partial< t_DataType, t_ReturnType>::GetNormalDimension() const
    {
      const size_t dim_setting( this->GetDimensionSetting());
      return util::IsDefined( dim_setting) ? dim_setting : m_Property->GetType().GetDimension();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer Partial< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "selects particular values (by index) of another descriptor");
      parameters.AddInitializer( "", "descriptor of interest", io::Serialization::GetAgent( &m_Property));
      parameters.AddInitializer
      (
        "indices",
        "desired indices (0-offset) of the descriptor to keep",
        io::Serialization::GetAgent( &m_Columns)
      );
      return parameters;
    } // GetParameters

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType, typename t_ReturnType>
    iterate::Generic< Base< t_DataType, t_ReturnType> > Partial< t_DataType, t_ReturnType>::GetInternalDescriptors()
    {
      if( m_Property.IsDefined())
      {
        return iterate::Generic< Base< t_DataType, t_ReturnType> >( &m_Property, &m_Property + 1);
      }
      return iterate::Generic< Base< t_DataType, t_ReturnType> >();
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType, typename t_ReturnType>
    void Partial< t_DataType, t_ReturnType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< t_ReturnType> &STORAGE
    )
    {
      const linal::VectorConstReference< t_ReturnType> full_property( m_Property->operator()( ITR));
      t_ReturnType *itr_storage( STORAGE.Begin());
      for
      (
        storage::Vector< size_t>::const_iterator itr( m_Columns.Begin()), itr_end( m_Columns.End());
        itr != itr_end;
        ++itr, ++itr_storage
      )
      {
        *itr_storage = full_property( *itr);
      }
    }

  } // namespace descriptor
} // namespace bcl
