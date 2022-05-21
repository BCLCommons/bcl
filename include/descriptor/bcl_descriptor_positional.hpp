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
#include "bcl_descriptor_positional.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor given position
    template< typename t_DataType, typename t_ReturnType>
    Positional< t_DataType, t_ReturnType>::Positional( const size_t &POSITION) :
      m_Column( POSITION)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType, typename t_ReturnType>
    Positional< t_DataType, t_ReturnType> *Positional< t_DataType, t_ReturnType>::Clone() const
    {
      return new Positional( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Positional< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Positional< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_names[] = { "1st", "2nd", "3rd", "4th", "5th"};
      return s_names[ m_Column];
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
    template< typename t_DataType, typename t_ReturnType>
    void Positional< t_DataType, t_ReturnType>::SetDimensionHook()
    {
      // Check that setting the dimension is valid in this case
      const size_t setting( this->GetDimensionSetting());
      BCL_Assert( setting, "Cannot use position descriptor " + GetAlias() + " for a " + this->GetObjectName() + " descriptor");
      BCL_Assert
      (
        setting - 1 >= m_Column,
        "Cannot get " + GetAlias() + " " + this->GetElementName() + " when only " + util::Format()( setting) + " is considered"
      );
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool Positional< t_DataType, t_ReturnType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      const size_t dimension( m_Property->GetType().GetDimension());
      if( dimension != size_t( 1))
      {
        ERR_STREAM << GetAlias() << " can only be used with element-wise descriptors, not " << m_Property.GetLabel();
        return false;
      }
      // ensure that if the dimensionality sub_itris increased, the feature size doubles as well
      const size_t n_features_dim_1( m_Property->GetSizeOfFeatures());
      m_Property->SetDimension( 2);
      if( n_features_dim_1 * 2 != m_Property->GetSizeOfFeatures())
      {
        ERR_STREAM << GetAlias() << " can only be used with simple element-wise descriptors, not " << m_Property.GetLabel();
        return false;
      }
      m_Property->SetDimension( 1);

      return true;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType, typename t_ReturnType>
    size_t Positional< t_DataType, t_ReturnType>::GetNormalSizeOfFeatures() const
    {
      return m_Property->GetSizeOfFeatures();
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_DataType, typename t_ReturnType>
    size_t Positional< t_DataType, t_ReturnType>::GetNormalDimension() const
    {
      const size_t dim_setting( this->GetDimensionSetting());
      return util::IsDefined( dim_setting) ? dim_setting : m_Column + 1;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer Positional< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "For pairwise or higher dimension descriptor generation, selects the result from the "
        + GetAlias()
        + " sub-object"
      );
      parameters.AddInitializer( "", "descriptor of interest", io::Serialization::GetAgent( &m_Property));
      return parameters;
    } // GetParameters

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType, typename t_ReturnType>
    void Positional< t_DataType, t_ReturnType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< t_ReturnType> &STORAGE
    )
    {
      STORAGE.CopyValues( m_Property->operator()( Iterator< t_DataType>( ITR( m_Column))));
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType, typename t_ReturnType>
    void Positional< t_DataType, t_ReturnType>::SetObjectHook()
    {
      m_Property->SetObject( *Base< t_DataType, t_ReturnType>::GetCurrentObject());
    }

  } // namespace descriptor
} // namespace bcl
