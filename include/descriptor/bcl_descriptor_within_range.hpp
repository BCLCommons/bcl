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
#include "bcl_descriptor_within_range.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_combine.h"
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_trigonometric_transition.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    WithinRange< t_DataType>::WithinRange() :
      m_RangeBegin( util::GetUndefined< float>()),
      m_RangeEnd( util::GetUndefined< float>()),
      m_InclusiveRange( true)
    {
    }

    //! @brief constructor
    //! @param RANGE_BEGIN where return values of 1 should start
    //! @param RANGE_END where return values of 1 should end
    //! @param INCLUSIVE whether points at the ends of the ranges are included 
    template< typename t_DataType>
    WithinRange< t_DataType>::WithinRange
    (
      const float &RANGE_BEGIN,
      const float &RANGE_END,
      const bool &INCLUSIVE
    ) :
      m_RangeBegin( RANGE_BEGIN),
      m_RangeEnd( RANGE_END),
      m_InclusiveRange( INCLUSIVE)
    {
      BCL_Assert( RANGE_BEGIN <= RANGE_END, "Beginning and end of ranges are reversed");
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    WithinRange< t_DataType> *WithinRange< t_DataType>::Clone() const
    {
      return new WithinRange( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &WithinRange< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &WithinRange< t_DataType>::GetAlias() const
    {
      static const std::string s_Name( "WithinRange");

      return s_Name;
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_DataType>
    size_t WithinRange< t_DataType>::GetNormalDimension() const
    {
      return m_Descriptor->GetType().GetDimension();
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t WithinRange< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

    //! @brief set descriptor to be queried
    //! @param DESCRIPTOR the descriptor to use
    template< typename t_DataType>
    void WithinRange< t_DataType>::SetDescriptor( const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR)
    {
      m_Descriptor = DESCRIPTOR;
    }

    //! @brief get the descriptor that is queried
    //! @return an implementation of the descriptor
    template< typename t_DataType>
    const util::Implementation< Base< t_DataType, float> > &WithinRange< t_DataType>::GetDescriptor() const
    {
      return m_Descriptor;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool WithinRange< t_DataType>::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      return m_RangeBegin <= m_RangeEnd;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer WithinRange< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "This descriptor takes on a value of 1 if a descriptor's mean value is within a range, or 0 otherwise"
      );

      parameters.AddInitializer
      (
        "descriptor",
        "the descriptor to consider",
        io::Serialization::GetAgent( &m_Descriptor)
      );

      parameters.AddInitializer
      (
        "begin",
        "beginning of the range",
        io::Serialization::GetAgent( &m_RangeBegin)
      );

      parameters.AddInitializer
      (
        "end",
        "end of the region where properties will be within range",
        io::Serialization::GetAgent( &m_RangeEnd)
      );

      parameters.AddInitializer
      (
        "inclusive",
        "whether to include the endpoints in the range",
        io::Serialization::GetAgent( &m_InclusiveRange),
        "true"
      );

      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > WithinRange< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    template< typename t_DataType>
    Type::Symmetry WithinRange< t_DataType>::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    template< typename t_DataType>
    bool WithinRange< t_DataType>::ConsiderRepeatedElements() const
    {
      return m_Descriptor->GetType().ConsiderRepeatedObjects();
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void WithinRange< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // instantiate the property as a vector with indices that correspond to atoms
      const linal::VectorConstReference< float> &result( m_Descriptor->operator ()( ITR));
      float mean( result.Sum() / result.GetSize());

      STORAGE( 0) = m_InclusiveRange ? ( mean >= m_RangeBegin && mean <= m_RangeEnd) : ( mean > m_RangeBegin && mean < m_RangeEnd);
    }

  } // namespace descriptor
} // namespace bcl
