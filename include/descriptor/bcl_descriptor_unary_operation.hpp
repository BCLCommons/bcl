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
#include "bcl_descriptor_unary_operation.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from an operation
    template< typename t_DataType>
    UnaryOperation< t_DataType>::UnaryOperation
    (
      const util::Implementation< math::AssignmentUnaryInterface> &OPERATION,
      const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR
    ) :
      m_Descriptor( DESCRIPTOR),
      m_Op( OPERATION)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    UnaryOperation< t_DataType> *UnaryOperation< t_DataType>::Clone() const
    {
      return new UnaryOperation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &UnaryOperation< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    template< typename t_DataType>
    const std::string &UnaryOperation< t_DataType>::GetAlias() const
    {
      return m_Op->GetAlias();
    }

    //! @brief Return the assignment operation for this descriptor
    //! @return the assignment operation for this descriptor
    template< typename t_DataType>
    const util::Implementation< math::AssignmentUnaryInterface> &
      UnaryOperation< t_DataType>::GetAssignmentOperation() const
    {
      return m_Op;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer UnaryOperation< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      const std::string description( m_Op->GetSerializer().GetClassDescription());
      parameters.SetClassDescription( description + " of a descriptor");

      parameters.AddInitializer
      (
        "",
        "descriptor for which to " + description,
        io::Serialization::GetAgent( &m_Descriptor)
      );

      return parameters;
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > UnaryOperation< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    CachePreference UnaryOperation< t_DataType>::GetNormalCachePreference() const
    {
      return m_Descriptor->GetCachePreference() == e_NeverCache ? e_NeverCache : e_PreferCache;
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void UnaryOperation< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // initialize the descriptor return vector with the first descriptor
      linal::VectorConstReference< float> result( m_Descriptor->operator ()( ITR));

      // check for vectors of incorrect size -> indicates descriptors that were not calculated
      if( result.GetSize() != m_Descriptor->GetSizeOfFeatures())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }

      STORAGE.CopyValues( result);

      // perform the operation
      m_Op->PerformOnEach( STORAGE);
    }

  } // namespace descriptor
} // namespace bcl
