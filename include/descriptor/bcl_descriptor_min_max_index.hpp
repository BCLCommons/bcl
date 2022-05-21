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
#include "bcl_descriptor_min_max_index.h"
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

    //! @brief default constructor, accepts bool of whether to auto-compute mean
    template< typename t_DataType>
    MinMaxIndex< t_DataType>::MinMaxIndex( const bool &MAX_INDEX) :
      m_ComputeMax( MAX_INDEX)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    MinMaxIndex< t_DataType> *MinMaxIndex< t_DataType>::Clone() const
    {
      return new MinMaxIndex( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &MinMaxIndex< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &MinMaxIndex< t_DataType>::GetAlias() const
    {
      static const std::string s_NameMax( "MaxIndex"), s_NameMin( "MinIndex");

      return m_ComputeMax ? s_NameMax : s_NameMin;
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_DataType>
    size_t MinMaxIndex< t_DataType>::GetNormalDimension() const
    {
      return m_Property->GetType().GetDimension();
    }

    //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t MinMaxIndex< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

    //! @brief set properties to be logarithmized
    //! @param PROPERTY SmallMoleculeProperty
    template< typename t_DataType>
    void MinMaxIndex< t_DataType>::SetProperty( const util::Implementation< Base< t_DataType, float> > &PROPERTY)
    {
      m_Property = PROPERTY;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    template< typename t_DataType>
    const util::Implementation< Base< t_DataType, float> > &MinMaxIndex< t_DataType>::GetProperty() const
    {
      return m_Property;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer MinMaxIndex< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the index of "
        + std::string( m_ComputeMax ? "maximum" : "minimum") + " value for the given descriptor"
      );

      parameters.AddInitializer
      (
        "",
        "descriptor for which to find the index with the " + std::string( m_ComputeMax ? "maximum" : "minimum")
        + " corresponding value",
        io::Serialization::GetAgent( &m_Property)
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
    iterate::Generic< Base< t_DataType, float> > MinMaxIndex< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Property, &m_Property + 1);
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    template< typename t_DataType>
    Type::Symmetry MinMaxIndex< t_DataType>::GetSymmetry() const
    {
      return m_Property->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    template< typename t_DataType>
    bool MinMaxIndex< t_DataType>::ConsiderRepeatedElements() const
    {
      return m_Property->GetType().ConsiderRepeatedObjects();
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void MinMaxIndex< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // get a reference to the underlying property
      linal::VectorConstReference< float> ref( m_Property->operator ()( ITR));

      // compute the max or min index, whichever was desired by the user
      if( m_ComputeMax)
      {
        STORAGE( 0) = long( std::max_element( ref.Begin(), ref.End()) - ref.Begin());
      }
      else
      {
        STORAGE( 0) = long( std::min_element( ref.Begin(), ref.End()) - ref.Begin());
      }
    }

  } // namespace descriptor
} // namespace bcl
