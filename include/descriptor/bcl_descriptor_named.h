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

#ifndef BCL_DESCRIPTOR_NAMED_H_
#define BCL_DESCRIPTOR_NAMED_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Named
    //! @brief Allows definition of a name for a particular instance of a descriptor
    //!
    //! @see @link example_descriptor_named.cpp @endlink
    //! @author mendenjl
    //! @date Jun 25, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Named :
      public Base< t_DataType, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, t_ReturnType> > m_Descriptor; //!< descriptor to calculate internally

      //! Alias for this descriptor
      std::string m_Alias;

      //! description for this descriptor
      std::string m_Description;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param REPLACE_MEAN_WITH_NON_NAN_VALUE true to ignore the replacement property and just use the mean non-nan
      //!        value
      Named() :
        m_Descriptor(),
        m_Alias(),
        m_Description()
      {
      }

      //! @brief constructor from descriptor, desired alias, and description
      Named
      (
        const Base< t_DataType, t_ReturnType> &DESCRIPTOR,
        const std::string &ALIAS,
        const std::string &DESCRIPTION
      ) :
        m_Descriptor( DESCRIPTOR),
        m_Alias( ALIAS),
        m_Description( DESCRIPTION)
      {
      }

      //! @brief copy constructor
      Named *Clone() const
      {
        return new Named( *this);
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

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const
      {
        return m_Alias;
      }

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const
      {
        return m_Descriptor->GetSymmetry();
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const
      {
        return m_Descriptor->ConsiderRepeatedElements();
      }

      //! @brief return the internal descriptor
      //! @return the internal descriptor
      const Base< t_DataType, t_ReturnType> &GetInternalDescriptor() const
      {
        return *m_Descriptor;
      }

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to forward calls from SetObject on to all internal implementations
      //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
      //! that function is called
      iterate::Generic< Base< t_DataType, t_ReturnType> > GetInternalDescriptors()
      {
        // if auto-replacing with the mean of the non-nan values, then do not return the replacement property
        return
          iterate::Generic< Base< t_DataType, t_ReturnType> >( &m_Descriptor, &m_Descriptor + 1);
      }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer parameters;
        parameters.SetClassDescription
        (
          m_Description + std::string( m_Description[ m_Description.size() - 1] != '.' ? "." : "")
          + " Alias for : " + m_Descriptor.GetString()
        );

        return parameters;
      }

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_Descriptor->GetSizeOfFeatures();
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        return m_Descriptor->GetType().GetDimension();
      }

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! @note users do not normally need to override this function
      bool InjectDimensions() const
      {
        return true;
      }

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return m_Descriptor->DimensionIsWellDefined();
      }

    private:

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< t_ReturnType> &STORAGE
      )
      {
        STORAGE.CopyValues( m_Descriptor->operator()( ITR));
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        return m_Descriptor.ReadInitializerSuccessHook( LABEL, ERR_STREAM);
      }

    }; // class Named

    BCL_EXPIMP_TEMPLATE template class BCL_API Named< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Named< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Named< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Named< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Named< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Named< chemistry::AtomConformationalInterface, char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_NAMED_H_
