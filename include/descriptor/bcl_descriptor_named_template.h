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

#ifndef BCL_DESCRIPTOR_NAMED_TEMPLATE_H_
#define BCL_DESCRIPTOR_NAMED_TEMPLATE_H_

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
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NamedTemplate
    //! @brief Allows definition of a name for a particular instance of a descriptor
    //!
    //! @see @link example_descriptor_template.cpp @endlink
    //! @author mendenjl
    //! @date Jan 30, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class NamedTemplate :
      public Base< t_DataType, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, t_ReturnType> > m_Descriptor; //!< descriptor to calculate internally

      util::ObjectDataLabel m_Definition; //!< Definition of the function to define
      storage::Vector< std::string> m_ArgumentNames; //!< Names the user provided for the arguments
      storage::Vector< util::ObjectDataLabel> m_ArgumentValues; //!< Values the user provided for the arguments

      //! Alias for this descriptor
      std::string m_Alias;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      NamedTemplate() :
        m_Descriptor(),
        m_ArgumentValues(),
        m_Alias()
      {
      }

      //! @brief constructor from member data
      NamedTemplate
      (
        const util::ObjectDataLabel &DEFINITION,
        const std::string &ALIAS,
        const storage::Vector< std::string> &ARGUMENTS
      ) :
        m_Descriptor(),
        m_Definition( DEFINITION),
        m_ArgumentNames( ARGUMENTS),
        m_ArgumentValues( ARGUMENTS.GetSize()),
        m_Alias( ALIAS)
      {
      }

      //! @brief copy constructor
      NamedTemplate *Clone() const
      {
        return new NamedTemplate( *this);
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
          " Equivalent to "
          + m_Definition.ToString()
          + " with positional arguments "
          + util::Join( ", ", m_ArgumentNames)
        );
        if( m_ArgumentNames.GetSize())
        {
          parameters.AddInitializer
          (
            "",
            "positional arguments",
            io::Serialization::GetAgentWithSizeLimits
            (
              &m_ArgumentValues,
              m_ArgumentValues.GetSize(),
              m_ArgumentValues.GetSize()
            )
          );
        }

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
        return m_Descriptor.IsDefined() && m_Descriptor->InjectDimensions();
      }

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return m_Descriptor.IsDefined() && m_Descriptor->DimensionIsWellDefined();
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
        util::ObjectDataLabel def_with_replacements( m_Definition);
        for
        (
          size_t argument_number( 0), n_args( m_ArgumentNames.GetSize());
          argument_number < n_args;
          ++argument_number
        )
        {
          def_with_replacements.ReplaceValue
          (
            m_ArgumentNames( argument_number),
            m_ArgumentValues( argument_number).ToString()
          );
        }
        return m_Descriptor.TryRead( def_with_replacements, ERR_STREAM);
      }

    }; // class NamedTemplate

    BCL_EXPIMP_TEMPLATE template class BCL_API NamedTemplate< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NamedTemplate< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NamedTemplate< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NamedTemplate< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NamedTemplate< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NamedTemplate< chemistry::AtomConformationalInterface, char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_NAMED_TEMPLATE_H_
