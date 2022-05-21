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

#ifndef BCL_DESCRIPTOR_TEMPLATE_H_
#define BCL_DESCRIPTOR_TEMPLATE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_named.h"
#include "bcl_descriptor_named_template.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Template
    //! @brief Allows run-time definition of templates. This is essentially Define for descriptors where there are still
    //!        be arguments
    //!
    //! @see @link example_descriptor_template.cpp @endlink
    //! @author mendenjl
    //! @date Jan 30, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Template :
      public BaseSequence< t_DataType, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      util::ObjectDataLabel m_Signature;  //!< Template Signature; what the user calls this template with
      util::ObjectDataLabel m_Definition; //!< Definition of the signature

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
      Template()
      {
      }

      //! @brief constructor from object data label
      explicit Template( const util::ObjectDataLabel &NAME_LABEL, const util::ObjectDataLabel &LABEL) :
        m_Signature( NAME_LABEL),
        m_Definition( LABEL)
      {
        ReadInitializerSuccessHook( LABEL, util::GetLogger());
      }

      //! @brief copy constructor
      Template *Clone() const
      {
        return new Template( *this);
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
        static const std::string s_alias( "Template");
        return s_alias;
      }

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return size_t( 0);
      }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    protected:

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
          "Create a descriptor template; which is a partial definition for a descriptor, "
          "with remaining user-defined arguments that can be set later"
        );
        parameters.AddInitializer
        (
          "signature",
          "Signature of the descriptor template, e.g. 3DA12(X). "
          "Parameters of the signature become (anonymous) arguments whenever the template is called, so the template can "
          " be called with 3DA12(Atom_Identity)",
          io::Serialization::GetAgent( &m_Signature)
        );
        parameters.AddInitializer
        (
          "",
          "Definition of the descriptor template. Each argument in the signature should normally appear in the "
          "definition, e.g. Template(signature=3DA12(X),3DA(property=X,steps=48,step size=0.25,temperature=100))",
          io::Serialization::GetAgent( &m_Definition)
        );
        return parameters;
      }

    private:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< t_ReturnType> &STORAGE)
      {
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        storage::Vector< std::string> args( m_Signature.GetNumberArguments());
        for( size_t arg_n( 0), n_args( args.GetSize()); arg_n < n_args; ++arg_n)
        {
          if( m_Signature.GetArgument( arg_n).GetNumberArguments())
          {
            ERR_STREAM << "Replacement strings in the signature cannot have arguments themselves";
            return false;
          }
          if( !m_Signature.GetArgument( arg_n).GetName().empty())
          {
            ERR_STREAM << "Replacement strings in the signature cannot have argument names";
            return false;
          }
          args( arg_n) = m_Signature.GetArgument( arg_n).GetValue();
        }
        // try to create a descriptor with the alias
        std::stringstream stream;
        NamedTemplate< t_DataType, t_ReturnType> named_instance( m_Definition, m_Signature.GetValue(), args);
        if( util::Enumerated< Base< t_DataType, t_ReturnType> >::HaveImplementationWithAlias( m_Signature.GetValue()))
        {
          typename util::Enumerated< Base< t_DataType, t_ReturnType> >::const_iterator
            itr_exists( util::Enumerated< Base< t_DataType, t_ReturnType> >::Find( m_Signature.GetValue()));
          if( itr_exists->second->GetSerializer().GetClassDescription() != named_instance.GetSerializer().GetClassDescription())
          {
            ERR_STREAM << "Conflicting template for " << m_Signature.GetValue()
                       << " already defined as: " << itr_exists->second->GetSerializer().GetClassDescription()
                       << " which conflicts with the new definition: "
                       << named_instance.GetSerializer().GetClassDescription();
            return false;
          }
          // implementation already present, skip
          return true;
        }
        util::Enumerated< Base< t_DataType, t_ReturnType> >::AddInstance( named_instance.Clone());
        return true;
      }

    }; // class Template

    BCL_EXPIMP_TEMPLATE template class BCL_API Template< char, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Template< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Template< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Template< chemistry::AtomConformationalInterface, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_TEMPLATE_H_
