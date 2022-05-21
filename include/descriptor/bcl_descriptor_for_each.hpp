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
#include "bcl_descriptor_for_each.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_iterator.h"
#include "bcl_descriptor_named.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"
#include "model/bcl_model_feature_label_set.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief virtual copy constructor
    template< typename t_DataType, typename t_ReturnType>
    ForEach< t_DataType, t_ReturnType> *ForEach< t_DataType, t_ReturnType>::Clone() const
    {
      return new ForEach( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &ForEach< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &ForEach< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_name( "ForEach");
      return s_name;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool ForEach< t_DataType, t_ReturnType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // this is calle really just to clear out the
      Combine< t_DataType, t_ReturnType>::operator =( Combine< t_DataType, t_ReturnType>());

      for
      (
        storage::Vector< util::ObjectDataLabel>::const_iterator itr( m_Values.Begin()), itr_end( m_Values.End());
        itr != itr_end;
        ++itr
      )
      {
        util::ObjectDataLabel templ_copy( m_TemplateLabel);
        templ_copy.ReplaceValue( m_Variable, itr->ToNamedString(), size_t( -1));
        util::Implementation< Base< t_DataType, t_ReturnType> > new_desc;
        if( !new_desc.TryRead( templ_copy, ERR_STREAM))
        {
          return false;
        }
        this->PushBack( new_desc);
      }
      for
      (
        storage::Vector< util::ObjectDataLabel>::const_iterator
          itr_desc( m_ValuesDescriptors.Begin()), itr_desc_end( m_ValuesDescriptors.End());
        itr_desc != itr_desc_end;
        ++itr_desc
      )
      {
        if( itr_desc->IsScalar())
        {
          util::Implementation< Base< t_DataType, t_ReturnType> > new_desc;
          typename util::Enumerated< Base< t_DataType, t_ReturnType> >::const_iterator
            itr_vd( util::Enumerated< Base< t_DataType, t_ReturnType> >::Find( itr_desc->GetValue()));
          if( itr_vd == util::Enumerated< Base< t_DataType, t_ReturnType> >::End())
          {
            ERR_STREAM << *itr_desc << " is not a previously-defined list or regular base descriptor. descriptors(X) may "
                       << "only be used if X is a basic descriptor (no arguments) or a list previously defined like Define(YZ=Combine(Y,Z)),";
            return false;
          }
          bool result( AddDescriptorArguments( *itr_vd->second, *itr_desc, ERR_STREAM));
          if( !result)
          {
            return false;
          }
          continue;
        }
        util::Implementation< Base< t_DataType, t_ReturnType> > new_desc;
        if( !new_desc.TryRead( *itr_desc, ERR_STREAM))
        {
          return false;
        }
        bool result( AddDescriptorArguments( *new_desc, *itr_desc, ERR_STREAM));
        if( !result)
        {
          return false;
        }

      }
      BCL_MessageDbg
      (
        "ForEach expanded to: " + ( Combine< t_DataType, t_ReturnType>::GetSerializer().GetLabel().ToString( 120, 0, 1))
      );
      Combine< t_DataType, t_ReturnType>::ReadInitializerSuccessHook( LABEL, ERR_STREAM);

      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer ForEach< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Allows creation of a set of descriptors by substituting a user-specified set of values in for a specific "
        "parameter"
      );
      parameters.AddInitializer
      (
        "template",
        "Template descriptor; should contain one or more variables to substitute e.g. "
        "template=3DA(property=X,step size=0.25,number steps=48)",
        io::Serialization::GetAgent( &m_TemplateLabel)
      );
      parameters.AddInitializer
      (
        "variable",
        "Variable name; this name should be unique from any other descriptor used inside this label",
        io::Serialization::GetAgent( &m_Variable)
      );
      parameters.AddOptionalInitializer
      (
        "values",
        "Values that the variable will take on for this descriptor. Can be omitted if descriptors parameter is given",
        io::Serialization::GetAgentWithSizeLimits( &m_Values, size_t( 0), size_t( -1))
      );
      parameters.AddOptionalInitializer
      (
        "descriptors",
        "If the replacement values are descriptors of " + this->GetObjectName()
        + ", then the descriptors can be taken from any previously-defined "
        "list, or given directly. E.g., if Define(YZ=Combine(Y,Z)), then ForEach(template=Foo(Bar),variable=Bar,descriptors(YZ,T)) expands "
        "to: Combine(Foo(Y),Foo(Z),Foo(T))",
        io::Serialization::GetAgent( &m_ValuesDescriptors)
      );

      return parameters;
    } // GetParameters

    //! @brief interpret the variables in the given combine as arguments to be used in the template
    //! @param BASE the base to use. If it is a defined list or list directly, the internal descriptors will be used instead
    //! @param LABEL label used to initialize the base
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool ForEach< t_DataType, t_ReturnType>::AddDescriptorArguments
    (
      const Base< t_DataType, t_ReturnType> &BASE,
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // test whether the type can be dynamic cast to the named descriptor
      // perform the dynamic cast here directly because if SiPtr performs the cast, it will issue an unproductive
      // warning
      const Named< t_DataType, t_ReturnType> *named_ptr =
        dynamic_cast< const Named< t_DataType, t_ReturnType> *>( &BASE);
      const Combine< t_DataType, t_ReturnType> *combine_ptr =
        dynamic_cast< const Combine< t_DataType, t_ReturnType> *>( named_ptr ? &named_ptr->GetInternalDescriptor() : &BASE);
      if( !combine_ptr)
      {
        // not a named descriptor. Use the descriptor directly
        util::ObjectDataLabel templ_copy( m_TemplateLabel);
        templ_copy.ReplaceValue( m_Variable, LABEL.ToString(), size_t( -1));
        util::Implementation< Base< t_DataType, t_ReturnType> > new_desc;
        if( !new_desc.TryRead( templ_copy, ERR_STREAM))
        {
          return false;
        }
        this->PushBack( new_desc);
        return true;
      }
      for
      (
        typename Combine< t_DataType, t_ReturnType>::const_iterator
          itr( combine_ptr->Begin()), itr_end( combine_ptr->End());
        itr != itr_end;
        ++itr
      )
      {
        util::ObjectDataLabel templ_copy( m_TemplateLabel);
        templ_copy.ReplaceValue( m_Variable, itr->GetString(), size_t( -1));
        util::Implementation< Base< t_DataType, t_ReturnType> > new_desc;
        if( !new_desc.TryRead( templ_copy, ERR_STREAM))
        {
          return false;
        }
        this->PushBack( new_desc);
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
