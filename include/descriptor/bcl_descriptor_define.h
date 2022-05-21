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

#ifndef BCL_DESCRIPTOR_DEFINE_H_
#define BCL_DESCRIPTOR_DEFINE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_named.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Define
    //! @brief Allows run-time definition of new aliases for complex properties
    //!
    //! @see @link example_descriptor_define.cpp @endlink
    //! @author mendenjl
    //! @date Jun 20, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Define :
      public BaseSequence< t_DataType, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      std::string m_NewAlias; //!< Parameter name; allows the serialization framework to properly serialize this class,
                              //!< for which the LHS of the only parameter is variable
      util::ObjectDataLabel m_Label; //!< Label; lhs is the new alias, rhs is the full descriptor

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
      Define()
      {
      }

      //! @brief constructor from object data label
      explicit Define( const util::ObjectDataLabel &LABEL) :
        m_NewAlias( LABEL.GetName()),
        m_Label( LABEL)
      {
        ReadInitializerSuccessHook( LABEL, util::GetLogger());
      }

      //! @brief copy constructor
      Define *Clone() const
      {
        return new Define( *this);
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
        static const std::string s_alias( "Define");
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
        parameters.SetClassDescription( "Define an alias for a given descriptor");
        parameters.AddInitializer
        (
          m_NewAlias,
          "The LHS of the = sign will become an alias that can be used later (in the same file or command) to refer to the RHS. "
          "Existing aliases cannot be overridden. Accepts ",
          io::Serialization::GetAgent( &m_Label)
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

      //! @brief a function that derived classes can override to perform some action on a class before any data members
      //!        are read, e.g. resetting certain data members so that a post-read command can update parameters that
      //!        were not set
      //! @param SERIALIZER The serializer that will be used
      //! @param ERR_STREAM stream to write any errors encountered to
      //! @return result of any validation performed internally
      io::ValidationResult PreReadHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        if( SERIALIZER.GetNumberArguments() != 1)
        {
          return io::ValidationResult( false);
        }
        if( SERIALIZER.GetArgument( 0).GetName().empty())
        {
          if( SERIALIZER.GetArgument( 0).GetValue() == "help")
          {
            return io::e_Help;
          }
          else
          {
            ERR_STREAM << "syntax: NewAlias=Some(Other(Complicated),Descriptor=5); lhs missing!";
            return io::ValidationResult( false);
          }
        }
        m_NewAlias = SERIALIZER.GetArgument( 0).GetName();
        return io::ValidationResult( true);
      }

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        util::Implementation< Base< t_DataType, t_ReturnType> > impl;
        if( !impl.TryRead( m_Label, ERR_STREAM))
        {
          ERR_STREAM << "Failed to read label for " << m_NewAlias;
          return false;
        }
        // try to create a descriptor with the alias
        util::Implementation< Base< t_DataType, t_ReturnType> > impl_exist;
        std::stringstream stream;
        Named< t_DataType, t_ReturnType> named_instance( *impl, m_NewAlias, "User-defined");
        if
        (
          util::Enumerated< Base< t_DataType, t_ReturnType> >::HaveImplementationWithAlias( m_NewAlias)
          && impl_exist.TryRead( util::ObjectDataLabel( m_NewAlias), stream)
        )
        {
          if( impl_exist.GetSerializer().GetClassDescription() != named_instance.GetSerializer().GetClassDescription())
          {
            ERR_STREAM << "Conflicting alias for " << m_NewAlias
                       << " already defined as: " << impl_exist.GetSerializer().GetClassDescription()
                       << " which conflicts with the new definition: "
                       << named_instance.GetSerializer().GetClassDescription();
            return false;
          }
          // implementation already present, skip
          return true;
        }
        else if( util::Enumerated< Base< t_DataType, t_ReturnType> >::HaveImplementationWithAlias( m_NewAlias))
        {
          ERR_STREAM << "Have impl for " << m_NewAlias << ", but it requires parameters";
          return true;
        }
        util::Enumerated< Base< t_DataType, t_ReturnType> >::AddInstance( named_instance.Clone());
        return true;
      }

    }; // class Define

    BCL_EXPIMP_TEMPLATE template class BCL_API Define< char, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Define< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Define< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Define< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Define< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Define< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Define< chemistry::AtomConformationalInterface, char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_DEFINE_H_
