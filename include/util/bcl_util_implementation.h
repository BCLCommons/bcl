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

#ifndef BCL_UTIL_IMPLEMENTATION_H_
#define BCL_UTIL_IMPLEMENTATION_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_enumerated.h"
#include "bcl_util_functional_type.h"
#include "bcl_util_implementation_interface.h"
#include "command/bcl_command_guesser.h"
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Implementation
    //! @brief Wrapper for an implementation of an interface that uses data labels
    //!
    //! @see @link example_util_implementation.cpp @endlink
    //! @author mendenjl
    //! @date Jan 5, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Interface>
    class Implementation :
      virtual public ImplementationInterface
    {

    private:

    //////////
    // data //
    //////////

      t_Interface *m_Implementation; //!< A pointer to the implementation of the interface

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Implementation() :
        m_Implementation( NULL)
      {
      }

      //! @brief copy constructor
      Implementation( const Implementation &ORIGINAL) :
        m_Implementation( ORIGINAL.m_Implementation ? ORIGINAL.m_Implementation->Clone() : NULL)
      {
      }

      //! @brief constructor from a particular implementation
      Implementation( const t_Interface &ORIGINAL) :
        m_Implementation( ORIGINAL.Clone())
      {
      }

      //! @brief constructor from a pointer a new implementation
      Implementation( t_Interface *const &ORIGINAL) :
        m_Implementation( ORIGINAL)
      {
      }

      //! @brief constructor from DATA_LABEL string, asserts on failure
      Implementation( const std::string &DATA_LABEL) :
        m_Implementation( NULL)
      {
        AssertRead( ObjectDataLabel( DATA_LABEL));
      }

      //! @brief constructor from DATA_LABEL, asserts on failure
      Implementation( const ObjectDataLabel &DATA_LABEL) :
        m_Implementation( NULL)
      {
        AssertRead( DATA_LABEL);
      }

      //! @brief constructor from DATA_LABEL, resets on failure
      //! @param DATA_LABEL label to use to construct the object
      //! @param ERR_STREAM stream to wrtie errors to
      Implementation( const ObjectDataLabel &DATA_LABEL, std::ostream &ERR_STREAM) :
        m_Implementation( NULL)
      {
        if( !TryRead( DATA_LABEL, ERR_STREAM))
        {
          Reset();
        }
      }

      //! @brief virtual destructor
      ~Implementation()
      {
        Reset();
      }

      //! @brief Clone the implementation
      //! @return pointer to new Implementation
      virtual Implementation *Clone() const
      {
        return new Implementation( *this);
      }

      //! @brief guaranteed deep-copy of the implementation
      virtual Implementation HardCopy() const
      {
        return Implementation( GetLabel( true));
      }

      //! @brief return an empty copy of the implementation
      //! @return pointer to an empty copy of the implementation
      virtual ImplementationInterface *Empty() const
      {
        return new Implementation();
      }

      //! assignment operator
      Implementation &operator =( const Implementation &ORIGINAL)
      {
        if( m_Implementation != ORIGINAL.m_Implementation)
        {
          Reset();
          if( ORIGINAL.m_Implementation != NULL)
          {
            m_Implementation = ORIGINAL.m_Implementation->Clone();
          }
        }
        return *this;
      }

      //! @brief assign to a different implementation via a string
      Implementation &operator =( const std::string &DATA_LABEL)
      {
        AssertRead( ObjectDataLabel( DATA_LABEL));
        return *this;
      }

      //! @brief assign to a different implementation via a data label
      Implementation &operator =( const ObjectDataLabel &DATA_LABEL)
      {
        AssertRead( DATA_LABEL);
        return *this;
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

      //! @brief get name of the property
      //! @return name of the property
      const std::string &GetAlias() const
      {
        return IsDefined() ? m_Implementation->GetAlias() : GetUndefinedString();
      }

      //! @brief get the label containing only initialization parameters
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const
      {
        return IsDefined() ? GetCompleteSerializer().GetLabel( WITH_DATA) : ObjectDataLabel( "", GetUndefinedString());
      }

      //! @brief determine the type of value that can be parsed
      DataType::Type GetSerializedType() const
      {
        return DataType::e_DynamicObject;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief remove any existing implementation
      void Reset()
      {
        if( IsDefined())
        {
          delete m_Implementation;
          m_Implementation = NULL;
        }
      }

      //! @brief test whether the implementation is defined
      //! @return true if implementation is defined
      bool IsDefined() const
      {
        return m_Implementation;
      }

    ///////////////
    // operators //
    ///////////////

      //! overload operator -> to return const pointer on object t_Interface
      const t_Interface *operator ->() const
      {
        BCL_Assert( IsDefined(), "Implementation was never given for " + GetStaticClassName< t_Interface>());
        return m_Implementation;
      }

      //! overload operator * to return const reference of object t_Interface
      const t_Interface &operator *() const
      {
        BCL_Assert( IsDefined(), "Implementation was never given for " + GetStaticClassName< t_Interface>());
        return *m_Implementation;
      }

      //! overload operator -> to return const pointer on object t_Interface
      t_Interface *operator ->()
      {
        BCL_Assert( IsDefined(), "Implementation was never given for " + GetStaticClassName< t_Interface>());
        return m_Implementation;
      }

      //! overload operator * to return const reference of object t_Interface
      t_Interface &operator *()
      {
        BCL_Assert( IsDefined(), "Implementation was never given for " + GetStaticClassName< t_Interface>());
        return *m_Implementation;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief a function that derived classes can override to perform some action on a class before any data members
      //!        are read, e.g. resetting certain data members so that a post-read command can update parameters that
      //!        were not set
      //! @param SERIALIZER The serializer that will be used
      //! @param ERR_STREAM stream to write any errors encountered to
      //! @return result of any validation performed internally
      io::ValidationResult PreReadHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM);

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        if( IsDefined())
        {
          return m_Implementation->ReadInitializerSuccessHook( SERIALIZER, ERR_STREAM);
        }
        return true;
      }

      //! @brief a function that derived classes can override to take additional actions whenever Read is called successfully
      //!        AND any data members are specified for the given class
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadDataSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        if( IsDefined())
        {
          return m_Implementation->ReadDataSuccessHook( SERIALIZER, ERR_STREAM);
        }
        return true;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        return IsDefined() ? m_Implementation->GetCompleteSerializer() : io::Serializer();
      }

      //! @brief Get a set of all class names used by the serializer. Useful for introspection
      //! @param TYPES set to insert type names into
      //! @param INCLUDE_OPTIONAL true to also count optional members
      //! @param INCLUDE_DATA true to also include data-containing members
      void InsertDataTypes
      (
        storage::Map< std::string, size_t> &TYPES,
        const bool &INCLUDE_OPTIONAL = true,
        const bool &INCLUDE_DATA = false,
        const size_t &MAX_DEPTH = size_t( -1)
      ) const
      {
        if( IsDefined())
        {
          m_Implementation->InsertDataTypes( TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA);
        }
        else
        {
          SerializableInterface::InsertDataTypes( TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA);
        }
      }

      //! @brief write this instance of the enumerated class to the stream with its default parameters
      //! @param STREAM stream to write to
      //! @param INDENT indentation to give each line
      std::ostream &WriteHelp( std::ostream &STREAM, const size_t INDENT = 0) const;

      //! @brief write this instance of the enumerated class to the stream with its default parameters
      //! @param STREAM stream to write to
      io::FixedLineWidthWriter &WriteHelp( io::FixedLineWidthWriter &STREAM) const;

      //! @brief write all the enumerated member's help to the stream
      //! @param STREAM stream to write to
      static io::FixedLineWidthWriter &WriteInstancesHelp( io::FixedLineWidthWriter &STREAM);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief static option to check whether the serializer could create a valid implementation
      //! @param SERIALIZER The serializer to test
      //! @return true if the label could be used to successfully initialize this object
      static std::pair< Implementation, bool> CreateIfPossible( const ObjectDataLabel &SERIALIZER)
      {
        std::stringstream err;
        std::pair< Implementation, bool> impl_success;
        impl_success.second = impl_success.first.TryRead( SERIALIZER, err);
        return impl_success;
      }

      //! @brief unset whether the help was already displayed for this implementation
      //! This is useful if, while testing an example, the help is written out, and other examples
      //! may have written out help for implementations used by this example.  In that case, it is useful to reset
      //! the help for all implementations that may be called
      static void ResetHaveDisplayedHelp()
      {
        GetHaveDisplayedHelp() = false;
      }

      //! @brief set whether the help was already displayed for this implementation
      //! This is useful for preventing an implementation with a high number of options from overwhelming the output
      //! of another implementation
      static void SetHaveDisplayedHelp()
      {
        GetHaveDisplayedHelp() = true;
      }

      //! @brief typedef for a function pointer to a help writer function
      typedef io::FixedLineWidthWriter &( *HelpWriterFunctionPtr)
      (
        io::FixedLineWidthWriter &,
        const storage::Map< std::string, OwnPtr< t_Interface> > &,
        bool //!< Indicates whether to display full instance help (true), or brief help (false)
      );

      //! @brief set the function used to write help for the enumerated instances in util::Implementation. This function
      //!        can be set externally to e.g. group different implementations by some criteria
      //! @param FUNCTION the function to set the help writer to
      //! @return true unless the help writer has already been set to a different function (then assert);
      static bool SetHelpWriter( HelpWriterFunctionPtr FUNCTION)
      {
        BCL_Assert
        (
          GetHelpWriter() == NULL || GetHelpWriter() == FUNCTION,
          "Help writer for " + GetStaticClassName< t_Interface>() + " was already set!"
        );
        GetHelpWriter() = FUNCTION;
        return true;
      }

      //! @brief set the name used to describe undefined instances of t_Interface
      //! @param FUNCTION the function to set the help writer to
      //! @return true unless the help writer has already been set to a different function (then assert);
      static bool SetUndefinedInstanceName( const std::string &NAME)
      {
        GetUndefinedString() = NAME;
        return true;
      }

      //! @brief get the name used to describe undefined instances of t_Interface
      //! @return undefined instance name
      static const std::string &GetUndefinedInstanceName()
      {
        return GetUndefinedString();
      }

    private:

      //! @brief get the function used to write help for the t_Interface
      //! @return the function used to write help for the t_Interface; null if the default implementation should be used
      static HelpWriterFunctionPtr &GetHelpWriter()
      {
        static HelpWriterFunctionPtr s_function( NULL);
        return s_function;
      }

      //! @brief return a static bool that tells whether or not help was already written for this implementation
      //! @return changable reference to bool indicating whether help was already displayed
      static bool &GetHaveDisplayedHelp()
      {
        static bool s_have_displayed_help( false);
        return s_have_displayed_help;
      }

      //! @brief get the string for an undefined object of this interface
      static std::string &GetUndefinedString()
      {
        // Get the undefined string, which is used for the name if no implementation is available for a class
        static std::string s_undefined( GetStaticClassName< t_Interface>());
        return s_undefined;
      }

      //! @brief a function to convert an implementation into a shared pointer without an expensive clone
      //! @param IMPL the implementation to convert; after this call, it will be undefined
      template< typename t_Iface>
      friend ShPtr< t_Iface> ConvertToShPtr( Implementation< t_Iface> &IMPL);

      //! @brief release ownership of the internal pointer, only used to transfer ownership to a ShPtr
      t_Interface *const ReleaseOwnership()
      {
        t_Interface *old_ptr( m_Implementation);
        m_Implementation = NULL;
        return old_ptr;
      }

    }; // template class Implementation

    //! @brief a function to convert an implementation into a shared pointer without an expensive clone
    //! @param IMPL the implementation to convert; after this call, it will be undefined
    template< typename t_Iface>
    ShPtr< t_Iface> ConvertToShPtr( Implementation< t_Iface> &IMPL)
    {
      return ShPtr< t_Iface>( IMPL.ReleaseOwnership());
    }

    //! @brief a function that derived classes can override to perform some action on a class before any data members
    //!        are read, e.g. resetting certain data members so that a post-read command can update parameters that
    //!        were not set
    //! @param SERIALIZER The serializer that will be used
    //! @param ERR_STREAM stream to write any errors encountered to
    //! @return result of any validation performed internally
    template< typename t_Interface>
    io::ValidationResult Implementation< t_Interface>::PreReadHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      // check whether the implementation is the same, if so, no need to reset
      if( IsDefined())
      {
        if( m_Implementation->GetAlias() == SERIALIZER.GetValue())
        {
          // already have the correct implementation
          return io::ValidationResult( io::e_Allowed);
        }
        Reset();
      }

      // case: undefined instance
      if( SERIALIZER.GetValue() == GetUndefinedString())
      {
        Reset();
        return io::ValidationResult( io::e_Complete);
      }

      // look for the alias in the map
      typename Enumerated< t_Interface>::const_iterator
        itr( Enumerated< t_Interface>::GetInstanceMap().Find( SERIALIZER.GetValue())),
        itr_end( Enumerated< t_Interface>::End());

      // case: found instance
      if( itr != itr_end)
      {
        m_Implementation = itr->second->Clone();
        return m_Implementation->PreReadHook( SERIALIZER, ERR_STREAM);
      }

      // case: requested help
      if( SERIALIZER.GetValue() == io::ValidationResult::GetHelpString())
      {
        ERR_STREAM << "Implementations of " << GetStaticClassName< t_Interface>() << ":\n";
        WriteHelp( ERR_STREAM, 0);
        return io::ValidationResult( io::e_Help);
      }

      // case: default implementation
      if( !SERIALIZER.IsEmpty())
      {
        itr = Enumerated< t_Interface>::GetDefaultImplementation();
        if( itr != itr_end)
        {
          // use the default implementation
          m_Implementation = itr->second->Clone();
          return m_Implementation->PreReadHook( SERIALIZER, ERR_STREAM);
        }
      }
      else
      {
        // empty serializer, presumably undefined
        Reset();
        return io::ValidationResult( io::e_Complete);
      }

      // if the alias was not found and no default has been given write out the available implementations to the user
      command::Guesser::GetDefaultGuesser().WriteGuesses
      (
        SERIALIZER.GetValue(),
        Enumerated< t_Interface>::GetInstanceMap().GetKeysAsVector(),
        ERR_STREAM,
        "implementation"
      );
      ERR_STREAM << "^^^ In serializer: " << SERIALIZER.ToString() << '\n';
      return io::ValidationResult( io::e_Invalid);
    }

    //! @brief write this instance of the enumerated class to the stream with its default parameters
    //! @param STREAM stream to write to
    //! @param INDENT indentation to give each line
    template< typename t_Interface>
    std::ostream &Implementation< t_Interface>::WriteHelp( std::ostream &STREAM, const size_t INDENT) const
    {
      // allow a little extra space at the end of lines; often this help is written by a class that does not call it
      // with the indent properly set
      io::FixedLineWidthWriter writer( 2, GetLogger().GetMaxLineWidth() - 10);
      writer.SetBclIndent( INDENT);
      this->WriteHelp( writer);
      STREAM << writer.String();
      return STREAM;
    }

    //! @brief write this instance of the enumerated class to the stream with its default parameters
    //! @param STREAM stream to write to
    template< typename t_Interface>
    io::FixedLineWidthWriter &Implementation< t_Interface>::WriteHelp( io::FixedLineWidthWriter &STREAM) const
    {
      if( IsDefined())
      {
        return m_Implementation->WriteHelp( STREAM);
      }

      if( !GetHaveDisplayedHelp())
      {
        // prevent recursive output
        GetHaveDisplayedHelp() = true;
      }
      else
      {
        STREAM << "choose any implementation of " << GetUndefinedString() << " (already listed)";
        return STREAM;
      }

      if( GetHelpWriter() != NULL)
      {
        return ( *GetHelpWriter())( STREAM, Enumerated< t_Interface>::GetInstanceMap(), false);
      }

      // write a short description
      STREAM << "choose any implementation of " << GetUndefinedString() << ":";
      STREAM.AddIndent( 2);

      storage::Vector< storage::Vector< SiPtr< const t_Interface> > > interfaces_by_type
      (
        ( size_t( FunctionalType::s_NumberTypes))
      );

      // write the default data label of each instance out to STREAM
      for
      (
        typename storage::Map< std::string, OwnPtr< t_Interface> >::const_iterator
          itr( Enumerated< t_Interface>::GetInstanceMap().Begin()),
          itr_end( Enumerated< t_Interface>::GetInstanceMap().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->first != "")
        {
          interfaces_by_type( itr->second->InferFunctionalType( this->GetClassIdentifier())).PushBack( *itr->second);
        }
      }
      for( size_t type_num( 0); type_num < FunctionalType::s_NumberTypes; ++type_num)
      {
        const storage::Vector< SiPtr< const t_Interface> > &impls_of_type( interfaces_by_type( type_num));
        if( impls_of_type.IsEmpty())
        {
          continue;
        }
        STREAM.NewLineIndent();
        STREAM.NewLineIndent();
        STREAM << FunctionalType::GetTypeName( FunctionalType::Type( type_num));
        STREAM.AddIndent( 2);
        STREAM.NewLineIndent();
        for
        (
          typename storage::Vector< SiPtr< const t_Interface> >::const_iterator
            itr( impls_of_type.Begin()), itr_end( impls_of_type.End());
          itr != itr_end;
          ++itr
        )
        {
          Implementation< t_Interface>( ( *itr)->Clone()).GetCompleteSerializer().WriteBriefHelp( STREAM);
          STREAM.NewLineIndent();
        }
        STREAM.PopIndent();
      }
      STREAM.PopIndent();
      if( Enumerated< t_Interface>::GetDefaultImplementation() != Enumerated< t_Interface>::GetInstanceMap().End())
      {
        STREAM << "\nOther strings will be interpreted as follows:\n";
        Enumerated< t_Interface>::GetDefaultImplementation()->second->GetCompleteSerializer().WriteBriefHelp( STREAM);
      }
      return STREAM;
    }

    //! @brief write all the enumerated member's help to the stream
    //! @param STREAM stream to write to
    template< typename t_Interface>
    io::FixedLineWidthWriter &Implementation< t_Interface>::WriteInstancesHelp( io::FixedLineWidthWriter &STREAM)
    {
      if( !GetHaveDisplayedHelp())
      {
        // prevent recursive output
        GetHaveDisplayedHelp() = true;
      }
      else
      {
        STREAM << "choose any implementation of " << GetUndefinedString() << " (already listed)";
        return STREAM;
      }

      if( GetHelpWriter() != NULL)
      {
        return ( *GetHelpWriter())( STREAM, Enumerated< t_Interface>::GetInstanceMap(), true);
      }

      // write a short description
      STREAM << "choose any implementation of " << GetUndefinedString() << ":\n";

      STREAM.AddIndent( 2);
      // write the default data label of each instance out to STREAM
      for
      (
        typename storage::Map< std::string, OwnPtr< t_Interface> >::const_iterator
          itr( Enumerated< t_Interface>::GetInstanceMap().Begin()),
          itr_end( Enumerated< t_Interface>::GetInstanceMap().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->first != "")
        {
          STREAM.WriteOnOneLine( "* " + itr->second->GetAlias() + " : ");
          itr->second->WriteHelp( STREAM);
        }
      }
      STREAM.PopIndent();

      if( Enumerated< t_Interface>::GetDefaultImplementation() != Enumerated< t_Interface>::GetInstanceMap().End())
      {
        STREAM << "Other strings will be interpreted as follows:\n";
        Enumerated< t_Interface>::GetDefaultImplementation()->second->GetCompleteSerializer().WriteBriefHelp( STREAM);
      }
      return STREAM;
    }

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_IMPLEMENTATION_H_

