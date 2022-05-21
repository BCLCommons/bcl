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

#ifndef BCL_COMMAND_FLAG_INTERFACE_H_
#define BCL_COMMAND_FLAG_INTERFACE_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_interface.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FlagInterface
    //! @brief The interface class for command line flags
    //!
    //! @remarks example unnecessary
    //! @author heinzes1, woetzen, karakam
    //! @date 10/12/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FlagInterface :
      public util::ObjectInterface
    {

    protected:

      //! true if the flag is unused now; maintained only to allow backwards compatibility with older scripts but has no
      //! influence on the outcome of the program
      bool m_IsUnused;

    public:

    //////////////
    // typedefs //
    //////////////

      //! typedef for the function type that will be called after setting or resetting the flag
      typedef void ( *t_Signal)();

      //! @brief default constructor
      FlagInterface( const bool &DEPRECATED = false) :
        m_IsUnused( DEPRECATED)
      {
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns name of flag
      //! @return name of flag as std::string
      virtual const std::string &GetName() const = 0;

      //! @brief set name - only possible if there was no name given yet
      //! @param NAME - what you want to name the flag
      virtual void SetName( const std::string &NAME) = 0;

      //! @brief returns description of flag
      //! @return description of flag as string
      virtual const std::string &GetDescription() const = 0;

      //! @brief set flag to true if it was false
      virtual void SetFlag() = 0;

      //! @brief set flag to false if it was true
      virtual void UnsetFlag() = 0;

      //! @brief reset the flag
      //! @detail resets all internal parameters, then removes all internally held parameters
      virtual void ResetFlag() = 0;

      //! @brief returns true if was set in command line
      //! @return whether this flag was set in the command line
      virtual bool GetFlag() const = 0;

      //! @brief returns the function to be called whenever this flag is updated
      //! @return the function to be called whenever this flag is updated
      virtual t_Signal GetSignal() const = 0;

      //! @brief set that the flag is unused
      void SetIsUnused()
      {
        m_IsUnused = true;
      }

      //! @brief returns true if the flag is unused and may be removed from a future release
      bool IsUnused() const
      {
        return m_IsUnused;
      }

      //! @brief returns the number of parameters
      //! @return number of parameters after flag
      size_t GetSize() const
      {
        return GetParameterList().GetSize();
      }

      //! @brief returns m_ParameterList
      //! @return the parameter list as a ShPtrVector
      virtual const util::ShPtrVector< ParameterInterface> &GetParameterList() const = 0;

      //! @brief returns m_ParameterList
      //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
      virtual util::ShPtrVector< ParameterInterface> &GetParameterList() = 0;

      //! @brief returns the first parameter
      //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
      virtual const util::ShPtr< ParameterInterface> &GetFirstParameter() const = 0;

      //! @brief returns the # of required parameters for this flag
      //! @return the effective # of required parameters for this flag
      virtual size_t GetNumberRequiredParameters() const = 0;

      //! @brief Get vector of strings for the parameter
      storage::Vector< std::string> GetStringList() const
      {
        // call the GetObjectList with correct template specialization
        return GetObjectList< std::string>();
      }

      //! @brief template function that converts the parameters string list into an object vector of specified type
      //! this function requires a constructor/conversion function from std::string to t_DataType to be defined
      //! @return an object vector of specified type constructed from parameter list
      template< typename t_DataType>
      storage::Vector< t_DataType> GetObjectList() const
      {
        // make the vector of objects & allocate memory for all of the parameters
        storage::Vector< t_DataType> object_list;
        object_list.AllocateMemory( GetParameterList().GetSize());

        // iterate over all parameters and add each to the string list
        for
        (
          util::ShPtrVector< ParameterInterface>::const_iterator
            param_itr( GetParameterList().Begin()), param_itr_end( GetParameterList().End());
          param_itr != param_itr_end;
          ++param_itr
        )
        {
          // convert the string into an object
          object_list.PushBack( t_DataType( ( *param_itr)->GetValue()));
        }

        // end
        return object_list;
      }

      //! @brief convert the strings after the parameter into a set of unique t_DataType objects if their is a constructor t_DataType( std::string)
      //! @tparam t_DataType the type of the objects
      //! @return set of unique objects constructed from parameter strings
      template< typename t_DataType>
      storage::Set< t_DataType> GetObjectSet() const
      {
        const storage::Vector< t_DataType> object_list( GetObjectList< t_DataType>());
        return storage::Set< t_DataType>( object_list.Begin(), object_list.End());
      }

      //! @brief template function that converts the parameters string list into an numerical vector of specified type
      //! this function requires a meaningful ConvertStringToNumeralValue function from std::string to t_DataType to be defined
      //! @return an object vector of specified type constructed from parameter list
      template< typename t_DataType>
      storage::Vector< t_DataType> GetNumericalList() const
      {
        // make the vector of t_DataType & allocate memory for all of the parameters
        storage::Vector< t_DataType> numeric_list;
        numeric_list.AllocateMemory( GetParameterList().GetSize());

        // iterate over all parameters and add each to the string list
        for
        (
          util::ShPtrVector< ParameterInterface>::const_iterator
            param_itr( GetParameterList().Begin()), param_itr_end( GetParameterList().End());
          param_itr != param_itr_end;
          ++param_itr
        )
        {
          // convert the string into a numerical value
          numeric_list.PushBack( util::ConvertStringToNumericalValue< t_DataType>( ( *param_itr)->GetValue()));
        }

        // end
        return numeric_list;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief fill the parameter list from storage vector of strings
      //! @param PARAMETER_LIST a storage vector of strings containing the parameters as strings
      //! @param ERROR_STREAM the stream to which errors should be written
      //! @return true if successful, false otherwise
      virtual bool ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM) = 0;

      //! @brief checks that for every parameter was given a default value or it was passed through commandline
      //! @brief or if it is a dynamic list, if it meets the size specifications
      //! @param ERROR_STREAM the stream to which errors should be written
      //! @return true if successful, false otherwise
      virtual bool IsValidList( std::ostream &ERROR_STREAM) const = 0;

      //! @brief checks that for every parameter was given a default value or it was passed through commandline
      //! @brief or if it is a dynamic list, if it meets the size specifications
      //! @return true if successful, false otherwise
      //! @see IsValidList
      bool IsValidList() const
      {
        std::stringstream errorstream;
        return IsValidList( errorstream);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM the stream to which the help is written to
      //! @return the given stream to which the help was written to
      virtual std::ostream &WriteHelp( std::ostream &OSTREAM) const = 0;

      //! @brief writes the user provided commandline
      //! @param OSTREAM the stream to which the commandline is written
      //! @return given ostream to which the commandline was written
      virtual std::ostream &WriteUserCommand( std::ostream &OSTREAM) const = 0;

      //! @brief writes the usage command, complete with required parameters and flags
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      virtual std::ostream &WriteRequiredUsage( std::ostream &OSTREAM) const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the global bool for whether signals should be emitted for flags at this time
      //! Generally, when performing a dry run for commands, signals should not be emitted
      //! @brief detect whether signals should be emitted for flags at this time.  Generally, when performing a dry run
      //!        for commands, signals should not be emitted
      static const bool &ShouldSignal()
      {
        return GetShouldSignalBoolean();
      }

    private:

      friend class CommandState;

      //! @brief Set the flags to not signal
      static void SetNoSignal()
      {
        GetShouldSignalBoolean() = false;
      }

      //! @brief Set the flags to not signal
      static void SetSignal()
      {
        GetShouldSignalBoolean() = true;
      }

      //! @brief detect whether signals should be emitted for flags at this time.  Generally, when performing a dry run
      //!        for commands, signals should not be emitted
      static bool &GetShouldSignalBoolean()
      {
        static bool s_should_signal( true);
        return s_should_signal;
      }

    }; // FlagInterface

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_FLAG_INTERFACE_H_
