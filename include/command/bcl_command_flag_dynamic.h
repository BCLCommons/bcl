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

#ifndef BCL_COMMAND_FLAG_DYNAMIC_H_
#define BCL_COMMAND_FLAG_DYNAMIC_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_flag_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FlagDynamic
    //! @brief This class is a command line helper class derived from FlagInterface
    //!
    //! @see @link example_command_flag_dynamic.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date 10/12/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FlagDynamic :
      public FlagInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string                            m_Name;                //!< name of the flag
      std::string                            m_Description;         //!< description of the flag
      bool                                   m_WasSetInCommandline; //!< stores whether the flag was set in the commandline
      util::ShPtrVector< ParameterInterface> m_ParameterList;       //!< contains a list of Parameters
      size_t                                 m_MinNumberParameters; //!< minimal number of parameters for a dynamic list
      size_t                                 m_MaxNumberParameters; //!< maximal number of parameters for a dynamic list
      util::ShPtrVector< ParameterInterface> m_TemplateParameters;  //!< template parameters - every parameter has to be like this one
      FlagInterface::t_Signal                m_Signal;              //!< Thunk to call on setting or resetting the flag

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
      FlagDynamic();

      //! @brief construct FlagDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETER,
      //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
      //! @param NAME - the name of the flag as a std::string
      //! @param DESCRIPTION - a description of the flag as a std::string
      //! @param TEMPLATE_PARAMETER
      //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
      //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
      //! @param SIGNAL optional function to call whenever calling set or reset on the flag
      FlagDynamic
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const ParameterInterface &TEMPLATE_PARAMETER,
        const size_t MIN_NUMBER_PARAMETERS = size_t( 0),
        const size_t MAX_NUMBER_PARAMETERS = util::GetUndefinedSize_t(),
        t_Signal SIGNAL = NULL
      );

      //! @brief construct FlagDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETERS,
      //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
      //! @param NAME - the name of the flag as a std::string
      //! @param DESCRIPTION - a description of the flag as a std::string
      //! @param TEMPLATE_PARAMETER
      //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
      //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
      //! @param SIGNAL optional function to call whenever calling set or reset on the flag
      FlagDynamic
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const storage::Vector< Parameter> &TEMPLATE_PARAMETERS,
        const size_t MIN_NUMBER_PARAMETERS = size_t( 0),
        const size_t MAX_NUMBER_PARAMETERS = util::GetUndefinedSize_t(),
        t_Signal SIGNAL = NULL
      );

      //! @brief copy constructor
      FlagDynamic *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief name of class
      //! @return the name of the class as a std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns name of flag
      //! @return name of flag as std::string
      const std::string &GetName() const;

      //! @brief set name - only possible if there was no name given yet
      //! @param NAME - what you want to name the flag
      void SetName( const std::string &NAME);

      //! @brief returns description of flag
      //! @return description of flag as string
      const std::string &GetDescription() const;

      //! @brief set flag to true if it was false
      void SetFlag();

      //! @brief set flag to false if it was true
      void UnsetFlag();

      //! @brief reset the flag
      //! @detail resets all internal parameters, then removes all internally held parameters
      void ResetFlag();

      //! @brief returns if m_IsSet in command line
      //! @return whether this flag was set in the command line
      bool GetFlag() const;

      //! @brief returns the function to be called whenever this flag is updated
      //! @return the function to be called whenever this flag is updated
      t_Signal GetSignal() const;

      //! @brief returns m_ParameterList
      //! @return the parameter list as a ShPtrVector
      const util::ShPtrVector< ParameterInterface> &GetParameterList() const;

      //! @brief returns m_ParameterList
      //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
      util::ShPtrVector< ParameterInterface> &GetParameterList();

      //! @brief returns the first parameter
      //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
      const util::ShPtr< ParameterInterface> &GetFirstParameter() const;

      //! @brief returns the # of required parameters for this flag
      //! @return the effective # of required parameters for this flag
      size_t GetNumberRequiredParameters() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief fill the parameter list from storage vector of strings
      //! @param PARAMETER_LIST - the parameter list as a storage vector of strings
      //! @param ERROR_STREAM - the stream to which errors should be written
      //! @return true on success, false on failure
      //! does not allow you to add more than m_MaxNumberParameters parameters to the list
      //! @see SetFlag
      //! @see PushBack
      bool ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM);

      //! @brief checks if there are too few items in your parameter list
      //! @param ERROR_STREAM - the stream to which errors should be written
      //! @return true if the parameter list is of a valid size, false otherwise
      //! @see MeetSizeSpecification
      bool IsValidList( std::ostream &ERROR_STREAM) const;

    private:

      //! @brief returns true if the number of parameters is between m_MinNumberParameters and m_MaxNumberParameters
      //! @return true if the parameter list meets the size specification, false otherwise
      bool MeetSizeSpecification() const;

      //! @brief adds PARAMETER to the m_ParameterList
      //! @param PARAMETER - the parameter to be added to the parameter list
      //! @param ERROR_STREAM - the stream to which errors should be written
      //! @return true on success, false otherwise
      //! @see FlagInterface::PushBack
      bool PushBack( const std::string &PARAMETER, std::ostream &ERROR_STREAM);

    public:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the command line
      //! @param OSTREAM - the stream to write to
      //! @return the stream after you wrote to it
      std::ostream &WriteHelp( std::ostream &OSTREAM) const;

      //! @brief writes the user provided commandline
      std::ostream &WriteUserCommand( std::ostream &OSTREAM) const;

      //! @brief writes the usage command, complete with required parameters and flags
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &WriteRequiredUsage( std::ostream &OSTREAM) const;

    protected:

      //! @brief write Flag with Params to ostream
      //! @param OSTREAM - the stream to write to
      //! @return the stream after you wrote to it
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read Flag with params from istream
      //! @param ISTREAM - the stream to read from
      //! @return the stream after you read from it
      std::istream &Read( std::istream &ISTREAM);

    }; //class FlagDynamic

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_FLAG_DYNAMIC_H_
