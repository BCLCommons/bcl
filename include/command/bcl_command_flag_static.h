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

#ifndef BCL_COMMAND_FLAG_STATIC_H_
#define BCL_COMMAND_FLAG_STATIC_H_

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
    //! @class FlagStatic
    //! @brief This class is a command line helper class derived from FlagInterface
    //!
    //! @see @link example_command_flag_static.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date 10/12/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FlagStatic :
      public FlagInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string m_Name;                                      //!< name of the flag
      std::string m_Description;                               //!< description of the flag
      bool m_WasSetInCommandline;                              //!< stores whether the flag was set in the commandline
      util::ShPtrVector< ParameterInterface> m_ParameterList;  //!< contains a list of Parameters
      FlagInterface::t_Signal m_Signal;                        //!< Thunk to call on setting or resetting the flag

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FlagStatic();

      //! @brief construct FlagStatic from NAME, DESCRIPTION
      //! @param NAME the name of the flag as a string
      //! @param DESCRIPTION a string description of the flag
      //! @param SIGNAL optional function to call whenever calling set or reset on the flag
      FlagStatic
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        FlagInterface::t_Signal SIGNAL = NULL
      );

      //! @brief construct FlagStatic from NAME, DESCRIPTION and single PARAMETER
      //! @param NAME the name of the flag as a string
      //! @param DESCRIPTION a string description of the flag
      //! @param PARAMETER the parameter associated with this flag
      //! @param SIGNAL optional function to call whenever calling set or reset on the flag
      FlagStatic
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const ParameterInterface &PARAMETER,
        FlagInterface::t_Signal SIGNAL = NULL
      );

      //! @brief construct FlagStatic from NAME, DESCRIPTION and single PARAMETER
      //! @param NAME the name of the flag as a string
      //! @param DESCRIPTION a string description of the flag
      //! @param PARAMETERS the parameters associated with this flag
      //! @param SIGNAL optional function to call whenever calling set or reset on the flag
      FlagStatic
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const util::ShPtrVector< ParameterInterface> &PARAMETERS,
        FlagInterface::t_Signal SIGNAL = NULL
      );

      //! @brief copy constructor
      //! @return pointer to FlagStatic object
      FlagStatic *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief name of class
      //! @return name of class as string
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
      //! @param PARAMETER_LIST a storage vector of strings containing the parameters as strings
      //! @param ERROR_STREAM the stream to which errors should be written
      //! @return true if successful, false otherwise
      bool ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM);

      //! @brief checks that for every parameter was given a default value or it was passed through commandline
      //! @brief or if it is a dynamic list, if it meets the size specifications
      //! @param ERROR_STREAM the stream to which errors should be written
      //! @return true if successful, false otherwise
      bool IsValidList( std::ostream &ERROR_STREAM) const;

      //! @brief adds PARAMETER to the m_ParameterList
      //! @param PARAMETER the parameter to be added
      void PushBack( const util::ShPtr< ParameterInterface> &PARAMETER);

    public:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the user provided commandline
      //! @param OSTREAM the stream to which the commandline is written
      //! @return given ostream to which the commandline was written
      std::ostream &WriteUserCommand( std::ostream &OSTREAM) const;

      //! @brief writes the help for the commandline
      //! @param OSTREAM the stream to which the help is written to
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM) const;

      //! @brief writes the usage command, complete with required parameters and flags
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &WriteRequiredUsage( std::ostream &OSTREAM) const;

    protected:

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

    }; // class FlagStatic

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_FLAG_STATIC_H_
