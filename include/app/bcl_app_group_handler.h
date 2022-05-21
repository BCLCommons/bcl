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

#ifndef BCL_APP_GROUP_HANDLER_H_
#define BCL_APP_GROUP_HANDLER_H_

// include the namespace header
#include "bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_app_interface.h"
#include "command/bcl_command_flag_interface.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GroupHandler
    //! @brief represents a group of applications
    //! @remarks every application should add itself to its corresponding groups
    //!
    //! @author mendenjl
    //! @date Feb 26, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GroupHandler :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string                                m_GroupName;         //!< name of the group
      std::string                                m_GroupDescription;  //!< description of the group
      storage::Map< std::string, std::string>    m_NameToDescription; //!< map of implementations
      storage::Map< std::string, std::string>    m_NameToAliases;     //!< map of implementations
      util::ShPtrVector< command::FlagInterface> m_Flags;             //!< flags to apply to all members of the group
      size_t                                     m_LongestNameSize;   //!< Length of the longest app name for this group
      static size_t                              s_LongestGroupSize;  //!< Length of the longest group name in any
      static size_t                              s_LongestNameSize;   //!< Length of the longest app name in any

    public:

      //! delimiter between the group name and the application name
      static const char s_GroupDelimiter = ':';

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default group handler
      GroupHandler();

      //! @brief constructor from string of group name
      GroupHandler( const std::string &NAME, const std::string &DESCRIPTION);

      //! @brief Clone function
      //! @return pointer to new GroupHandler
      GroupHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the group
      const std::string &GetName() const;

      //! @brief get the description of the group
      const std::string &GetDescription() const;

      //! @brief get the length of the longest application name for this group
      size_t GetLengthLongestAppNameThisGroup() const
      {
        return m_LongestNameSize;
      }

      //! @brief get the length of the longest application name for any group
      static size_t GetLengthLongestAppName()
      {
        return s_LongestNameSize;
      }

      //! @brief get the length of the longest application name for any group
      static size_t GetLengthLongestAppGroupName()
      {
        return s_LongestGroupSize;
      }

      //! @brief test whether a particular application is in this group
      //! @param APP the application of interst
      bool Contains( const Interface &APP) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add group flags to a given command
      //! @param COMMAND the command to add group flags to
      void AddGroupFlags( command::Command &COMMAND);

      //! @brief add a flag to the application group
      //! @param FLAG the new flag to add
      void AddFlagToGroup( const util::ShPtr< command::FlagInterface> &FLAG);

      //! @brief add an application to this group
      ApplicationType &AddInstance( const util::ShPtr< Interface> &APP);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent to write before the help
      //! @param APP_INDENT the amount of indent to write before each application
      //! @param INCLUDE_DESCRIPTION description for each application
      //! @param INCLUDE_MEMBERS whether to write the help for every group member
      //! @param INCLUDE_MEMBER_DESCRIPTIONS whether to write a brief description for each group member
      //! @param INCLUDE_MEMBER_ALIASES whether to write aliases (even if they are deprecated) for each group member
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp
      (
        std::ostream &OSTREAM,
        const size_t INDENT = 0,
        const size_t &APP_INDENT = s_LongestNameSize,
        const bool &INCLUDE_DESCRIPTION = true,
        const bool &INCLUDE_MEMBERS = true,
        const bool &INCLUDE_MEMBER_DESCRIPTIONS = true,
        const bool &INCLUDE_MEMBER_ALIASES = true
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      //! @see ParameterCheckInterface::Write
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class GroupHandler

  } // namespace app
} // namespace bcl

#endif // BCL_APP_GROUP_HANDLER_H_
