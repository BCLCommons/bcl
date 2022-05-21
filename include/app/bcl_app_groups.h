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

#ifndef BCL_APP_GROUPS_H_
#define BCL_APP_GROUPS_H_

// include the namespace header
#include "bcl_app.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"
#include "bcl_app_group_handler.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Groups
    //! @brief Class collecting all application_groups
    //! @detail Application groups are a set of applications that either:
    //! @detail   1. are tightly coupled to a particular workflow
    //! @detail   2. primarily operate off of inputs that other members of the group output
    //!
    //! @author mendenjl
    //! @see @link example_app_groups.cpp @endlink
    //! @date Mar 4, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Groups :
      public util::Enumerate< GroupHandler, Groups>
    {
      friend class util::Enumerate< GroupHandler, Groups>;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Groups();

    public:

    //////////
    // data //
    //////////

      ApplicationGroup e_Bcl;        //!< Apps whose input is usually absent or is a BCL-serialized storage object
      ApplicationGroup e_Protein;    //!< Apps whose input (or output) is usually proteins in PDB format
      ApplicationGroup e_Sequence;   //!< Apps whose input or output is usually a sequence / FASTA file
      ApplicationGroup e_Molecule;   //!< Apps whose input is usually molecules in SDF format
      ApplicationGroup e_Descriptor; //!< Apps whose input is just a set of descriptors
      ApplicationGroup e_Model;      //!< Apps whose input or output is usually a machine learning model
      ApplicationGroup e_ChemInfo;   //!< Apps that typically use both machine learning models and molecules
      ApplicationGroup e_BioInfo;    //!< Apps that typically use both machine learning models and proteins or sequences
      ApplicationGroup e_Restraint;  //!< Apps that mainly produce or work with restraints
      ApplicationGroup e_Density;    //!< Apps that use density maps and protein models or sequences

      //! Applications that should probably never be released due to e.g.
      //! 1. lack of general utility to external users
      //! 2. complicated, undocumented, input file requirements
      //! 3. Plans to merge functionality into other applications
      ApplicationGroup e_InternalBiol; //!< Internal biology applications
      ApplicationGroup e_InternalChem; //!< Internal chemistry applications
      ApplicationGroup e_Utility;      //!< Internal utility applications
      ApplicationGroup e_All;          //!< Group used for listing descriptions of all applications

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add a particular group to the app
      //! @param GROUP application group of interest
      //! @return the created enum
      ApplicationGroup &AddGroup( const GroupHandler &GROUP);

      //! @brief add new Interface derived application to a single application group
      //! @param APPLICATION the newly-created application
      //! @param GROUP the group to add it to
      //! @return the newly created application enum
      ApplicationType &AddAppToGroup( Interface *const &APPLICATION, ApplicationGroup &GROUPS);

      //! @brief add new Interface derived application to several groups
      //! @param APPLICATION the newly-created application
      //! @param GROUPS the set of groups to add it to
      //! @return the newly created application enum
      ApplicationType &AddAppToGroups( Interface *const &APPLICATION, const storage::Vector< ApplicationGroup> &GROUPS);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief writes the list of enums
      //! @param OSTREAM the stream to which the help is written to
      //! @return the given stream to which the list was written to
      //! Virtual to allow derived classes alter how the help is displayed without overriding Enum
      virtual std::ostream &WriteList( std::ostream &OSTREAM) const;

      //! @brief find all the groups that a particular app belongs to
      //! @param APP the application of interest
      //! @return all application group names that the given app belongs to
      storage::Vector< std::string> GetApplicationGroupsForApp( const Interface &APP);

    }; // class Apps

    BCL_API Groups &GetAppGroups();

  } // namespace app

  namespace util
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< app::GroupHandler, app::Groups>;

  } // namespace util

} // namespace bcl

#endif // BCL_APP_GROUPS_H_
