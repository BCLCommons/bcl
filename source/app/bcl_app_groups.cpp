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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "app/bcl_app_groups.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "io/bcl_io_fixed_line_width_writer.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Groups::Groups() :
      util::Enumerate< GroupHandler, Groups>( false),
      e_Bcl(          AddGroup( GroupHandler( "bcl", "applications that use only BCL-generated input"))),
      e_Protein   (   AddGroup( GroupHandler( "protein", "applications that primarily use or create proteins (usually PDB files)"))),
      e_Sequence  (   AddGroup( GroupHandler( "sequence", "applications that primarily use AA sequences (usually FASTA files)"))),
      e_Molecule  (   AddGroup( GroupHandler( "molecule", "applications that primarily use or create molecules (usually SDF files)"))),
      e_Descriptor(   AddGroup( GroupHandler( "descriptor", "applications that primarily use or create descriptor objects"))),
      e_Model     (   AddGroup( GroupHandler( "model", "applications that primarily use or create machine learning models"))),
      e_ChemInfo  (   AddGroup( GroupHandler( "cheminfo", "applications that primarily use or create molecules using machine learning models"))),
      e_BioInfo(      AddGroup( GroupHandler( "bioinfo", "applications that primarily use AA sequences and machine learning models"))),
      e_Restraint (   AddGroup( GroupHandler( "restraint", "applications that primarily use or create restraints that guide protein folding"))),
      e_Density   (   AddGroup( GroupHandler( "density", "applications that primarily use or create density maps from sequence or protein models"))),
      e_InternalBiol( AddGroup( GroupHandler( "bioutil", "bcl development and internal analysis tools for proteins or sequences"))),
      e_InternalChem( AddGroup( GroupHandler( "molutil", "bcl development and internal analysis tools for molecules"))),
      e_Utility(      AddGroup( GroupHandler( "util", "utilities used primarily for bcl development and internal analysis"))),
      e_All(          AddGroup( GroupHandler( "all", "")))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Groups::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a particular group to the app
    //! @param GROUP application group of interest
    //! @return the created enum
    ApplicationGroup &Groups::AddGroup( const GroupHandler &GROUP)
    {
      return AddEnum( GROUP.GetName(), GROUP);
    }

    //! @brief add new Interface derived application to a single application group
    //! @param APPLICATION the newly-created application
    //! @param GROUP the group to add it to
    //! @return the newly created application enum
    ApplicationType &Groups::AddAppToGroup( Interface *const &APPLICATION, ApplicationGroup &GROUP)
    {
      if( GetVersion().IsLicense() && !APPLICATION->IsReleaseApplication())
      {
        static ApplicationType s_undefined_app;
        // do not add non-release apps to the enum; this keeps internal apps invisible to external users
        return s_undefined_app;
      }
      return GROUP->AddInstance( util::ShPtr< Interface>( APPLICATION));
    }

    //! @brief add new Interface derived application to several groups
    //! @param APPLICATION the newly-created application
    //! @param GROUPS the set of groups to add it to
    //! @return the newly created application enum
    ApplicationType &Groups::AddAppToGroups
    (
      Interface *const &APPLICATION,
      const storage::Vector< ApplicationGroup> &GROUPS
    )
    {
      // handle if no groups were provided
      BCL_Assert( !GROUPS.IsEmpty(), "Tried to add an application without specifying any groups");

      // create a shared pointer with the new app
      util::ShPtr< Interface> sp_app( APPLICATION);

      // handle all but the last group
      for
      (
        storage::Vector< ApplicationGroup>::const_iterator itr( GROUPS.Begin()), itr_end( GROUPS.End() - 1);
        itr != itr_end;
        ++itr
      )
      {
        ApplicationGroup( *itr)->AddInstance( sp_app);
      }

      // add it to the last-specified group
      return ApplicationGroup( GROUPS.LastElement())->AddInstance( sp_app);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes the list of enums
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the list was written to
    //! Virtual to allow derived classes alter how the help is displayed without overriding Enum
    std::ostream &Groups::WriteList( std::ostream &OSTREAM) const
    {
      OSTREAM << '\n';
      io::FixedLineWidthWriter writer( 0, util::GetLogger().GetMaxLineWidth() - 4);
      const size_t description_indent( GroupHandler::GetLengthLongestAppGroupName() + 2);
      writer << "<group>:";
      writer << std::string( description_indent - std::min( writer.GetLinePosition(), description_indent), ' ') << "Description" << '\n';
      writer << "<group>:Help";
      writer << std::string( description_indent - std::min( writer.GetLinePosition(), description_indent), ' ');
      writer.SetIndent( description_indent);
      writer << "prints all <application>'s and descriptions for the group";
      OSTREAM << writer.String() << std::flush;

      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        ( *itr)->WriteHelp( OSTREAM, 0, 2, true, true, false);
        OSTREAM << std::flush;
      }
      e_Undefined->WriteHelp( OSTREAM, 0, 2, true, true, false);
      return OSTREAM;
    }

    //! @brief find all the groups that a particular app belongs to
    //! @param APP the application of interest
    //! @return all application group names that the given app belongs to
    storage::Vector< std::string> Groups::GetApplicationGroupsForApp( const Interface &APP)
    {
      storage::Vector< std::string> app_groups;
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( ( *itr)->Contains( APP) && !( *itr)->GetName().empty())
        {
          app_groups.PushBack( ( *itr)->GetName());
        }
      }
      return app_groups;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    Groups &GetAppGroups()
    {
      return Groups::GetEnums();
    }

  } // namespace app

  namespace util
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< app::GroupHandler, app::Groups>;

  } // namespace util
} // namespace bcl
