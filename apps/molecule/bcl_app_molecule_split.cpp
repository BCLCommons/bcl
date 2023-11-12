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
#include "command/bcl_command_parameter_check_ranged.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "molecule/bcl_app_molecule_split.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    const ApplicationType MoleculeSplit::MoleculeSplit_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeSplit(), GetAppGroups().e_Molecule)
    );

    //! @brief standard constructor
    MoleculeSplit::MoleculeSplit() :
      m_ImplementationFlag
      (
        new command::FlagStatic
        (
          "implementation",
          "method to split molecules",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::FragmentSplitInterface>())
          )
        )
      ),
      m_MinFragSizeFlag
      (
        new command::FlagStatic
        (
          "min_frag_size",
          "only output a split fragment if it contains at least this many heavy atoms",
          command::Parameter
          (
            "minimum number of atoms",
            "default: output fragments of all sizes",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()), "0"
          )
        )
      ),
      m_RecenterFlag
      (
        new command::FlagStatic
        (
          "recenter",
          "translate molecules such that the middle of the molecule is at the origin"
        )
      ),
      m_PreserveMDLPropertiesFlag
      (
        new command::FlagStatic
        (
          "preserve_mdl_properties",
          "keep the MDL properties of original molecule on all split components"
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "file to write split molecules into",
          command::Parameter
          (
            "output",
            "file to write split molecules into"
          )
        )
      )
    {
    }

    //! copy constructor; skips i/o streams
    MoleculeSplit::MoleculeSplit( const MoleculeSplit &PARENT) :
          m_ImplementationFlag( PARENT.m_ImplementationFlag),
          m_MinFragSizeFlag( PARENT.m_MinFragSizeFlag),
          m_RecenterFlag( PARENT.m_RecenterFlag),
          m_PreserveMDLPropertiesFlag( PARENT.m_PreserveMDLPropertiesFlag),
          m_OutputFilenameFlag( PARENT.m_OutputFilenameFlag),
          m_MoleculeIndex( 0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeSplit *MoleculeSplit::Clone() const
    {
      return new MoleculeSplit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeSplit::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeSplit::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // split molecular complexes into separate small molecules
      sp_cmd->AddFlag( m_ImplementationFlag);

      // minimum fragment output size
      sp_cmd->AddFlag( m_MinFragSizeFlag);

      //! whether to recenter the molecules
      sp_cmd->AddFlag( m_RecenterFlag);

      //! keep the MDL properties of original molecule on all split components
      sp_cmd->AddFlag( m_PreserveMDLPropertiesFlag);

      //! output filename
      sp_cmd->AddFlag( m_OutputFilenameFlag);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags
        (
          *sp_cmd,
          storage::Set< command::FlagTypeEnum>( command::e_Pthread)
        );

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeSplit::Main() const
    {
      // open the output file
      io::File::MustOpenOFStream( m_OutputFile, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      m_MoleculeIndex = 0;

      chemistry::FragmentFeed itr_fragments;

      m_Split = m_ImplementationFlag->GetFirstParameter()->GetValue();

      if( m_ImplementationFlag->GetFlag())
      {
        for( ; itr_fragments.NotAtEnd(); ++itr_fragments, ++m_MoleculeIndex)
        {
          const storage::Map< std::string, std::string> &properties( itr_fragments->GetStoredProperties().GetMDLProperties());
          chemistry::FragmentEnsemble isolated_objects( ( *m_Split)( *itr_fragments));
          for
          (
            storage::List< chemistry::FragmentComplete>::iterator
            itr( isolated_objects.Begin()), itr_end( isolated_objects.End());
            itr != itr_end;
            ++itr
          )
          {
            if( m_PreserveMDLPropertiesFlag->GetFlag())
            {
              for
              (
                  auto property_itr( properties.Begin()), property_itr_end( properties.End());
                  property_itr != property_itr_end;
                  ++property_itr
              )
              {
                itr->StoreProperty( property_itr->first, property_itr->second);
              }
            }
            Write( *itr);
          }
        }
      }
      else
      {
        for( ; itr_fragments.NotAtEnd(); ++itr_fragments, ++m_MoleculeIndex)
        {
          chemistry::FragmentComplete fragment( *itr_fragments);
          Write( fragment);
        }
      }

      io::File::CloseClearFStream( m_OutputFile);

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeSplit::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoleculeSplit::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes a single molecule out
    //! @param MOLECULE the molecule or fragment to write
    void MoleculeSplit::Write( chemistry::FragmentComplete &MOLECULE) const
    {
      if( m_MinFragSizeFlag->GetFlag())
      {
        if( MOLECULE.GetSize() - MOLECULE.GetNumberHydrogens() >= m_MinFragSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>())
        {
          if( m_RecenterFlag->GetFlag())
          {
            MOLECULE.Translate( -MOLECULE.GetCenter());
          }
          MOLECULE.WriteMDL( m_OutputFile);
        }
      }
      else
      {
        if( m_RecenterFlag->GetFlag())
        {
          MOLECULE.Translate( -MOLECULE.GetCenter());
        }
        MOLECULE.WriteMDL( m_OutputFile);
      }
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeSplit::GetReadMe() const
    {
      static std::string s_read_me =
        "MoleculeSplit splits all molecules in an ensemble by one of several methods:\n"
        " All methods add a property (ParentIndex) indicating the parent molecule index\n"
        " -isolates splits molecular complexes (an sdf entry containing disconnected molecules) into constituent molecules\n"
        " -largest keeps the largest component of each sdf entry (an sdf entry containing a connected molecule remains unchanged\n"
        " -rings splits molecules into constituent ring systems\n"
        " -chains splits molecules into simple chains\n"
        " -rings can be combined with -chains to split molecules into constituent chains and rings\n";

      return s_read_me;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeSplit::GetDescription() const
    {
      return "Generate molecule fragments via multiple different splitting schemes";
    }

  } // namespace app
} // namespace bcl
