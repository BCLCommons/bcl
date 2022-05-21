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
#include "bcl_app_molecule_unique.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_set.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeUnique::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // molecule comparison level - constitution, configuration, exact
      sp_cmd->AddFlag( m_CompareFlag);

      // File to write molecules into
      sp_cmd->AddFlag( m_OutputFilenameFlag);

      // file to write duplicates into
      sp_cmd->AddFlag( m_OutputDuplicatesFilenameFlag);

      // bin size to use while comparing conformations
      sp_cmd->AddFlag( m_ConformerComparerFlag);

      // whether to merge descriptors
      sp_cmd->AddFlag( m_MergeDescriptorsFlag);

      // whether to merge descriptors
      sp_cmd->AddFlag( m_OverwriteDescriptorsFlag);

      // whether to merge descriptors
      sp_cmd->AddFlag( m_SameMoleculeAndSameNameFlag);

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

    //! @brief update the duplicates/unique count and write the molecules out
    //! @param WAS_UNIQUE whether the molecule was unique
    //! @param FRAGMENT the fragment to write
    void MoleculeUnique::PostProcess( const bool &WAS_UNIQUE) const
    {
      m_RemovalTimer.Stop();
      chemistry::FragmentFeed &feed( *m_SpFeed);
      m_OutputTimer.Start();
      if( WAS_UNIQUE)
      {
        // update the count
        ++m_UniqueCount;

        // write out the molecule, if the user gave an output path for unique molecules
        if( m_OutputFilenameFlag->GetFlag())
        {
          feed->WriteMDL( m_OutputUnique);
        }
      }
      else
      {
        // update the count
        ++m_DuplicatesCount;

        // write out the molecule, if the user gave an output path for duplicates
        if( m_OutputDuplicatesFilenameFlag->GetFlag())
        {
          feed->WriteMDL( m_OutputDuplicates);
        }
      }
      m_OutputTimer.Stop();
      m_InputTimer.Start();
      ++feed;
      m_InputTimer.Stop();
      m_RemovalTimer.Start();
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeUnique::Main() const
    {
      if( !m_OutputFilenameFlag->GetFlag() && !m_OutputDuplicatesFilenameFlag->GetFlag())
      {
        BCL_MessageCrt( "No output filename given! Only duplicate counts will be given");
      }
      BCL_Assert
      (
        !m_OutputDuplicatesFilenameFlag->GetFlag() || ( !m_MergeDescriptorsFlag->GetFlag() && !m_OverwriteDescriptorsFlag->GetFlag()),
        "Duplicates are not written when merging/overwriting descriptors"
      );
      BCL_Assert
      (
        !m_OverwriteDescriptorsFlag->GetFlag() || !m_MergeDescriptorsFlag->GetFlag(),
        "Cannot overwrite and merge descriptors"
      );

      // initialize the feed
      m_SpFeed = util::ShPtr< chemistry::FragmentFeed>( new chemistry::FragmentFeed);

      // get a reference to the feed
      chemistry::FragmentFeed &feed( *m_SpFeed);

      // reset counters
      m_DuplicatesCount = 0;
      m_UniqueCount = 0;
      m_OutputTimer.Reset();
      m_RemovalTimer.Reset();
      m_InputTimer.Reset();

      // where to output files to
      if( m_OutputFilenameFlag->GetFlag())
      {
        io::File::MustOpenOFStream( m_OutputUnique, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      }
      if( m_OutputDuplicatesFilenameFlag->GetFlag())
      {
        io::File::MustOpenOFStream( m_OutputDuplicates, m_OutputDuplicatesFilenameFlag->GetFirstParameter()->GetValue());
      }

      m_RemovalTimer.Start();
      if( !m_MergeDescriptorsFlag->GetFlag() && !m_OverwriteDescriptorsFlag->GetFlag())
      {
        // if not merging descriptors, molecules can be read in incrementally and then hashed.
        if( m_CompareFlag->GetFirstParameter()->GetValue() == "Exact")
        {
          storage::Map< std::string, storage::Set< std::string> > molecule_hashes;
          while( feed.NotAtEnd())
          {
            // create a hash
            const std::string molecule_hash
            (
              sdf::MdlHandler::CreateConformationalHashString( feed->GetAtomInfo(), feed->GetBondInfo())
            );
            PostProcess
            (
              molecule_hashes[ m_SameMoleculeAndSameNameFlag->GetFlag() ? feed->GetName() : ""].Insert( molecule_hash).second
            );
          }
        }
        else if( m_CompareFlag->GetFirstParameter()->GetValue() == "Configurations")
        {
          storage::Map< std::string, chemistry::ConfigurationSet> configurations;
          while( feed.NotAtEnd())
          {
            PostProcess
            (
              configurations[ m_SameMoleculeAndSameNameFlag->GetFlag() ? feed->GetName() : ""].Insert( chemistry::FragmentConfigurationShared( *feed)).second
            );
          }
        }
        else if( m_CompareFlag->GetFirstParameter()->GetValue() == "Conformations")
        {
          const double tolerance( m_ConformerComparerFlag->GetParameterList().LastElement()->GetNumericalValue< double>());
          util::Implementation< chemistry::ConformationComparisonInterface> comparison( m_ConformerComparerFlag->GetFirstParameter()->GetValue());
          chemistry::ConformationSet conformations( comparison, tolerance);

          if( !m_SameMoleculeAndSameNameFlag->GetFlag())
          {
            while( feed.NotAtEnd())
            {
              PostProcess( conformations.Insert( chemistry::FragmentConformationShared( *feed)).second);
            }
          }
          else
          {
            storage::Map< std::string, chemistry::ConformationSet> name_to_set;
            while( feed.NotAtEnd())
            {
              chemistry::ConformationSet &set( name_to_set[ feed->GetName()]);
              if( set.GetConformationsMap().IsEmpty())
              {
                set.Setup( comparison, tolerance);
              }
              PostProcess( set.Insert( chemistry::FragmentConformationShared( *feed)).second);
            }
          }
        }
        else
        {
          storage::Map< std::string, chemistry::ConstitutionSet> constitutions;
          while( feed.NotAtEnd())
          {
            PostProcess
            (
              constitutions[ m_SameMoleculeAndSameNameFlag->GetFlag() ? feed->GetName() : ""].Insert( chemistry::FragmentConstitutionShared( *feed)).second
            );
          }
        }
      }
      else
      {
        BCL_Assert( !m_SameMoleculeAndSameNameFlag->GetFlag(), "merge_descriptors and " + m_SameMoleculeAndSameNameFlag->GetName() + " are not currently supported together");
        // read all molecules from the feed
        chemistry::FragmentEnsemble ensemble_all;

        // keep track of whether each molecule is unique
        if( m_CompareFlag->GetFirstParameter()->GetValue() == "Exact")
        {
          std::map< std::string, chemistry::FragmentEnsemble::iterator> molecule_hashes;
          for( ; feed.NotAtEnd(); ++feed)
          {
            // create a hash
            const std::string molecule_hash
            (
              sdf::MdlHandler::CreateConformationalHashString( feed->GetAtomInfo(), feed->GetBondInfo())
            );
            std::pair< std::map< std::string, chemistry::FragmentEnsemble::iterator>::iterator, bool> itr_map
            (
              molecule_hashes.insert( std::make_pair( molecule_hash, ensemble_all.End()))
            );
            if( !itr_map.second)
            {
              if( m_MergeDescriptorsFlag->GetFlag())
              {
                itr_map.first->second->GetStoredPropertiesNonConst().Merge( feed->GetStoredProperties());
              }
              else
              {
                *itr_map.first->second = *feed;
              }
              ++m_DuplicatesCount;
            }
            else
            {
              ensemble_all.PushBack( *feed);
              itr_map.first->second = ensemble_all.GetMolecules().Last();
              ++m_UniqueCount;
            }
          }
        }
        else if( m_CompareFlag->GetFirstParameter()->GetValue() == "Configurations")
        {
          chemistry::ConfigurationSet configurations;
          std::map
          <
            util::SiPtr< const chemistry::FragmentConfigurationShared>,
            chemistry::FragmentEnsemble::iterator
          > iterator_map;
          for( ; feed.NotAtEnd(); ++feed)
          {
            // insert the iterator into the configurational map
            std::pair< chemistry::ConfigurationSet::const_iterator, bool> itr_conf_set
            (
              configurations.Insert( chemistry::FragmentConfigurationShared( *feed))
            );
            util::SiPtr< const chemistry::FragmentConfigurationShared> ptr_configuration( &**itr_conf_set.first);

            if( !itr_conf_set.second)
            {
              if( m_MergeDescriptorsFlag->GetFlag())
              {
                iterator_map[ ptr_configuration]->GetStoredPropertiesNonConst().Merge( feed->GetStoredProperties());
              }
              else
              {
                *iterator_map[ ptr_configuration] = *feed;
              }

              ++m_DuplicatesCount;
            }
            else
            {
              ensemble_all.PushBack( *feed);
              iterator_map[ ptr_configuration] = ensemble_all.GetMolecules().Last();
              ++m_UniqueCount;
            }
          }
        }
        else if( m_CompareFlag->GetFirstParameter()->GetValue() == "Conformations")
        {
          const double tolerance( m_ConformerComparerFlag->GetParameterList().LastElement()->GetNumericalValue< double>());
          util::Implementation< chemistry::ConformationComparisonInterface> comparison( m_ConformerComparerFlag->GetFirstParameter()->GetValue());
          chemistry::ConformationSet conformations( comparison, tolerance);

          std::map
          <
            util::SiPtr< const chemistry::FragmentConformationShared>,
            chemistry::FragmentEnsemble::iterator
          > iterator_map;
          for( ; feed.NotAtEnd(); ++feed)
          {
            // insert the iterator into the configurational map
            std::pair< chemistry::ConformationSet::const_iterator, bool> itr_conf_set
            (
              conformations.Insert( chemistry::FragmentConformationShared( *feed))
            );
            util::SiPtr< const chemistry::FragmentConformationShared> ptr_conformation( &**itr_conf_set.first);

            if( !itr_conf_set.second)
            {
              if( m_MergeDescriptorsFlag->GetFlag())
              {
                iterator_map[ ptr_conformation]->GetStoredPropertiesNonConst().Merge( feed->GetStoredProperties());
              }
              else
              {
                *iterator_map[ ptr_conformation] = *feed;
              }
              ++m_DuplicatesCount;
            }
            else
            {
              ensemble_all.PushBack( *feed);
              iterator_map[ ptr_conformation] = ensemble_all.GetMolecules().Last();
              ++m_UniqueCount;
            }
          }
        }
        else
        {
          chemistry::ConstitutionSet constitutions;
          std::map
          <
            util::SiPtr< const chemistry::FragmentConstitutionShared>,
            chemistry::FragmentEnsemble::iterator
          > iterator_map;
          for( ; feed.NotAtEnd(); ++feed)
          {
            // insert the iterator into the configurational map
            std::pair< chemistry::ConstitutionSet::const_iterator, bool> itr_set
            (
              constitutions.Insert( chemistry::FragmentConstitutionShared( *feed))
            );
            util::SiPtr< const chemistry::FragmentConstitutionShared> ptr_constitution( &**itr_set.first);

            if( !itr_set.second)
            {
              if( m_MergeDescriptorsFlag->GetFlag())
              {
                iterator_map[ ptr_constitution]->GetStoredPropertiesNonConst().Merge( feed->GetStoredProperties());
              }
              else
              {
                *iterator_map[ ptr_constitution] = *feed;
              }

              ++m_DuplicatesCount;
            }
            else
            {
              ensemble_all.PushBack( *feed);
              iterator_map[ ptr_constitution] = ensemble_all.GetMolecules().Last();
              ++m_UniqueCount;
            }
          }
        }
        if( m_OutputUnique.is_open())
        {
          ensemble_all.WriteMDL( m_OutputUnique);
        }
      }
      m_RemovalTimer.Stop();

      BCL_MessageStd
      (
        "Loaded " + util::Format()( m_DuplicatesCount + m_UniqueCount)
        + " molecules in " + m_InputTimer.GetTotalTime().GetTimeAsHourMinuteSecond()
      );

      BCL_MessageStd
      (
        "Found " + util::Format()( m_DuplicatesCount)
        + " duplicates (Level: " + m_CompareFlag->GetFirstParameter()->GetValue()
        + ") in " + m_RemovalTimer.GetTotalTime().GetTimeAsHourMinuteSecond()
      );

      // print output times
      if( m_OutputFilenameFlag->GetFlag() || m_OutputDuplicatesFilenameFlag->GetFlag())
      {
        BCL_MessageStd
        (
          "Total time writing molecules: "
          + m_OutputTimer.GetTotalTime().GetTimeAsHourMinuteSecond()
        );
      }

      // where to output files to
      if( m_OutputFilenameFlag->GetFlag())
      {
        BCL_MessageStd
        (
          "Wrote ensemble of " + util::Format()( m_UniqueCount)
          + " unique molecules to " + m_OutputFilenameFlag->GetFirstParameter()->GetValue()
        );
        io::File::CloseClearFStream( m_OutputUnique);
      }
      // where to output files to
      if( m_OutputDuplicatesFilenameFlag->GetFlag())
      {
        BCL_MessageStd
        (
          "Wrote ensemble of " + util::Format()( m_DuplicatesCount)
          + " duplicate molecules to " + m_OutputDuplicatesFilenameFlag->GetFirstParameter()->GetValue()
        );
        io::File::CloseClearFStream( m_OutputDuplicates);
      }

      //end
      return 0;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeUnique::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoleculeUnique::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! @brief standard constructor
    MoleculeUnique::MoleculeUnique() :
      m_CompareFlag
      (
        new command::FlagStatic
        (
          "compare",
          "Choice of comparison levels to determine duplicate molecules:\n"
          "        * Constitutions: Same Atoms & Connectivity (default)\n"
          "        * Configurations: Same Atoms, Connectivity, and Stereochemistry\n"
          "        * Conformations: Same Configuration but different 3D conformation\n"
          "        * Exact: Same Atoms in same order, connectivity, bond orders, and 3D coordinates",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create( "Exact", "Configurations", "Conformations", "Constitutions")
            ),
            "Constitutions"
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "sdf filename for where to write out molecules",
          command::Parameter
          (
            "filename",
            "sdf filename for where to write out molecules",
            ""
          )
        )
      ),
      m_OutputDuplicatesFilenameFlag
      (
        new command::FlagStatic
        (
          "output_dupes",
          "sdf filename for where to write out duplicate molecules",
          command::Parameter
          (
            "filename",
            "sdf filename for where to write out molecules",
            ""
          )
        )
      ),
      m_ConformerComparerFlag
      (
        new command::FlagStatic
        (
          "conformation_comparer",
          "if conformations are being compared, this flag sets the method for comparison",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "rmsd_type",
                "rmsd method to be used for comparing conformations",
                command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>()),
                "DihedralBins(bin size=30)"
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "tolerance",
                "amount of tolerance allowed between two conformers",
                command::ParameterCheckRanged< double>( 0.0, std::numeric_limits< double>::max()),
                "1.0 "
              )
            )
          )
        )
      ),
      m_MergeDescriptorsFlag
      (
        new command::FlagStatic
        (
          "merge_descriptors",
          "If true, merge descriptors of the duplicates. Note that this significantly increases memory requirements"
        )
      ),
      m_OverwriteDescriptorsFlag
      (
        new command::FlagStatic
        (
          "overwrite_descriptors",
          "If true, overwrite descriptors of the duplicates. Note that this significantly increases memory requirements and "
          "is mutually exclusive with merging the descriptors"
        )
      ),
      m_SameMoleculeAndSameNameFlag
      (
        new command::FlagStatic
        (
          "same_molecule_same_name",
          "If set, only match molecules if they have identical names and otherwise match based on the comparer"
        )
      ),
      m_InputTimer( false),
      m_OutputTimer( false),
      m_RemovalTimer( false)
    {
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeUnique::GetDescription() const
    {
      return "MoleculeUnique removes duplicate molecules at constitution, configuration or conformation levels, defined by user";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeUnique::GetReadMe() const
    {

      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::MoleculeUnique, terms of use, "
        "appropriate citation, installation procedures, BCL::MoleculeUnique execution, "
        "technical support, and future research directions.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::MoleculeUnique?\n"
        "BCL::MoleculeUnique is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::MoleculeUnique is a utility identifying unique "
        "molecules at constitution, configuration or conformation level from an ensemble of molecules. A variety of " ""
        "comparison algorithms are available such as RMSD.\n"
        "\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::MoleculeUnique.\n"
        "\n"
        "Butkiewicz M, Lowe EW, Mueller R, Mendenhall JL, Teixeira PL, Weaver CD, Meiler J. "
        "Benchmarking Ligand-Based Virtual High-Throughput Screening with the PubChem Database "
        "Molecules, 18, (1), 735-756. ; 2013\n"
        "Link:  www.http://meilerlab.org/index.php/publications/show/2013\nJournal link: http://www.mdpi.com/1420-3049/18/1/735\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::MoleculeUnique.\n"
        "Running BCL::MoleculeUnique requires an sdf file containing molecules that need to be uniquely identified.\n"
        "\n"
        "2) Run BCL::MoleculeUnique to identify unique molecules\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::MoleculeUnique executable program. The syntax for"
        "running the application looks like the following"
        "\n"
        "bcl.exe MoleculeUnique -input_filenames <filename.sdf> -output <filename>"
        "\n\nFLAGS:\n\n"
        "-input_filenames <filename> -> file containing ensemble of molecules that need to be identified uniquely\n"
        "-output <filename> -> file to which unique molecules will be written out\n"
        "-compare <Comparison Options> -> the level at which molecules need to be compared\n"
        "\tOptions for comparing molecules include : \n"
        "\t\t1. Constitutions - Duplicate means a molecule with the same connectivity (same atom types and bonds, reordering allowed)\n"
        "\t\t2. Configurations - In addition to connectivity, chirality and stereoisometry of the bonds must also be equal\n"
        "\t\t3. Conformations - In addition to configuration, all dihedral angles are within same bins ( bin of sizes 30 degrees centered around 0, 30, 60 degrees and so on)\n"
        "\t\t4. Exact - Duplicate means the SDF entries are identical (e.g. same 3D coordinates, order of atoms, etc.)"
        "\n-bin_size <double> -> dihedral bin size to use when comparing conformations\n"
        "-output_dupes <filename> -> file to which duplicate molecules will be written out \n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe MoleculeUnique -help\n"
        "\n"
        "For more general information about the product, type bcl.exe MoleculeUnique -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::MoleculeUnique.\n"
        "BCL::MoleculeUnique is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );
      return s_read_me;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &MoleculeUnique::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::MoleculeUnique is a utility identifying unique or highly similar molecules from a file"
        "Features of BCL::MoleculeProperties\n"
        "<ul>"
        "  <li>Find duplicate molecules in an sdf file by various criteria such as\n"
        "    <ul>\n"
        "      <li>Same constitution - identical gasteiger atom types, bond types, and connectictivity</li>"
        "      <li>Same configuration - like constitution, but also considers chirality and double bond isometry</li>"
        "      <li>Similar conformation - Find a rotamer set based on RMSD tolerance, dihedral angle tolerance, tanimoto similarity, and more</li>"
        "      <li>Exact - Duplicate entries in a molecule file i.e. same 3D coordinates and atom ordering</li>"
        "    </ul>"
        "  </li>"
        "  <li>Duplicates and unique molecules can be written out to separate files</li>"
        "  <li>Compressed molecule files (bz2, gzip) are supported</li>"
        "</ul>\n\n"
        "!molecule_unique.png!"
      );

      return s_web_text;
    }

    const ApplicationType MoleculeUnique::MoleculesUnique_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeUnique(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
