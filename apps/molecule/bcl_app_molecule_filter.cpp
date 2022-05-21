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
#include "bcl_app_molecule_filter.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_rmsd.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "storage/bcl_storage_triplet.h"

namespace bcl
{
  namespace app
  {

    //! @brief get a description for the given type
    //! @param TYPE the desired type
    //! @return a string description for that type
    const std::string &FilterStatistics::GetTypeDescription( const Type &TYPE)
    {
      static const std::string s_type_names[ size_t( s_NumberFilterStatistics) + 1] =
      {
        "with defined gasteiger atom types ",
        "with reasonable 3D coordinates ",
        "with planar amide bonds",
        "with planar unbridged aromatic ring systems",
        "with simple connectivity -> not a molecular complex "
        "containing at least one fragment from ",
        "matching at least one molecule from ",
        "with property ",
        "with property ",
        GetStaticClassName< Type>()
      };
      return s_type_names[ TYPE];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FilterStatistics::FilterStatistics() :
      m_Timer( false),
      m_Type( s_NumberFilterStatistics),
      m_Description(),
      m_Count( 0)
    {
      m_Timer.Reset();
    }

    //! @brief constructor from description
    FilterStatistics::FilterStatistics( const FilterStatistics::Type &TYPE, const std::string &DESCRIPTION) :
      m_Timer( false),
      m_Type( TYPE),
      m_Description( GetTypeDescription( TYPE) + " " + DESCRIPTION),
      m_Count( 0)
    {
      m_Timer.Reset();
    }

    //! @brief register the start of a task
    void FilterStatistics::StartFilter()
    {
      m_Timer.Start();
    }

    //! @brief stop the task, increment the count
    //! @param MATCHED did the molecule match the filter
    void FilterStatistics::StopFilter( const bool &MATCHED)
    {
      m_Timer.Stop();
      if( MATCHED)
      {
        ++m_Count;
      }
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    MoleculeFilter::MoleculeFilter() :
      m_OutputFilenameForMoleculesSatisfyingCriteriaFlag
      (
        new command::FlagStatic
        (
          "output_matched",
          "sdf filename for where to write out molecules that met all the criteria",
          command::Parameter
          (
            "output_matched",
            "sdf filename for where to write out molecules that met all the criteria",
            "matched.sdf.bz2"
          )
        )
      ),
      m_OutputFilenameForMoleculesFailingCriteriaFlag
      (
        new command::FlagStatic
        (
          "output_unmatched",
          "sdf filename for where to write out molecules that failed any of the criteria",
          command::Parameter
          (
            "output_unmatched",
            "sdf filename for where to write out molecules that that failed any of the criteria",
            "unmatched.sdf.bz2"
          )
        )
      ),
      m_Match3DCoordinatesFlag
      (
        new command::FlagStatic
        (
          "3d",
          "match molecules that have 3d coordinates on at least one neighbor of every atom with 4 neighbors"
        )
      ),
      m_MatchDefinedAtomTypesFlag
      (
        new command::FlagStatic
        (
          "defined_atom_types",
          "match molecules with no undefined atom types"
        )
      ),
      m_MatchPlanarAmideBondsFlag
      (
        new command::FlagStatic
        (
          "only_planar_amide_bonds",
          "match molecules for which any amide or thioamide bonds outside rings are within a specified number of degrees of planar",
          command::Parameter
          (
            "tolerance",
            "number of degrees out of plane to allow amide bonds to be",
            command::ParameterCheckRanged< double>( 0.0, 180.0),
            "20.0"
          )
        )
      ),
      m_MatchPlanarAromaticRingsFlag
      (
        new command::FlagStatic
        (
          "only_planar_aromatic_rings",
          "match molecules for which any simple, non-bridged aromatic rings are planar"
        )
      ),
      m_MatchSimpleMoleculesFlag
      (
        new command::FlagStatic
        (
          "simple",
          "matches molecules which are simple (e.g. not molecular complexes)"
        )
      ),
      m_MoleculesContainingFragmentsFilenameFlag
      (
        new command::FlagStatic
        (
          "contains_fragments_from",
          "filename for input sdf of fragments to find (constitutional level). Molecules with any of the fragments satisfy this criterion",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "contains_fragments_from",
                "filename for input sdf of fragments to find. Molecules with any of the fragments satisfy this criterion",
                command::ParameterCheckFileExistence(),
                ""
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "atom_comparison_type",
                "string value of the atomtype info for substructure matching calculation",
                command::ParameterCheckSerializable
                (
                  chemistry::ConformationGraphConverter::AtomComparisonTypeEnum()
                ),
                "ElementType"
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "bond_comparison_type",
                "string value of the bondtype info for substructure matching calculation",
                command::ParameterCheckSerializable
                (
                  chemistry::ConfigurationalBondTypeData::DataEnum()
                ),
                "BondOrderAmideOrAromaticWithRingness"
              )
            )
          )
        )
      ),
      m_MoleculesContainingConformersFilenameFlag
      (
        new command::FlagStatic
        (
          "contains_conformers",
          "filename for input sdf of conformers to find (conformational level). Molecules within tolerance satisfy this criterion",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "contains_conformers",
                "filename for input sdf of conformers to find. Molecules within tolerance satisfy this criterion",
                command::ParameterCheckFileExistence(),
                ""
              )
            )
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
                "comparison",
                "comparison to perform between input molecules and conformers",
                command::ParameterCheckSerializable( math::Comparisons< float>::EnumType()),
                "less_equal"
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
      m_MoleculesToFindFilenameFlag
      (
        new command::FlagStatic
        (
          "contains",
          "filename for input sdf of molecules to find (constitutional level). Molecules that match any of the molecules in the ensemble satisfy this criterion",
          command::Parameter
          (
            "contains",
            "filename for input sdf of molecules to find. Molecules that match any of the molecules in the ensemble satisfy this criterion",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      ),
      m_HasPropertiesFlag
      (
        new command::FlagDynamic
        (
          "has_properties",
          "miscellaneous properties that matching molecules will already have",
          command::Parameter
          (
            "property",
            "miscellaneous property that matching molecules will already have"
          )
        )
      ),
      m_HasPropertyWithStringFlag
      (
        new command::FlagDynamic
        (
          "property_has_string",
          "",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "property",
              "miscellaneous property that matching molecules will have"
            ),
            command::Parameter
            (
              "target_string",
              "string value that the property should have for matched molecules"
            )
          )
        )
      ),
      m_HasPropertyContainingStringFlag
      (
        new command::FlagDynamic
        (
          "property_contains_string",
          "",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "property",
              "miscellaneous property that matching molecules will have"
            ),
            command::Parameter
            (
              "target_string",
              "string that the property should contain for matched molecules"
            )
          )
        )
      ),
      m_ComparePropertyValuesFlag
      (
        new command::FlagDynamic
        (
          "compare_property_values",
          "Match molecules with property values that satisfy the comparison",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "property",
              "miscellaneous property that matching molecules will have"
            ),
            command::Parameter
            (
              "comparison",
              "",
              command::ParameterCheckSerializable( math::Comparisons< float>::EnumType())
            ),
            command::Parameter
            (
              "target_string",
              "string value that the property should have for matched molecules"
            )
          )
        )
      ),
      m_AnyFlag
      (
        new command::FlagStatic
        (
          "any",
          "if any of the criteria specified is match, the molecule will be matched (by default, all criteria must match)"
        )
      ),
      m_AmideBondTolerance( 20.0)
    {
    }

    //! copy constructor; ignores everything but flags
    MoleculeFilter::MoleculeFilter( const MoleculeFilter &PARENT) :
      m_OutputFilenameForMoleculesSatisfyingCriteriaFlag( PARENT.m_OutputFilenameForMoleculesSatisfyingCriteriaFlag),
      m_OutputFilenameForMoleculesFailingCriteriaFlag( PARENT.m_OutputFilenameForMoleculesFailingCriteriaFlag),
      m_Match3DCoordinatesFlag( PARENT.m_Match3DCoordinatesFlag),
      m_MatchDefinedAtomTypesFlag( PARENT.m_MatchDefinedAtomTypesFlag),
      m_MatchPlanarAmideBondsFlag( PARENT.m_MatchPlanarAmideBondsFlag),
      m_MatchPlanarAromaticRingsFlag( PARENT.m_MatchPlanarAromaticRingsFlag),
      m_MatchSimpleMoleculesFlag( PARENT.m_MatchSimpleMoleculesFlag),
      m_MoleculesContainingFragmentsFilenameFlag( PARENT.m_MoleculesContainingFragmentsFilenameFlag),
      m_MoleculesContainingConformersFilenameFlag( PARENT.m_MoleculesContainingConformersFilenameFlag),
      m_ConformerComparerFlag( PARENT.m_ConformerComparerFlag),
      m_MoleculesToFindFilenameFlag( PARENT.m_MoleculesToFindFilenameFlag),
      m_HasPropertiesFlag( PARENT.m_HasPropertiesFlag),
      m_HasPropertyWithStringFlag( PARENT.m_HasPropertyWithStringFlag),
      m_HasPropertyContainingStringFlag( PARENT.m_HasPropertyContainingStringFlag),
      m_ComparePropertyValuesFlag( PARENT.m_ComparePropertyValuesFlag),
      m_AnyFlag( PARENT.m_AnyFlag),
      m_AmideBondTolerance( PARENT.m_AmideBondTolerance)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeFilter *MoleculeFilter::Clone() const
    {
      return new MoleculeFilter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeFilter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeFilter::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be filtered
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // add AtomMdlLine to molecule
      sdf::MdlHandler::AddAtomMdlLineFlag( *sp_cmd);

      // Where to place the molecules that satisfy all the criteria given beloe
      sp_cmd->AddFlag( m_OutputFilenameForMoleculesSatisfyingCriteriaFlag);

      // Where to place the molecules that fail any of the criteria given below
      sp_cmd->AddFlag( m_OutputFilenameForMoleculesFailingCriteriaFlag);

      // Basic criteria

      // Check that all tetrahedrally-hybridized atoms with 4 neighbors has at least one neighbor w/a z-coordinate
      sp_cmd->AddFlag( m_Match3DCoordinatesFlag);

      // Check that all atom types are defined
      sp_cmd->AddFlag( m_MatchDefinedAtomTypesFlag);

      // whether to match only simple molecules, e.g. not molecular complexes
      sp_cmd->AddFlag( m_MatchSimpleMoleculesFlag);

      // Check that amide bonds are planar
      sp_cmd->AddFlag( m_MatchPlanarAmideBondsFlag);

      // Check that aromatic rings that should be planar (e.g. not cyclophanes) are planar
      sp_cmd->AddFlag( m_MatchPlanarAromaticRingsFlag);

      // input file containing sdfs of fragments
      // Molecules with at least 1 of the fragments satisfy the criteria
      sp_cmd->AddFlag( m_MoleculesContainingFragmentsFilenameFlag);

      // input file containing sdfs of conformers
      sp_cmd->AddFlag( m_MoleculesContainingConformersFilenameFlag);

      // input method of comparison for conformers
      // Molecules within tolerance satisfy the criteria
      sp_cmd->AddFlag( m_ConformerComparerFlag);

      // Molecule criteria

      // input file containing sdfs of molecules
      // Molecules in the ensemble that are also in filename satisfy this criteria
      sp_cmd->AddFlag( m_MoleculesToFindFilenameFlag);

      // Property-based criteria

      // Molecules with these misc properties satisfy this criteria
      sp_cmd->AddFlag( m_HasPropertiesFlag);

      // 2-parameter: property-name(name of a misc property), string(the value of the property)
      // Molecules with a specified misc property (property-name) with a value string exactly equal to string satisfy this criteria
      // Example: -has_property_with_string MiscProperty("mGlur5 Category",1) "Inhibitor"
      sp_cmd->AddFlag( m_HasPropertyWithStringFlag);

      // 2-parameter: property-name(name of a misc property), string(the value of the property)
      // Molecules with a specified misc property (property-name) with a value string containing string satisfy this property
      // Example: -has_property_with_value MiscProperty("Atom_SigmaCharge",1) nan
      // finds all the molecules whose sigma charge (which must already be stored in the sdf file) had a nan in it
      sp_cmd->AddFlag( m_HasPropertyContainingStringFlag);

      // numeric value of property with only one value

      // 3 parameter: property-name (name of a misc property already in the sdf file),
      //              comparison operator to use (one of >,<,>=,<=,==,!=)
      //              value (any numeric value)
      // if the molecule has a misc property called property-name and its value satisfies the comparison operator when
      // compared to value, it satisfies this property
      sp_cmd->AddFlag( m_ComparePropertyValuesFlag);

      //! Any flag
      sp_cmd->AddFlag( m_AnyFlag);

    ///////////////////
    // default flags //
    ///////////////////

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeFilter::GetDescription() const
    {
      return "MoleculeFilter filters an ensemble of molecules by a user defined criteria";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeFilter::GetReadMe() const
    {
      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::MoleculeFilter, terms of use, "
        "appropriate citation, installation procedures, BCL::MoleculeFilter execution, "
        "technical support, and future research directions.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::MoleculeFilter?\n"
        "BCL::MoleculeFilter is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::MoleculeFilter filters an ensemble of molecules "
        "by a variety of criteria. A variety of comparison algorithms are available such as substructure search.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::MoleculeFilter.\n"
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
        "VI. RUNNING BCL::MoleculeFilter.\n"
        "Running BCL::MoleculeFilter requires an sdf file containing molecules that need to be filtered according to some properties.\n"
        "\n"
        "2) Run BCL::MoleculeFilter to filter molecules based on some properties\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::MoleculeFilter executable program. The syntax for"
        "running the application looks like the following"
        "\n"
        "bcl.exe MoleculeFilter -input_filenames <filename.sdf> -output_matched <filename>"
        "\n\nFLAGS:\n\n"
        "-input_filenames <filename> -> file containing ensemble of molecules that need to be identified uniquely\n"
        "-output_matched <filename> -> file to write out molecules that match desired criterion \n"
        "-output_unmatched <filename> -> file to write out molecules that do not match desired criterion \n"
        "-contains_fragments_from <filename> -> Molecule contains at least one of a fragment from a given file\n"
        "-contains_conformers <filename> <tolerance> -> Molecule matches conformer within tolerance of a given file\n"
        "-contains <filename> -> Molecule is identical to another molecule in a particular file\n"
        "-has_properties <property> -> Molecule contains a particular misc. property (e.g. mGluR5 Category)\n"
        "-property_has_string <property> <string> -> A particular misc property value compares equal to a given string (e.g. Inhibitor)\n"
        "-property_contains_string <property> <string> -> A particular misc property value contains a particular string (e.g. nan)\n"
        "-compare_property_values <property> <comparison> <string> -> A particular misc property value satisfies a comparison (e.g. EC50 < 4.0)\n"
        "-simple -> Filter molecules that are molecular complexes\n"
        "-3d -> Filter molecules that have 3d coordinates (if it needs them)\n"
        "-defined_atom_types -> filter molecules that have defined atom types\n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe MoleculeFilter -help\n"
        "\n"
        "For more general information about the product, type bcl.exe MoleculeFilter -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::MoleculeFilter.\n"
        "BCL::MoleculeFilter is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );

      return s_read_me;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &MoleculeFilter::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::MoleculeFilter filters an ensemble of molecules by a variety of criteria including \n"
        "<ul>"
        "<li>valid 3D structures</li>"
        "<li>valid gasteiger atom types</li>"
        "<li>descriptors criteria</li>"
        "<li>sub-structure search</li>"
        "<li>whole molecule search</li>"
        "</ul>\n\n"
        "Features of BCL::MoleculeFilter\n"
        "<ul>"
        "<li>multiple filters can be used at once</li>"
        "<li>filtered-out molecules can also be written out to a separate file</li>"
        "<li>ability to use novel BCL descriptors for filtering</li>"
        "<li>fast and accurate atom and bond-type based structure searches</li>"
        "<li>Compressed molecule files (bz2, gzip) are supported</li>"
        "</ul>\n\n"
        "!moleculefilter.png!\n"
      );

      return s_web_text;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeFilter::Main() const
    {
      Initialize();

      // create filter statistics for loading
      FilterStatistics loading_stats;
      FilterStatistics writing_satisfied_stats;
      FilterStatistics writing_unsatisfied_stats;

      loading_stats.StartFilter();
      chemistry::FragmentFeed feed;
      loading_stats.StopFilter( false);
      while( feed.NotAtEnd())
      {
        chemistry::FragmentComplete fragment( *feed);
        if( MatchesFilters( fragment))
        {
          if( m_OutputMatched.is_open())
          {
            writing_satisfied_stats.StartFilter();
            fragment.WriteMDL( m_OutputMatched);
            writing_satisfied_stats.StopFilter( true);
          }
        }
        else if( m_OutputUnmatched.is_open())
        {
          writing_unsatisfied_stats.StartFilter();
          fragment.WriteMDL( m_OutputUnmatched);
          writing_unsatisfied_stats.StopFilter( true);
        }
        loading_stats.StartFilter();
        ++feed;
        loading_stats.StopFilter( true);
      }

      // output statistics like:
      BCL_MessageStd
      (
        "Loaded " + util::Format()( loading_stats.GetMatchedCount())
        + " molecules in " + loading_stats.GetTimer().GetTotalTime().GetTimeAsHourMinuteSecond()
      );

      size_t last_filter_match_number( loading_stats.GetMatchedCount());
      for( size_t filter_id( 0); filter_id < FilterStatistics::s_NumberFilterStatistics; ++filter_id)
      {
        for
        (
          std::vector< FilterStatistics>::const_iterator
            itr_stat( m_FilterStatistics[ filter_id].begin()), itr_stat_end( m_FilterStatistics[ filter_id].end());
          itr_stat != itr_stat_end;
          ++itr_stat
        )
        {
          if( !m_AnyFlag->GetFlag())
          {
            BCL_MessageStd
            (
              util::Format()( itr_stat->GetMatchedCount()) + " molecules passed filter "
              + itr_stat->GetDescription()
              + " in " + itr_stat->GetTimer().GetTotalTime().GetTimeAsHourMinuteSecond()
              + ", " + util::Format()( last_filter_match_number - itr_stat->GetMatchedCount())
              + " additional molecules filtered out"
            );
          }
          else
          {
            BCL_MessageStd
            (
              util::Format()( itr_stat->GetMatchedCount()) + " molecules matched at filter "
              + itr_stat->GetDescription()
              + " in " + itr_stat->GetTimer().GetTotalTime().GetTimeAsHourMinuteSecond()
            );
          }
          last_filter_match_number = itr_stat->GetMatchedCount();
        }
      }

      if( writing_satisfied_stats.GetMatchedCount())
      {
        BCL_MessageStd
        (
          "Wrote " + util::Format()( writing_satisfied_stats.GetMatchedCount())
          + " molecules that matched all criteria in "
          + writing_satisfied_stats.GetTimer().GetTotalTime().GetTimeAsHourMinuteSecond()
        );
      }
      if( writing_unsatisfied_stats.GetMatchedCount())
      {
        BCL_MessageStd
        (
          "Wrote " + util::Format()( writing_unsatisfied_stats.GetMatchedCount())
          + " molecules that failed any criteria in "
          + writing_unsatisfied_stats.GetTimer().GetTotalTime().GetTimeAsHourMinuteSecond()
        );
      }

      // close the files
      io::File::CloseClearFStream( m_OutputMatched);
      io::File::CloseClearFStream( m_OutputUnmatched);

      // end
      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief prepare this object to handle the flags that were passed
    void MoleculeFilter::Initialize() const
    {
      // open output files
      if( m_OutputFilenameForMoleculesSatisfyingCriteriaFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          m_OutputMatched,
          m_OutputFilenameForMoleculesSatisfyingCriteriaFlag->GetFirstParameter()->GetValue()
        );
      }
      if( m_OutputFilenameForMoleculesFailingCriteriaFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          m_OutputUnmatched,
          m_OutputFilenameForMoleculesFailingCriteriaFlag->GetFirstParameter()->GetValue()
        );
      }
      // alert the user if no output file was given
      else if( !m_OutputFilenameForMoleculesSatisfyingCriteriaFlag->GetFlag())
      {
        BCL_MessageCrt( "No output filename given, only messages will be output");
      }

      for( size_t stat_number( 0); stat_number < FilterStatistics::s_NumberFilterStatistics; ++stat_number)
      {
        m_FilterStatistics[ stat_number].clear();
      }

      // add other filter statistics

      if( m_Match3DCoordinatesFlag->GetFlag())
      {
        m_FilterStatistics[ FilterStatistics::e_3d].push_back( FilterStatistics( FilterStatistics::e_3d));
      }
      if( m_MatchDefinedAtomTypesFlag->GetFlag())
      {
        m_FilterStatistics[ FilterStatistics::e_DefinedAtomTypes].push_back
        (
          FilterStatistics( FilterStatistics::e_DefinedAtomTypes)
        );
      }
      if( m_MatchPlanarAmideBondsFlag->GetFlag())
      {
        m_AmideBondTolerance = m_MatchPlanarAmideBondsFlag->GetNumericalList< double>()( 0);
        m_FilterStatistics[ FilterStatistics::e_PlanarAmideBonds].push_back
        (
          FilterStatistics( FilterStatistics::e_PlanarAmideBonds, " to within " + util::Format()( m_AmideBondTolerance) + " degrees")
        );
      }
      if( m_MatchPlanarAromaticRingsFlag->GetFlag())
      {
        m_FilterStatistics[ FilterStatistics::e_PlanarUnbridgedAromaticRings].push_back
        (
          FilterStatistics( FilterStatistics::e_PlanarUnbridgedAromaticRings)
        );
      }
      if( m_MatchSimpleMoleculesFlag->GetFlag())
      {
        m_FilterStatistics[ FilterStatistics::e_Simple].push_back
        (
          FilterStatistics( FilterStatistics::e_Simple)
        );
      }
      if( m_MoleculesContainingFragmentsFilenameFlag->GetFlag())
      {
        m_FilterStatistics[ FilterStatistics::e_ContainingFragments].push_back
        (
          FilterStatistics
          (
            FilterStatistics::e_ContainingFragments,
            m_MoleculesContainingFragmentsFilenameFlag->GetFirstParameter()->GetValue()
          )
        );
      }
      if( m_MoleculesToFindFilenameFlag->GetFlag())
      {
        m_FilterStatistics[ FilterStatistics::e_Matches].push_back
        (
          FilterStatistics
          (
            FilterStatistics::e_Matches,
            m_MoleculesToFindFilenameFlag->GetFirstParameter()->GetValue()
          )
        );
      }

      // split by the presence / absence of a property
      for
      (
        size_t property_number( 0), number_properties( m_HasPropertiesFlag->GetSize());
        property_number < number_properties;
        ++property_number
      )
      {
        const std::string property( m_HasPropertiesFlag->GetParameterList()( property_number)->GetValue());
        m_PropertyConstraints.PushBack
        (
          std::make_pair( property, util::StringReplacement( util::StringReplacement::e_Any, ""))
        );
        m_FilterStatistics[ FilterStatistics::e_PropertiesConstraint].push_back
        (
          FilterStatistics( FilterStatistics::e_PropertiesConstraint, property)
        );
      }

      if( m_HasPropertyWithStringFlag->GetSize() > 1)
      {
        // the argument list should be tiled, e.g. SomeProperty ItsValue SomeOtherProperty ItsValue, etc.
        for
        (
          size_t property_number( 0), number_properties( m_HasPropertyWithStringFlag->GetSize());
          property_number < number_properties - 1;
          property_number += 2
        )
        {
          const std::string property( m_HasPropertyWithStringFlag->GetParameterList()( property_number)->GetValue());
          const std::string value( m_HasPropertyWithStringFlag->GetParameterList()( property_number + 1)->GetValue());

          m_PropertyConstraints.PushBack
          (
            std::make_pair( property, util::StringReplacement( util::StringReplacement::e_Exact, value))
          );
          m_FilterStatistics[ FilterStatistics::e_PropertiesConstraint].push_back
          (
            FilterStatistics( FilterStatistics::e_PropertiesConstraint, property + " equal to " + value)
          );
        }
      }
      if( m_HasPropertyContainingStringFlag->GetSize() > 1)
      {
        // the argument list should be tiled, e.g. SomeProperty ItsValue SomeOtherProperty ItsValue, etc.
        for
        (
          size_t property_number( 0), number_properties( m_HasPropertyContainingStringFlag->GetSize());
          property_number < number_properties - 1;
          property_number += 2
        )
        {
          const std::string property( m_HasPropertyContainingStringFlag->GetParameterList()( property_number)->GetValue());
          const std::string value( m_HasPropertyContainingStringFlag->GetParameterList()( property_number + 1)->GetValue());

          m_PropertyConstraints.PushBack
          (
            std::make_pair( property, util::StringReplacement( util::StringReplacement::e_Any, value))
          );
          m_FilterStatistics[ FilterStatistics::e_PropertiesConstraint].push_back
          (
            FilterStatistics( FilterStatistics::e_PropertiesConstraint, property + " containing \"" + value + "\"")
          );
        }
      }
      if( m_ComparePropertyValuesFlag->GetSize() > 2)
      {
        // the argument list should be tiled, e.g. SomeProperty Comparison Value SomeOtherProperty Comparison ItsValue, etc.
        for
        (
          size_t property_number( 0), number_properties( m_ComparePropertyValuesFlag->GetSize());
          property_number < number_properties - 2;
          property_number += 3
        )
        {
          const std::string lhs( m_ComparePropertyValuesFlag->GetParameterList()( property_number)->GetValue());
          const std::string cmp( m_ComparePropertyValuesFlag->GetParameterList()( property_number + 1)->GetValue());
          const std::string rhs( m_ComparePropertyValuesFlag->GetParameterList()( property_number + 2)->GetValue());
          const descriptor::CheminfoProperty lhs_value_getter( lhs);
          const descriptor::CheminfoProperty rhs_value_getter( rhs);

          BCL_Assert
          (
            lhs_value_getter->GetSizeOfFeatures() == 1
            && rhs_value_getter->GetSizeOfFeatures() == 1,
            "one of the properties given to compare does not return just a single value"
          );

          m_PropertyComparisons.PushBack
          (
            storage::Triplet
            <
              descriptor::CheminfoProperty,
              math::Comparisons< float>::Comparison,
              descriptor::CheminfoProperty
            >( lhs_value_getter, math::Comparisons< float>::Comparison( cmp), rhs_value_getter)
          );
          m_FilterStatistics[ FilterStatistics::e_ComparePropertyValues].push_back
          (
            FilterStatistics( FilterStatistics::e_ComparePropertyValues, lhs + " " + cmp + " " + rhs)
          );
        }
      }

      // reset substructure graph arrays
      m_ExactMatches.Reset();
      m_SearchSubstructures.Reset();

      // substructure searches
      if( m_MoleculesContainingFragmentsFilenameFlag->GetFlag())
      {
        io::IFStream input; // stream for reading in fragments
        io::File::MustOpenIFStream
        (
          input,
          m_MoleculesContainingFragmentsFilenameFlag->GetFirstParameter()->GetValue()
        );
        // read in molecules
        chemistry::FragmentEnsemble ensemble
        (
          input,
          sdf::GetCommandLineHydrogensPref() == sdf::e_Remove ? sdf::e_Remove : sdf::e_Maintain
        );
        // close stream
        io::File::CloseClearFStream( input);

        // make graphs for each molecule
        m_SearchSubstructures.AllocateMemory( ensemble.GetSize());
        // make graphs for each molecule
        chemistry::ConformationGraphConverter::AtomComparisonTypeEnum atom_type
        (
          m_MoleculesContainingFragmentsFilenameFlag->GetParameterList()( 1)->GetValue()
        );
        chemistry::ConfigurationalBondTypeData::DataEnum bond_type
        (
          m_MoleculesContainingFragmentsFilenameFlag->GetParameterList()( 2)->GetValue()
        );
        chemistry::ConformationGraphConverter graph_maker
        (
          atom_type,
          bond_type
        );
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
          itr != itr_end;
          ++itr
        )
        {
          m_SearchSubstructures.PushBack( graph_maker( *itr)); // Substructures to look for
        }
        // close stream
        io::File::CloseClearFStream( input);
      }

      // conformer searches
      if( m_MoleculesContainingConformersFilenameFlag->GetFlag())
      {
        io::IFStream input; // stream for reading in conformers
        io::File::MustOpenIFStream
        (
          input,
          m_MoleculesContainingConformersFilenameFlag->GetFirstParameter()->GetValue()
        );
        // read in molecules
        m_ConformersToMatch.ReadMoreFromMdl
        (
          input,
          sdf::GetCommandLineHydrogensPref() == sdf::e_Remove ? sdf::e_Remove : sdf::e_Maintain
        );
        // close stream
        io::File::CloseClearFStream( input);
      }

      // exactly matching molecules to look for
      if( m_MoleculesToFindFilenameFlag->GetFlag())
      {
        io::IFStream input; // stream for reading in fragments
        io::File::MustOpenIFStream
        (
          input,
          m_MoleculesToFindFilenameFlag->GetFirstParameter()->GetValue()
        );

        // read in molecules
        chemistry::FragmentEnsemble ensemble( input, sdf::GetCommandLineHydrogensPref());
        // close stream
        io::File::CloseClearFStream( input);

        // make graphs for each molecule
        chemistry::ConformationGraphConverter graph_maker;
        size_t index( 0);
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
          itr != itr_end;
          ++itr, ++index
        )
        {
          // create a graph from the molecule
          const graph::ConstGraph< size_t, size_t> graph( graph_maker( *itr));

          // add the molecule to the exact matches map
          m_ExactMatches[ MakeHashStringFromGraph( graph)].PushBack
          (
            storage::Pair< size_t, graph::ConstGraph< size_t, size_t> >( index, graph)
          );
        }
      }
    }

    //! @brief make a hash string from a map of size-t's to size-t's
    //! @param MAP the map to make a hashable string from
    //! @return a string containing key,value,key,value pairs
    std::string MoleculeFilter::MakeHashStringFromMap( const storage::Map< size_t, size_t> &MAP)
    {
      if( MAP.GetSize() == 0)
      {
        return std::string();
      }

      // separator will go between each size_t printed in the map
      // since these are size_ts, it is fine to just use a space
      const char separator( ' ');

      // make a hash beginning with the first element
      std::string hash( util::Format()( MAP.Begin()->first) + separator + util::Format()( MAP.Begin()->second));

      // add all the rest of the elements to the map
      for
      (
        storage::Map< size_t, size_t>::const_iterator itr( ++MAP.Begin()), itr_end( MAP.End());
        itr != itr_end;
        ++itr
      )
      {
        hash += separator + util::Format()( itr->first) + separator + util::Format()( itr->second);
      }

      return hash;
    }

    //! @brief make a hash string from a graph
    //! @param GRAPH a graph
    //! This string can be used as a key to a map that holds graphs with identical hash strings, thus narrowing the
    //! number of graphs that must be searched to determine whether a new scaffold is unique
    //! @return a string that specifies something about the graph that is vertex-invariant
    //!         e.g. does not depend on the ordering of the vertices
    std::string MoleculeFilter::MakeHashStringFromGraph( const graph::ConstGraph< size_t, size_t> &GRAPH)
    {
      // make maps containing counts of each vertex and edge type
      const storage::Map< size_t, size_t> &vertex_counts( GRAPH.GetVertexTypeCountMap());
      const storage::Map< size_t, size_t> &edge_counts( GRAPH.GetEdgeTypeCountMap());

      // make hash strings from each map, separated by "$"
      return MakeHashStringFromMap( vertex_counts) + "$" + MakeHashStringFromMap( edge_counts);
    }

    //! @brief apply all the filters
    //! @return true if the molecule matched all the filters
    bool MoleculeFilter::MatchesFilters( chemistry::FragmentComplete &MOLECULE) const
    {
      // whether to require all criteria to match
      bool need_all_match( !m_AnyFlag->GetFlag());
      if( m_MatchDefinedAtomTypesFlag->GetFlag())
      {
        FilterStatistics &stats( m_FilterStatistics[ FilterStatistics::e_DefinedAtomTypes][ 0]);
        stats.StartFilter();
        bool had_undefined( MOLECULE.HasNonGasteigerAtomTypes());
        stats.StopFilter( !had_undefined);
        if( had_undefined == need_all_match)
        {
          return !had_undefined;
        }
      }

      if( m_Match3DCoordinatesFlag->GetFlag())
      {
        FilterStatistics &stats( m_FilterStatistics[ FilterStatistics::e_3d][ 0]);
        stats.StartFilter();
        bool matched( !MOLECULE.HasBadGeometry());
        stats.StopFilter( matched);
        if( matched == !need_all_match)
        {
          return matched;
        }
      }

      // planar amide bonds
      if( m_MatchPlanarAmideBondsFlag->GetFlag())
      {
        FilterStatistics &stats( m_FilterStatistics[ FilterStatistics::e_PlanarAmideBonds][ 0]);
        stats.StartFilter();
        bool had_undefined( !MOLECULE.AreAmideBondsPlaner( m_AmideBondTolerance));
        stats.StopFilter( !had_undefined);
        if( had_undefined == need_all_match)
        {
          return !had_undefined;
        }
      }

      // planar aromatic nonbridged rings
      if( m_MatchPlanarAromaticRingsFlag->GetFlag())
      {
        FilterStatistics &stats( m_FilterStatistics[ FilterStatistics::e_PlanarUnbridgedAromaticRings][ 0]);
        stats.StartFilter();
        bool had_undefined( !MOLECULE.AreAromaticRingsPlaner());
        stats.StopFilter( !had_undefined);
        if( had_undefined == need_all_match)
        {
          return !had_undefined;
        }
      }

      // check for molecular complexes
      if( m_MatchSimpleMoleculesFlag->GetFlag())
      {
        FilterStatistics &stats( m_FilterStatistics[ FilterStatistics::e_Simple][ 0]);
        stats.StartFilter();
        chemistry::ConformationGraphConverter graph_maker;
        bool matched( graph::Connectivity::IsConnected( graph_maker( MOLECULE)));
        stats.StopFilter( matched);
        if( !matched == need_all_match)
        {
          return matched;
        }
      }

      // filter by substructure
      if( m_MoleculesContainingFragmentsFilenameFlag->GetFlag())
      {
        FilterStatistics &stats( m_FilterStatistics[ FilterStatistics::e_ContainingFragments][ 0]);

        stats.StartFilter();
        chemistry::ConformationGraphConverter::AtomComparisonTypeEnum atom_type
        (
          m_MoleculesContainingFragmentsFilenameFlag->GetParameterList()( 1)->GetValue()
        );
        chemistry::ConfigurationalBondTypeData::DataEnum bond_type
        (
          m_MoleculesContainingFragmentsFilenameFlag->GetParameterList()( 2)->GetValue()
        );
        chemistry::ConformationGraphConverter conformation_graph_maker
        (
          atom_type,
          bond_type
        );

        // create a graph for this molecule
        graph::ConstGraph< size_t, size_t> mol_graph( conformation_graph_maker( MOLECULE));

        storage::Vector< double> contained_fragment_indices; // track indices of matching fragments
        size_t index( 0);

        graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
        isomorphism.SetGraphExternalOwnership( mol_graph);

        for
        (
          storage::Vector< graph::ConstGraph< size_t, size_t> >::iterator
            itr_frags( m_SearchSubstructures.Begin()), itr_frags_end( m_SearchSubstructures.End());
          itr_frags != itr_frags_end;
          ++itr_frags, ++index
        )
        {
          isomorphism.SetSubgraphExternalOwnership( *itr_frags);
          if( isomorphism.FindIsomorphism())
          {
            contained_fragment_indices.PushBack( double( index));
          }
        }
        const bool match( !contained_fragment_indices.IsEmpty());
        stats.StopFilter( match);
        if( match)
        {
          MOLECULE.StoreProperty( "MatchedFragmentIndices", contained_fragment_indices);
        }
        if( need_all_match != match)
        {
          return match;
        }
      }

      // filter by conformer
      if( m_MoleculesContainingConformersFilenameFlag->GetFlag())
      {
        bool match( false);
        util::Implementation< chemistry::ConformationComparisonInterface> comparer( m_ConformerComparerFlag->GetFirstParameter()->GetValue());
        math::Comparisons< float>::Comparison conf_comparison_type( m_ConformerComparerFlag->GetParameterList()( 1)->GetValue());
        const double tolerance( m_ConformerComparerFlag->GetParameterList().LastElement()->GetNumericalValue< double>());
        for
        (
            chemistry::FragmentEnsemble::const_iterator ens_itr( m_ConformersToMatch.Begin());
            ens_itr != m_ConformersToMatch.End();
            ++ens_itr
        )
        {
          // compute comparison
          float comparison_value( ( *comparer)( *ens_itr, MOLECULE));

          if( ( **conf_comparison_type)( comparison_value, tolerance))
          {
            match = true;
            break;
          }
        }
        if( match != need_all_match)
        {
          return match;
        }
      }

      // filter by exact structure
      if( m_MoleculesToFindFilenameFlag->GetFlag())
      {
        FilterStatistics &stats( m_FilterStatistics[ FilterStatistics::e_Matches][ 0]);

        stats.StartFilter();
        chemistry::ConformationGraphConverter conformation_graph_maker;

        // create a graph for this molecule
        graph::ConstGraph< size_t, size_t> mol_graph( conformation_graph_maker( MOLECULE));

        // create a hash string for this molecule
        const std::string hash( MakeHashStringFromGraph( mol_graph));

        storage::Map
        <
          std::string,
          storage::List< storage::Pair< size_t, graph::ConstGraph< size_t, size_t> > >
        >::iterator itr_hash_map( m_ExactMatches.Find( hash));

        if( itr_hash_map == m_ExactMatches.End())
        {
          stats.StopFilter( false);
          if( need_all_match)
          {
            return false;
          }
        }
        else
        {
          graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
          isomorphism.SetSubgraphExternalOwnership( mol_graph);

          // track indices of fragments that were matched
          storage::Vector< double> matched_indices;

          for
          (
            storage::List< storage::Pair< size_t, graph::ConstGraph< size_t, size_t> > >::iterator
              itr_frags( itr_hash_map->second.Begin()), itr_frags_end( itr_hash_map->second.End());
            itr_frags != itr_frags_end;
            ++itr_frags
          )
          {
            isomorphism.SetGraphExternalOwnership( itr_frags->Second());
            if( isomorphism.FindIsomorphism())
            {
              matched_indices.PushBack( double( itr_frags->First()));
            }
          }

          const bool match( !matched_indices.IsEmpty());
          if( match)
          {
            MOLECULE.StoreProperty( "MatchIndices", matched_indices);
          }
          stats.StopFilter( match);
          if( match != need_all_match)
          {
            return match;
          }
        }
      }

      if( m_PropertyConstraints.GetSize())
      {
        std::vector< FilterStatistics>::iterator itr_stats
        (
          m_FilterStatistics[ FilterStatistics::e_PropertiesConstraint].begin()
        );
        for
        (
          storage::List< std::pair< std::string, util::StringReplacement> >::const_iterator
            itr( m_PropertyConstraints.Begin()), itr_end( m_PropertyConstraints.End());
          itr != itr_end;
          ++itr, ++itr_stats
        )
        {
          itr_stats->StartFilter();

          if
          (
            !MOLECULE.IsPropertyStored( itr->first)
            ||
            itr->second.FindNextMatch( util::TrimString( MOLECULE.GetMDLProperty( itr->first)), 0)
              == std::string::npos
          )
          {
            // molecule did not match the constraint
            itr_stats->StopFilter( false);
            if( need_all_match)
            {
              return false;
            }
          }
          else
          {
            itr_stats->StopFilter( true);
            if( !need_all_match)
            {
              return true;
            }
          }
        }
      }

      if( m_PropertyComparisons.GetSize())
      {
        std::vector< FilterStatistics>::iterator itr_stats
        (
          m_FilterStatistics[ FilterStatistics::e_ComparePropertyValues].begin()
        );
        for
        (
          storage::List
          <
            storage::Triplet
            <
              descriptor::CheminfoProperty,
              math::Comparisons< float>::Comparison,
              descriptor::CheminfoProperty
            >
          >::iterator itr( m_PropertyComparisons.Begin()), itr_end( m_PropertyComparisons.End());
          itr != itr_end;
          ++itr, ++itr_stats
        )
        {
          itr_stats->StartFilter();

          // get the property on the left hand side of the comparison
          const linal::Vector< float> value_lhs( itr->First()->SumOverObject( MOLECULE));
          const linal::Vector< float> value_rhs( itr->Third()->SumOverObject( MOLECULE));

          if
          (
            value_lhs.GetSize() != size_t( 1) || value_rhs.GetSize() != size_t( 1)
            || !( **itr->Second())( value_lhs( 0), value_rhs( 0))
          )
          {
            // do not match if the properties failed or the comparison failed
            itr_stats->StopFilter( false);
            if( need_all_match)
            {
              return false;
            }
          }
          else
          {
            itr_stats->StopFilter( true);
            if( !need_all_match)
            {
              return true;
            }
          }
        }
      }

      // all filters passed (if need all match was set) or none of them passed (if -any was set)
      return need_all_match;
    }

    // static instance
    const ApplicationType MoleculeFilter::MoleculeFilter_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeFilter(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
