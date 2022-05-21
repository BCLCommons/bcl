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
#include "bcl_app_conformer_generator.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_rmsd.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_symmetry_rmsd.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_ligand_design_helper.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "math/bcl_math_statistics.h"
#include "mc/bcl_mc_temperature_default.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_template_instantiations.h"

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    ConformerGenerator::ConformerGenerator() :
      m_RotamerLibraryFlag
      (
        new command::FlagStatic
        (
          "rotamer_library",
          "flag for rotamer library to use",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::RotamerLibraryInterface>()),
            "cod"
          )
        )
      ),
      m_NativeEnsembleFileFlag
      (
        new command::FlagDynamic
        (
          "natives",
          "sdf file containing the conformations that need to be compared to generated conformers (saame order as ensemble)",
          command::Parameter
          (
            "native",
            "an sdf input file"
          ),
          0,
          1
        )
      ),
      m_ConformerComparerFlag
      (
        new command::FlagStatic
        (
          "conformation_comparer",
          "method to compare conformers with",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "rmsd_type",
                "method to use to find rmsd between conformers. If it takes too long use RMSD or Diheral(method=Max)",
                command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>()),
                "SymmetryRMSD"
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "tolerance",
                "amount of tolerance allowed between two conformers. If set to 0 and -cluster is not set, skip checking for redundancy (faster)",
                command::ParameterCheckRanged< double>( 0.0, std::numeric_limits< double>::max()),
                "0.0"
              )
            )
          )
        )
      ),
      m_MaxIterations
      (
        new command::FlagStatic
        (
          "max_iterations",
          "maximum number of iterations for searching conformations",
          command::Parameter
          (
            "max_iterations",
            "maximum number of iterations for searching conformations",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            "8000"
          )
        )
      ),
      m_OutputLigandData
      (
        new command::FlagDynamic
        (
          "output_ligand_data",
          "For debugging",
          command::Parameter
          (
            "output_ligand_data",
            "filename base to output fragments and conformations for every ligand whose conformations need to be sampled"
          ),
          0,
          1
        )
      ),
      m_TopModels
      (
        new command::FlagStatic
        (
          "top_models",
          "number of conformers required",
          command::Parameter
          (
            "top_models",
            "number of conformers required",
            command::ParameterCheckRanged< size_t>( size_t( 1), math::GetHighestBoundedValue< size_t>()),
            "100"
          )
        )
      ),
      m_OutputFailed3D
      (
        new command::FlagDynamic
        (
          "failed_3D",
          "filename base to output molecules whose 3D structures could not be generated",
          command::Parameter
          (
            "failed_3D",
            "filename base to output molecules whose 3D structures could not be generated"
          ),
          0,
          1
        )
      ),
      m_RmsdScorePrefix
      (
        new command::FlagDynamic
        (
          "output_rmsd_score",
          "filename base to output a file containing rmsd to native and score of conformers",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "output_rmsd_score",
              "filename to output sdf containing rmsd to native and score of conformers"
            ),
            command::Parameter
            (
              "rmsd_type",
              "method to use to find rmsd between conformers",
              command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>()),
              "SymmetryRMSD"
            )
          ),
          0,
          2
        )
      ),
      m_ConformersSingleFile
      (
        new command::FlagDynamic
        (
          "conformers_single_file",
          "output all conformations of ensemble to one output file",
          command::Parameter
          (
            "conformers_single_file",
            "output all conformations of ensemble to one output file"
          ),
          0,
          1
        )
      ),
      m_ConformersSeparateFile
      (
        new command::FlagDynamic
        (
          "conformers_separate_files",
          "this flag will generate separate sdf file for each ligand",
          command::Parameter
          (
            "conformers_separate_files",
            "this flag will generate separate sdf file for each ligand"
          ),
          0,
          1
        )
      ),
      m_FixedNumberConformers
      (
        new command::FlagDynamic
        (
          "minimum_number_conformations",
          "if one wants to output a fixed number of conformers",
          command::Parameter
          (
            "minimum_number_conformations",
            "if one wants to output a fixed number of conformers",
            "5"
          )
        )
      ),
      m_ChangeChirality
      (
        new command::FlagStatic
        (
          "change_chirality",
          "flag to specify if chirality and isometry should be sampled as well during conformation sampling"
          "By default, the provided isometry/chirality in the input structure is maintained"
        )
      ),
      m_Generate3DCoordinates
      (
        new command::FlagStatic
        (
          "generate_3D",
          "Use rotamer library instead to generate initial conformation, ignore input structure and bond lengths / angles. "
        )
      ),
      m_Cluster
      (
        new command::FlagStatic
        (
          "cluster",
          "(ignored; present for backwards compatibility) To cluster conformations - so as to improve the top_models coverage "
          "conformational space, set this flag. By default"
          "the most likely conformations are instead selected. Clustering is always beneficial unless the application is to "
          "generate a boltzmann-like distribution of conformer states. Rather than setting -skip_cluster, it is recommended to "
          "reduce the number of iterations performed, which will as a side effect speed up the clustering"
        )
      ),
      m_NoCluster
      (
        new command::FlagStatic
        (
          "skip_cluster",
          "Disable clustering, appropriate when a Boltzmann-like conformer ensemble is desired rather than one that beest covers "
          "likely ligand conformational space"
        )
      ), 
      m_RandomDihedralMutateWeight
      (
        new command::FlagStatic
        (
          "rnd_dihedral_mutate_weight",
          "Relative weight for random mutation of dihedral angles (the bulk of the weight is made up of fragment-based mutates)",
          command::Parameter
          (
            "mutate_weight",
            "weight of random dihedral mutates relative to fragment-propensity directed mutates",
            command::ParameterCheckRanged< double>( 0.0, 1000.0),
            "0.0"
          )
        )
      ),
      m_SkipDihedralSamplingFlag
      (
        new command::FlagStatic
        (
          "skip_rotamer_dihedral_sampling",
          "Skip dihedral sampling according to rotamers, but still allow wiggling dihedrals up to +/- 30 degrees around the "
          "input conformation"
        )
      ),
      m_SkipBondAnglesFlag
      (
        new command::FlagStatic
        (
          "skip_bond_angle_sampling",
          "Skip bond angle and bond length sampling. On average (with clustering), this speeds up sampling ~10% at a "
          "minor cost in sampling coverage"
        )
      ),
      m_SkipRingSamplingFlag
      (
        new command::FlagStatic
        (
          "skip_ring_sampling",
          "Skip ring conformation sampling"
        )
      ),
      m_MaxClashResolutionIterationsFlag
      (
        new command::FlagStatic
        (
          "clash_resolution",
          "maximum number of tries  (times number of dihedral angles plus number of bond angles) to resolve clashes per"
          "molecule, before clustering). Larger numbers slow the calculation and make it more likely to find kinetically "
          "inaccessible conformational space, but a value of 0 makes it hard to sample conformations of strained molecules, so "
          "a value around 0.05-0.5 is usually best",
          command::Parameter
          (
            "max_clash_resolution_iterations",
            "maximum number of iterations for resolving clashes",
            command::ParameterCheckRanged< double>( 0, 10),
            "0.1"
          )
        )
      ),
      m_MaxClashToleranceFlag
      (
        new command::FlagStatic
        (
          "max_clash_tolerance",
          "maximum average angstroms clash present across all atoms in the molecule. A clash is defined as atoms more than "
          "three atoms separated that are closer than the sum of their vdw radii - 0.7 Angstroms (see paper)",
          command::Parameter
          (
            "max_clash_tolerance",
            "maximum average angstroms clash present across all atoms in the molecule",
            command::ParameterCheckRanged< double>( 0.0, 1.0),
            "0.1"
          )
        )
      )
    {
    }

      //! copy constructor
      ConformerGenerator::ConformerGenerator( const ConformerGenerator &APP) :
        m_RotamerLibraryFlag( APP.m_RotamerLibraryFlag),
        m_NativeEnsembleFileFlag( APP.m_NativeEnsembleFileFlag),
        m_ConformerComparerFlag( APP.m_ConformerComparerFlag),
        m_MaxIterations( APP.m_MaxIterations),
        m_OutputLigandData( APP.m_OutputLigandData),
        m_TopModels( APP.m_TopModels),
        m_OutputFailed3D( APP.m_OutputFailed3D),
        m_RmsdScorePrefix( APP.m_RmsdScorePrefix),
        m_ConformersSingleFile( APP.m_ConformersSingleFile),
        m_ConformersSeparateFile( APP.m_ConformersSeparateFile),
        m_FixedNumberConformers( APP.m_FixedNumberConformers),
        m_ChangeChirality( APP.m_ChangeChirality),
        m_Generate3DCoordinates( APP.m_Generate3DCoordinates),
        m_Cluster( APP.m_Cluster),
        m_NoCluster( APP.m_NoCluster),
        m_RandomDihedralMutateWeight( APP.m_RandomDihedralMutateWeight),
        m_SkipDihedralSamplingFlag( APP.m_SkipDihedralSamplingFlag),
        m_SkipBondAnglesFlag( APP.m_SkipBondAnglesFlag),
        m_SkipRingSamplingFlag( APP.m_SkipRingSamplingFlag),
        m_MaxClashResolutionIterationsFlag( APP.m_MaxClashResolutionIterationsFlag),
        m_MaxClashToleranceFlag( APP.m_MaxClashToleranceFlag)
      {
      }

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      ConformerGenerator *ConformerGenerator::Clone() const
      {
        return new ConformerGenerator( *this);
      }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConformerGenerator::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ConformerGenerator::GetDescription() const
    {
      return "Conformer generator generates small molecule conformations for ensemble of molecules that are provided";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ConformerGenerator::GetReadMe() const
    {
      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::ConformerGenerator, terms of use, "
        "appropriate citation, installation procedures, BCL::ConformerGenerator execution, "
        "technical support, and future research directions.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::ConformerGenerator?\n"
        "BCL::ConformerGenerator is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::ConformerGenerator is a utility that generates "
        "3D conformations for molecules of interest."
        "\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::ConformerGenerator.\n"
        "\n"
        ""
        ""
        ""
        ""
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::ConformerGenerator.\n"
        "Running BCL::ConformerGenerator requires an sdf file containing molecules whose conformations need to be sampled.\n"
        ""
        "\n"
        "2) Run BCL::ConformerGenerator to fragment molecules\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::ConformerGenerator executable program. The syntax for"
        "running the application looks like the following"
        "\n"
        "bcl.exe ConformerGenerator -ensemble_filename <filename.sdf> -conformers_single_file <filename> -rotamer_library <file or db>"
        "\n\nFLAGS:\n\n"
        "-ensemble_filenames <filename> -> file containing ensemble of molecules whose conformations need to be sampled\n"
        "-conformers_single_file <filename> -> file to which conformers will be written out\n"
        "-rotamer_library <file or path>-> specify file or db containing fragment rotamers"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe ConformerGenerator -help\n"
        "\n"
        "For more general information about the product, type bcl.exe ConformerGenerator -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::ConformerGenerator.\n"
        "BCL::ConformerGenerator is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );
      return s_read_me;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &ConformerGenerator::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::ConformerGenerator is a utility for conformation sampling of small molecules"
        "Features of BCL::ConformerGenerator\n"
        "<ul>"
        "  <li>Conformation sampling using rotamer library derived from strucure database\n</li>"
        "  <li>Output conformations for an ensemble of molecules in a single file or separate file for each molecule\n</li>"
        "  <li>Compressed molecule files (bz2, gzip) are supported\n</li>"
        "</ul>\n\n"
        "!conformer_generator.png!"
      );

      return s_web_text;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ConformerGenerator::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd, "ensemble");

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // insert all the flags and params
      sp_cmd->AddFlag( m_NativeEnsembleFileFlag);
      sp_cmd->AddFlag( m_ConformerComparerFlag);
      sp_cmd->AddFlag( m_MaxIterations);
      sp_cmd->AddFlag( m_OutputLigandData);
      sp_cmd->AddFlag( m_TopModels);
      sp_cmd->AddFlag( m_OutputFailed3D);
      sp_cmd->AddFlag( m_RmsdScorePrefix);
      sp_cmd->AddFlag( m_ConformersSingleFile);
      sp_cmd->AddFlag( m_ConformersSeparateFile);
      sp_cmd->AddFlag( m_FixedNumberConformers);
      sp_cmd->AddFlag( m_ChangeChirality);
      sp_cmd->AddFlag( m_Generate3DCoordinates);
      sp_cmd->AddFlag( m_Cluster);
      sp_cmd->AddFlag( m_NoCluster);
      sp_cmd->AddFlag( m_RandomDihedralMutateWeight);
      sp_cmd->AddFlag( m_SkipDihedralSamplingFlag);
      sp_cmd->AddFlag( m_SkipBondAnglesFlag);
      sp_cmd->AddFlag( m_SkipRingSamplingFlag);
      sp_cmd->AddFlag( m_MaxClashResolutionIterationsFlag);
      sp_cmd->AddFlag( m_MaxClashToleranceFlag);
      sp_cmd->AddFlag( m_RotamerLibraryFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief function that writes rmsd and score data to require file stream
    //! @param MOLECULE molecule of interest whose conformations need to be sampled
    //! @param MOLECULE_INDEX index of molecule of interest in ensemble
    //! @param ENSEMBLE ensemble that contains conformations of molecule of interest
    void ConformerGenerator::RmsdScoreFile
    (
      const chemistry::FragmentComplete &MOLECULE,
      const size_t &MOLECULE_INDEX,
      const chemistry::FragmentEnsemble &ENSEMBLE
    ) const
    {

      // create streams to store rmsd and score
      std::stringstream rmsd_string;
      std::stringstream bin_string;
      std::stringstream score_string;
      std::stringstream score_rmsd;

      chemistry::ConformationComparisonByDihedralBins comparison_bin( 30.0);
      float number_rotatable_bonds( descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( MOLECULE)( 0));
      rmsd_string << util::Format()( MOLECULE_INDEX)<< '-' << util::Format()( number_rotatable_bonds) << ' ';
      bin_string << util::Format()( MOLECULE_INDEX) << '-' << util::Format()( number_rotatable_bonds) << ' ';
      score_string << util::Format()( MOLECULE_INDEX) << '-' << util::Format()( number_rotatable_bonds) << ' ';
      score_rmsd << " Molecule " + util::Format()( MOLECULE_INDEX) << '\n';

      util::Implementation< chemistry::ConformationComparisonInterface> rmsd_calculator( m_RmsdScorePrefix->GetParameterList()( 1)->GetValue());

      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End());
        itr != itr_end;
        ++itr
      )
      {
        // find rmsd between molecule conformation provided and conformation. Update stream
        std::string rmsd_value;
        std::string bin_value;

        // get score for conformation and update the score stream
        std::string score_value( util::Format()( itr->GetStoredProperties().GetMDLProperty( "ConfScore")));

        if( itr->GetNumberHydrogens())
        {
          chemistry::FragmentComplete noh( *itr);
          noh.RemoveH();
          rmsd_value = util::Format()( ( *rmsd_calculator)( m_NativeEnsemble( MOLECULE_INDEX), noh));
          bin_value = util::Format()( comparison_bin( m_NativeEnsemble( MOLECULE_INDEX), noh));
          score_rmsd << util::Format()( comparison_bin( m_NativeEnsemble( MOLECULE_INDEX), noh)) << ' ' << rmsd_value << ' ' << score_value << '\n';
        }
        else
        {
          rmsd_value = util::Format()( ( *rmsd_calculator)( m_NativeEnsemble( MOLECULE_INDEX), *itr));
          bin_value = util::Format()( comparison_bin( m_NativeEnsemble( MOLECULE_INDEX), *itr));
          score_rmsd << util::Format()( comparison_bin( m_NativeEnsemble( MOLECULE_INDEX), *itr)) << ' ' << rmsd_value << ' ' << score_value << '\n';
          itr->GetChangeSignal().Emit( *itr);
        }

        rmsd_string << rmsd_value << ' ';
        bin_string << bin_value << ' ';
        score_string << score_value << ' ';
      }

      m_RmsdFile << rmsd_string.str() << '\n';
      m_DihedralBinFile << bin_string.str() << '\n';
      m_RmsdScore << score_rmsd.str() << '\n';
    }

    //! @brief function that creates pymol representation of fragments and conformations of molecule of interest
    //! @param MOLECULE molecule of interest whose conformations need to be sampled
    //! @param MOLECULE_INDEX index of molecule of interest in ensemble
    //! @param ENSEMBLE ensemble that contains conformations of molecule of interest
    //! @param FRAGMENT_ISO fragments that are part of molecule of interest
    void ConformerGenerator::PymolRepresentation
    (
      const chemistry::FragmentComplete &MOLECULE,
      const size_t &MOLECULE_INDEX,
      const chemistry::FragmentEnsemble &ENSEMBLE,
      const chemistry::FragmentEnsemble &FRAGMENT_ENSEMBLE
    ) const
    {
      // create directory where conformations, fragments and their pymol script files will be written
      std::string file_prefix( m_OutputLigandData->GetFirstParameter()->GetValue());
      io::Directory new_directory( file_prefix + "_" + util::Format()( MOLECULE_INDEX));
      new_directory.Make();

      // open the output file for writing out conformations of the molecule
      io::OFStream conformers;
      io::File::MustOpenOFStream( conformers, new_directory.GetPath() + '/' + "conformers" + util::Format()( MOLECULE_INDEX) + ".sdf");
      ENSEMBLE.WriteMDL( conformers);
      io::File::CloseClearFStream( conformers);

      io::OFStream molecule;
      io::File::MustOpenOFStream( molecule, new_directory.GetPath() + '/' + "molecule" + util::Format()( MOLECULE_INDEX) + ".sdf");
      MOLECULE.WriteMDL( molecule);
      io::File::CloseClearFStream( molecule);

      io::OFStream output_fragments;
      io::File::MustOpenOFStream( output_fragments, new_directory.GetPath() + '/' + "all_fragments" + util::Format()( MOLECULE_INDEX) + ".sdf");
      FRAGMENT_ENSEMBLE.WriteMDL( output_fragments);
      io::File::CloseClearFStream( output_fragments);

      io::OFStream pymol_output;
      io::File::MustOpenOFStream( pymol_output, new_directory.GetPath() + "_" + util::Format()( MOLECULE_INDEX) + ".pml");
      pymol_output << "load " + new_directory.GetPath() + '/' + "molecule" + util::Format()( MOLECULE_INDEX) + ".sdf" << '\n';
      pymol_output << "set valence, 1" << '\n';

      // get priority dihedral angles for molecule of interest
      const storage::Vector< storage::VectorND< 4, size_t> > dihedral_bonds_priority
      (
        chemistry::PriorityDihedralAngles()( MOLECULE).Second()
      );

      // highlight non ring priority dihedral bonds in pymol script
      size_t dihedral_bond( 0);
      for
      (
        storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
          itr( dihedral_bonds_priority.Begin()), itr_end( dihedral_bonds_priority.End());
        itr != itr_end;
        ++itr, ++dihedral_bond
      )
      {
        const size_t atom_a( itr->First() + 1);
        const size_t atom_b( itr->Second()+ 1);
        const size_t atom_c( itr->Third()+ 1);
        const size_t atom_d( itr->Fourth() + 1);
        pymol_output << "dihedral dihedral_angle_" << dihedral_bond;
        pymol_output << ", id " << util::Format()( atom_a) << ", id " << util::Format()( atom_b) << ", id " << util::Format()( atom_c) << ", id " << util::Format()( atom_d) << '\n';
      }
      pymol_output << "save " << "conformers" << "_" << util::Format()( MOLECULE_INDEX) << ".pse" << "\n";
      io::File::CloseClearFStream( pymol_output);

    }

    //! @brief controls output from the app of interest
    //! @param MOLECULES the molecule of interest whose conformations have to be sampled
    //! @param ENSEMBLE ensemble of conformations that were sampled for the molecule of interest
    //! @param NUMBER_OF_CONFORMATIONS number of conformations desired
    //! @param MOLECULE_INDEX index of the molecule in the ensemble
    void ConformerGenerator::Output
    (
      const chemistry::FragmentComplete MOLECULE,
      const chemistry::FragmentEnsemble &ENSEMBLE,
      const size_t NUMBER_OF_CONFORMATIONS,
      const size_t MOLECULE_INDEX
    ) const
    {
      size_t number_conformations( NUMBER_OF_CONFORMATIONS);

      util::ShPtrList< chemistry::FragmentComplete> final_list;
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
          itr_mols( ENSEMBLE.Begin()), itr_mols_end( ENSEMBLE.End());
        itr_mols != itr_mols_end;
        ++itr_mols
      )
      {
        final_list.PushBack( util::ShPtr< chemistry::FragmentComplete>( itr_mols->Clone()));
      }

      if( m_FixedNumberConformers->GetFlag())
      {
        number_conformations = util::ConvertStringToNumericalValue< size_t>( m_FixedNumberConformers->GetFirstParameter()->GetValue());
      }

      if( m_FixedNumberConformers->GetFlag())
      {
        if( final_list.GetSize() < number_conformations)
        {
          size_t current_size( final_list.GetSize());
          for( ; current_size != number_conformations; ++current_size)
          {
            final_list.PushBack( *final_list.Begin());
          }
        }
      }
      number_conformations = std::min( number_conformations, final_list.GetSize());

      size_t count( 0);

      // generate rmsd and score file for conformations
      if( m_ConformersSingleFile->GetFlag())
      {
        for
        (
          util::ShPtrList< chemistry::FragmentComplete>::const_iterator itr( final_list.Begin());
          count < number_conformations;
          ++itr, ++count
        )
        {
          ( *itr)->WriteMDL( m_OutputConformers);
        }
      }
      else
      {
        std::string file_prefix;
        if( m_ConformersSeparateFile->GetFlag())
        {
          file_prefix = m_ConformersSeparateFile->GetFirstParameter()->GetValue();
        }
        else
        {
          file_prefix = "top_models";
        }
        io::OFStream output;
        io::File::MustOpenOFStream( output, file_prefix + "_" + util::Format()( MOLECULE_INDEX) + ".sdf.gz");
        for
        (
          util::ShPtrList< chemistry::FragmentComplete>::const_iterator itr( final_list.Begin());
            count < number_conformations;
          ++itr, ++count
        )
        {
          ( *itr)->WriteMDL( output);
        }
        io::File::CloseClearFStream( output);
      }
    }

    //! @brief generates conformations for a molecule of interest
    //! @param MOLECULES the molecule of interest whose conformations have to be sampled
    //! @param MOLECULE_INDEX index of the molecule in the ensemble
    void ConformerGenerator::ConformationSampling
    (
      const chemistry::FragmentComplete &MOLECULE,
      const size_t &MOLECULE_INDEX
    ) const
    {
      if( m_OutputLigandData->GetFlag())
      {
        std::string file_prefix( m_OutputLigandData->GetFirstParameter()->GetValue());
        io::Directory new_directory( file_prefix + "_" + util::Format()( MOLECULE_INDEX));
        new_directory.Make();
        m_SampleConformations->SetOutputScoreFilename( new_directory.GetPath() + "/scored_fragments.sdf");
      }
      storage::Pair< chemistry::FragmentEnsemble, chemistry::FragmentEnsemble> conformation_result( m_SampleConformations->operator ()( MOLECULE));

      chemistry::FragmentEnsemble &ensemble( conformation_result.First());
      chemistry::FragmentEnsemble &fragment_ensemble( conformation_result.Second());

      if( ( !m_RmsdScorePrefix->GetFlag() && m_TopModels->GetFlag()) || m_FixedNumberConformers->GetFlag() || m_ConformersSingleFile->GetFlag())
      {
        Output( MOLECULE, ensemble, ensemble.GetSize(), MOLECULE_INDEX);
      }
      // generate rmsd and score file for conformations
      if( m_RmsdScorePrefix->GetFlag())
      {
        RmsdScoreFile( MOLECULE, MOLECULE_INDEX, ensemble);
      }
      // write out pymol representations
      if( m_OutputLigandData->GetFlag())
      {
        PymolRepresentation( MOLECULE, MOLECULE_INDEX, ensemble, fragment_ensemble);
      }
      // if conformers needed to be written for each molecule separately
      if( m_ConformersSeparateFile->GetFlag())
      {
        std::string file_prefix( m_ConformersSeparateFile->GetFirstParameter()->GetValue());

        // open the output file for writing out conformations of the molecule
        io::OFStream output;
        io::File::MustOpenOFStream( output, file_prefix + "_" + util::Format()( MOLECULE_INDEX) + ".sdf.gz");
        ensemble.WriteMDL( output);
        io::File::CloseClearFStream( output);
      }
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ConformerGenerator::Main() const
    {
      BCL_Assert( !m_Cluster->GetFlag() || !m_NoCluster->GetFlag(), "-skip_cluster and -cluster are mutually exclusive!");
      // get a feed to get molecules that need to be fragmented
      chemistry::FragmentFeed feed( "ensemble");

      util::Implementation< chemistry::RotamerLibraryInterface> rotamer_library( m_RotamerLibraryFlag->GetFirstParameter()->GetValue());

      io::IFStream input;

      // get native conformations from flag
      if( m_NativeEnsembleFileFlag->GetFlag())
      {
        io::File::MustOpenIFStream( input, m_NativeEnsembleFileFlag->GetFirstParameter()->GetValue());
        chemistry::FragmentEnsemble native_ensemble( input, sdf::e_Remove);
        m_NativeEnsemble = storage::Vector< chemistry::FragmentComplete>( native_ensemble.Begin(), native_ensemble.End());
      }
      io::File::CloseClearFStream( input);
      size_t number_conformations( util::ConvertStringToNumericalValue< size_t>( m_TopModels->GetFirstParameter()->GetValue()));
      if( m_FixedNumberConformers->GetFlag())
      {
        number_conformations = util::ConvertStringToNumericalValue< size_t>( m_FixedNumberConformers->GetFirstParameter()->GetValue());
      }

      // get comparer for comparing conformations
      m_ConformationComparison = m_ConformerComparerFlag->GetParameterList().LastElement()->GetNumericalValue< double>();
      chemistry::RotamerLibraryFile empty_rotlib;
      // initialize object for searching fragments of desired molecules
      m_SampleConformations = util::ShPtr< chemistry::SampleConformations>
      (
        new chemistry::SampleConformations
        (
          rotamer_library.IsDefined() ? *rotamer_library : empty_rotlib,
          m_ConformerComparerFlag->GetFirstParameter()->GetValue(),
          m_ConformationComparison,
          number_conformations,
          util::ConvertStringToNumericalValue< size_t>( m_MaxIterations->GetFirstParameter()->GetValue()),
          m_ChangeChirality->GetFlag(),
          util::ConvertStringToNumericalValue< double>( m_RandomDihedralMutateWeight->GetFirstParameter()->GetValue()),
          m_Generate3DCoordinates->GetFlag(),
          m_MaxClashToleranceFlag->GetFirstParameter()->GetNumericalValue< double>(), // clash tolerance
          !m_NoCluster->GetFlag(),
          m_MaxClashResolutionIterationsFlag->GetFirstParameter()->GetNumericalValue< double>(), // max clash resolution cycles
          "" // filename for all scores to be written out to
        )
      );
      m_SampleConformations->SetSamplingPreferences
      (
        !m_SkipDihedralSamplingFlag->GetFlag(),
        !m_SkipRingSamplingFlag->GetFlag(),
        !m_SkipBondAnglesFlag->GetFlag(),
        m_ChangeChirality->GetFlag()
      );
      // if rmsd and score is desired then initialize the required objects
      if( m_RmsdScorePrefix->GetFlag())
      {
        if( !m_NativeEnsembleFileFlag->GetFlag())
        {
          BCL_Assert( m_NativeEnsemble.IsEmpty() == 0, "Native file not provided for comparison");
        }

        std::string file_prefix( m_RmsdScorePrefix->GetFirstParameter()->GetValue());
        io::File::MustOpenOFStream( m_RmsdFile, file_prefix + "_R.txt");
        io::File::MustOpenOFStream( m_DihedralBinFile, file_prefix + "_DB.txt");
        io::File::MustOpenOFStream( m_RmsdScore, file_prefix + "_RS.txt");
      }

      // output conformations of all molecules in a single file
      if( m_ConformersSingleFile->GetFlag())
      {
        std::string file_prefix( m_ConformersSingleFile->GetFirstParameter()->GetValue());

        // open the output file for writing out conformations of the molecule
        io::File::MustOpenOFStream( m_OutputConformers, file_prefix);
      }

      if( m_OutputFailed3D->GetFlag())
      {
        std::string file_prefix( m_OutputFailed3D->GetFirstParameter()->GetValue());

        // open the output file for writing out conformations of the molecule
        io::File::MustOpenOFStream( m_Failed3DStream, file_prefix + ".sdf");
      }

      // keep track of which molecule we are on
      for( size_t ensemble_number( 0); feed.NotAtEnd(); ++feed, ++ensemble_number)
      {
        //TODO use status
        // sample conformations
        BCL_Message( util::Message::e_Critical, " iterating through ensemble in ANC - " + util::Format()( ensemble_number));
        ConformationSampling( *feed, ensemble_number);
      }

      // close the output file streams
      if( m_OutputLigandData->GetFlag())
      {
        io::File::CloseClearFStream( m_OutputSatisfied);
      }

      if( m_RmsdScorePrefix->GetFlag())
      {
        io::File::CloseClearFStream( m_RmsdFile);
        io::File::CloseClearFStream( m_DihedralBinFile);
        io::File::CloseClearFStream( m_ScoreFile);
        io::File::CloseClearFStream( m_RmsdScore);
      }

      if( m_ConformersSingleFile->GetFlag())
      {
        io::File::CloseClearFStream( m_OutputConformers);
      }

      if( m_OutputFailed3D->GetFlag())
      {
        io::File::CloseClearFStream( m_Failed3DStream);
      }
      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const ApplicationType ConformerGenerator::ConformerGenerator_Instance
    (
      GetAppGroups().AddAppToGroup( new ConformerGenerator(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl

