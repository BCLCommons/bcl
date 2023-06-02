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

// App header
#include "app/bcl_app_apps.h"

// include header for this application
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_molecule_feature_mapper.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_default.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "molecule/bcl_app_molecule_features.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"


namespace bcl
{
  namespace app
  {

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeFeatures::GetDescription() const
    {
      return "Generate absolute and relative pharmacophore/feature maps";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeFeatures::GetReadMe() const
    {
      static std::string s_read_me =
          "MoleculeFeatures is an application used to map the contribution of each atom to global molecular features. "
          "Specifically, MoleculeFeatures can be used to map the per-atom input sensitivities from QSAR "
          "neural networks to chemical structures. This helps to interpret which part of the molecule is driving "
          "the ANN signal. In addition, MoleculeFeatures can be used to create relative pharmacophore maps, which "
          "compare the per-atom contribution to some property (e.g. neural network QSAR score) of corresponding "
          "atoms in two or more molecules.";
      return s_read_me;
    }

    //! @brief initializes the command object for this application
    //! @return a ShPtr to a Command containing all of this applications parameters
    util::ShPtr< command::Command> MoleculeFeatures::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // add AtomMdlLine to molecule
      sdf::MdlHandler::AddAtomMdlLineFlag( *sp_cmd);

      // flags for input/output
      sp_cmd->AddFlag( m_MoleculesFlag);
      sp_cmd->AddFlag( m_ReferenceFlag);
      sp_cmd->AddFlag( m_MutuallyMatchingAtomsFlag);
      sp_cmd->AddFlag( m_ScorerFlag);
      sp_cmd->AddFlag( m_ElementPerturberFlag);
      sp_cmd->AddFlag( m_OutputFilenameFlag);
      sp_cmd->AddFlag( m_IgnoreHFlag);
      sp_cmd->AddFlag( m_NormalizeFlag);
      sp_cmd->AddFlag( m_StatisticFlag);
      sp_cmd->AddFlag( m_SplitLargestFlag);
      sp_cmd->AddFlag( m_PerturbTypeFlag);
      sp_cmd->AddFlag( m_AtomTypeFlag);
      sp_cmd->AddFlag( m_BondTypeFlag);
      sp_cmd->AddFlag( m_ScoreDistributionTypeFlag);
      sp_cmd->AddFlag( m_ColorFlag);
      sp_cmd->AddFlag( m_ColorSpectrumFlag);
      sp_cmd->AddFlag( m_OutputIntermediatesFlag);

      // default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    //! @brief the Main function
    //! @return 0 for success
    int MoleculeFeatures::Main() const
    {

      chemistry::FragmentEnsemble molecules;
      chemistry::FragmentFeed feed( m_MoleculesFlag->GetStringList(), sdf::GetCommandLineHydrogensPref());
      for( ; feed.NotAtEnd(); ++feed)
      {
        if( feed->GetNumberAtoms() > 0)
        {
          molecules.PushBack( *feed);
        }
      }

      // read in reference molecules
      chemistry::FragmentEnsemble ref_molecules;
      if( m_ReferenceFlag->GetFlag())
      {
        chemistry::FragmentFeed ref_feed( m_ReferenceFlag->GetStringList(), sdf::GetCommandLineHydrogensPref());
        for( ; ref_feed.NotAtEnd(); ++ref_feed)
        {
          if( ref_feed->GetNumberAtoms() > 0)
          {
            ref_molecules.PushBack( *ref_feed);
          }
        }
      }

      util::ObjectDataLabel score_label( m_ScorerFlag->GetFirstParameter()->GetValue());

      chemistry::MoleculeFeatureMapper mapper( score_label);

      std::string out_filename( m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      io::OFStream out;
      io::File::MustOpenOFStream( out, out_filename);
      out << "molecule,atom,element,effect" << std::endl;

      // open the information file
//      io::OFStream out1;
//      io::File::MustOpenOFStream( out1, out_filename + ".frags.sdf");

      // Write preliminary parts of the PML
      std::string pml_filename( out_filename + ".pml");
      io::OFStream pml;
      io::File::MustOpenOFStream( pml, pml_filename);
      pml << "load " << m_MoleculesFlag->GetFirstParameter()->GetValue() << ", molecule" << std::endl;
      pml << "load " << m_MoleculesFlag->GetFirstParameter()->GetValue() << ", molecule2" << std::endl;
      pml << "alter molecule and state 1, b=0" << std::endl;
      pml << "set label_size, -0.4" << std::endl;
      pml << "show sticks, molecule2 " << std::endl;
      pml << "hide sticks, molecule" << std::endl;
      pml << "show surface, molecule" << std::endl;
      pml << "set valence, 1" << std::endl;
      pml << "set ray_trace_mode, 1" << std::endl;
      pml << "set ray_opaque_background, 1" << std::endl;
      pml << "set antialias, 2" << std::endl;
      pml << "set ray_shadows, off" << std::endl;
      pml << "set ray_shadow, off" << std::endl;

      size_t mol_no( 0);
      for
      (
        chemistry::FragmentEnsemble::const_iterator itr_mol( molecules.Begin()), itr_mol_end( molecules.End());
        itr_mol != itr_mol_end;
        ++itr_mol, ++mol_no
      )
      {
        BCL_MessageStd( "Molecule " + util::Format()( mol_no));
        std::vector< std::map< size_t, float> > result;
        if( m_PerturbTypeFlag->GetFirstParameter()->GetValue() == "Atoms")
        {
          // check if there are elements that are passed for perturbation
          storage::Set< chemistry::ElementType> elements_to_perturb;
          if( m_ElementPerturberFlag->GetFlag())
          {
            storage::Vector< std::string> ele_strings( m_ElementPerturberFlag->GetStringList());
            for( auto ele_itr( ele_strings.Begin()), ele_itr_end( ele_strings.End()); ele_itr != ele_itr_end; ++ele_itr)
            {
              elements_to_perturb.InsertElement( chemistry::GetElementTypes().ElementTypeLookup( *ele_itr));
            }
          }

          result = mapper.Perturb
              (
                *itr_mol,
                m_IgnoreHFlag->GetFlag(),
                m_SplitLargestFlag->GetFlag(),
                storage::Vector< size_t>(),
                storage::Vector< chemistry::ElementType>( elements_to_perturb.Begin(), elements_to_perturb.End())
              );
        }
        else if( m_PerturbTypeFlag->GetFirstParameter()->GetValue() == "Molecules" && m_ScoreDistributionTypeFlag->GetFirstParameter()->GetValue() == "Naive")
        {
          result = mapper.CompareSubstructuresNaive
              (
                *itr_mol,
                m_ReferenceFlag->GetFlag() ? ref_molecules : molecules,
                chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomTypeFlag->GetFirstParameter()->GetValue()),
                chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeFlag->GetFirstParameter()->GetValue()),
                m_NormalizeFlag->GetFlag(),
                m_IgnoreHFlag->GetFlag(),
                m_StatisticFlag->GetFlag(),
                m_ReferenceFlag->GetFlag() ? util::GetUndefinedSize_t() : mol_no,
                m_OutputIntermediatesFlag->GetFlag()
              );
        }
        else if( m_PerturbTypeFlag->GetFirstParameter()->GetValue() == "Molecules" && m_ScoreDistributionTypeFlag->GetFirstParameter()->GetValue() == "Rigorous")
        {
          result = mapper.CompareSubstructuresRigorous
              (
                *itr_mol,
                m_ReferenceFlag->GetFlag() ? ref_molecules : molecules,
                chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomTypeFlag->GetFirstParameter()->GetValue()),
                chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeFlag->GetFirstParameter()->GetValue()),
                m_MutuallyMatchingAtomsFlag->GetFlag(),
                m_NormalizeFlag->GetFlag(),
                m_StatisticFlag->GetFlag(),
                m_ReferenceFlag->GetFlag() ? util::GetUndefinedSize_t() : mol_no,
                m_OutputIntermediatesFlag->GetFlag()
              );
        }
        else
        {
          util::Implementation< chemistry::FragmentSplitInterface> splitter( m_PerturbTypeFlag->GetFirstParameter()->GetValue());
          result = mapper.PerturbByFragment
              (
                *itr_mol,
                splitter,
                chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomTypeFlag->GetFirstParameter()->GetValue()),
                chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeFlag->GetFirstParameter()->GetValue()),
                m_MutuallyMatchingAtomsFlag->GetFlag(),
                m_NormalizeFlag->GetFlag(),
                m_IgnoreHFlag->GetFlag(),
                m_StatisticFlag->GetFlag(),
                m_ReferenceFlag->GetFlag() ? util::GetUndefinedSize_t() : mol_no,
                m_OutputIntermediatesFlag->GetFlag()
              );
        }
        if( result.empty())
        {
          continue;
        }

        storage::Vector< sdf::AtomInfo> atom_infos( itr_mol->GetAtomInfo());

        // Find the appropriate range for coloring in the PML
        float max_rm( std::max_element( result[ 0].begin(), result[ 0].end())->second);
        float min_rm( std::min_element( result[ 0].begin(), result[ 0].end())->second);
        float range( std::max( std::fabs( max_rm), std::fabs( min_rm)));

        for
        (
          std::map< size_t, float>::const_iterator itr_atom( result[ 0].begin()), itr_atom_end( result[ 0].end());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          if( !util::IsDefined( itr_atom->second))
          {
            continue;
          }
          const sdf::AtomInfo &atom( atom_infos( itr_atom->first));

          // molecule,atom,element,effect
          out << mol_no << ",";
          out << itr_atom->first << ",";
          out << atom.GetAtomType()->GetElementType()->GetChemicalSymbol() << ",";
          out << itr_atom->second << std::endl;

          pml << "alter \"molecule\" and state " << mol_no + 1 << " and id " << itr_atom->first + 1
              << ", b=" << itr_atom->second << std::endl;
        }

        // Color the atoms in the PML
        pml << "spectrum b, " << m_ColorSpectrumFlag->GetFirstParameter()->GetValue() << ", selection=molecule and state " << mol_no + 1
            << ", minimum=" << float( -1 * range) << ", maximum=" << range << std::endl;
        pml << "unset surface_color" << std::endl;
        pml << "set ray_shadow, off" << std::endl;
        pml << "hide everything, hydrogen" << std::endl;
        pml << "util.cbaw molecule2" << std::endl;
        pml << "orient" << std::endl;
        pml << "bg_color white" << std::endl;
        pml << "set surface_mode, 3" << std::endl;
        pml << "set transparency_mode, 1" << std::endl;
        pml << "set transparency, 0.5" << std::endl;

      }
      // set global color options
      if( m_ColorFlag->GetFlag())
      {
        // Color the atoms in the PML
        pml << "spectrum b, " << m_ColorSpectrumFlag->GetFirstParameter()->GetValue() << ", selection=molecule"
            << ", minimum=" << float( m_ColorFlag->GetParameterList().FirstElement()->GetNumericalValue< float>())
            << ", maximum=" << float( m_ColorFlag->GetParameterList().LastElement()->GetNumericalValue< float>()) << std::endl;
      }
//      io::File::CloseClearFStream( out1);
      io::File::CloseClearFStream( out);
      io::File::CloseClearFStream( pml);

      return 0;
    }
     // Main

    //! @brief standard constructor
    MoleculeFeatures::MoleculeFeatures() :
      m_MoleculesFlag
      (
        new command::FlagStatic
        (
          "input_filenames",
          "filename containing the molecules to analyze",
          command::Parameter
          (
            "filename",
            "the name of an SDF file containing the molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ReferenceFlag
      (
        new command::FlagDynamic
        (
          "reference_mols",
          "filename containing the reference molecule(s) for substructure comparisons the perturbation type 'Molecules'; "
          "by default all molecules from 'input_filenames' are compared against one another",
          command::Parameter
          (
            "reference_mols",
            "the name of an SDF file containing the reference molecule(s)"
          )
        )
      ),
      m_MutuallyMatchingAtomsFlag
      (
        new command::FlagDynamic
        (
          "mutually_matching_atoms",
          "whether to perturb substructures based on spatially matched atoms instead of topological substructure",
          command::Parameter
          (
            "mutually_matching_atoms",
            "if set, will use spatially aligned mutually matched atom pairs for substructure perturbation in the 'Molecules' perturbation routine; "
            "recommended if the molecules being compared are not from a congeneric series of small molecules; "
            "this is faster than a substructure search, but requires that the molecules are already superimposed via small molecule alignment "
            "(e.g. with BCL MolAlign) or from docking (e.g. RosettaLigand or BCL cheminfo:MoleculeFit)"
          )
        )
      ),
      m_ScorerFlag
      (
        new command::FlagStatic
        (
          "model",
          "models used to compute the activity or score of a molecule",
          command::Parameter
          (
            "descriptor",
            "the descriptor to use to calculate activities (probably should be PredictionMean())",
            command::ParameterCheckSerializable
            (
              descriptor::CheminfoProperty()
            )
          )
        )
      ),
      m_ElementPerturberFlag
      (
        new command::FlagDynamic
        (
          "perturb_element_types",
          "sequentially mutate each atom only if it is of one of the specified element types.",
          command::Parameter
          (
            "element",
            "if provided, will use element list to focus per-atom feature contribution calculations",
            command::ParameterCheckDefault(),
            ""
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "filename where output information will be written",
          command::Parameter
          (
            "filename",
            "file to write output to"
          )
        )
      ),
      m_IgnoreHFlag
      (
        new command::FlagDynamic
        (
          "ignore_h",
          "whether to ignore hydrogens",
          command::Parameter
          (
            "ignore_h",
            "if set, will ignore hydrogens during substructure comparisons and feature mapping"
          )
        )
      ),
      m_NormalizeFlag
      (
        new command::FlagDynamic
        (
          "normalize",
          "if performing fragment-based feature mapping (one of the non-'Atoms' perturbation types), "
          "normalize per atom scores by the number of atoms in the fragment; "
          "if used simultaneously with 'ignore_h', will only normalize by heavy atom counts",
          command::Parameter
          (
            "normalize",
            "if set, will normalize the per atom score from each fragment by the total number of atoms in that fragment"
          )
        )
      ),
      m_StatisticFlag
      (
        new command::FlagDynamic
        (
          "average",
          "if performing molecule comparisons type feature mapping, average per atom results across all comparisons; otherwise "
          "report cumulative sums",
          command::Parameter
          (
            "average",
            "if set, will average the per atom score from each comparison by the total number of molecule comparisons"
          )
        )
      ),
      m_SplitLargestFlag
      (
        new command::FlagDynamic
        (
          "split_largest",
          "if performing atoms-based feature mapping, remove the small fragments created when atoms are removed",
          command::Parameter
          (
            "split_largest",
            "if set, will remove small fragments formed when atoms are removed prior to scoring"
          )
        )
      ),
      m_PerturbTypeFlag
      (
        new command::FlagStatic
        (
          "perturb_type",
          "type of perturbation to perform on the molecule to build the feature map; "
          "'Atoms' will sequentially remove individual atoms, 'Molecules' will sequentially "
          "compare the other input molecules, and all of the other types are based on the "
          "molecule:Split interface implementations",
          command::Parameter
          (
            "perturb_type",
            "perturbation type used to build the feature map",
            "Rigid"
          )
        )
      ),
      m_AtomTypeFlag
      (
        new command::FlagStatic
        (
          "atom_type_comparison",
          "level of detail with which to compare vertices (atoms) in substructure graph searches during perturbation type "
          "'Molecules' at the 'Rigorous' level.",
          command::Parameter
          (
            "atom_type_data",
            "Resolution of atom types during substructure searches",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create
              (
                "Identity", "ElementType", "AtomType", "AtomTypeAndChirality", "AtomTypeAndComplexRingChirality", "AtomTypeAndSymmetry", "AtomTypeAndHasSymmetry", "AtomTypeAndNumberHydrogens",
                "AtomTypeAndNumberHydrogensOnRingAtoms", "AtomTypeAndNumberHydrogensOnRingsAndDistinguishHydrogens", "AtomTypeAndDistinguishHydrogens", "CIPPriorityHighToLow"
              )
            ),
            "ElementType"
          )
        )
      ),
      m_BondTypeFlag
      (
        new command::FlagStatic
        (
          "bond_type_comparison",
          "level of detail with which to compare edges (bonds) in substructure graph searches during perturbation type "
          "'Molecules' at the 'Rigorous' level.",
          command::Parameter
          (
            "bond_type_data",
            "Resolution of bond types during substructure searches",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create
              (
                "Identity", "BondOrder", "NumberOfElectrons", "Conjugation", "IsConjugated", "IsAromatic", "IsAmide", "IsInRing", "BondOrderInRingOrAromatic", "BondOrderOrAromatic",
                "BondOrderAmideOrAromatic", "BondOrderOrAromaticWithRingness", "BondOrderAmideOrAromaticWithRingness", "FuzzyBondOrderWithIsometryOrAromaticWithRingness",
                "FuzzyBondOrderAmideOrAromaticWithRingness", "ConstitutionalBondType", "BondOrderWithIsometry", "Isometry", "IsIsometric", "BondOrderWithIsometryOrAromatic",
                "BondOrderAmideWithIsometryOrAromaticWithRingness", "ConfigurationalBondType"
              )
            ),
            "BondOrderAmideOrAromaticWithRingness"
          )
        )
      ),
      m_ScoreDistributionTypeFlag
      (
        new command::FlagStatic
        (
          "atom_score_distribution_type",
          "method with which to distribute the per-atom contributions when using the 'Molecules' perturbation type. \n"
          "'Naive' will evenly distribute the total score difference between the two molecules to the atoms of the "
          "differing substructures. This can be modified by passing the 'normalize', 'ignore_h', and/or 'average' flags, "
          "but it is possible that score differences will be attributed to atoms that are inconsequential. The benefit "
          "is that this approach is faster than the alternative. \n"
          "'Rigorous' will individually remove the sets of connected atoms (i.e. fragment) that differ between the parent molecule "
          "and the comparison molecule(s) and then rescore the perturbed structure. The score differences will then be attributed to "
          "the atoms of the specific fragment removed at that time. Note that this cannot be done with cached properties (e.g. experimental "
          "values saved as MDL properties on the SDF file). ",
          command::Parameter
          (
            "distribution_type",
            "score distribution type used to build the feature map",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "Naive", "Rigorous")),
            "Naive"
          )
        )
      ),
      m_ColorFlag
      (
        new command::FlagDynamic
        (
          "color_minmax",
          "Set global min and max values for all molecules during PyMol visualization; "
          "defaults to individual min and max values for each molecule",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "color_min",
              "minimum boundary on spectrum coloring of feature map molecule surface",
              command::ParameterCheckRanged< float>( -100.0, 100.0),
              "0.0"
            ),
            command::Parameter
            (
              "color_max",
              "minimum boundary on spectrum coloring of feature map molecule surface",
              command::ParameterCheckRanged< float>( -100.0, 100.0),
              "1.0"
            )
          )
        )
      ),
      m_ColorSpectrumFlag
      (
        new command::FlagStatic
        (
          "color_spectrum",
          "PyMol color spectrum to use for feature mapping",
          command::Parameter
          (
            "color_spectrum",
            "valid PyMol color spectrum",
            "blue_white_red"
          )
        )
      ),
      m_OutputIntermediatesFlag
      (
        new command::FlagDynamic
        (
          "output_intermediates",
          "if performing fragment-based feature mapping (one of the non-'Atoms' perturbation types), "
          "output intermediate fragments to SDF files",
          command::Parameter
          (
            "output_intermediates",
            "if set, will output intermediate fragments from perturbations"
          )
        )
      )
    {
    }

    // Construct the static instance of this application, and add it to the ChemInfo group
    const ApplicationType MoleculeFeatures::MoleculeFeatures_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeFeatures(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
