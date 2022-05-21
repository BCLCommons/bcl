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
#include "bcl_app_align_to_scaffold.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

    // Static instance initialization
    const ApplicationType AlignToScaffold::AlignToScaffold_Instance
    (
      GetAppGroups().AddAppToGroup( new AlignToScaffold(), GetAppGroups().e_Molecule)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    AlignToScaffold::AlignToScaffold() :
          m_InputScaffold
          (
            new command::Parameter
            (
              "scaffold_filename",
              "filename for input sdf of scaffold",
              command::ParameterCheckFileExistence()
            )
          ),
          m_EnsembleFileName
          (
            new command::Parameter
            (
              "ensemble_filename",
              "filename for input sdf of molecules, which will be aligned to scaffold",
              command::ParameterCheckFileExistence()
            )
          ),
          m_OutputFilename
          (
            new command::Parameter( "output", "filename for output sdf of aligned ensemble")
          ),
          m_AlignScaffoldAtoms
          (
            new command::FlagDynamic
            (
              "align_scaffold_atoms",
              "if you know which atoms you need to align, give a list of 3 atoms from the scaffold you want to align"
              " (if this is specified you must also specify which atoms to align in the ensemble)",
              command::Parameter
              (
                "align_atoms",
                "indices (0-indexed) of atoms you want to align"
              ),
              0,
              3
            )
          ),
          m_AlignEnsembleAtoms
          (
            new command::FlagDynamic
            (
              "align_ensemble_atoms",
              "if you know which atoms you need to align, give a list of 3 atoms from the ensemble that you want to align"
              " (if this is specified you must also specify which atoms to align in the scaffold)",
              command::Parameter
              (
                "align_atoms",
                "indices (0-indexed) of atoms you want to align"
              ),
              0,
              3
            )
          ),
          m_AlignRigid
          (
            new command::FlagStatic
            (
              "align_rigid",
              "align the largest rigid component of the scaffold"
            )
          ),
          m_AtomComparisonType
          (
            new command::FlagStatic
            (
              "atom_type",
              "atom type comparison",
              command::Parameter
              (
                "atom_comparison_type",
                "string value of the atom type info for substructure matching calculation",
                command::ParameterCheckSerializable
                (
                  chemistry::ConformationGraphConverter::AtomComparisonTypeEnum()
                ),
                "ElementType"
              )
            )
          ),
          m_BondComparisonType
          (
            new command::FlagStatic
            (
              "bond_type",
              "bond type comparison",
              command::Parameter
              (
                "bond_comparison_type",
                "string value of the bond type info for substructure matching calculation",
                command::ParameterCheckSerializable
                (
                  chemistry::ConfigurationalBondTypeData::DataEnum()
                ),
                "BondOrderOrAromaticWithRingness"
              )
            )
          ),
          m_SolutionTypeFlag
          (
            new command::FlagStatic
            (
              "solution_type",
              "graph search solution type",
              command::Parameter
              (
                "",
                "",
                command::ParameterCheckAllowed(
                  storage::Vector< std::string>::Create
                  (
                    "Connected",
                    "Unconnected",
                    "GreedyUnconnected"
                  )
                ),
                "Connected"
              )
            )
          )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    AlignToScaffold *AlignToScaffold::Clone() const
    {
      return new AlignToScaffold( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AlignToScaffold::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string AlignToScaffold::GetDescription() const
    {
      return "Align molecules by substructure";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &AlignToScaffold::GetReadMe() const
    {
      static std::string s_read_me =
          "AlignToScaffold is an application for aligning small molecules based on their common substructure atoms.";
      return s_read_me;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> AlignToScaffold::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // scaffold sdf that ensemble will be aligned to
      sp_cmd->AddParameter( m_InputScaffold);

      // hydrogen preferences
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // ensemble containing the molecules to be aligned to scaffold
      sp_cmd->AddParameter( m_EnsembleFileName);
      // ensemble that will be written out
      sp_cmd->AddParameter( m_OutputFilename);
      // atoms in the scaffold to align
      sp_cmd->AddFlag( m_AlignScaffoldAtoms);
      // atoms in the ensemble to align
      sp_cmd->AddFlag( m_AlignEnsembleAtoms);
      // atoms in the ensemble to align
      sp_cmd->AddFlag( m_AlignRigid);
      //! bond comparison
      sp_cmd->AddFlag( m_BondComparisonType);
      //! atom comparison
      sp_cmd->AddFlag( m_AtomComparisonType);
      //! isomorphism solution type
      sp_cmd->AddFlag( m_SolutionTypeFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int AlignToScaffold::Main() const
    {
      //read in scaffold sdf and ensemble sdf
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_InputScaffold->GetValue());
      chemistry::FragmentEnsemble scaffold( input, sdf::GetCommandLineHydrogensPref());
      io::File::MustOpenIFStream( input, m_EnsembleFileName->GetValue());
      chemistry::FragmentEnsemble ensemble( input, sdf::GetCommandLineHydrogensPref());
      io::File::CloseClearFStream( input);

      // Make sure align_scaffold_atoms and align_ensemble_atoms are both given or not given
      BCL_Assert
      (
        m_AlignScaffoldAtoms->GetFlag() == m_AlignEnsembleAtoms->GetFlag(),
        m_AlignScaffoldAtoms->GetName() + " and " + m_AlignEnsembleAtoms->GetName() + " must be given together or"
        "not at all"
      );

      // Make sure these flags have multiples of 3 given
      if( m_AlignScaffoldAtoms->GetFlag())
      {
        BCL_Assert
        (
          ( m_AlignScaffoldAtoms->GetParameterList().GetSize() == 3)
          && ( m_AlignEnsembleAtoms->GetParameterList().GetSize() == 3),
          m_AlignScaffoldAtoms->GetName() + " and " + m_AlignEnsembleAtoms->GetName() + " must have 3 specified atoms each"
        );
      }

      // Convert scaffold to graph in order to find isomorphism. Keep element type rather than atom type
      // because parts of scaffold, though aligned in the small molecules by element type, may be of a
      // different atom type (non-cyclic scaffold may extend into a ring in one or more of the molecules)
      graph::ConstGraph< size_t, size_t> scaffold_graph
      (
        chemistry::ConformationGraphConverter
        (
          chemistry::ConformationGraphConverter::e_ElementType,
          chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
        )( scaffold.GetMolecules().FirstElement())
      );

      //Initialize output file
      chemistry::FragmentEnsemble output_ensemble;

      // Reorder vectors, used to align specific atoms
      storage::Vector< size_t> reorder_scaff_vec;
      storage::Vector< size_t> reorder_ensemble_vec;

      // If align_*_atoms was specified populate the reorder vectors
      if( m_AlignScaffoldAtoms->GetFlag())
      {
        size_t indx_0( util::ConvertStringToNumericalValue< size_t>( m_AlignScaffoldAtoms->GetParameterList()( 0)->GetValue()));
        size_t indx_1( util::ConvertStringToNumericalValue< size_t>( m_AlignScaffoldAtoms->GetParameterList()( 1)->GetValue()));
        size_t indx_2( util::ConvertStringToNumericalValue< size_t>( m_AlignScaffoldAtoms->GetParameterList()( 2)->GetValue()));

        reorder_scaff_vec = storage::Vector< size_t>::Create( indx_0, indx_1, indx_2);

        indx_0 = util::ConvertStringToNumericalValue< size_t>( m_AlignEnsembleAtoms->GetParameterList()( 0)->GetValue());
        indx_1 = util::ConvertStringToNumericalValue< size_t>( m_AlignEnsembleAtoms->GetParameterList()( 1)->GetValue());
        indx_2 = util::ConvertStringToNumericalValue< size_t>( m_AlignEnsembleAtoms->GetParameterList()( 2)->GetValue());

        reorder_ensemble_vec = storage::Vector< size_t>::Create( indx_0, indx_1, indx_2);
      }
      else if( m_AlignRigid->GetFlag())
      {
        chemistry::FragmentSplitRigid rigid( 3);
        auto scaff_graph
        (
          chemistry::ConformationGraphConverter::CreateGraphWithAtoms( scaffold.GetMolecules().FirstElement())
        );
        auto rigid_components( rigid.GetComponentVertices( scaffold.GetMolecules().FirstElement(), scaff_graph));
        storage::List< storage::Vector< size_t> >::const_iterator itr_largest( rigid_components.Begin());
        for( auto itr( rigid_components.Begin()), itr_end( rigid_components.End()); itr != itr_end; ++itr)
        {
          if( itr->GetSize() > itr_largest->GetSize())
          {
            itr_largest = itr;
          }
        }
        if( itr_largest != rigid_components.End())
        {
          reorder_scaff_vec = *itr_largest;
          reorder_ensemble_vec = *itr_largest;
        }
      }

      // Iterate through all molecules in the input ensemble
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
        itr_compared( ensemble.Begin()), itr_compared_end( ensemble.End());
        itr_compared != itr_compared_end;
        ++itr_compared
      )
      {
        //Store the coordinates of the scaffold atoms for alignment
        util::SiPtrVector< const linal::Vector3D> scaffold_coords
        (
          scaffold.GetMolecules().FirstElement().GetAtomCoordinates()
        );

        //Store the atom coordinates of the current model
        util::SiPtrVector< const linal::Vector3D> molecule_coords( itr_compared->GetAtomCoordinates());

        if( !reorder_scaff_vec.IsEmpty())
        {
          scaffold_coords.Reorder( reorder_scaff_vec);
          molecule_coords.Reorder( reorder_ensemble_vec);
        }
        else
        {
          //Find the isomorphism between the scaffold and the current small molecule
          storage::Map< size_t, size_t> isomorphism( FindIsomorphism( *itr_compared, scaffold_graph));

          BCL_MessageDbg( "isomorphism between scaffold and molecule" + util::Format()( isomorphism));

          scaffold_coords.Reorder( isomorphism.GetKeysAsVector());

          molecule_coords.Reorder( isomorphism.GetMappedValues());
        }

        //Generate transformation matrix based on common isomorphism between scaffold and current small molecule
        math::TransformationMatrix3D transform
        (
          quality::RMSD::SuperimposeCoordinates( scaffold_coords, molecule_coords)
        );
        //Store the atom information of the small molecule
        storage::Vector< sdf::AtomInfo> atom_vector( itr_compared->GetAtomInfo());

        //Transform the coordinates of the small molecule atoms based on the transformation matrix
        for
        (
          storage::Vector< sdf::AtomInfo>::iterator itr( atom_vector.Begin()), itr_end( atom_vector.End());
          itr != itr_end;
          ++itr
        )
        {
          linal::Vector3D temp( itr->GetCoordinates());
          itr->SetCoordinates( temp.Transform( transform));
        }

        //Store transformed coordinates and atom information in a new molecule and add it to output ensemble
        chemistry::FragmentComplete new_molecule
        (
          chemistry::AtomVector< chemistry::AtomComplete>( atom_vector, itr_compared->GetBondInfo()),
          itr_compared->GetName()
        );

        // Copy the properties from the original molecule to the new molecule
        new_molecule.StoreProperties( *itr_compared);
        output_ensemble.PushBack( new_molecule);
      }

      // write out the transformed molecule ensemble
      io::OFStream output;
      io::File::MustOpenOFStream( output, m_OutputFilename->GetValue());
      output_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      // end
      return 0;
    }

    // Find the isomorphism (overlap) between the scaffold and small molecule. Will be called for each molecule in the
    // input ensemble.
    storage::Map< size_t, size_t> AlignToScaffold::FindIsomorphism
    (
      const chemistry::FragmentComplete &MOLECULE,
      graph::ConstGraph< size_t, size_t> &SCAFFOLD_GRAPH
    ) const
    {
      // empty initialization
      graph::CommonSubgraphIsomorphism< size_t, size_t> csi_substructure;

      // check solution type
      // TODO: can this be directly initialized from string instead?
      if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "Unconnected")
      {
        csi_substructure = graph::CommonSubgraphIsomorphism< size_t, size_t>( graph::CommonSubgraphIsomorphismBase::SolutionType::e_Unconnected);
      }
      else if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "GreedyUnconnected")
      {
        csi_substructure = graph::CommonSubgraphIsomorphism< size_t, size_t>( graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected);
      }
      else
      {
        csi_substructure = graph::CommonSubgraphIsomorphism< size_t, size_t>( graph::CommonSubgraphIsomorphismBase::SolutionType::e_Connected);
      }

      graph::SubgraphIsomorphism< size_t, size_t> si_substructure;

      // Convert current small molecule to a graph keeping element type due to common overlap of elements that may
      // exist in different atom types between scaffold and small molecule
      graph::ConstGraph< size_t, size_t> molecule_graph
      (
        chemistry::ConformationGraphConverter
        (
          chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomComparisonType->GetFirstParameter()->GetValue()),
          chemistry::ConfigurationalBondTypeData::DataEnum( m_BondComparisonType->GetFirstParameter()->GetValue())
        )
        ( MOLECULE)
      );

      si_substructure.SetGraphExternalOwnership( molecule_graph);
      si_substructure.SetSubgraphExternalOwnership( SCAFFOLD_GRAPH);
      if( si_substructure.FindIsomorphism())
      {
        const storage::Vector< size_t> iso( si_substructure.GetIsomorphism());
        storage::Map< size_t, size_t> iso_map;
        for( size_t i( 0), sz( iso.GetSize()); i < sz; ++i)
        {
          iso_map[ i] = iso( i);
        }
        return iso_map;
      }

      csi_substructure.SetGraphA
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &SCAFFOLD_GRAPH, false)
      );
      csi_substructure.SetGraphB
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &molecule_graph, false)
      );

      //Find isomorphism between the scaffold graph and the small molecule graph
      csi_substructure.FindIsomorphism( SCAFFOLD_GRAPH.GetSize());
      storage::Map< size_t, size_t> isomorphism( csi_substructure.GetIsomorphism());

      // Send message if scaffold is not entirely contained within the small molecule. This may or may not be
      // important depending on the context of this application's use.
      if( isomorphism.GetSize() != SCAFFOLD_GRAPH.GetSize())
      {
        BCL_MessageCrt( "the isomorphism is not equal to the size of scaffold");
      }
      return isomorphism;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AlignToScaffold::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AlignToScaffold::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
