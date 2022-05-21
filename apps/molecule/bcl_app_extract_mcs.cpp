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
#include "bcl_app_extract_mcs.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_csi_substructure.h"
#include "io/bcl_io_file.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

    const ApplicationType ExtractMCS::ExtractMCS_Instance
    (
      GetAppGroups().AddAppToGroup( new ExtractMCS(), GetAppGroups().e_Molecule)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! default constructor
    //! @brief standard constructor
    ExtractMCS::ExtractMCS() :
        m_InputFileFlag
        (
          new command::FlagStatic
          (
            "input_filenames",
            "sdf file containing the ensembles from which to find scaffolds",
            command::Parameter
            (
              "input",
              "an sdf input file",
              command::ParameterCheckFileExistence()
            )
          )
        ),
        m_MinSizeFlag
        (
          new command::FlagStatic
          (
            "min_size",
            "minimum size of the fragments to keep it in the fragments list",
            command::Parameter
            (
              "min_size",
              "minimum size of the fragments to keep it in the fragments list",
              command::ParameterCheckRanged< size_t>(),
              "2"
            )
          )
        ),
        m_BondTypeData
        (
          new command::FlagStatic
          (
            "bond_coloring_scheme",
            "how to color bond data",
            command::Parameter
            (
              "scheme",
              "used to compare bonds",
              command::ParameterCheckSerializable( chemistry::ConfigurationalBondTypeData::DataEnum()),
              "BondOrderInRingOrAromatic"
            )
          )
        ),
        m_OutputFileFlag
        (
          new command::FlagStatic
          (
            "output",
            "filename to output sdf containing scaffolds",
            command::Parameter
            (
              "output_filename",
              "filename to output sdf containing scaffolds"
            )
          )
        )
    {
    }

    //! copy constructor, only copy the flags
    ExtractMCS::ExtractMCS( const ExtractMCS &PARENT) :
          m_InputFileFlag( PARENT.m_InputFileFlag),
          m_MinSizeFlag( PARENT.m_MinSizeFlag),
          m_BondTypeData( PARENT.m_BondTypeData),
          m_OutputFileFlag( PARENT.m_OutputFileFlag)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ExtractMCS
    ExtractMCS *ExtractMCS::Clone() const
    {
      return new ExtractMCS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ExtractMCS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ExtractMCS::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // insert all the flags and params
      sp_cmd->AddFlag( m_InputFileFlag);
      sp_cmd->AddFlag( m_MinSizeFlag);
      sp_cmd->AddFlag( m_BondTypeData);
      sp_cmd->AddFlag( m_OutputFileFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ExtractMCS::Main() const
    {
      // read in ensemble
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_InputFileFlag->GetFirstParameter()->GetValue());
      // remove any hydrogens because they slow down the scaffold search and are unnecessary
      m_Ensemble = chemistry::FragmentEnsemble( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);

      // store the size of the ensemble
      const size_t ensemble_size( m_Ensemble.GetSize());

      if( ensemble_size <= 1)
      {
        BCL_MessageCrt( "Cannot find common scaffolds without at least two molecules");
        return 1;
      }

      // Get the coloring scheme to be used for bonds
      m_BondColoringScheme =
        chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeData->GetFirstParameter()->GetValue());

      // For each small molecule, instantiate a graphical representation based on atom type and bond order/aromatic
      // store the graphs in vectors
      m_EnsembleSimpleGraphs.Reset();
      m_EnsembleSimpleGraphs.AllocateMemory( ensemble_size);
      m_EnsembleGraphs.Reset();
      m_EnsembleGraphs.AllocateMemory( ensemble_size);
      chemistry::ConformationGraphConverter graph_converter
      (
        chemistry::ConformationGraphConverter::e_AtomType,
        m_BondColoringScheme
      );
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( m_Ensemble.Begin()), itr_end( m_Ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        m_EnsembleSimpleGraphs.PushBack( graph_converter( *itr));
        m_EnsembleGraphs.PushBack( chemistry::ConformationGraphConverter::CreateGraphWithAtoms( *itr));
      }

      // open the output file so that we can write out scaffolds as we go
      io::OFStream output;
      io::File::MustOpenOFStream( output, m_OutputFileFlag->GetFirstParameter()->GetValue());

      chemistry::FragmentComplete mcs;
      FindMCS( mcs);

      BCL_MessageStd
      (
        "Saved the scaffold to " + m_OutputFileFlag->GetFirstParameter()->GetValue()
      );
      mcs.WriteMDL( output);

      // close the file stream
      io::File::CloseClearFStream( output);

      return 0;
    }

    //! @brief test whether an isomorphism is large enough
    //! @param SUBGRAPH_ISOMORPHISM_SIZE the size of the isomorphism
    //! @return true if the subgraph isomorphism size is larger than or equal to the m_MinSizeFlag
    bool ExtractMCS::IsIsomorphismLargerThanMinSize( const size_t &SUBGRAPH_ISOMORPHISM_SIZE) const
    {
      static const size_t s_MinIsomorphismSize( m_MinSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>());
      return SUBGRAPH_ISOMORPHISM_SIZE >= s_MinIsomorphismSize;
    }

    //! @brief compare the molecules given by the indices in a vector
    void ExtractMCS::FindMCS( chemistry::FragmentComplete &MCS) const
    {
      // this object will hold the isomorphism between any graphs that are compared (vertices colored by atom type index)
      graph::CommonSubgraphIsomorphism< size_t, size_t>
        isomorphism( graph::CommonSubgraphIsomorphismBase::e_Connected);

      // store the minimum size for an isomorphism to be interesting
      const size_t min_size( m_MinSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>());

      // store the size of the ensemble
      const size_t ensemble_size( m_EnsembleGraphs.GetSize());
      if( !ensemble_size)
      {
        return;
      }

      chemistry::FragmentComplete mcs( *m_Ensemble.Begin());
      chemistry::ConformationGraphConverter graph_converter
      (
        chemistry::ConformationGraphConverter::e_AtomType,
        m_BondColoringScheme
      );

      size_t mol_index( 1);
      for
      (
        chemistry::FragmentEnsemble::const_iterator itr_mol( ++m_Ensemble.Begin()), itr_mol_end( m_Ensemble.End());
        itr_mol != itr_mol_end;
        ++itr_mol, ++mol_index
      )
      {
        graph::ConstGraph< size_t, size_t> mcs_graph( graph_converter( mcs));
        graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> mcs_atom_graph( chemistry::ConformationGraphConverter::CreateGraphWithAtoms( mcs));

        // set up the isomorphisms with pointers to the graphs
        // copying the graphs is fairly expensive, so we avoid it whenever possible
        isomorphism.SetGraphA
        (
          util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &mcs_graph, false)
        );
        isomorphism.SetGraphB
        (
          util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &m_EnsembleSimpleGraphs( mol_index), false)
        );

        // get the actual isomorphism, using estimated upper bounds on its size
        isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds(), min_size);

        // if the isomorphism is too small, ignore it
        if( !IsIsomorphismLargerThanMinSize( isomorphism.GetIsomorphism().GetSize()))
        {
          continue;
        }

        // get the actual subgraph from graph a (this choice is arbitrary)
        graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> atom_subgraph
        (
          mcs_atom_graph.GetSubgraph
          (
            isomorphism.GetSubgraphIsomorphismsOfGraphA().FirstElement().GetVertexIndices()
          )
        );

        // convert the scaffold graph back into a small molecule and add properties of interest
        std::string new_name;
        chemistry::FragmentComplete new_fc
        (
          chemistry::ConformationGraphConverter::CreateAtomsFromGraph( atom_subgraph),
          new_name
        );
        mcs = new_fc;
      }
      MCS = mcs;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ExtractMCS::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ExtractMCS::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
