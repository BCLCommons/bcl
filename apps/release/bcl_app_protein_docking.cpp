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

// include header of this class
#include "bcl_app_protein_docking.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_exposure_prediction.h"
#include "biol/bcl_biol_membrane.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "mc/bcl_mc_optimization_docking.h"
#include "opti/bcl_opti_optimization_interface.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_printer_quality_docking.h"
#include "pdb/bcl_pdb_printer_score.h"
#include "quality/bcl_quality_rmsd.h"
#include "score/bcl_score_protein_model_aa_neighborhood_docking.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    //! @brief default constructor
    ProteinDocking::ProteinDocking() :
      m_FlagReceptor
      (
        new command::FlagStatic
        (
          "receptor",
          "\tPDB file of the receptor",
          command::Parameter( "receptor PDB file", "")
        )
      ),
      m_FlagLigand
      (
        new command::FlagStatic
        (
          "ligand",
          "\tPDB file of the ligand",
          command::Parameter( "ligand PDB file", "")
        )
      ),
      m_FlagNative
      (
        new command::FlagStatic
        (
          "native",
          "\tPDB file of the native structure of the complex",
          command::Parameter( "native complex PDB file", "")
        )
      ),
      m_FlagMembrane
      (
        new command::FlagStatic
        (
          "membrane",
          "\tMembrane for the complex",
          command::Parameter( "PDBTM XML file", "")
        )
      ),
      m_FlagOptimizer
      (
        new command::FlagStatic
        (
          "optimizer",
          "\tfile containing specifications of the optimizer to be used",
          command::Parameter( "optimizer", "")
        )
      ),
      m_FlagNumberModels
      (
        new command::FlagStatic
        (
          "number_models",
          "\tnumber of models to generate",
          command::Parameter( "number of models", "")
        )
      ),
      m_FlagPackingDensity
      (
        new command::FlagStatic
        (
          "packing_density",
          "\tfile that contains residue packing densities",
          command::Parameter( "residue packing density file", "", "")
        )
      ),
      m_FlagOutputPrefix
      (
        new command::FlagStatic
        (
          "output_prefix",
          "\tprefix to the output filename",
          command::Parameter( "output prefix", "")
        )
      )
//      m_FlagSymmetry
//      (
//        new command::FlagStatic
//        (
//          "symmetry",
//          "\tsymmetry of the system",
//          command::Parameter( "symmetry", "")
//        )
//      )
    {
      // nothing else to do
    }

    //! @brief returns a pointer to a copy of this object
    ProteinDocking *ProteinDocking::Clone() const
    {
      return new ProteinDocking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &ProteinDocking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief a detailed description of this application
    std::string ProteinDocking::GetDescription() const
    {
      std::string description
      (
        "application for protein protein docking, optionally with "
        "biochemical or biophysical restraints."
      );
      return description;
    }

    //! @brief initializes the command object for that executable
    //! @return initialized command object
    util::ShPtr< command::Command> ProteinDocking::InitializeCommand() const
    {
      // create a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add flags or parameters
      sp_cmd->AddFlag( m_FlagReceptor);
      sp_cmd->AddFlag( m_FlagLigand);
      sp_cmd->AddFlag( m_FlagNative);
      sp_cmd->AddFlag( m_FlagMembrane);
      sp_cmd->AddFlag( m_FlagOptimizer);
      sp_cmd->AddFlag( m_FlagNumberModels);
      sp_cmd->AddFlag( m_FlagPackingDensity),
      sp_cmd->AddFlag( m_FlagOutputPrefix);

      // add all default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return the command object
      return sp_cmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the main function of this application
    //! @return exit code - 0 for success
    int ProteinDocking::Main() const
    {
      // PDB factory for generating the receptor model and ligand model
      pdb::Factory pdb_factory;

      // create the receptor model
      const std::string receptor_pdb_file( m_FlagReceptor->GetFirstParameter()->GetValue());
      assemble::ProteinModel receptor( pdb_factory.ProteinModelFromPDBFilename( receptor_pdb_file));
      std::string receptor_chain_ids( receptor.GetChainIDs());

      // create the ligand model
      std::string ligand_pdb_file( m_FlagLigand->GetFirstParameter()->GetValue());
      assemble::ProteinModel ligand( pdb_factory.ProteinModelFromPDBFilename( ligand_pdb_file));
      std::string ligand_chain_ids( ligand.GetChainIDs());

      // both the receptor and the ligand are assumed to sit at the center of the membrane initially
      linal::Vector3D translation( 1000.0, 1000.0, 0.0);
      math::RotationMatrix3D rotation3d( coord::GetAxes().e_Z, random::GetGlobalRandom().Random( 2.0 * math::g_Pi));
      translation.Rotate( rotation3d);
      ligand.Translate( translation);
      PreDock( receptor, ligand);

      // combine the ligand and the receptor to create a start model
      util::ShPtrVector< assemble::Chain> chains;
      chains.Append( receptor.GetChains());
      chains.Append( ligand.GetChains());
      assemble::ProteinModel starting_model( chains);

      // preprocess the starting model, such as add native structure, membrane, etc
      PreProcess( starting_model);

      // write out starting model
      io::OFStream write_pdb;
      const std::string preprocessed_model_filename
      (
        m_FlagOutputPrefix->GetFirstParameter()->GetValue() + "_starting_model_preprocessed.pdb"
      );
      io::File::MustOpenOFStream( write_pdb, preprocessed_model_filename);
      BCL_MessageStd( "Writing the preprocessed starting conformation to " + preprocessed_model_filename);
      pdb_factory.WriteModelToPDB( starting_model, write_pdb);
      io::File::CloseClearFStream( write_pdb);

      // create the optimizer
      mc::OptimizationDocking optimizer;
      optimizer.AssertRead( util::ObjectDataLabel( m_FlagOptimizer->GetFirstParameter()->GetValue()));

      // modifies the factory so that the score table and quality measures are written to the PDB file
      pdb_factory.AppendPrinter
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
        (
          new pdb::PrinterScore
          (
            util::ShPtr< score::ProteinModelScoreSum>( new score::ProteinModelScoreSum( optimizer.GetScoreFunction()))
          )
        )
      );

      // get ligand chain id
      pdb_factory.AppendPrinter
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
        (
          new pdb::PrinterQualityDocking
          (
            optimizer.GetQualityMeasures(),
            ligand_chain_ids
          )
        )
      );

      // generate requested number of models for the current complex
      size_t number_models( m_FlagNumberModels->GetFirstParameter()->GetNumericalValue< size_t>());
      for( size_t i( 0); i < number_models; ++i)
      {
        // use the same starting complex for each iteration
        util::ShPtr< assemble::ProteinModel> sp_starting_model( starting_model.HardCopy());
        BCL_MessageStd
        (
          "Start optimizating model #" + util::Format()( i) + "/" + util::Format()( number_models)
        );

        // model index to be appended to the filename when writing the optimized model to pdb file
        optimizer.Optimize( *sp_starting_model);

        // writes the MODEL into a pdb file
        io::OFStream write;

        // output prefix
        const std::string output_prefix( m_FlagOutputPrefix->GetFirstParameter()->GetValue());
        const std::string target_pdb_file
        (
          output_prefix + "_" + util::Format()( i) + ".pdb"
        );
        io::File::MustOpenOFStream( write, target_pdb_file);
        BCL_MessageStd( "Writing optimized model #" + util::Format()( i) + " to " + target_pdb_file);
        pdb_factory.WriteModelToPDB( *sp_starting_model, write);
        io::File::CloseClearFStream( write);
      }

      //
      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief reads in the native protein structure and add it to ProteinModelData of MODEL
    //! @param MODEL protein model to be pre-processed
    //! @return MODEL pre-processed protein model
    void ProteinDocking::PreProcess( assemble::ProteinModel &MODEL) const
    {

      // if it's a membrane protein complex or native structure is given
      if( m_FlagNative->GetFlag() || m_FlagMembrane->GetFlag())
      {
        // protein model data of the starting model
        util::ShPtr< assemble::ProteinModelData> sp_data( MODEL.GetProteinModelData());

        // pdb factory for creating a protein model from the native structure
        if( m_FlagNative->GetFlag())
        {
          // get PDB file name of the native structure
          const std::string native_pdb_file( m_FlagNative->GetFirstParameter()->GetValue());

          // construct a native model from the given PDB file
          const pdb::Factory pdb_factory;
          assemble::ProteinModel native_model( pdb_factory.ProteinModelFromPDBFilename( native_pdb_file));

          // add sse information for native model
          for
          (
            auto chain_itr( native_model.GetChains().Begin()), chain_itr_end( native_model.GetChains().End());
            chain_itr != chain_itr_end;
            ++chain_itr
          )
          {
            *chain_itr = util::ShPtr< assemble::Chain>
            (
              assemble::ConstructChainWithSSEsFromConformation( ( *chain_itr)->GetSequence()).Clone()
            );
          }

          //
          const std::string receptor_pdb_file( m_FlagReceptor->GetFirstParameter()->GetValue());
          assemble::ProteinModel receptor( pdb_factory.ProteinModelFromPDBFilename( receptor_pdb_file));
          std::string receptor_chain_ids( receptor.GetChainIDs());
          assemble::ProteinModel native_receptor( native_model.GetChains( receptor_chain_ids));

          // superimpose MODEL over native_model by receptor
          math::TransformationMatrix3D transformation
          (
            quality::RMSD().SuperimposeCoordinates
            (
              receptor.GetDefinedAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)),
              native_receptor.GetDefinedAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA))
            )
          );

          // transform the native model
          native_model.Transform( transformation);

          // inserts the native structure into ProteinModelData
          sp_data->Insert
          (
            assemble::ProteinModelData::e_NativeModel,
            util::ShPtr< assemble::ProteinModel>( new assemble::ProteinModel( native_model))
          );

          // write transformed native model to pdb file
          const std::string transformed_pdb_file( io::File::RemoveLastExtension( native_pdb_file) + "_transformed.pdb");
          io::OFStream write;
          io::File::MustOpenOFStream( write, transformed_pdb_file);
          pdb_factory.WriteModelToPDB( native_model, write);
          io::File::CloseClearFStream( write);
        }

        // if membrane protein complex, insert membrane object
        if( m_FlagMembrane->GetFlag())
        {
          // read in the PDBTM XML file
          io::IFStream pdbtm_xml;
          io::File::MustOpenIFStream( pdbtm_xml, m_FlagMembrane->GetFirstParameter()->GetValue());

          // construct a membrane from the given PDBTM XML file
          storage::Pair< biol::Membrane, math::TransformationMatrix3D> membrane_mat3d
          (
            biol::Membrane::MembraneAndTransformationFromPDBTMXML
              (
                pdbtm_xml,
                biol::Membrane::GetParameterTransitionThickness()->GetNumericalValue< double>(),
                biol::Membrane::GetParameterGapThickness()->GetNumericalValue< double>()
              )
          );
          util::ShPtr< biol::Membrane> sp_membrane( util::CloneToShPtr( membrane_mat3d.First()));

          // insert the membrane to protein model data
          sp_data->Insert
          (
            assemble::ProteinModelData::e_Membrane,
            sp_membrane
          );
        }

        // modifies ProteinModelData of MODEL
        MODEL.SetProteinModelData( sp_data);

        // assign residue packing density
        if( m_FlagPackingDensity->GetFlag())
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_FlagPackingDensity->GetFirstParameter()->GetValue());
          score::ProteinModelAANeighborhoodDocking::ReadPredictions
          (
            read,
            MODEL
          );
          io::File::CloseClearFStream( read);
        }
      }

      // add sse information for MODEL
//      for
//      (
//        auto chain_itr( MODEL.GetChains().Begin()), chain_itr_end(MODEL.GetChains().End());
//          chain_itr != chain_itr_end;
//        ++chain_itr
//      )
//      {
//        // center of the current chain
//        linal::Vector3D chain_center( ( *chain_itr)->GetCenter());
//        *chain_itr = util::ShPtr< assemble::Chain>
//        (
//          assemble::ConstructChainWithSSEsFromConformation( ( *chain_itr)->GetSequence()).Clone()
//        );
//        ( *chain_itr)->Translate( chain_center - ( *chain_itr)->GetCenter());
//      }
    }

    //! @brief Bring the ligand in contact with the receptor, the position of the ligand gets modified.
    //! @param RECEPTOR The receptor model
    //! @param LIGAND The ligand model
    void ProteinDocking::PreDock
    (
      const assemble::ProteinModel &RECEPTOR, assemble::ProteinModel &LIGAND
    ) const
    {
      // get all ligand atoms and all receptor atoms
      util::SiPtrVector< const biol::Atom> lig_atoms( LIGAND.GetAtoms());
      util::SiPtrVector< const biol::Atom> rec_atoms( RECEPTOR.GetAtoms());

      // maximum possible double values
      double shortest_distance( std::numeric_limits< double>::max());
      linal::Vector3D shortest_distance_vector3d
      (
        coord::CenterOfMass( RECEPTOR.GetAtomCoordinates(), true)
          - coord::CenterOfMass( LIGAND.GetAtomCoordinates(), true)
      );

      // iterate over ligand atoms
      for
      (
        auto lig_atom_itr( lig_atoms.Begin()), lig_atom_itr_end( lig_atoms.End());
          lig_atom_itr != lig_atom_itr_end;
        ++lig_atom_itr
      )
      {
        // iterate over atoms in the receptor
        for
        (
          auto rec_atom_itr( rec_atoms.Begin()), rec_atom_itr_end( rec_atoms.End());
            rec_atom_itr != rec_atom_itr_end;
          ++rec_atom_itr
        )
        {
          // current distance
          double current_distance
          (
            linal::Distance( ( *lig_atom_itr)->GetCoordinates(), ( *rec_atom_itr)->GetCoordinates())
          );

          // update shortest distance
          if( util::IsDefined( current_distance) && current_distance < shortest_distance)
          {
            shortest_distance = current_distance;
            shortest_distance_vector3d = ( *rec_atom_itr)->GetCoordinates() - ( *lig_atom_itr)->GetCoordinates();
          }
        }
      }

      // bring the ligand to the receptor so that they are in contact
      LIGAND.Translate( shortest_distance_vector3d);
    }

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &ProteinDocking::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &ProteinDocking::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! single instance of this class
    const ApplicationType ProteinDocking::ProteinDocking_Instance
    (
      GetAppGroups().AddAppToGroup( new ProteinDocking(), GetAppGroups().e_Protein)
    );

  } // namespace app
} // namespace bcl
