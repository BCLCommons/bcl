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
#include "bcl_app_analyze_spin_label_parameters.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const ApplicationType AnalyzeSpinLabelParameters::AnalyzeSpinLabelParameters_Instance
    (
      GetAppGroups().AddAppToGroup( new AnalyzeSpinLabelParameters(), GetAppGroups().e_Restraint)
    );

    //! header of the statistics table
    const storage::TableHeader AnalyzeSpinLabelParameters::s_ResultHeader
    (
      storage::TableHeader::Create( "seq_id", "exposure", "r", "theta", "phi")
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeSpinLabelParameters::AnalyzeSpinLabelParameters() :
      m_PDBList
      (
        new command::FlagStatic
        (
          "pdb_list",
          "\tlist containing the paths to the pdbs that should be analyzed",
          command::Parameter( "pdb_list_filename", "\tlist containing the paths to the pdbs")
        )
      ),
      m_Native
      (
        new command::FlagStatic
        (
          "native",
          "\tpath to the native the spin labels should be compared to",
          command::Parameter( "native_filename", "\tpath to the native")
        )
      ),
      m_OutputPrefix
      (
        new command::FlagDynamic
        (
          "output_prefix",
          "\tprefix for the output files",
          command::Parameter( "prefix", "\tprefix for the output files", "")
        )
      )
    {
    }

    //! @brief clone function
    //! @return pointer to a new AnalyzeSpinLabelParameters
    AnalyzeSpinLabelParameters *AnalyzeSpinLabelParameters::Clone() const
    {
      return new AnalyzeSpinLabelParameters( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &AnalyzeSpinLabelParameters::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns a shared pointer to the command object
    //! @return shared pointer to the command object
    util::ShPtr< command::Command> AnalyzeSpinLabelParameters::InitializeCommand() const
    {
      // add the BCL and the application specific flags to the command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());
      sp_cmd->AddFlag( m_PDBList);
      sp_cmd->AddFlag( m_Native);
      sp_cmd->AddFlag( m_OutputPrefix);
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      return sp_cmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the main function of this application
    //! @return exit code - 0 for success
    int AnalyzeSpinLabelParameters::Main() const
    {
      // read in the native structure
      io::IFStream native_file;
      io::File::MustOpenIFStream( native_file, m_Native->GetFirstParameter()->GetValue());
      const pdb::Handler native_handler( native_file, true);
      const pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
      const assemble::ProteinModel native( factory.ProteinModelFromPDB( native_handler));
      io::File::CloseClearFStream( native_file);

      // compute the exposure of the residues in the native structure
      const util::ShPtr< storage::Vector< double> > sp_exposure( ComputeExposure( native));

      // gather the statistics for the spin label conformations in the models provided through the pdb list
      io::IFStream pdb_list;
      io::File::MustOpenIFStream( pdb_list, m_PDBList->GetFirstParameter()->GetValue());
      const storage::Vector< std::string> pdb_paths( util::StringLineListFromIStream( pdb_list));
      io::File::CloseClearFStream( pdb_list);
      util::ShPtr< storage::Table< double> > sp_statistics( new storage::Table< double>( s_ResultHeader));
      typedef storage::Vector< std::string>::const_iterator const_list_it;
      for( const_list_it pdb_it( pdb_paths.Begin()), pdb_it_end( pdb_paths.End()); pdb_it != pdb_it_end; ++pdb_it)
      {
        // read in the current model
        BCL_MessageTop( "analyzing " + *pdb_it);
        io::IFStream pdb_file;
        io::File::MustOpenIFStream( pdb_file, *pdb_it);
        const pdb::Handler pdb_handler( pdb_file, true);
        const assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb_handler));

        // find the spin labels in the sequence
        util::ShPtrList< biol::AABase> spin_labels;
        storage::List< std::string> names;
        const util::SiPtrVector< const biol::AASequence> sequences( model.GetSequences());
        typedef util::SiPtrVector< const biol::AASequence>::const_iterator const_seq_it;
        for( const_seq_it seq_it( sequences.Begin()), seq_it_end( sequences.End()); seq_it != seq_it_end; ++seq_it)
        {
          spin_labels.Append( ( **seq_it).GetData( biol::GetAATypes().R1A));
          names.Append( *pdb_it);
        }

        // compute the statistics for the spin label conformations
        const util::ShPtr< storage::Table< double> > sp_statistics_tmp( ComputeStatistics( spin_labels, sp_exposure, names));
        sp_statistics->Append( *sp_statistics_tmp, true);
      }

      // write out the statistics
      io::OFStream statistics_table;
      const std::string table_file( m_OutputPrefix->GetFirstParameter()->GetValue() + "_statistics.tbl");
      io::File::MustOpenOFStream( statistics_table, table_file);
      sp_statistics->WriteFormatted( statistics_table);
      io::File::CloseClearFStream( statistics_table);

      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &AnalyzeSpinLabelParameters::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &AnalyzeSpinLabelParameters::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief computes the statistics for the given spin labels
    //! @param SPIN_LABELS spin labels for which to compute the statistics for
    //! @param SP_EXPOSURE exposure values for the residues of the protein
    //! @return statistics regarding the given spin labels
    util::ShPtr< storage::Table< double> > AnalyzeSpinLabelParameters::ComputeStatistics
    (
      const util::ShPtrList< biol::AABase> &SPIN_LABELS,
      const util::ShPtr< storage::Vector< double> > &SP_EXPOSURE,
      const storage::List< std::string> &MODEL_PATHS
    )
    {
      // compute statistics for each spin label
      util::ShPtr< storage::Table< double> > sp_statistics( new storage::Table< double>( s_ResultHeader));
      typedef util::ShPtrList< biol::AABase>::const_iterator const_sl_it;
      storage::List< std::string>::const_iterator name_it( MODEL_PATHS.Begin());
      for( const_sl_it sl_it( SPIN_LABELS.Begin()), sl_it_end( SPIN_LABELS.End()); sl_it != sl_it_end; ++sl_it)
      {
        // get the relevant atom coordinates
        const biol::AABase &spin_label( **sl_it);
        const linal::Vector3D &n( spin_label.GetAtom( biol::GetAtomTypes().N).GetCoordinates());
        const linal::Vector3D &ca( spin_label.GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
        const linal::Vector3D &cb( spin_label.GetAtom( biol::GetAtomTypes().CB).GetCoordinates());
        const linal::Vector3D &o1( spin_label.GetAtom( biol::GetAtomTypes().O1).GetCoordinates());

        // compute the reference plane
        const linal::Vector3D ca_n( ( n - ca).Normalize());
        const linal::Vector3D ca_cb( ( cb - ca).Normalize());
        const linal::Vector3D ca_cb_ortho( ( ca_cb - ( ( ca_cb * ca_n) * ca_n)).Normalize());
        const linal::Vector3D plane_normal( linal::CrossProduct( ca_n, ca_cb).Normalize());

        // compute the spherical coordinates
        const double r( ( o1 - cb).Norm());
        const linal::Vector3D cb_o1( ( o1 - cb).Normalize());
        const double theta( linal::ProjAngle( cb_o1, plane_normal));
        const linal::Vector3D proj_1( ( cb_o1 * ca_n) * ca_n);
        const linal::Vector3D proj_2( ( cb_o1 * ca_cb_ortho) * ca_cb_ortho);
        const linal::Vector3D proj_res( ( proj_1 + proj_2).Normalize());
        const double phi( linal::ProjAngle( proj_res, ca_cb));

        // add the coordinates to the statistics table
        const double exposure( ( *SP_EXPOSURE)( spin_label.GetSeqID() - 1));
        const int seq_id( spin_label.GetSeqID());
        storage::Vector< double> statistics( storage::Vector< double>::Create( seq_id, exposure, r, theta, phi));
        sp_statistics->InsertRow( *name_it, statistics, true);
        ++name_it;
      }

      return sp_statistics;
    }

    //! @brief computes the exposure of the residues in the given protein model
    //! @param MODEL protein model to calculate the exposure for
    //! @return shared pointer to the exposure of the residues in the given protein model
    util::ShPtr< storage::Vector< double> > AnalyzeSpinLabelParameters::ComputeExposure
    (
      const assemble::ProteinModel &MODEL
    )
    {
      // compute the neighbor lists for each residue in this model
      const util::SiPtrVector< const biol::AABase> amino_acids( MODEL.GetAminoAcids());
      util::ShPtr< storage::Vector< double> > sp_exposure( new storage::Vector< double>);
      const assemble::AANeighborCount neighbor_count;
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer>
      > sp_generator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          neighbor_count.GetThresholdRange().GetMax(),
          neighbor_count.GetMinimalSequenceSeparation(),
          false,
          false
        )
      );
      const assemble::AANeighborListContainer neighbor_container( sp_generator->operator()( MODEL));

      // calculate the exposure for each residue in this model
      typedef assemble::AANeighborListContainer::const_iterator const_nl_it;
      for( const_nl_it it( neighbor_container.Begin()), it_end( neighbor_container.End()); it != it_end; ++it)
      {
        sp_exposure->PushBack( neighbor_count( it->second));
      }

      return sp_exposure;
    }

  } // namespace app
} // namespace bcl
