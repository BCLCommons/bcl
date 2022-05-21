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
#include "bcl_app_loop_template.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const ApplicationType LoopTemplate::LoopTemplates_Instance
    (
      GetAppGroups().AddAppToGroup( new LoopTemplate(), GetAppGroups().e_InternalBiol)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopTemplate::LoopTemplate() :
      m_PDBList
      (
        new command::FlagStatic
        (
          "pdb_list",
          "\tlist containing the paths to the pdbs from which to create the loop library",
          command::Parameter( "pdb_list_filename", "\tlist containing the paths to the pdbs")
        )
      ),
      m_TemplateLibraryFile
      (
        new command::FlagStatic
        (
          "library_output_path",
          "\tname of the template library",
          command::Parameter( "file name", "\tname of the template library", "loops.ls")
        )
       ),
      m_StatisticsFile
      (
        new command::FlagStatic
        (
          "statistics_output_path",
          "\tname of the statistics file",
          command::Parameter( "file name", "\tname of the statistics file", "loop_statistics.tbl")
        )
      )
    {
    }

    //! @brief clone function
    //! @return pointer to a new LoopTemplate
    LoopTemplate *LoopTemplate::Clone() const
    {
      return new LoopTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &LoopTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns a shared pointer to the command object
    //! @return shared pointer to the command object
    util::ShPtr< command::Command> LoopTemplate::InitializeCommand() const
    {
      // add the BCL and the application specific flags to the command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());
      sp_cmd->AddFlag( m_PDBList);
      sp_cmd->AddFlag( m_TemplateLibraryFile);
      sp_cmd->AddFlag( m_StatisticsFile);
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      return sp_cmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the main function of this application
    //! @return exit code - 0 for success
    int LoopTemplate::Main() const
    {
      // gather the loop conformations from the given pdb list
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      io::IFStream pdb_list;
      io::File::MustOpenIFStream( pdb_list, m_PDBList->GetFirstParameter()->GetValue());
      const storage::Vector< std::string> pdb_paths( util::StringLineListFromIStream( pdb_list));
      io::File::CloseClearFStream( pdb_list);
      storage::Vector< fold::LoopParameters> all_loops;
      for( auto pdb_it( pdb_paths.Begin()), pdb_it_end( pdb_paths.End()); pdb_it != pdb_it_end; ++pdb_it)
      {
        // read in the current model
        BCL_MessageTop( "Getting loop templates from " + *pdb_it);
        const assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( *pdb_it));

        // get the loops in the model
        all_loops.Append( GetLoops( model));
      }

      // write the results to the output file
      WriteResults( all_loops, m_TemplateLibraryFile->GetFirstParameter()->GetValue());

      // write out statistics about the loops
      WriteStatistics( all_loops, m_StatisticsFile->GetFirstParameter()->GetValue());

      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes out the results file
    //! @param PARAMETERS parameters of the gathered loop conformations
    //! @param FILE_NAME name of the output file
    void LoopTemplate::WriteResults
    (
      const storage::Vector< fold::LoopParameters> &PARAMETERS, const std::string &FILE_NAME
    )
    {
      // open a stream to the file
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILE_NAME);

      // write the values into the file
      write << PARAMETERS.GetSize() << std::endl;
      for( auto it( PARAMETERS.Begin()), it_end( PARAMETERS.End()); it != it_end; ++it)
      {
        const fold::LoopParameters &param( *it);
        write << param << std::endl;
      }
      io::File::CloseClearFStream( write);
    }

    //! @brief writes out statistics about the loops
    //! @param PARAMETERS parameters of the gathered loop conformations
    //! @param FILE_NAME name of the output file
    void LoopTemplate::WriteStatistics
    (
      const storage::Vector< fold::LoopParameters> &PARAMETERS, const std::string &FILE_NAME
    )
    {
      // open a stream to the file
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILE_NAME);

      // write the header of the statistics table
      write << "seq_length\tt_x\tt_y\tt_z\trot_x\trot_y\trot_z" << std::endl;

      // write the loop statistics
      for( auto it( PARAMETERS.Begin()), it_end( PARAMETERS.End()); it != it_end; ++it)
      {
        // get the statistics for this loop
        const fold::LoopParameters &param( *it);
        const size_t seq_length( param.GetSequenceDistance());
        const linal::Vector3D &translation( param.GetTranslation());
        const linal::Vector3D &rotation( param.GetRotation());

        // write the statistics into the table
        write << seq_length << "\t" << translation( 0) << "\t" << translation( 1) << "\t";
        write << translation( 2) << "\t" << rotation( 0) << "\t" << rotation( 1) << "\t";
        write << rotation( 2) << std::endl;
      }
      io::File::CloseClearFStream( write);
    }

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &LoopTemplate::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &LoopTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the loop conformations in the given protein model
    //! @detail loop conformations are parameterized through their dihedral angles
    //! @param MODEL protein model to get loop conformations from
    //! @return the loop conformations in the given protein model
    storage::Vector< fold::LoopParameters> LoopTemplate::GetLoops( const assemble::ProteinModel &MODEL)
    {
      // get the anchor points of the loops in the model
      storage::Vector< fold::LoopParameters> loop_parameters;
      storage::Vector< storage::Pair< storage::VectorND< 2, int>, char> > anchor_points;
      const util::SiPtrVector< const assemble::SSE> loops( MODEL.GetSSEs( biol::GetSSTypes().COIL));
      for( auto loop_it( loops.Begin()), loop_it_end( loops.End()); loop_it != loop_it_end; ++loop_it)
      {
        storage::VectorND< 2, int> anchors_ids
        (
          ( **loop_it).GetFirstAA()->GetSeqID() - 1, ( **loop_it).GetLastAA()->GetSeqID() + 1
        );
        if
        (
          anchors_ids( 0) < 1 ||
          anchors_ids( 1) > ( int) MODEL.GetChain( ( **loop_it).GetChainID())->GetNumberAAs()
        )
        {
          continue;
        }

        // anchor points of the loop
        storage::Pair< storage::VectorND< 2, int>, char> anchors
        (
          anchors_ids, ( **loop_it).GetFirstAA()->GetChainID()
        );

        anchor_points.PushBack( anchors);
      }

      // compute the loop parameters
      for( auto an_it( anchor_points.Begin()), an_end( anchor_points.End()); an_it != an_end; ++an_it)
      {
        // get the anchor residues
        const util::SiPtrVector< const biol::AABase> aas( MODEL.GetChain( an_it->Second())->GetAminoAcids());
        const int anchor_first_id( ( *an_it).First().First());
        const int anchor_second_id( ( *an_it).First().Second());
        if( anchor_first_id <= 1 || anchor_second_id >= ( int) aas.GetSize() - 1)
        {
          continue;
        }
        const biol::AABase &aa_first( *aas( ( *an_it).First().First() - 1));
        const biol::AABase &aa_second( *aas( ( *an_it).First().Second() - 1));

        // store the dihedral angles
        storage::Vector< double> dihedral_angles;

        // add psi angle of n-terminal anchor
        const biol::AABase &following_aa( *aas( ( *an_it).First().First()));
        const double n_psi( aa_first.CalculatePsi( following_aa.GetAtom( biol::GetAtomTypes().N)));
        dihedral_angles.Append( n_psi);

        // add phi and psi angles of residues in the loop
        for( int seq_id( anchor_first_id + 1); seq_id < anchor_second_id; ++seq_id)
        {
          const biol::AABase &aa_prev( *aas( seq_id - 2));
          const biol::AABase &aa_curr( *aas( seq_id - 1));
          const biol::AABase &aa_next( *aas( seq_id));
          storage::Pair< double, double> angles
          (
            aa_curr.CalculatePhiPsi
            (
              aa_prev.GetAtom( biol::GetAtomTypes().C), aa_next.GetAtom( biol::GetAtomTypes().N)
            )
          );
          dihedral_angles.Append( angles.First());
          dihedral_angles.Append( angles.Second());
        }

        // add phi angle of c-terminal anchor
        const biol::AABase &previous_aa( *aas( ( *an_it).First().Second() - 2));
        const double c_phi( aa_second.CalculatePhi( previous_aa.GetAtom( biol::GetAtomTypes().C)));
        dihedral_angles.Append( c_phi);

        // add the parameters to the list if they are sane
        util::ShPtr< fold::LoopParameters> new_loop
        (
          fold::LoopParameters::Create( aa_first, aa_second, dihedral_angles)
        );
        if( IsValid( *new_loop))
        {
          loop_parameters.PushBack( *fold::LoopParameters::Create( aa_first, aa_second, dihedral_angles));
        }
        else
        {
          BCL_MessageCrt( "template is invalid");
        }
      }

      return loop_parameters;
    }

    //! @brief checks if the given loop parameters are valid
    //! @param LOOP loop parameters to be checked
    //! @return true, if the given loop parameters are valid
    bool LoopTemplate::IsValid( const fold::LoopParameters &LOOP)
    {
      // check validity of dihedral angles
      const storage::Vector< double> &angles( LOOP.GetAngles());
      for( auto da_it( angles.Begin()), da_it_end( angles.End()); da_it != da_it_end; ++da_it)
      {
        // check if angle is defined
        if( !util::IsDefined( *da_it))
        {
          return false;
        }
      }

      // check validity of translation vector
      if( !LOOP.GetTranslation().IsDefined())
      {
        return false;
      }

      // check validity of Euler angles
      const linal::Vector3D &rotation( LOOP.GetRotation());
      for( size_t i( 0); i < 3; ++i)
      {
        // check if Euler angle is defined
        if( !util::IsDefined( rotation( i)))
        {
          return false;
        }
      }

      return true;
    }

  } // namespace app
} // namespace bcl
