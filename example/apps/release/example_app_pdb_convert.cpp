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

// include example header
#include "example.h"
// include the header of the class which this example is for
//#include "release/bcl_app_pdb_convert.cpp"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_pdb_convert.cpp
  //!
  //! @author woetzen
  //! @date November 17, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppPDBConvert :
    public ExampleInterface
  {
  public:

    ExampleAppPDBConvert *Clone() const
    {
      return new ExampleAppPDBConvert( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      const app::ApplicationType app_enum_pdb_convert( "PDBConvert");
      BCL_ExampleAssert( app_enum_pdb_convert.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // pdb convert can be called without any flags
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // check that flags are needed
        BCL_ExampleCheck( pdb_convert_helper.CheckCommandString( false), true);
      }

      // chain A as fasta
      // PDBConvert 1C1D.pdb -chains A -fasta -output_prefix 1C1D
      // output: 1C1DA.fasta
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1C1D")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "chains", "A");
        pdb_convert_helper.SetFlag( "fasta");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // read in the fasta
        const std::string fasta_filename( output_prefix + "A.fasta");
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, fasta_filename);

        // create a sequence and read in the fasta sequence
        const biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read, biol::GetAAClasses().e_AA, 'A'));
        io::File::CloseClearFStream( read);

        // check length
        BCL_ExampleCheck( seq.GetSize(), 355);
      }

      // chain B as fasta
      // PDBConvert 1C1D.pdb -chains B -fasta -output_prefix 1C1D
      // output: 1C1DB.fasta
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1C1D")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "chains", "B");
        pdb_convert_helper.SetFlag( "fasta");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // read in the fasta
        const std::string fasta_filename( output_prefix + "B.fasta");
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, fasta_filename);

        // create a sequence and read in the fasta sequence
        const biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read, biol::GetAAClasses().e_AA, 'B'));
        io::File::CloseClearFStream( read);

        // check length
        BCL_ExampleCheck( seq.GetSize(), 355);
      }

      // chain A and B as fasta
      // PDBConvert 1C1D.pdb -chains A B -fasta -output_prefix 1C1D
      // output: 1C1DA.fasta 1C1DB.fasta
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1C1D")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "chains", storage::Vector< std::string>::Create( "A", "B"));
        pdb_convert_helper.SetFlag( "fasta");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);
      }

      // 02 SEQRES from fasta and protein params
      // PDBConvert -pdb_from_fasta 1C1DA.fasta 1C1DB.fasta -output_prefix 1C1D
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const storage::Vector< std::string> fasta_file_names
        (
          storage::Vector< std::string>::Create
          (
            AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1C1DA.fasta"),
            AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1C1DB.fasta")
          )
        );
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1C1D")
        );

        // flags
        pdb_convert_helper.SetFlag( "protein_params");
        pdb_convert_helper.SetFlag( "pdb_from_fasta", fasta_file_names);
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name( output_prefix + "fasta.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name, biol::GetAAClasses().e_AA));

        // empty model, but two chains A and B with sequence (no sses or coordinates)
        BCL_ExampleAssert( model.GetChain( 'A').IsDefined(), true);
        BCL_ExampleCheck( model.GetChain( 'A')->GetSequence()->GetSize(), 355);
        BCL_ExampleAssert( model.GetChain( 'B').IsDefined(), true);
        BCL_ExampleCheck( model.GetChain( 'B')->GetSequence()->GetSize(), 355);
      }

      // renumbering of sequence id and atom id
      // PDBConvert 1IE9.pdb -bcl_pdb Split -chains A -output_prefix 1IE9
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1IE9")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Split");
        pdb_convert_helper.SetFlag( "chains", "A");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "aaclass", "AABackBone");

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name_out( output_prefix + "A.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name_out, biol::GetAAClasses().e_AABackBone));

        // locate amino acid by pdb id or seq id
        const assemble::LocatorAA locate_by_pdbid( 'A', 10, true);
        const assemble::LocatorAA locate_by_seqid( 'A', 10, false);

        BCL_ExampleCheck( locate_by_seqid.Locate( model)->GetSeqID(), 10);
        BCL_ExampleCheck( locate_by_pdbid.Locate( model)->GetPdbID(), 10);
        BCL_ExampleCheck( locate_by_pdbid.Locate( model)->GetSeqID() == locate_by_seqid.Locate( model)->GetSeqID(), true);
        BCL_ExampleCheck( locate_by_seqid.Locate( model)->GetCA().GetPdbID(), 47);
      }

      // prevent renumbering of sequence id and atom id
      // PDBConvert 1IE9.pdb -bcl_pdb Split -chains A -output_prefix 1IE9 -write_pdb_res_ids -write_pdb_atom_ids
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1IE9_no_renumber")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Split");
        pdb_convert_helper.SetFlag( "chains", "A");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "aaclass", "AABackBone");
        pdb_convert_helper.SetFlag( "write_pdb_res_ids");
        pdb_convert_helper.SetFlag( "write_pdb_atom_ids");

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name_out( output_prefix + "A.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name_out, biol::GetAAClasses().e_AABackBone));

        // locate amino acid by pdb id or seq id
        const assemble::LocatorAA locate_by_pdbid( 'A', 10, true);
        const assemble::LocatorAA locate_by_seqid( 'A', 10, false);

        BCL_ExampleCheck( locate_by_seqid.Locate( model)->GetSeqID(), 10);
        BCL_ExampleCheck( locate_by_seqid.Locate( model)->GetCA().GetPdbID(), 74);
        BCL_ExampleCheck( locate_by_pdbid.Locate( model).IsDefined(), false);
      }

      // select sses based on size
      // PDBConvert 1ubi.pdb -bcl_pdb Split -chains A -output_prefix 1ubi_select -min_sse_size 5 5 999
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1ubi_select")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Split");
        pdb_convert_helper.SetFlag( "chains", "A");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "min_sse_size", storage::Vector< std::string>::Create( "5", "5", "999"));

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name_out( output_prefix + "A.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name_out, biol::GetAAClasses().e_AABackBone));

        // number of sses (including coils)
        BCL_ExampleCheck( model.GetNumberSSEs(), 10);
      }

      // renaming chain ids flag A to ' ' and B to C
      // PDBConvert 1C1D.pdb -bcl_pdb -rename_chain_id A, B,C -output_prefix 1C1D
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "1C1D.pdb")
        );
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1C1D")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Full");
        pdb_convert_helper.SetFlag( "rename_chain_id", storage::Vector< std::string>::Create( "A,", "B,C"));

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name_out( output_prefix + "bcl.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name_out, biol::GetAAClasses().e_AA));

        // empty model, but two chains '  and C with sequence (no sses or coordinates)
        BCL_ExampleAssert( model.GetChain( ' ').IsDefined(), true);
        BCL_ExampleAssert( model.GetChain( 'C').IsDefined(), true);
      }

      // superimpose on other pdb
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string pdb_file_name_superimpose( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1ubi_select")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Full");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix + "superimposed");
        pdb_convert_helper.SetFlag( "superimpose", storage::Vector< std::string>::Create( pdb_file_name_superimpose, "RMSD", "CA"));

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name_out( output_prefix + "superimposedbcl.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name, biol::GetAAClasses().e_AABackBone));
        const assemble::ProteinModel model_ideal( Proteins::GetModel( pdb_file_name_superimpose, biol::GetAAClasses().e_AABackBone));
        const assemble::ProteinModel model_superimposed( Proteins::GetModel( pdb_file_name_out, biol::GetAAClasses().e_AABackBone));

        // calculate the rmsd
        const double rmsd_before( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD_NoSuperimposition, model_ideal, model, storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)));
        const double rmsd_after( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD_NoSuperimposition, model_superimposed, model, storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)));

        BCL_ExampleCheck( rmsd_before > rmsd_after, true);
      }

      // generate biomolecule
      // PDBConvert 1lgh.pdb -bcl_pdb Full -biomolecule 1 -output_prefix 1lgh
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1lgh.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1lgh")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Full");
        pdb_convert_helper.SetFlag( "biomolecule", "1");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // read in the fasta
        const std::string pdb_filename_out( output_prefix + "bcl.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_filename_out, biol::GetAAClasses().e_AA));

        // check number of chains
        BCL_ExampleCheck( model.GetChains().GetSize(), 16);
      }

      // generate and translate biomolecule from pdbtm xml file
      // PDBConvert 1LGH.pdb -bcl_pdb Full -biomolecule 1 -pdbtm_xml 1lgh.xml -output_prefix 1LGHpdbtm
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1lgh.pdb"));
        const std::string xml_file_name( AddExampleInputPathToFilename( e_Biology, "1lgh.xml"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1lghpdbtm")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Full");
        pdb_convert_helper.SetFlag( "biomolecule", "1");
        pdb_convert_helper.SetFlag( "pdbtm_xml", xml_file_name);
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // read in the fasta
        const std::string pdb_filename_out( output_prefix + "bcl.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_filename_out, biol::GetAAClasses().e_AA));

        // check number of chains
        BCL_ExampleCheck( model.GetChains().GetSize(), 16);
      }

      // convert unnatural amino acids to their natural "parent" they were derived from
      // PDBConvert 2YV8.pdb -bcl_pdb Full -output_prefix 2YV8 -convert_to_natural_aa_type
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "2yv8.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "2yv8")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Full");
        pdb_convert_helper.SetFlag( "convert_to_natural_aa_type");
        pdb_convert_helper.SetFlag( "write_pdb_res_ids");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // read in the fasta
        const std::string pdb_filename_out( output_prefix + "bcl.pdb");
        const assemble::ProteinModel orig_model( Proteins::GetModel( pdb_file_name, biol::GetAAClasses().e_AA));
        const assemble::ProteinModel model( Proteins::GetModel( pdb_filename_out, biol::GetAAClasses().e_AA));

        // check amino acid type
        const assemble::LocatorAA loc_aa( 'A', 115, true);
        BCL_ExampleCheck( loc_aa.Locate( orig_model)->GetType(), biol::GetAATypes().MSE);
        BCL_ExampleCheck( loc_aa.Locate( orig_model)->GetType()->GetParentType(), biol::GetAATypes().MET);
        BCL_ExampleCheck( loc_aa.Locate( model)->GetType(), biol::GetAATypes().MET);
      }

// stripping sidecahins
// PDBConvert 1UBI.pdb -bcl_pdb -aaclass AABackBone -output_prefix 1UBI

// adding sidechains
// PDBConvert 1UBI.pdb -bcl_pdb -aaclass AABackBone -output_prefix 1UBI -side_chains

      // adding zero coordinates for undefined atom coordinates and generating rosetta loop file http://www.rosettacommons.org/manuals/rosetta3_user_guide/app_loop.html
      // PDBConvert 1ubi.pdb -bcl_pdb Split -output_prefix 1ubi_nosplit -min_sse_size 0 0 999 -write_zero_coordinates -loop_file_rosetta KIC
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1ubi_noloop")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "bcl_pdb", "Split");
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "min_sse_size", storage::Vector< std::string>::Create( "0", "0", "999"));
        pdb_convert_helper.SetFlag( "write_zero_coordinates");
        pdb_convert_helper.SetFlag( "loop_file_rosetta", "KIC");

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name_out( output_prefix + "A.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name_out, biol::GetAAClasses().e_AABackBone));

        // number of sses (including coils)
        const assemble::LocatorAA loc_aa( 'A', 8, true);
        BCL_ExampleCheck( loc_aa.Locate( model)->GetCA().GetCoordinates().Sum(), 0.0);
      }

      // write topology file
      // PDBConvert 1ubi.pdb -output_prefix 1ubi -topology
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1ubi")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "topology");

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);
      }

      // write dssp file and reassign secondary structure
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( pdb::GetNamespaceIdentifier(), "1ubi_dssp")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "dssp");
        pdb_convert_helper.SetFlag( "bcl_pdb");

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // create a protein model
        const std::string pdb_file_name_out( output_prefix + "bcl.pdb");
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file_name_out, biol::GetAAClasses().e_AABackBone));

        BCL_ExampleAssert( model.GetNumberSSE( biol::GetSSTypes().STRAND), 5);
        BCL_ExampleAssert( model.GetNumberSSE( biol::GetSSTypes().HELIX), 1);
      }

      // test sdf flag
      {
        ApplicationExampleHelper pdb_convert_helper( app_enum_pdb_convert);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "all_aas_connected.corina.pdb"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "testallaas")
        );

        // flags
        pdb_convert_helper.AddParameter( pdb_file_name);
        pdb_convert_helper.SetFlag( "output_prefix", output_prefix);
        pdb_convert_helper.SetFlag( "sdf");

        // check the command line
        BCL_ExampleAssert( pdb_convert_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( pdb_convert_helper.RunCommand(), 0);

        // check that the produced sdf is correct
        const std::string sdf_file_name_out( output_prefix + "_.sdf");
        const std::string sdf_file_name_out_correct( output_prefix + "_.correct.sdf");

        if
        (
          BCL_ExampleCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance( sdf_file_name_out, sdf_file_name_out_correct, 0.0001),
            true
          )
        )
        {
          remove( sdf_file_name_out.c_str());
        }
      }

      // reset all pdb factory flags, since this application changes them
      pdb::Factory::ResetFlagDefaults();

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppPDBConvert

  const ExampleClass::EnumType ExampleAppPDBConvert::s_Instance
  (
    GetExamples().AddEnum( ExampleAppPDBConvert())
  );

} // namespace bcl
