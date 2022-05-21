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
#include "biol/bcl_biol_aa_back_bone_completer.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_back_bone_completer.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAABackBoneCompleter :
    public ExampleInterface
  {
  public:

    ExampleBiolAABackBoneCompleter *Clone() const
    {
      return new ExampleBiolAABackBoneCompleter( *this);
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

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // get the protein model
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone)
      );
      // create pointers on the chain and the sequence
      util::ShPtr< assemble::Chain> sp_chain( model.GetChain( 'A'));
      util::ShPtr< biol::AASequence> sp_sequence( sp_chain->GetSequence());

      // locate strand_1_7
      util::SiPtr< const assemble::SSE> strand_1_7( assemble::LocatorSSE( 'A', 1, 7).Locate( model));

      // make a copy of the sequence of this strand
      util::ShPtrVector< biol::AABase> new_strand_sequence( strand_1_7->GetData().HardCopy());

      // initialize an oxygen atom with undefined coordinates
      biol::Atom oxygen_undefined( biol::GetAtomTypes().O);

      // set the coordinates for oxygens of residues 2 to 4 to undefined
      BCL_MessageStd
      (
        "Removing oxygens for following residues:" +
        new_strand_sequence( 1)->GetIdentification() + " and " +
        new_strand_sequence( 2)->GetIdentification() + " and " +
        new_strand_sequence( 3)->GetIdentification()
      );
      new_strand_sequence( 1)->SetAtom( oxygen_undefined);
      new_strand_sequence( 2)->SetAtom( oxygen_undefined);
      new_strand_sequence( 3)->SetAtom( oxygen_undefined);

      // create a new strand and insert it into the model
      util::ShPtr< assemble::SSE> new_strand_1_7( new assemble::SSE( new_strand_sequence, biol::GetSSTypes().STRAND));
      model.Replace( new_strand_1_7);

      // locate helix_23_34
      util::SiPtr< const assemble::SSE> helix_23_34( assemble::LocatorSSE( 'A', 23, 34).Locate( model));

      // make a copy of the sequence of this strand
      util::ShPtrVector< biol::AABase> new_helix_sequence( helix_23_34->GetData());

      // set the coordinates for oxygens of residues 26 to 28
      BCL_MessageStd
      (
        "Removing oxygens for following residues:" +
        new_helix_sequence( 1)->GetIdentification() + " and " +
        new_helix_sequence( 2)->GetIdentification() + " and " +
        new_helix_sequence( 3)->GetIdentification()
      );
      new_helix_sequence( 3)->SetAtom( oxygen_undefined);
      new_helix_sequence( 4)->SetAtom( oxygen_undefined);
      new_helix_sequence( 5)->SetAtom( oxygen_undefined);

      // create a new strand and insert it into the model
      util::ShPtr< assemble::SSE> new_helix_23_34( new assemble::SSE( new_helix_sequence, biol::GetSSTypes().HELIX));
      model.Replace( new_helix_23_34);

      // set the flag for writing hydrogens to true
      pdb::Factory::GetFlagWriteHydrogens()->SetFlag();

    /////////////////
    // data access //
    /////////////////

      // default constructor
      BCL_MessageStd( "test default constructor");
      biol::AABackBoneCompleter completer_def;

      // constructor from add_hydrogen and add_oxygen booleans
      BCL_MessageStd( "test constructor from add_hydrogen and add_oxygen booleans");
      biol::AABackBoneCompleter completer( true, true, false);

      // constructor from add_hydrogen, add_oxygen and add_ha booleans
      BCL_MessageStd( "test constructor from add_hydrogen, add_oxygen and add_ha booleans");
      biol::AABackBoneCompleter completer_ha( true, true, true);

      BCL_MessageStd( "test clone constructor");
      util::ShPtr< biol::AABackBoneCompleter> sp_completer( completer.Clone());
      BCL_Example_Check
      (
        sp_completer->GetAddAmideHydrogens() == completer.GetAddAmideHydrogens() &&
        sp_completer->GetAddCarbonylOxygens() == completer.GetAddCarbonylOxygens(),
        "The cloned completer is different!\n" + util::Format()( *sp_completer) + "\nvs\n" + util::Format()( completer)
      );

    /////////////////
    // data access //
    /////////////////

       // test the GetClassIdentifier
       BCL_ExampleCheck( completer.GetClassIdentifier(), GetStaticClassName( completer));

       // test GetAddAmideHydrogens()
       BCL_ExampleCheck( completer.GetAddAmideHydrogens(), true);

       // test GetAddCarbonylOxygens()
       BCL_MessageStd( "test GetAddCarbonylOxygens()");
       BCL_ExampleIndirectCheck( completer.GetAddCarbonylOxygens(), true, "Constructor");

       // test GetAddHAHydrogens()
       BCL_MessageStd( "test GetAddHAHydrogens()");
       BCL_ExampleIndirectCheck( completer_ha.GetAddHAHydrogens(), true, "Constructor");

    ////////////////
    // operations //
    ////////////////

       // create a new model with hydrogens and oxygen
       BCL_MessageStd( "testing CompleterProteinModel()");
       util::ShPtr< assemble::ProteinModel> new_model( completer.CompleteProteinModel( model));

       // write the model to pdb
       BCL_MessageStd( "writing model with hydrogens to pdb");
       Proteins::WriteModelToPDB( *new_model, AddExampleOutputPathToFilename( completer, "1ubi_with_hydrogens.pdb"));

       // locate the residues with missing oxygen and make sure they now have defined carbonyl oxygen coordinates
       BCL_MessageStd( "checking residues that had undefined oxygen coordinates previously");
       BCL_MessageStd( "checking residue 3");

       // locate the residue and the atom for residue 3
       util::SiPtr< const biol::AABase> residue_3( assemble::LocatorAA( 'A', 3).Locate( *new_model));
       const biol::Atom oxygen_3( residue_3->GetAtom( biol::GetAtomTypes().O));
       BCL_MessageStd( "coordinates of oxygen 3: " + util::Format()( oxygen_3.GetCoordinates()));
       // make sure the oxygen is defined
       BCL_Example_Check
       (
         oxygen_3.GetType().IsDefined() && oxygen_3.GetCoordinates().IsDefined(),
         "Carbonyl oxygen on residue 3 is not defined!"
       );

       // locate the residue and the atom for residue 27
       BCL_MessageStd( "checking residues 27");
       util::SiPtr< const biol::AABase> residue_27( assemble::LocatorAA( 'A', 27).Locate( *new_model));
       const biol::Atom oxygen_27( residue_27->GetAtom( biol::GetAtomTypes().O));
       BCL_MessageStd( "coordinates of oxygen 27: " + util::Format()( oxygen_27.GetCoordinates()));
       // make sure the oxygen is defined
       BCL_Example_Check
       (
         oxygen_27.GetType().IsDefined() && oxygen_27.GetCoordinates().IsDefined(),
         "Carbonyl oxygen or residue 27 is not defined!"
       );

       // create a new model with H, HA, and O atoms
       BCL_MessageStd( "testing CompleterProteinModel()");
       util::ShPtr< assemble::ProteinModel> new_model_ha( completer_ha.CompleteProteinModel( model));

       // write the model to pdb
       BCL_MessageStd( "writing model with H, HA, and O to pdb");
       Proteins::WriteModelToPDB( *new_model_ha, AddExampleOutputPathToFilename( completer, "1ubi_with_h_ha_o.pdb"));

       // locate the residue and the atom for residue 3
       util::SiPtr< const biol::AABase> residue_3_ha( assemble::LocatorAA( 'A', 3).Locate( *new_model_ha));
       const biol::Atom ha_3( residue_3_ha->GetAtom( biol::GetAtomTypes().HA));
       BCL_MessageStd( "coordinates of HA 3: " + util::Format()( ha_3.GetCoordinates()));
       // make sure the HA is defined
       BCL_Example_Check
       (
         ha_3.GetType().IsDefined() && ha_3.GetCoordinates().IsDefined(),
         "HA on residue 3 is not defined!"
       );

       // add H N CA
       {
         const biol::AABase &aa8( *model.GetChain( 'A')->GetSequence()->GetAA( 8));
         const storage::Map< biol::AtomType, biol::Atom> pseudo_atoms( biol::AABackBoneCompleter::GenerateHNCA( *strand_1_7->GetLastAA()));
         // write strand with additional atoms
         io::OFStream write;
         const std::string filename( AddExampleOutputPathToFilename( biol::AABackBoneCompleter(), "strand_HNCA.pdb"));
         BCL_ExampleMustOpenOutputFile( write, filename);
         size_t atom_id( 1);
         pdb::Handler new_handler;
         new_handler.AppendLines( pdb::Factory::WriteAASequenceToLines( *strand_1_7, atom_id));
         // iterate over atoms
         for( storage::Map< biol::AtomType, biol::Atom>::const_iterator itr( pseudo_atoms.Begin()), itr_end( pseudo_atoms.End()); itr != itr_end; ++itr)
         {
           new_handler.PushBack( pdb::Factory::WriteAtomToLine( itr->second, aa8, aa8.GetChainID(), atom_id));
         }
         new_handler.WriteLines( write);
         io::File::CloseClearFStream( write);
       }

    //////////////////////
    // input and output //
    //////////////////////

       // test read write functions
       BCL_MessageStd( "testing read write functions");
       WriteBCLObject( completer);
       biol::AABackBoneCompleter completer_read;
       ReadBCLObject( completer_read);
       BCL_Example_Check
       (
         completer_read.GetAddAmideHydrogens() == completer.GetAddAmideHydrogens() &&
           completer_read.GetAddCarbonylOxygens() == completer.GetAddCarbonylOxygens(),
         "The completer read is different!\n" + util::Format()( completer_read) + "\nvs\n" + util::Format()( completer)
       );

    //////////////////////
    // helper functions //
    //////////////////////

      // set the flag for writing hydrogens to true
      pdb::Factory::GetFlagWriteHydrogens()->UnsetFlag();

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAABackBoneCompleter

  const ExampleClass::EnumType ExampleBiolAABackBoneCompleter::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAABackBoneCompleter())
  );

} // namespace bcl
