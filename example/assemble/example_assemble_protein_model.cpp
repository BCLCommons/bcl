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
#include "assemble/bcl_assemble_protein_model.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_model.cpp
  //!
  //! @author woetzen, fischea
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinModel :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinModel *Clone() const
    { return new ExampleAssembleProteinModel( *this);}

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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      //read pdbfile
      io::IFStream read;
      BCL_MessageStd( "read pdb file: " + pdb_filename);
      BCL_ExampleMustOpenInputFile( read, pdb_filename);

      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      pdb::Factory factory( biol::GetAAClasses().e_AABackBone);

      //instantiate sequences
      BCL_MessageStd( "building sequences from pdb chains");
      util::ShPtrVector< biol::AASequence> sequences( factory.AASequencesFromPDB( pdb));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      assemble::ProteinModel protein_model( factory.ProteinModelFromPDB( pdb));

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( protein_model.GetClassIdentifier(), ( GetStaticClassName< assemble::ProteinModel>()));

      // check access to SSEs
      const util::SiPtrVector< const biol::AABase> residues( protein_model.GetAminoAcids());
      const biol::AABase &aa_first( *residues( 6));
      const biol::AABase &aa_second( *residues( 363));
      const assemble::SSE &sse_first( *protein_model.GetSSE( aa_first));
      const assemble::SSE &sse_second( *protein_model.GetSSE( aa_second));
      BCL_ExampleCheck( sse_first.GetIdentification(), protein_model.GetSSEs()( 0)->GetIdentification());
      BCL_ExampleCheck( sse_second.GetIdentification(), protein_model.GetSSEs()( 50)->GetIdentification());

    ////////////////
    // operations //
    ////////////////

      // set all secondary structure element types in protein_model to ideal conformation
      protein_model.SetToIdealConformation();

      // write idealized protein models to an example pdb
      BCL_MessageStd( "write ideal_model2.pdb");
      Proteins::WriteModelToPDB( protein_model, AddExampleOutputPathToFilename( protein_model, "ideal_model2.pdb"));
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel chopped_model( factory.ProteinModelFromPDB( pdb, ssetype_min_size));

      //chop all secondary structure element types of chopped_model to a size of 5 / 9 AAs for strand / helix
      protein_model.ChopSSEs( storage::VectorND< 3, size_t>( 9, 5, 5));

      // write chopped protein models to an example pdb
      BCL_MessageStd( "write chopped_model.pdb");
      Proteins::WriteModelToPDB( chopped_model, AddExampleOutputPathToFilename( chopped_model, "chopped_model.pdb"));

      // test AddLoops function
      // create protein model without loops
      assemble::ProteinModel add_loops_model( factory.ProteinModelFromPDB( pdb, ssetype_min_size));

      // loops with zero coordinates are added
      add_loops_model.AddLoops( true, false);
      // write "add_loops_model" to pdb
      BCL_MessageStd( "write add_loops_model.pdb");
      Proteins::WriteModelToPDB
      (
        add_loops_model, AddExampleOutputPathToFilename( add_loops_model, "add_loops_model.pdb")
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinModel

  const ExampleClass::EnumType ExampleAssembleProteinModel::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinModel())
  );

} // namespace bcl
