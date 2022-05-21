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
#include "assemble/bcl_assemble_protein_model_with_cache.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_model_with_cache.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinModelWithCache :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinModelWithCache *Clone() const
    { return new ExampleAssembleProteinModelWithCache( *this);}

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

      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      assemble::ProteinModelWithCache protein_model( factory.ProteinModelFromPDB( pdb), false);

      // check that the size is correct
      BCL_ExampleCheck( protein_model.GetSize(), 710);

      // check that the returned iterator begins at the first AA
      BCL_ExampleCheck( protein_model.GetIterator()->GetPdbID(), 1);

      // check that the returned iterator ends at the last AA
      BCL_ExampleCheck( ( --protein_model.GetIterator().GotoEnd())->GetPdbID(), 355);
      BCL_ExampleCheck( protein_model.GetIterator().GetSize(), 710);

      // try building a model that has defined coordinates
      assemble::ProteinModelWithCache protein_model_with_coords( factory.ProteinModelFromPDB( pdb), true);

      // Requiring coordinates should not change the sequence; it should only force removal of any sspred::Methods that
      // rely on 3D coordinates for residues that lack coordinates
      BCL_ExampleCheck( protein_model_with_coords.GetSize(), 710);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinModelWithCache

  const ExampleClass::EnumType ExampleAssembleProteinModelWithCache::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinModelWithCache())
  );

} // namespace bcl
