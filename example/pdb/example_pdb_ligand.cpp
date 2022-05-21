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
#include "pdb/bcl_pdb_ligand.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_ligand.cpp
  //!
  //! @author woetzen
  //! @date Feb 1, 2012
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbLigand :
    public ExampleInterface
  {
  public:

    ExamplePdbLigand *Clone() const
    { return new ExamplePdbLigand( *this);}

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
      const std::string pdb_str( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, pdb_str);
      pdb::Handler pdb_reader( read);
      io::File::CloseClearFStream( read);

    ////////////////
    // operations //
    ////////////////

      // get all ligands from the pdb
      const util::ShPtrList< pdb::Ligand> ligands( pdb_reader.GetLigands());

      BCL_MessageStd( "ligands found\n" + util::Format()( ligands));
      BCL_MessageStd( "het full name\n" + util::Format()( pdb_reader.GetHead().GetHetFullname()));
      BCL_MessageStd( "het formula\n" + util::Format()( pdb_reader.GetHead().GetHetFormula()));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbLigand

  const ExampleClass::EnumType ExamplePdbLigand::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbLigand())
  );

} // namespace bcl
