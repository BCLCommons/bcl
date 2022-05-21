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
#include "assemble/bcl_assemble_printer_protein_model.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_storage_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_printer_protein_model.cpp
  //!
  //! @author weinerbe, fischea
  //! @date Nov 22, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssemblePrinterProteinModel :
    public ExampleInterface
  {

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

    ExampleAssemblePrinterProteinModel *Clone() const
    {
      return new ExampleAssemblePrinterProteinModel( *this);
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
      // get a protein model
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

      // create a storage object
      util::ShPtr< assemble::ProteinStorageFile> sp_storage
      (
        new assemble::ProteinStorageFile
        (
          AddExampleOutputPathToFilename( assemble::PrinterProteinModel(), ""),
          assemble::ProteinStorageFile::e_Overwrite
        )
      );

      // get superimposition
      const quality::SuperimposeMeasure &superimpose( quality::GetSuperimposeMeasures().e_RMSD);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::PrinterProteinModel def_construct;
      BCL_ExampleCheck( def_construct.GetPrefix(), "");

      // construct from members
      const std::string prefix( "printer_");
      assemble::PrinterProteinModel storage_construct( prefix, sp_storage, superimpose);

    /////////////////
    // data access //
    /////////////////

      // test GetPrefix
      BCL_ExampleCheck( storage_construct.GetPrefix(), prefix);

      // test SetPrefix
      def_construct.SetPrefix( prefix);
      BCL_ExampleIndirectCheck( def_construct.GetPrefix(), prefix, "SetPrefix");

    ////////////////
    // operations //
    ////////////////

      // initialize
      storage_construct.Initialize( 13, 2);

      // test print function
      opti::Tracker< assemble::ProteinModel, double> tracker;
      tracker.Track
      (
        util::ShPtr< storage::Pair< assemble::ProteinModel, double> >
        (
          new storage::Pair< assemble::ProteinModel, double>( protein_model, 0.0)
        )
      );
      storage_construct.Print( tracker);

      return 0;
    }

  }; //end ExampleAssemblePrinterProteinModel

  const ExampleClass::EnumType ExampleAssemblePrinterProteinModel::s_Instance
  (
    GetExamples().AddEnum( ExampleAssemblePrinterProteinModel())
  );

} // namespace bcl
