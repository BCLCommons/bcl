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
#include "pdb/bcl_pdb_printer_biomatrix.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_printer_biomatrix.cpp
  //!
  //! @author weinerbe
  //! @date Feb 17, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbPrinterBiomatrix :
    public ExampleInterface
  {
  public:

    ExamplePdbPrinterBiomatrix *Clone() const
    {
      return new ExamplePdbPrinterBiomatrix( *this);
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
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1j4n_multimer.pdb")).GetChain( 'A')
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const pdb::PrinterBiomatrix def_construct;
      BCL_ExampleIndirectCheck( def_construct( protein_model).IsEmpty(), true, "default construction");

    ///////////////
    // operators //
    ///////////////

      // create the multiplier
      util::ShPtr< assemble::ProteinModelData> sp_protein_model_data( new assemble::ProteinModelData());
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        new assemble::ProteinModelMultiplier( coord::GetAxes().e_Z, 4, protein_model)
      );
      sp_protein_model_data->Insert( assemble::ProteinModelData::e_Multiplier, sp_multiplier);
      protein_model.SetProteinModelData( sp_protein_model_data);

      // test () operator
      const util::ShPtrList< pdb::Line> lines( def_construct( protein_model));
      BCL_ExampleCheck( lines.IsEmpty(), false);

      // write out the lines
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( def_construct, "biomatrix_printer.pdb"));
      pdb::Handler handler;
      handler.AppendLines( lines);
      write << handler;
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbPrinterBiomatrix

  const ExampleClass::EnumType ExamplePdbPrinterBiomatrix::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbPrinterBiomatrix())
  );

} // namespace bcl
