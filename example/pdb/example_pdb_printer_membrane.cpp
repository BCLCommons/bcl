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
#include "pdb/bcl_pdb_printer_membrane.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_printer_membrane.cpp
  //!
  //! @author weinerbe
  //! @date Jan 13, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbPrinterMembrane :
    public ExampleInterface
  {
  public:

    ExamplePdbPrinterMembrane *Clone() const
    {
      return new ExamplePdbPrinterMembrane( *this);
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
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"))
      );

      // add a membrane
      const util::ShPtr< biol::Membrane> sp_membrane( new biol::Membrane());
      util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Membrane, sp_membrane);
      protein_model.SetProteinModelData( sp_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const pdb::PrinterMembrane def_construct;

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const util::ShPtrList< pdb::Line> lines( def_construct( protein_model));
      BCL_ExampleCheck( lines.IsEmpty(), false);

      // write out the lines
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( def_construct, "membrane_printer.pdb"));
      pdb::Handler pdb_handler;
      pdb_handler.AppendLines( lines);
      pdb_handler.WriteLines( write);
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbPrinterMembrane

  const ExampleClass::EnumType ExamplePdbPrinterMembrane::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbPrinterMembrane())
  );

} // namespace bcl
