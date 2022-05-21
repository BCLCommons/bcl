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
#include "pdb/bcl_pdb_printer_loop_closure.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "fold/bcl_fold_handler_locator_loop_domain.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_protocol_loop_close.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_printer_loop_closure.cpp
  //!
  //! @author weinerbe
  //! @date Dec 12, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbPrinterLoopClosure :
    public ExampleInterface
  {
  public:

    ExamplePdbPrinterLoopClosure *Clone() const
    {
      return new ExamplePdbPrinterLoopClosure( *this);
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
      // get a protein model with missing loops
      assemble::ProteinModel model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

      // add hydrogens, which are required for the rms
      model = *fold::ProtocolLoopClose::AddNitrogenHydrogens( model);

      // create a handler
      fold::HandlerLocatorLoopDomain handler( false);

      // get the loop domain locators created from "model"
      util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > sp_locators
      (
        handler.CreateLocatorLoopDomainsForInteriorCoil( model).Clone()
      );

      // construct model data
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());
      sp_model_data->Insert( assemble::ProteinModelData::e_LoopDomainLocators, sp_locators);
      model.SetProteinModelData( sp_model_data);

      // the threshold for loop closure
      const double closure_threshold( 0.8);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      pdb::PrinterLoopClosure def_construct;
      BCL_ExampleIndirectCheck( def_construct( model).IsEmpty(), true, "default construction");

      // test constructor from threshold
      const pdb::PrinterLoopClosure theshold_construct( closure_threshold);

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const util::ShPtrList< pdb::Line> lines( theshold_construct( model));
      BCL_ExampleCheck( lines.IsEmpty(), false);

      // write out the lines
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( def_construct, "loop_printer.pdb"));
      pdb::Handler pdb_handler;
      pdb_handler.AppendLines( lines);
      pdb_handler.WriteLines( write);
      io::File::CloseClearFStream( write);

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( theshold_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleIndirectCheck( def_construct( model).IsEmpty(), false, "read/write");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbPrinterLoopClosure

  const ExampleClass::EnumType ExamplePdbPrinterLoopClosure::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbPrinterLoopClosure())
  );

} // namespace bcl
