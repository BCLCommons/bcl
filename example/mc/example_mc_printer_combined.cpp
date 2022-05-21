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
#include "mc/bcl_mc_printer_combined.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_printer_protein_model.h"
#include "assemble/bcl_assemble_printer_protein_model_multimer.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_printer_combined.cpp
  //!
  //! @author linders, fischea
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMcPrinterCombined :
    public ExampleInterface
    {
  public:

    ExampleMcPrinterCombined *Clone() const
    {
      return new ExampleMcPrinterCombined( *this);
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

    /////////////////
    // preparation //
    /////////////////

      // store the native multimer
      util::ShPtr< assemble::ProteinModel> sp_native_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1j4n_multimer.pdb")).Clone()
      );

      // get a protein model
      assemble::ProteinModel protein_model( sp_native_model->GetChain( 'A').HardCopy());

      // create the multiplier
      util::ShPtr< assemble::ProteinModelData> sp_protein_model_data( new assemble::ProteinModelData());
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        new assemble::ProteinModelMultiplier( coord::GetAxes().e_Z, 4, protein_model)
      );
      sp_protein_model_data->Insert( assemble::ProteinModelData::e_Multiplier, sp_multiplier);
      sp_protein_model_data->Insert
      (
        assemble::ProteinModelData::e_Identification,
        util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( "1J4N_mult"))
      );
      protein_model.SetProteinModelData( sp_protein_model_data);

      // move a helix
      util::ShPtr< assemble::SSE> sp_helix( protein_model.GetSSEs( biol::GetSSTypes().HELIX).FirstElement()->Clone());
      sp_helix->Translate( linal::Vector3D( 0.5, 0.5, 0.5));
      protein_model.Replace( sp_helix);

      const quality::SuperimposeMeasure &rmsd( quality::GetSuperimposeMeasures().e_RMSD);

      // create the score printer
      util::ShPtr< assemble::ProteinStorageFile> sp_storage
      (
        new assemble::ProteinStorageFile
        (
          AddExampleOutputPathToFilename( mc::PrinterCombined< assemble::ProteinModel, double>(), ""),
          assemble::ProteinStorageFile::e_Overwrite
        )
      );
      const std::string prefix( "dummy");
      const util::ShPtr< assemble::PrinterProteinModel> model_printer
      (
        new assemble::PrinterProteinModel( prefix, sp_storage, rmsd)
      );

      // create multimer printer
      util::ShPtr< assemble::PrinterProteinModelMultimer> multimer_printer
      (
        new assemble::PrinterProteinModelMultimer
        (
          prefix,
          sp_native_model,
          sp_storage,
          rmsd
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      mc::PrinterCombined< assemble::ProteinModel, double> printer_combined_default;

      BCL_Example_Check
      (
        printer_combined_default.GetPrinters().GetSize() == 0,
        "printer_combined_default should not contain any printers!"
      );

      // test copy constructor
      mc::PrinterCombined< assemble::ProteinModel, double> printer_combined( printer_combined_default);

      // test clone function
      util::ShPtr< mc::PrinterCombined< assemble::ProteinModel, double> > printer_combined_clone
      (
        printer_combined.Clone()
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::mc::PrinterCombined<bcl::assemble::ProteinModel,double>");
      const std::string static_class_name( GetStaticClassName< mc::PrinterCombined< assemble::ProteinModel, double> >());

      BCL_Example_Check
      (
        static_class_name == correct_static_class_name,
        "GetStaticClassName gives " + static_class_name + " but should give " + correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
       (
        static_class_name == printer_combined_clone->GetClassIdentifier(),
        "GetClassIdentifier gives " + printer_combined_clone->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ////////////////
    // operations //
    ////////////////

      // insert score printer into the printer combined
      printer_combined.Insert( model_printer);

      BCL_Example_Check
      (
        printer_combined.GetPrinters().GetSize() == 1,
        "printer_combined_default should contain exactly one printer, but contains " +
          util::Format()( printer_combined.GetPrinters().GetSize()) + " printers!"
      );

      printer_combined.Insert( multimer_printer);

      BCL_Example_Check
      (
        printer_combined.GetPrinters().GetSize() == 2,
        "printer_combined_default should contain exactly two printers, but contains " +
          util::Format()( printer_combined.GetPrinters().GetSize()) + " printers!"
      );

      // create a tracker
      opti::Tracker< assemble::ProteinModel, double> tracker;
      tracker.Track
      (
        util::ShPtr< storage::Pair< assemble::ProteinModel, double> >
        (
          new storage::Pair< assemble::ProteinModel, double>( protein_model, 2.25)
        )
      );

      // test PrintEnd function
      printer_combined.Print( tracker);

    //////////////////////
    // input and output //
    //////////////////////

      //implicitly test Write by using <<
      WriteBCLObject( printer_combined);

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

    }; //end ExampleMcPrinterCombined

    const ExampleClass::EnumType ExampleMcPrinterCombined::s_Instance
    (
      GetExamples().AddEnum( ExampleMcPrinterCombined())
    );

} // namespace bcl
