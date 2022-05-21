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
#include "assemble/bcl_assemble_printer_protein_model_multimer.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "io/bcl_io_directory_entry.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_printer_protein_model_multimer.cpp
  //!
  //! @author weinerbe
  //! @date May 10, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssemblePrinterProteinModelMultimer :
    public ExampleInterface
  {
  public:

    ExampleAssemblePrinterProteinModelMultimer *Clone() const
    {
      return new ExampleAssemblePrinterProteinModelMultimer( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::PrinterProteinModelMultimer def_construct;
      BCL_ExampleIndirectCheck( def_construct.GetPrefix(), "", "default constructor");

      // test constructor with native model
      const std::string prefix( "1j4n_");
      util::ShPtr< assemble::ProteinStorageFile> sp_storage
      (
        new assemble::ProteinStorageFile
        (
          AddExampleOutputPathToFilename( def_construct, ""),
          assemble::ProteinStorageFile::e_Overwrite
        )
      );
      assemble::PrinterProteinModelMultimer native_construct
      (
        prefix,
        sp_native_model,
        sp_storage,
        quality::GetSuperimposeMeasures().e_RMSD
      );

    /////////////////
    // data access //
    /////////////////

      // test GetPrefix
      BCL_ExampleCheck( native_construct.GetPrefix(), prefix);

    ////////////////
    // operations //
    ////////////////

      // remove any previously generated file since the printer appends
      io::DirectoryEntry( AddExampleOutputPathToFilename( def_construct, "1j4n_0000_final_model.pdb")).Remove();

      // test writing
      native_construct.Initialize( 0);
      opti::Tracker< assemble::ProteinModel, double> tracker;
      tracker.Track
      (
        util::ShPtr< storage::Pair< assemble::ProteinModel, double> >
        (
          new storage::Pair< assemble::ProteinModel, double>( protein_model, 0.0)
        )
      );
      native_construct.Print( tracker);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssemblePrinterProteinModelMultimer

  const ExampleClass::EnumType ExampleAssemblePrinterProteinModelMultimer::s_Instance
  (
    GetExamples().AddEnum( ExampleAssemblePrinterProteinModelMultimer())
  );

} // namespace bcl
