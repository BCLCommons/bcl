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
#include "pdb/bcl_pdb_printer_quality_docking.h"

// includes from bcl - sorted alphabetically
#include "example_interface.h"
#include "example_proteins.h"
#include "biol/bcl_biol_aa_classes.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "pdb/bcl_pdb_handler.h"
#include "quality/bcl_quality_measures.h"

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_printer_quality_docking.cpp
  //!
  //! @author lib14
  //! @date Sept 25, 2017
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePrinterQualityDocking :
    public ExampleInterface
  {

  public:

    //! @brief Clone function
    //! @return new Pointer to a copy of the actual object behind the pointer
    ExamplePrinterQualityDocking *Clone() const
    {
      return new ExamplePrinterQualityDocking( *this);
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

      // get a model for a protein complex
      const std::string native_filename( AddExampleInputPathToFilename( e_Biology, "4o6y.pdb"));
      const std::string docked_filename( AddExampleInputPathToFilename( e_Biology, "4o6y_docked.pdb"));
      assemble::ProteinModel docked_model
      (
        Proteins::GetModel( docked_filename, biol::GetAAClasses().e_AAComplete)
      );
      assemble::ProteinModel native_model
      (
        Proteins::GetModel( native_filename, biol::GetAAClasses().e_AAComplete)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const pdb::PrinterQualityDocking def_quality_printer_docking;

      // test constructor from quality measures and chain IDs
      const pdb::PrinterQualityDocking quality_printer_docking
      (
        storage::Set< quality::Measure>( quality::GetMeasures().e_RMSD_NoSuperimposition),
        "B"
      );

    ///////////////
    // operators //
    ///////////////

      // test () operator
      util::ShPtrList< pdb::Line> quality_lines_no_native( quality_printer_docking( docked_model));
      BCL_ExampleCheck( quality_lines_no_native.IsEmpty(), false);

      // write out the lines
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile
      (
        write,
        AddExampleOutputPathToFilename( quality_printer_docking, "printer_quality_docking_no_native.pdb")
      );
      pdb::Handler handler;
      handler.AppendLines( quality_lines_no_native);
      handler.WriteLines( write);
      io::File::CloseClearFStream( write);

      // set the native model
      util::ShPtr< assemble::ProteinModelData> sp_data
      (
        docked_model.GetProteinModelData()
      );
      sp_data->Insert
      (
        assemble::ProteinModelData::e_NativeModel,
        util::ShPtr< assemble::ProteinModel>( native_model.Clone())
      );
      docked_model.SetProteinModelData( sp_data);

      // write out lines after adding native model
      util::ShPtrList< pdb::Line> quality_lines_with_native( quality_printer_docking( docked_model));
      BCL_ExampleCheck( quality_lines_with_native.IsEmpty(), false);
      BCL_ExampleMustOpenOutputFile
      (
        write,
        AddExampleOutputPathToFilename( quality_printer_docking, "printer_quality_docking_with_native.pdb")
      );
      handler.AppendLines( quality_lines_with_native);
      handler.WriteLines( write);
      io::File::CloseClearFStream( write);

      return 0;

    }

    // single instance of this class
    static const ExampleClass::EnumType s_Instance;

  }; // end of ExamplePrinterQualityDocking

  const ExampleClass::EnumType ExamplePrinterQualityDocking::s_Instance
  (
    GetExamples().AddEnum( ExamplePrinterQualityDocking())
  );
} // namespace bcl

