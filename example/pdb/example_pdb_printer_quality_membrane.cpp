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
#include "pdb/bcl_pdb_printer_quality_membrane.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_printer_quality_membrane.cpp
  //!
  //! @author weinerbe
  //! @date Jun 1, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbPrinterQualityMembrane :
    public ExampleInterface
  {
  public:

    ExamplePdbPrinterQualityMembrane *Clone() const
    {
      return new ExamplePdbPrinterQualityMembrane( *this);
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
      // store the native
      storage::Map< biol::SSType, size_t> min_sse_size;
      min_sse_size[ biol::GetSSTypes().HELIX] = 5;
      util::ShPtr< assemble::ProteinModel> sp_native_model
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"),
          biol::GetAAClasses().e_AABackBone,
          min_sse_size
        ).Clone()
      );

      // hardcopy the model
      util::ShPtr< assemble::ProteinModel> sp_protein_model( sp_native_model->HardCopy());

      // add a membrane
      util::ShPtr< biol::Membrane> sp_membrane( new biol::Membrane());
      util::ShPtr< assemble::ProteinModelData> sp_data( sp_protein_model->GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Membrane, sp_membrane);

      // translate a helix
      const assemble::LocatorSSE helix_41_62_loc( 'A', 41, 62);
      util::ShPtr< assemble::SSE> sp_helix_41_62( helix_41_62_loc.Locate( *sp_protein_model)->HardCopy());
      sp_helix_41_62->Translate( linal::Vector3D( 0.0, 0.0, 1.0));
      sp_protein_model->Replace( sp_helix_41_62);

      // initialize measures and environments
      const storage::Set< quality::Measure> qualities( quality::GetMeasures().e_RMSD);
      const storage::Set< biol::EnvironmentType> environments( biol::GetEnvironmentTypes().e_MembraneCore);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const pdb::PrinterQualityMembrane def_construct;
      BCL_ExampleIndirectCheck( def_construct( *sp_protein_model).IsEmpty(), true, "default construction");

      // test constructor from members
      const pdb::PrinterQualityMembrane rmsd_construct( qualities, environments, sp_native_model);

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const util::ShPtrList< pdb::Line> lines( rmsd_construct( *sp_protein_model));
      BCL_ExampleCheck( lines.IsEmpty(), false);

      // write out the lines
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( def_construct, "membrane_quality_printer.pdb"));
      pdb::Handler handler;
      handler.AppendLines( lines);
      handler.WriteLines( write);
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbPrinterQualityMembrane

  const ExampleClass::EnumType ExamplePdbPrinterQualityMembrane::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbPrinterQualityMembrane())
  );

} // namespace bcl
