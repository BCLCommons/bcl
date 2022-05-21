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
#include "pdb/bcl_pdb_printer_body_assignment.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_body.h"
#include "restraint/bcl_restraint_contains_body_origin.h"
#include "restraint/bcl_restraint_handler_body.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_printer_body_assignment.cpp
  //!
  //! @author linders, weinerbe
  //! @date Nov 21, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbPrinterBodyAssignment :
    public ExampleInterface
  {
  public:

    ExamplePdbPrinterBodyAssignment *Clone() const
    {
      return new ExamplePdbPrinterBodyAssignment( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create the sse pool
      BCL_MessageStd( "Creating the SSEPool");
      const util::ShPtr< assemble::SSEPool> sp_sse_pool( new assemble::SSEPool( protein_model.GetSSEs()));
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());
      sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sp_sse_pool);
      protein_model.SetProteinModelData( sp_model_data);

      // create body occupancy determiner (method that determines if a body is occupied or not)
      const util::ShPtr
      <
        util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool>
      > occupancy( new restraint::ContainsBodyOrigin());

      // create the restraint handler
      restraint::HandlerBody restraint_handler( occupancy);

      io::IFStream read;
      // open the stream again
      BCL_ExampleMustOpenInputFile( read, pdb_filename);

      // create the restraints
      const util::ShPtrVector< restraint::Body> restraints( restraint_handler.CreateRestraintsBody( read));

      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const pdb::PrinterBodyAssignment def_construct;
      BCL_ExampleIndirectCheck( def_construct( protein_model).IsEmpty(), true, "default construction");

      // restraint constructor
      const pdb::PrinterBodyAssignment restraint_construct( util::CloneToShPtr( restraints));

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const util::ShPtrList< pdb::Line> lines( restraint_construct( protein_model));
      BCL_ExampleCheck( lines.IsEmpty(), false);

      // write out the lines
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( def_construct, "body_assignment_printer.pdb"));
      pdb::Handler handler;
      handler.AppendLines( lines);
      handler.WriteLines( write);
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbPrinterBodyAssignment

  const ExampleClass::EnumType ExamplePdbPrinterBodyAssignment::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbPrinterBodyAssignment())
  );

} // namespace bcl
