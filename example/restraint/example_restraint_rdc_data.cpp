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
#include "restraint/bcl_restraint_rdc_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_file.h"
#include "nmr/bcl_nmr_star_rdc_handler.h"
#include "score/bcl_score_restraint_residual_dipolar_coupling.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_rdc_data.cpp
  //!
  //! @author weinerbe
  //! @date Aug 9, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintRDCData :
    public ExampleInterface
  {
  public:

    ExampleRestraintRDCData *Clone() const
    {
      return new ExampleRestraintRDCData( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor from extension and handler
      restraint::RDCData data_construct( ( nmr::StarRDCHandler( ".txt")));

      // test ModifyScoreWeightSet
      fold::ScoreWeightSet weights;
      const std::string &score_scheme( score::RestraintResidualDipolarCoupling::GetDefaultScheme());
      const double score_weight( 5);
      restraint::GetFlagRestraintsFilePrefix()->ReadFromList
      (
        storage::Vector< std::string>::Create( AddExampleInputPathToFilename( e_Biology, "nmr_star_restraints")),
        util::GetLogger()
      );
      data_construct.InitializeScores();
      restraint::GetFlagRestraintsFilePrefix()->ResetFlag();
      data_construct.ModifyScoreWeightSet( weights);
      BCL_ExampleCheck( weights.GetWeight( fold::GetScores().GetEnumFromName( score_scheme)), score_weight);

    /////////////////
    // data access //
    /////////////////

      // test GetDataConstruct
      const size_t nr_restraints( 5);
      BCL_ExampleCheck( data_construct.GetRDC().GetData().GetSize(), nr_restraints);

      // test GetAllFlags
      const size_t nr_flags( 4);
      BCL_ExampleCheck( data_construct.GetLabel().GetNumberArguments(), nr_flags);

    ////////////////
    // operations //
    ////////////////

      // ModifyStartModel not implemented

      // ModifyPrinter not implemented

      // ModifyMutateTree not implemented

    //////////////////////
    // input and output //
    //////////////////////

      // ReadRestraints tested above

      // test WriteRestraints
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( data_construct, "nmr_star_angles.txt"));
      data_construct.GetDefaultHandler().WriteRestraints( write, data_construct.GetRDC());
      io::File::CloseClearFStream( write);

      // test read and write
      restraint::RDCData def_construct;
      WriteBCLObject( data_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleIndirectCheck
      (
        def_construct.GetRDC().GetData().GetSize() == data_construct.GetRDC().GetData().GetSize(),
        true,
        "read and write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintRDCData

  const ExampleClass::EnumType ExampleRestraintRDCData::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintRDCData())
  );

} // namespace bcl
