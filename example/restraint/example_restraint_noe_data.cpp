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
#include "restraint/bcl_restraint_noe_data.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_file.h"
#include "nmr/bcl_nmr_star_noe_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_noe_data.cpp
  //!
  //! @author weinerbe
  //! @date Jan 16, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintNOEData :
    public ExampleInterface
  {
  public:

    ExampleRestraintNOEData *Clone() const
    {
      return new ExampleRestraintNOEData( *this);
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
      const std::string extension( ".noe_star");
      restraint::NOEData data_construct( ( nmr::StarNOEHandler()));

      // read in restraints
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "nmr_star_restraints.txt"));
      data_construct.ReadRestraints( read);
      io::File::CloseClearFStream( read);

    /////////////////
    // data access //
    /////////////////

      // test GetDataConstruct
      const size_t nr_restraints( 5);
      BCL_ExampleCheck( data_construct.GetAtomDistanceRestraints().GetSize(), nr_restraints);

    ////////////////
    // operations //
    ////////////////

      // test ModifyScoreWeightSet
      fold::ScoreWeightSet weights;
      const double score_weight( 5);
      data_construct.InitializeScores();
      data_construct.ModifyScoreWeightSet( weights);
      BCL_ExampleCheck( weights.GetWeight( restraint::NOEData::e_ScoreNOERestraint), score_weight);

    //////////////////////
    // input and output //
    //////////////////////

      // test WriteRestraints
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( nmr::StarNOEHandler(), "nmr_star_noe.txt"));
      data_construct.GetDefaultHandler().WriteRestraints( write, data_construct.GetAtomDistanceRestraints());
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintNOEData

  const ExampleClass::EnumType ExampleRestraintNOEData::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintNOEData())
  );

} // namespace bcl
