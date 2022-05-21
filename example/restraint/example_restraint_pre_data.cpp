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
#include "restraint/bcl_restraint_pre_data.h"

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
  //! @example example_restraint_pre_data.cpp
  //!
  //! @author weinerbe
  //! @date Jan 16, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintPREData :
    public ExampleInterface
  {
  public:

    ExampleRestraintPREData *Clone() const
    {
      return new ExampleRestraintPREData( *this);
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
      restraint::PREData data_construct( ( nmr::StarNOEHandler( ".txt")));

      // test ModifyScoreWeightSet
      fold::ScoreWeightSet weights;
      const double score_weight( 5);
      restraint::GetFlagRestraintsFilePrefix()->ReadFromList
      (
        storage::Vector< std::string>::Create( AddExampleInputPathToFilename( e_Biology, "nmr_star_restraints")),
        util::GetLogger()
      );
      data_construct.InitializeScores();
      restraint::GetFlagRestraintsFilePrefix()->ResetFlag();
      data_construct.ModifyScoreWeightSet( weights);
      BCL_ExampleCheck( weights.GetWeight( restraint::PREData::e_ScorePRERestraint), score_weight);

    /////////////////
    // data access //
    /////////////////

      // test GetDataConstruct
      const size_t nr_restraints( 5);
      BCL_ExampleCheck( data_construct.GetAtomDistanceRestraints().GetSize(), nr_restraints);

      // test GetAllFlags
      const size_t n_flags( data_construct.GetLabel().GetNumberArguments());
      BCL_ExampleCheck( n_flags, 2);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test WriteRestraints
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( nmr::StarNOEHandler(), "nmr_star_pre.txt"));
      data_construct.GetDefaultHandler().WriteRestraints( write, data_construct.GetAtomDistanceRestraints());
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintPREData

  const ExampleClass::EnumType ExampleRestraintPREData::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintPREData())
  );

} // namespace bcl
