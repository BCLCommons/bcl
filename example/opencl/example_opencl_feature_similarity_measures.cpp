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

// include example header,
#include "example.h"
// include the header of the class which this example is for
#include "opencl/bcl_opencl_feature_similarity_measures.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_feature_similarity_measures.cpp
  //!
  //! @author loweew
  //! @date Oct 3, 2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclFeatureSimilarityMeasures :
    public ExampleInterface
  {
  public:

    ExampleOpenclFeatureSimilarityMeasures *Clone() const
    {
      return new ExampleOpenclFeatureSimilarityMeasures( *this);
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
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

      // constructing matrices of feature vectors
      float data_a[ 9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
      linal::Matrix< float> data_mat_a( 3, 3, data_a);

      float exp_cos [ 9] =
      {
          1          , 0.974631846, 0.959411946,
          0.974631846, 1          , 0.998190893,
          0.959411946, 0.998190893, 1
      };

      float exp_tan [ 9] =
      {
          1          , 0.542372881, 0.316455696,
          0.542372881, 1          , 0.818791946,
          0.316455696, 0.818791946, 1
      };

      float exp_dice[ 9] =
      {
          1          , 0.703296703, 0.480769231,
          0.703296703, 1          , 0.900369004,
          0.480769231, 0.900369004, 1
      };

      float exp_euc [ 9] =
      {
          1          , 0.161390478, 0.08777855 ,
          0.161390478, 1          , 0.161390478,
          0.08777855 , 0.161390478, 1
      };

      float exp_man [ 9] =
      {
          1          , 0.1, 0.052631579,
          0.1        , 1  , 0.1        ,
          0.052631579, 0.1, 1
      };

      linal::Matrix< float> exp_mat_tan ( 3, 3, exp_tan);
      linal::Matrix< float> exp_mat_cos ( 3, 3, exp_cos);
      linal::Matrix< float> exp_mat_dice( 3, 3, exp_dice);
      linal::Matrix< float> exp_mat_euc ( 3, 3, exp_euc);
      linal::Matrix< float> exp_mat_man ( 3, 3, exp_man);
      std::string s_cos( "Cosine"), s_tan( "Tanimoto"), s_dice( "Dice"), s_euc( "Euclidean"), s_man( "Manhattan");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor
      opencl::FeatureSimilarityMeasures< float> fsm     ( s_cos, opencl::GetTools().GetFirstCommandQueue());
      opencl::FeatureSimilarityMeasures< float> fsm_cos ( s_cos, opencl::GetTools().GetFirstCommandQueue());
      opencl::FeatureSimilarityMeasures< float> fsm_tan ( s_tan, opencl::GetTools().GetFirstCommandQueue());
      opencl::FeatureSimilarityMeasures< float> fsm_dice( s_dice, opencl::GetTools().GetFirstCommandQueue());
      opencl::FeatureSimilarityMeasures< float> fsm_euc ( s_euc, opencl::GetTools().GetFirstCommandQueue());
      opencl::FeatureSimilarityMeasures< float> fsm_man ( s_man, opencl::GetTools().GetFirstCommandQueue());

      // clone
      util::ShPtr< opencl::FeatureSimilarityMeasures< float> > sp_fsm( fsm.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      linal::Matrix< float> cos ( fsm_cos( data_mat_a));
      linal::Matrix< float> tan ( fsm_tan( data_mat_a));
      linal::Matrix< float> dice( fsm_dice( data_mat_a));
      linal::Matrix< float> euc ( fsm_euc( data_mat_a));
      linal::Matrix< float> man ( fsm_man( data_mat_a));

    ////////////////
    // operations //
    ////////////////

      BCL_MessageDbg( "cosine: " + util::Format()( cos));
      BCL_MessageDbg( "tanimoto: " + util::Format()( tan));
      BCL_MessageDbg( "dice: " + util::Format()( dice));
      BCL_MessageDbg( "euclidean: " + util::Format()( euc));
      BCL_MessageDbg( "manhattan: " + util::Format()( man));

      BCL_ExampleCheckWithinTolerance( cos,  exp_mat_cos , 0.001);
      BCL_ExampleCheckWithinTolerance( tan,  exp_mat_tan , 0.001);
      BCL_ExampleCheckWithinTolerance( dice, exp_mat_dice, 0.001);
      BCL_ExampleCheckWithinTolerance( euc,  exp_mat_euc , 0.001);
      BCL_ExampleCheckWithinTolerance( man,  exp_mat_man , 0.001);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclFeatureSimilarityMeasures

  const ExampleClass::EnumType ExampleOpenclFeatureSimilarityMeasures::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclFeatureSimilarityMeasures())
  );

} // namespace bcl
