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

//// (c) Copyright BCL @ Vanderbilt University 2014
//// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
//// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
//// (c) This file is part of the BCL software suite and is made available under license.
//// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
//// (c) For academic and non-profit users:
//// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
//// (c) For commercial users:
//// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
//// (c) For BCL developers at Vanderbilt University:
//// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
//// (c)
//// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the
//// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
//// (c) Vanderbilt University.
//// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
//// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org
//// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)
//
//// include example header
//#include "example.h"
//
//// include the header of the class which this example is for
//#include "descriptor/bcl_descriptor_central_2da_sign.h"
//
//// includes from bcl - sorted alphabetically
//#include "chemistry/bcl_chemistry_fragment_ensemble.h"
//#include "io/bcl_io_file.h"
//#include "sdf/bcl_sdf_fragment_factory.h"
//
//// external includes - sorted alphabetically
//
//namespace bcl
//{
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  //!
//  //! @example example_descriptor_central_2da_sign.cpp
//  //!
//  //! @author brownbp1, mendenjl
//  //! @date Jun 21, 2020
//  //! @remarks status complete
//  //! @remarks reviewed by nobody on
//  //!
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  class ExampleDescriptorCentral2DASign :
//    public ExampleInterface
//  {
//  public:
//
//    ExampleDescriptorCentral2DASign *Clone() const
//    {
//      return new ExampleDescriptorCentral2DASign( *this);
//    }
//
//  /////////////////
//  // data access //
//  /////////////////
//
//    //! @brief returns class name
//    //! @return the class name as const ref std::string
//    const std::string &GetClassIdentifier() const
//    {
//      return GetStaticClassName( *this);
//    }
//
//    int Run() const
//    {
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // test default constructor
//      descriptor::Central2DASign stereocenters_a;
//
//      // test Clone
//      descriptor::Central2DASign stereocenters_cp( stereocenters_a);
//
//    ////////////////////////////////////////
//    //  read Small molecule Stereocenters //
//    ////////////////////////////////////////
//
//      // use SmallMoleculeFactory to load molecule information into SmallMolecule
//      // this requires SmallMolecule and input stream
//      io::IFStream input_sdf;
//      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "corina_taxol.sdf"));
//      chemistry::FragmentComplete taxol = sdf::FragmentFactory::MakeFragment( input_sdf);
//      io::File::CloseClearFStream( input_sdf);
//
//      //calculate stereocenters for taxol
//      linal::Vector< float> calculated_stereocenters( stereocenters_a.CollectValuesOnEachElementOfObject( taxol));
//
//      float expected_vec_array [] =
//      {
//          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,-1,1,-1,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,
//          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//        0,0,0,0,0
//      };
//      linal::Vector< float> expected_vector( 113, expected_vec_array);
//      BCL_ExampleCheckWithinAbsTolerance( calculated_stereocenters, expected_vector, 0.001);
//
//      return 0;
//      } // Run
//
//    static const ExampleClass::EnumType s_Instance;
//
//  }; //end ExampleDescriptorCentral2DASign
//
//  const ExampleClass::EnumType ExampleDescriptorCentral2DASign::s_Instance
//  (
//    GetExamples().AddEnum( ExampleDescriptorCentral2DASign())
//  );
//
//} // namespace bcl
