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
#include "descriptor/bcl_descriptor_buser_metric.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_molecule_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "sdf/bcl_sdf_fragment_factory.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molprint_2d.cpp
  //!
  //! @author vuot2
  //! @date Feb 8, 2017
  //! @remarks status complete
  //! @remarks reviewed by
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorBuserMetric : public ExampleInterface
  {
  public:

    ExampleDescriptorBuserMetric *Clone() const
    {
      return new ExampleDescriptorBuserMetric( *this);
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

    //! @brief get the first molecule out of FILENAME, run PROPERTY on it and return the single result
    //! @param FILENAME the filename containing a molecule
    //! @param PROPERTY the property to run
    linal::Vector< float> TestMoleculeWithDescriptor
    (
      const std::string &FILENAME, descriptor::CheminfoProperty &PROPERTY
    ) const
    {
      io::IFStream input_sdf;
      io::File::TryOpenIFStream( input_sdf, FILENAME);
      chemistry::FragmentComplete mol( sdf::FragmentFactory::MakeFragment( input_sdf));
      io::File::CloseClearFStream( input_sdf);
      return PROPERTY->SumOverObject( mol);
    }

    int Run() const
    {
    /////////////////
    // preparation //
    /////////////////

      const std::string cyclohexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      const std::string mix_filename( AddExampleInputPathToFilename
        (
          e_Chemistry, "taxol_cylhexane_cyhexane1_2diol.sdf"
        ));
      const std::string cyclohexaneol_filename
      (
        AddExampleInputPathToFilename( e_Chemistry, "cyclohexane1_2diol.sdf")
      );

      // prepare atom properties vector
      descriptor::CheminfoProperty cyclohexane
      (
        "BuserMetric(atom hashing type=Element, filename=" + cyclohexane_filename + ")"
      );

      // prepare atom properties vector
      descriptor::CheminfoProperty mix
      (
        "BuserMetric(atom hashing type=Element, filename=" + mix_filename + ")"
      );

      // prepare atom properties vector
      descriptor::CheminfoProperty cyclohexaneol
      (
        "BuserMetric(atom hashing type=Element, filename=" + cyclohexaneol_filename + ")"
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test empty constructor
      descriptor::BuserMetric buser;

      BCL_ExampleCheck( buser.GetAlias(), "BuserMetric");

      // check if the default size of the descriptor is 1
      BCL_ExampleCheck( buser.GetNormalSizeOfFeatures(), 1);

      // BuserScore of cyclohexane and itself -> 0.0
      BCL_ExampleCheck
      (
        TestMoleculeWithDescriptor( cyclohexane_filename, cyclohexane)( 0),
        float( 0.0)
      );

      // cyclohexane and cyclohexan1_2diol
      // Check the Buser score of cyclohexane and cyclohexanediol
      BCL_ExampleCheckWithinTolerance
      (
        TestMoleculeWithDescriptor( cyclohexane_filename, cyclohexaneol)( 0),
        0.607625,
        0.00001
      );

      // cyclohexanediol and cyclohexane, cyclohexandiol and taxol
      BCL_ExampleCheckWithinTolerance
      (
        TestMoleculeWithDescriptor( cyclohexaneol_filename, mix)( 0),
        0.261065,
        0.00001
      );
      BCL_ExampleCheckWithinTolerance
      (
        TestMoleculeWithDescriptor( cyclohexaneol_filename, mix)( 1),
        0.607625,
        0.00001
      );
      BCL_ExampleCheckWithinTolerance
      (
        TestMoleculeWithDescriptor( cyclohexaneol_filename, mix)( 2),
        ( 0.607625 + 0.261065) / 2,
        0.00001
      );

      return 0;
    }
    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleDescriptorBuserMetric

  const ExampleClass::EnumType ExampleDescriptorBuserMetric::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorBuserMetric())
  );

} // namespace bcl
