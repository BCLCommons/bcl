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
#include "chemistry/bcl_chemistry_sub_fragment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_sub_fragment.cpp
  //! @details Tests ChemistrySubFragment class which contains small molecule configuration data
  //!
  //! @author kothiwsk
  //! @date
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySubFragment :
    public ExampleInterface
  {
  public:

    ExampleChemistrySubFragment *Clone() const
    {
      return new ExampleChemistrySubFragment( *this);
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

      // read in molecule and construct conformation from conformation interface
      io::IFStream input_sdf;
      const std::string diazepam_filename( AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, diazepam_filename);

      // load information into small_mol_conformation
      chemistry::FragmentComplete small_mol_conformation;
      small_mol_conformation
      = sdf::FragmentFactory::MakeFragment( input_sdf, sdf::e_Saturate);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

      chemistry::SubFragment sub_fragment_a( small_mol_conformation);

      size_t mol_size( small_mol_conformation.GetNumberAtoms());

      // test GetThisToNode
      BCL_ExampleCheck( sub_fragment_a.GetThisToNode().GetSize(), mol_size);

      // test GetThisToParent
      BCL_ExampleCheck( sub_fragment_a.GetThisToParent().GetSize(), mol_size);

      graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> molecule_graph_a
      (
        chemistry::ConformationGraphConverter::CreateGraphWithAtoms( small_mol_conformation)
      );

      storage::Vector< size_t> sub_vector_a( storage::Vector< size_t>::Create( 3, 4, 5, 6, 7, 8, 9, 10, 11, 12));
      chemistry::SubFragment sub_fragment_b( sub_fragment_a, sub_vector_a);

      // test GetThisToNode
      BCL_ExampleCheck( sub_fragment_b.GetThisToNode(), storage::Vector< size_t>::Create( 3, 4, 5, 6, 7, 8, 9, 10, 11, 12));

      // test GetThisToParent
      BCL_ExampleCheck( sub_fragment_b.GetThisToParent(), storage::Vector< size_t>::Create( 3, 4, 5, 6, 7, 8, 9, 10, 11, 12));

      storage::Vector< size_t> sub_vector_b( storage::Vector< size_t>::Create( 3, 4, 5, 6, 7));
      chemistry::SubFragment sub_fragment_c( sub_fragment_b, sub_vector_b);

      // test GetThisToNode
      BCL_ExampleCheck( sub_fragment_c.GetThisToNode(), storage::Vector< size_t>::Create( 6, 7, 8, 9, 10));

      // test GetThisToParent
      BCL_ExampleCheck( sub_fragment_c.GetThisToParent(), storage::Vector< size_t>::Create( 3, 4, 5, 6, 7));

      storage::Vector< size_t> sub_vector_c( storage::Vector< size_t>::Create( 2, 3, 4));
      chemistry::SubFragment sub_fragment_d( sub_fragment_c, sub_vector_c);

      // test GetThisToNode
      BCL_ExampleCheck( sub_fragment_d.GetThisToNode(), storage::Vector< size_t>::Create( 8, 9, 10));

      // test GetThisToParent
      BCL_ExampleCheck( sub_fragment_d.GetThisToParent(), storage::Vector< size_t>::Create( 2, 3, 4));

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistrySubFragment::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySubFragment())
  );

} // namespace bcl
