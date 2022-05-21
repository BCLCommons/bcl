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
#include "descriptor/bcl_descriptor_structure_count.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_vcharge.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_structure_count.cpp
  //!
  //! @author geanesar
  //! @date Jul 31, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorStructureCount :
    public ExampleInterface
  {
  public:

    ExampleDescriptorStructureCount *Clone() const
    {
      return new ExampleDescriptorStructureCount( *this);
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
    float TestMoleculeWithDescriptor( const std::string &FILENAME, descriptor::CheminfoProperty &PROPERTY) const
    {
      io::IFStream input_sdf;
      BCL_ExampleMustOpenInputFile( input_sdf, FILENAME);
      chemistry::FragmentComplete mol( sdf::FragmentFactory::MakeFragment( input_sdf));
      io::File::CloseClearFStream( input_sdf);
      return PROPERTY->SumOverObject( mol)( 0);
    }

    int Run() const
    {
    /////////////////
    // preparation //
    /////////////////

      const std::string cyclohexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      const std::string benzene_filename( AddExampleInputPathToFilename( e_Chemistry, "benzene_prepared.sdf"));
      const std::string taxol_filename( AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      const std::string diazepam_filename( AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf"));

      // Read in small molecule
      io::IFStream input_sdf;
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));

      // load information into small_mol
      chemistry::FragmentComplete taxol( sdf::FragmentFactory::MakeFragment( input_sdf));

      io::File::CloseClearFStream( input_sdf);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // prepare atom properties vector
      descriptor::CheminfoProperty has_six_membered_carbon_ring
      (
        "StructureCount(ignore h=True,atom comparison=ElementType,bond comparison=Identity,"
        "filename=" + cyclohexane_filename + ")"
      );
      descriptor::CheminfoProperty has_cyclohexane
      (
        "StructureCount(ignore h=True,atom comparison=AtomType,bond comparison=BondOrderOrAromaticWithRingness,"
        "filename=" + cyclohexane_filename + ")"
      );
      descriptor::CheminfoProperty has_benzene
      (
        "StructureCount(ignore h=True,atom comparison=AtomType,bond comparison=BondOrderOrAromaticWithRingness,"
        "filename=" + benzene_filename + ")"
      );

      // test empty constructor
      descriptor::StructureCount structure_search;

      BCL_ExampleCheck( structure_search.GetAlias(), "StructureCount");

      BCL_ExampleCheck( TestMoleculeWithDescriptor( cyclohexane_filename, has_six_membered_carbon_ring), float( 1.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( cyclohexane_filename, has_cyclohexane), float( 1.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( cyclohexane_filename, has_benzene), float( 0.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( benzene_filename, has_six_membered_carbon_ring), float( 1.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( benzene_filename, has_cyclohexane), float( 0.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( benzene_filename, has_benzene), float( 1.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( taxol_filename, has_six_membered_carbon_ring), float( 5.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( taxol_filename, has_cyclohexane), float( 1.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( taxol_filename, has_benzene), float( 3.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( diazepam_filename, has_six_membered_carbon_ring), float( 2.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( diazepam_filename, has_cyclohexane), float( 0.0));
      BCL_ExampleCheck( TestMoleculeWithDescriptor( diazepam_filename, has_benzene), float( 2.0));

    /////////////////
    // data access //
    /////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( has_six_membered_carbon_ring, descriptor::CheminfoProperty()),
        true,
        "StructureCount I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorStructureCount

  const ExampleClass::EnumType ExampleDescriptorStructureCount::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorStructureCount())
  );

} // namespace bcl
