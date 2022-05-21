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
#include "descriptor/bcl_descriptor_atom_volume.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_volume.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 18, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomVolume :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomVolume *Clone() const
    {
      return new ExampleDescriptorAtomVolume( *this);
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

      // test default constructor
      descriptor::AtomVolume volume_a;

      // test copy constructor
      descriptor::AtomVolume volume_b( volume_a);

      // test Clone
      util::ShPtr< descriptor::AtomVolume> sp_volume_b( volume_b.Clone());

    ////////////////////////////////////////
    //  read Small molecule Surface Areas //
    ////////////////////////////////////////

      // let's get information for small_mol_conformation
      io::IFStream input_sdf;
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "benzene_moe_out.sdf"));

      // load information into small_mol
      chemistry::FragmentComplete small_mol = sdf::FragmentFactory::MakeFragment( input_sdf);

      // close the input stream
      io::File::CloseClearFStream( input_sdf);

      // test constructor from atom property
      descriptor::AtomVolume volume_c
      (
        descriptor::GetCheminfoProperties().calc_VDWaalsRadius,
        descriptor::GetCheminfoProperties().calc_CovalentRadius
      );

      // calculate vdwsurfacearea for small_mol
      linal::Vector< float> vdw_sa( volume_c.CollectValuesOnEachElementOfObject( small_mol));

      // print small molecule data
      BCL_MessageStd( "number of atoms: " + util::Format()( small_mol.GetNumberAtoms()));
      BCL_MessageStd
      (
        "vdw surface area vector: "
        + util::Format()( vdw_sa)
      );

      //Check surface area of first atom in vector for accuracy
      BCL_ExampleIndirectCheckWithinTolerance
      (
        volume_c.CollectValuesOnEachElementOfObject( small_mol),
        descriptor::GetCheminfoProperties().calc_VdwVolume->CollectValuesOnEachElementOfObject( small_mol),
        0.10,
        "Test that the estimated surface area is within 10% of the real surface area"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  };

  const ExampleClass::EnumType ExampleDescriptorAtomVolume::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomVolume())
  );
} // namespace bcl
