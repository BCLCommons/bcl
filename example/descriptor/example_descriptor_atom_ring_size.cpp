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
#include "descriptor/bcl_descriptor_atom_ring_size.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_ring_size.cpp
  //!
  //! @author mendenjl
  //! @date Apr 17, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomRingSize :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomRingSize *Clone() const
    {
      return new ExampleDescriptorAtomRingSize( *this);
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

      // default constructor
      descriptor::AtomRingSize max_ring_size( true), min_ring_size( false);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( max_ring_size.GetAlias(), "AtomMaxRingSize");
      BCL_ExampleCheck( min_ring_size.GetAlias(), "AtomMinRingSize");

    ///////////////
    // operators //
    ///////////////

      // map from filename to # of expected rings (first) and expected aromatic rings (second)
      storage::Map< std::string, storage::Pair< size_t, size_t> > filenames_to_expected_min_max_ring_size;

      // diazepam has 3 rings, 2 of which are aromatic
      filenames_to_expected_min_max_ring_size[ "diazepam.sdf"] = storage::Pair< size_t, size_t>( 7, 6);
      // taxol has 7 rings, 3 of which are aromatic
      filenames_to_expected_min_max_ring_size[ "taxol.sdf"] = storage::Pair< size_t, size_t>( 8, 4);
      // hexane has no rings
      filenames_to_expected_min_max_ring_size[ "hexane.sdf"]
        = storage::Pair< size_t, size_t>
          (
            size_t( descriptor::AtomRingSize::s_ChainAtomsMaxRingSize),
            size_t( descriptor::AtomRingSize::s_ChainAtomsMinRingSize)
          );

      for
      (
        storage::Map< std::string, storage::Pair< size_t, size_t> >::const_iterator
          itr_files( filenames_to_expected_min_max_ring_size.Begin()),
          itr_files_end( filenames_to_expected_min_max_ring_size.End());
        itr_files != itr_files_end;
        ++itr_files
      )
      {
        // create input stream for reading a smallmolecule ensemble
        io::IFStream input;
        BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, itr_files->first));
        // read in ensemble
        chemistry::FragmentEnsemble ensemble( input);
        // close stream
        io::File::CloseClearFStream( input);

        const chemistry::FragmentComplete &mol( ensemble.GetMolecules().FirstElement());

        BCL_ExampleIndirectCheck
        (
          max_ring_size.CollectValuesOnEachElementOfObject( mol).Max(),
          float( itr_files->second.First()),
          "Max ring size in " + itr_files->first
        );
        BCL_ExampleIndirectCheck
        (
          min_ring_size.CollectValuesOnEachElementOfObject( mol).Min(),
          float( itr_files->second.Second()),
          "Min ring size in " + itr_files->first
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomRingSize

  const ExampleClass::EnumType ExampleDescriptorAtomRingSize::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomRingSize())
  );

} // namespace bcl
