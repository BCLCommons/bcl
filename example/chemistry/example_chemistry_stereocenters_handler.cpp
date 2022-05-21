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
#include "chemistry/bcl_chemistry_stereocenters_handler.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_stereocenters_handler.cpp
  //!
  //! @author sliwosgr
  //! @date   9/08/2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //! @remarks reviewed by mendenjl, sliwosgr on Nov 08, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryStereocentersHandler :
    public ExampleInterface
  {
  public:

    ExampleChemistryStereocentersHandler *Clone() const
    {
      return new ExampleChemistryStereocentersHandler( *this);
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
      chemistry::StereocentersHandler stereocenters_a;

      // test Clone
      util::ShPtr< chemistry::StereocentersHandler> stereocenters_b( stereocenters_a.Clone());

    ////////////////////////////////////////
    //  read Small molecule Stereocenters //
    ////////////////////////////////////////

      // use SmallMoleculeFactory to load molecule information into SmallMolecule
      // this requires SmallMolecule and input stream
      io::IFStream input_sdf;
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "corina_taxol.sdf"));
      chemistry::FragmentComplete taxol = sdf::FragmentFactory::MakeFragment( input_sdf);
      io::File::CloseClearFStream( input_sdf);

      chemistry::StereocentersHandler at_properties;

      //Test GetConnectedByPriority with all atoms of molecule
      storage::Map< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> indices_map;
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr( taxol.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        indices_map[ *itr] = itr.GetPosition();
      }

      storage::Vector< float> indexes;
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
          itr_atoms( taxol.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        util::SiPtrVector< const chemistry::AtomConformationalInterface> atoms_by_priority
        (
          at_properties.GetUniqueConnectedSubstituentsByPriority( *itr_atoms)
        );
        for
        (
          util::SiPtrVector< const chemistry::AtomConformationalInterface>::const_iterator
            itr_abp( atoms_by_priority.Begin()),
            itr_abp_end( atoms_by_priority.End());
          itr_abp != itr_abp_end;
          ++itr_abp
        )
        {
          indexes.PushBack( indices_map[ *itr_abp]);
        }
      }

      std::string stored_priority_string( "BclCalculatedPriority");
      linal::Vector< float> calculated_priorities( indexes);
      linal::Vector< float> stored_priorities( taxol.GetMDLPropertyAsVector( stored_priority_string));
      BCL_ExampleCheck( calculated_priorities, stored_priorities);

      //calculate stereocenters for taxol
      linal::Vector< float> calculated_stereocenters( at_properties.CalculateFromConformation( taxol.GetAtomsIterator()));

      float expected_vec_array [] =
      {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,-1,1,-1,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0
      };
      linal::Vector< float> expected_vector( 113, expected_vec_array);
      BCL_ExampleCheckWithinAbsTolerance( calculated_stereocenters, expected_vector, 0.001);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  };

  const ExampleClass::EnumType ExampleChemistryStereocentersHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryStereocentersHandler())
  );
} // namespace bcl
